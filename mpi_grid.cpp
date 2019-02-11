/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and others.*/
/*Bug fixed and modified for the purpose by Alankar Dutta*/

#include "cddefines.h"
#include "cddrive.h"
#include "version.h"

/* This flag indicates where we are multi-processing.\
 * set with -DMPI on compiler command line
 * if this is true then all output happens at the end,
 * if false then each set of results is printed and flushed when it happens.
 * must be set true on parallel machines, false on pc */
#ifdef MPI 
#	include <mpi.h>
#endif
/*lint -e506  constant value Boolean */


#ifdef _MSC_VER
	/* disable conditional expression is constant */
#	pragma warning( disable : 4127 )
#endif

/* ========================================================================= */

#define TITLE "OIII spectra (quantities in log10)" 

/* this is name of resulting output files, FILENAME.lin will be main results */
#define FILENAME "OIII_spectra" 

/* path to directory where we will write, must end in / or \ for current system
 * final location will be PATHNAME + FILENAME */
/*#define PATHNAME "c:\\projects\\cloudy\\tests\\"*/
#define PATHNAME ""

/* these following pairs make 8x8 square grid */
/* lower and upper bounds on hydrogen density - log cm^-3*/
#define HDENINIT 2.0
#define HDENLIMIT 5.0

/* low and up bounds on Temperature K */
#define TEMPINIT 3.7
#define TEMPLIMIT 5.0

/* the increment on the two axes */
#define INCREMENT_T 0.02
#define INCREMENT_H 0.05

/* this is the number of iterations to do -1 means to convergence, 1 will do 1 iteration and so on*/ 
#define ITERATIONS  -1 

/* if true the only call cdNoExec */
#define VERY_QUICK_MODEL false

/* if true execution is fast */
#define QUICK_MET true

/* says whether to make the large cloudy output, if true then output created,
 * if false then not output, only the emission line intensities. */
#define DO_PRINT false

/* print last iteration */
#define PRINT_LAST false

/* set no buffering on output - if true then print as it happens*/
#define NO_BUFFERING false

/*print cloudy commands or not*/
#define PRINT_COMMAND false

/* ============================================================================ */

/* this structure will hold results for final printout */
#define NLINTOT 1500 /* the maximum number of lines */

typedef struct //create a custom datatype GRIDTAG (dictionary)
{
	/* save these parameters for each grid point,
	 * density, temperature, and spectrum */
	double hden , temperature , pred[NLINTOT] ;
	/* number of cautions and warnings for this calc, all should be zero */
	int nWarnings , nCautions , nTFail;
	/* flag returned by Cloudy, 0 == OK, 1 for failure */
	int exit_status;
	/* the execution time for this model */
	double etime;
	/* two extra parameters */
	double eden , cool;
} GRIDTAG;

/* this is where we will save results */
#ifdef MPI
/*This creates a custom MPI datatype handler corresponding to GRIDTAG that can be communicated across processes*/
void Build_derived_type(GRIDTAG gridtag, MPI_Datatype* message_type_ptr)
{
  int block_lengths[10];
  MPI_Aint displacements[10];
  MPI_Aint addresses[11];
  MPI_Datatype typelist[10];

  /*First specify the types */
  typelist[0] = MPI_DOUBLE;    /*hden*/
  typelist[1] = MPI_DOUBLE;    /*temperature*/
  typelist[2] = MPI_DOUBLE;    /*lines*/
  typelist[3] = MPI_INT;       /*nWarnings*/
  typelist[4] = MPI_INT;       /*nCautions*/
  typelist[5] = MPI_INT;       /*nTFail*/
  typelist[6] = MPI_INT;       /*exit_status*/
  typelist[7] = MPI_DOUBLE;    /*etime*/
  typelist[8] = MPI_DOUBLE;    /*eden*/
  typelist[9] = MPI_DOUBLE;    /*cool*/

  /* Specify the number of elements of each type */
  block_lengths[0] = block_lengths[1] = block_lengths[3] = 1;
  block_lengths[4] = block_lengths[5] = block_lengths[6] = 1;
  block_lengths[7] = 1;
  /* this one is special since it is NLINTOT long */
  block_lengths[2] = NLINTOT;
  /* two new parameters are each one element long */
  block_lengths[8] = block_lengths[9] = 1;

  /* Calculate the displacement of the members relative to indata */
  MPI_Address( &gridtag,             &addresses[0]);
  MPI_Address(&(gridtag.hden),       &addresses[1]);
  MPI_Address(&(gridtag.temperature),&addresses[2]);
  MPI_Address(&(gridtag.pred),       &addresses[3]);
  MPI_Address(&(gridtag.nWarnings),  &addresses[4]);
  MPI_Address(&(gridtag.nCautions),  &addresses[5]);
  MPI_Address(&(gridtag.nTFail),     &addresses[6]);
  MPI_Address(&(gridtag.exit_status),&addresses[7]);
  MPI_Address(&(gridtag.etime),      &addresses[8]);
  MPI_Address(&(gridtag.eden),       &addresses[9]);
  MPI_Address(&(gridtag.cool),       &addresses[10]);
 
  /* now far each is from the beginning of the structure */
  displacements[0] = addresses[1] - addresses[0];
  displacements[1] = addresses[2] - addresses[0];
  displacements[2] = addresses[3] - addresses[0];
  displacements[3] = addresses[4] - addresses[0];
  displacements[4] = addresses[5] - addresses[0];
  displacements[5] = addresses[6] - addresses[0];
  displacements[6] = addresses[7] - addresses[0];
  displacements[7] = addresses[8] - addresses[0];
  displacements[8] = addresses[9] - addresses[0];
  displacements[9] = addresses[10]- addresses[0];

  /*Create the derived type */
  MPI_Type_struct(10, block_lengths, displacements, typelist,message_type_ptr);
  
  /*Commit it so that it can be used */
  MPI_Type_commit(message_type_ptr);

} /*Build_derived_type */

# endif

int main( int argc, char *argv[] )
{
	DEBUG_ENTRY( "main()" );

	exit_type exit_status = ES_SUCCESS;

	try {
		bool lgAbort;
		int IntegerMod;
		int  myrank = 0; //will be changed if MPI is used
		double hden , absolute , relative , temperature;
		char chVer[10] , chPath[1000] , chFilename[1000] , chFile[2000] ;

		double *xpar , *ypar; //, *eden, *cool; Let hden be the x parameter (xpar) on the grid and temperature be the y parameter (ypar)

		/* following will become array of line wavelengths and the number of lines 
		 * in this array */
		long int mod, nLines, LimModels, nModels;

		/* these will be passed to cdGetLineList and will become arrays of
		 * labels and wavelengths */
		vector<string> chLabelString;
		vector<realnum> wl;

		long int 
			NumberWarnings1, NumberCautions1, NumberNotes1, NumberSurprises1, NumberTempFailures1,
			NumberPresFailures1, NumberIonFailures1, NumberNeFailures1 ;
		/* number of time Cloudy returned with error condition */
		int nErrorExits;

		FILE *ioDATA ;
		char chLine[100]; //This will store the commands to be executed

		long int n;
		/* number of processors, =1 for single non-MPI */
		int numprocs = 1;

		/* this will hold a single grid point */
		GRIDTAG grid;
		grid.exit_status = ES_SUCCESS;

		/* start MPI if -DMPI on command line */
#		ifdef MPI
		GRIDTAG gridtag; //element of custom datatype GRIDTAG
		GRIDTAG *grids; //pointer to custom datatype GRIDTAG
		
		MPI_Datatype message_type; /* This is the custom dictionary that will be communicated */

		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

		Build_derived_type(gridtag, &message_type); //create the MPI communicator for the dictionary

#		endif

		/* the following is the path for the directory where the output should appear,
		 * it needs to end with the OS standard directory delimiter, a "/" on unix */
		strcpy( chPath , PATHNAME );

		/* this is the first part of the name of the resulting data file.  The final
		 * file will have this name with a ".lin" after it. */
		strcpy( chFilename , FILENAME );

		/* ====================================================================== */
		/*
		 * make arrays of parameters for this grid */

		/* total number of models */
		LimModels = (long)((TEMPLIMIT - TEMPINIT)/INCREMENT_T + 1.1); //gets round up to nearest integer by adding 1.1
		LimModels *= (long)((HDENLIMIT - HDENINIT)/INCREMENT_H+1.1); //gets round up to nearest integer by adding 1.1
		/* add number of processors since this is most we would possibly have to 
		 * add on to make integer number of mpi calls */
		LimModels += numprocs; //-1

		xpar = (double*)malloc(LimModels*sizeof(double) );
		ypar = (double*)malloc(LimModels*sizeof(double) );
		//eden = (double*)malloc(LimModels*sizeof(double) );
		//cool = (double*)malloc(LimModels*sizeof(double) );

		if( xpar==NULL || ypar==NULL )
		{
			fprintf(stderr,"malloc failed! Couldn't allocate space required!\n");
#			ifdef MPI 
			MPI_Finalize(); 
#			endif
			cdEXIT(EXIT_FAILURE);
		}
		/* now define set of parameters */
		hden = HDENINIT;
		temperature = TEMPINIT;
		nModels = 0;
		nErrorExits = 0;

		/* on scalar machines do models in xy order */
		for( hden = HDENINIT; hden<= 1.00001 * HDENLIMIT; hden+=INCREMENT_H )
		{
			for( temperature=TEMPINIT; temperature<=TEMPLIMIT;  temperature+= INCREMENT_T )
			{
				xpar[nModels] = hden;
				ypar[nModels] = temperature;
				//par1[nModels] = COLDEN;
				//par2[nModels] = 0.;
				++nModels;
			}

		}

		/* increase total number to an integer multiple of number of processors */
		if( nModels%numprocs!=0 )
		{
			/* add extra bit */
			IntegerMod = nModels + (numprocs - nModels%numprocs);
		}
		else
		{
			IntegerMod = nModels;
		}

		/* make up data for the "extra" models */
		for( mod=nModels; mod<IntegerMod; ++mod)
		{
			xpar[mod] = xpar[0]; //redundant repetitions
			ypar[mod] = ypar[0]; //redundant repetitions
			//eden[mod] = eden[0];
			//cool[mod] = cool[0];
		}

		/* ====================================================================== */

#		ifdef MPI
		/* allocate memory for grids */	
		if( (grids = ((GRIDTAG *)malloc( (size_t)(numprocs)*sizeof(GRIDTAG )))) == NULL )
		{
			fprintf(stderr," main could not malloc grid\n");
			MPI_Finalize(); 
			cdEXIT(EXIT_FAILURE);
		}
#		endif

		/* initialize the code so that cdGetLineList can be called */
		cdInit();
		/* get list of lines from standard data file.  With no arguments this reads
		 * the file BLRLineList, part of the standard data distribution.  There could
		 * have been another file name within the quotes - that file would have been
		 * used instead */
		nLines = cdGetLineList("LineList.dat", chLabelString, wl);
		//Convert to vector<char*> from vector<string> due to legacy issue
		std::vector<char*> chLabel;
		for (const auto &str : chLabelString) 
		{
			char *charStr = new char[str.size() + 1];
			std::strcpy(charStr, str.c_str());
			chLabel.push_back(charStr); 
		}

		/* now copy version into string to print */
		cdVersion(chVer);

		/* open file for large output if this is desired */
		if( DO_PRINT )
		{
			/* now create final filename */
			strcpy( chFile , chPath );
			strcat( chFile , chFilename );
			strcat( chFile , ".out" );
			cdOutput( chFile );
		}
		else
		{
			/* there should be nothing coming here, but something might */
			ioQQQ = stderr;
		}

		/* calculation's results, all the lines in file ending .lin */
		strcpy( chFile , chPath );
		strcat( chFile , chFilename );
		strcat( chFile , ".lin" );
		ioDATA = open_data(chFile,"w");

		/* print the header information before we begin the calculation */
		/* print version number on both headers */
		fprintf(ioDATA,"Cloudy version %s, %s\n",chVer,TITLE);

		for( mod=myrank; mod<IntegerMod; mod+=numprocs )
		{
			hden = xpar[mod];
			temperature = ypar[mod];
			/* zero out grid, which will hold all info for a single grid point */
			memset(&grid, 0, sizeof(GRIDTAG));

			/* initialize the code for this run */
			cdInit();

			/* send flag saying whether to create large output */
			cdTalk( DO_PRINT );

			/* title for grid */
			sprintf(chLine,"title %s", TITLE);
			n = cdRead( chLine  );

			/* print only last iteration? */
#			if PRINT_LAST
			n = cdRead("print last iteration ");
#			endif

#			if (VERY_QUICK_MODEL ||(ITERATIONS==0))
			(if myrank==0) printf("Nothing executed!");
			cdNoExec();
#			endif

			/* iterate if not fast model */
#			if ( (!VERY_QUICK_MODEL) &&(ITERATIONS!=0) )

			/* Temperature and density */
			sprintf(chLine,"constant temperature %f", temperature);
			n = cdRead( chLine  );

			sprintf(chLine,"hden %f", hden);
			n = cdRead( chLine  );

			n = cdRead( "stop zone 1 "  );

#			if QUICK_MET
			n = cdRead( "init file \"fast.ini\" "  );
#			endif

#      			if ITERATIONS < 0
			n = cdRead( "iterate convergence "  );/**/
#			else 
			sprintf(chLine,"iterate %i", ITERATIONS);
			n = cdRead( chLine  );
#			endif
#			endif

#			ifdef MPI 
			/* set flag so that exit handler will call MPI_Finalize */
			cpu.set_MPI();
#			endif

			/* option to turn off file buffering if code might crash */
#			ifdef NO_BUFFERING
			n = cdRead("no buffering ");
#			endif

			/* if we have run off end of array of models say not to execute,
			 * still must do something on this processor so that the gather can occur
			 * at the bottom of the loop */
			if( mod >= nModels )
				cdNoExec();

			n = cdRead( "normalize to \"O  3\" 5006.84"  );

			/* this is the continuum we used in the large apjs paper */
			//n = cdRead( "agn 6.00 -1.40 -0.50 -1.0  "  );/**/

			/* a bare powerlaw with a break at 1 microns */
			/*n = cdRead( "table powerlaw -1.0 1 micron break "  );*/

			/* the much maligned Mathews&Ferland continuum */
			/*n = cdRead( "table agn  "  );*/

			/*n = cdRead( "blackbody 142,000  "  );*/

			/* broken power law used in Fred's paper */
			/*n = cdRead( "interpolate (0.00000001, -14.899), (0.0091126732, 0.000)  "  );
			  n = cdRead( "continue (1.00, -1.2242), (73.54, -4.2106)  "  );
			  n = cdRead( "continue (3675, -5.7395), (7354000, -11.2526)  "  );*/

			/* broken power law used in Mark's paper */
			/*n = cdRead( "interpolate (10.517, -30.850) (13.383, -23.684)   "  );
			  n = cdRead( "continue (16.684, -26.986) (17.383, -27.9996)    "  );
			  n = cdRead( "continue (18.383, -28.8899) (19.124, -29.268)    "  );
			  n = cdRead( "continue (22.383, -35.788)  "  );*/

			/* add cosmic IR background as well */
			//n = cdRead( "background z=0  "  );
			
			n = cdRead( "blackbody 40000"  );
			n = cdRead( "ionization parameter -4"  );
			
			//sprintf(chLine,"stop column density %f ", par1[mod]);
			//n = cdRead( chLine  );/**/

			//n = cdRead( "failures 5 "  );

			/* options to change abundances */
#			if 0
			n = cdRead( "metals 10 "  );/**/
			n = cdRead( "element nitrogen scale 10 "  );/**/
			n = cdRead( "element helium scale 1.66 "  );/**/
			n = cdRead( "metals 20 "  );/**/
			n = cdRead( "element nitrogen scale 20 "  );/**/
			n = cdRead( "element helium scale 2.41 "  );/**/
#			endif

			/*sprintf(chLine,"turbulence %f km/sec", vturb);
			  n = cdRead( chLine  );*/

			/* all line intensities will be relative to incident continuum near Lya */
			//n = cdRead( "normalize to label \"inci\" 1215 "  );

			/* actually call the code, true indicates error exit 
			 * this is always an abort */
			if( cdDrive() ) //execute cloudy
			{
				grid.exit_status = ES_FAILURE;
				nErrorExits++;
			}

			/* print header information for very first model only */
			if( mod==0 )
			{

				/* for very first model (mod==0) send input stream to data file, followed
				 * by keyword showing end and finally set of tables */
#				if PRINT_COMMAND				
				cdPrintCommands( ioDATA );
#				endif
				/* all the printout happens here */
				/* print header for the data file */
				//fprintf(ioDATA,"abort\twarn\tExecT\tdensity\ttemperature\teden\tcool\t");
				fprintf(ioDATA,"density\ttemperature\teden\tcool\t");
				for( n=0; n<nLines; ++n )
				{
					fprintf(ioDATA,"%4s (",chLabel[n]);
					fprintf(ioDATA,"%f A)",wl[n]);
					fprintf(ioDATA,"\t");
				}
				fprintf(ioDATA,"\n");
			}
			//MPI_Barrier(MPI_COMM_WORLD);
			/* flush this output */
			if( DO_PRINT )
				fflush(ioQQQ);

			/* keep track of all comments about the calculation */
			cdNwcns( 
				&lgAbort,
				&NumberWarnings1 ,   &NumberCautions1,     &NumberNotes1, 
				&NumberSurprises1,   &NumberTempFailures1, &NumberPresFailures1,
				&NumberIonFailures1, &NumberNeFailures1 );
			/* save abort status in exit_status */
			if( lgAbort )
				grid.exit_status = ES_FAILURE;
			grid.nWarnings = NumberWarnings1;
			grid.nCautions = NumberCautions1;
			grid.nTFail = NumberTempFailures1;

			/* print the exec time */
			grid.etime = cdExecTime();

			/* for very quick MPI run, remember processor number */
#			if (VERY_QUICK_MODEL && MPI)
			grid.nWarnings = myrank;
			grid.nCautions = mod;
#			endif

			/* save grid parameters for this model */
			grid.hden = hden;
			grid.temperature = temperature;
			grid.eden = log10(cdEDEN_last());
			grid.cool = log10(cdCooling_last());

			/* if we have run off end of array of models do not try to pull out
			 * results since they do not exist*/
			if( mod < nModels )
			{
				/* print line intensiies */
#				if !VERY_QUICK_MODEL
				for( n=0; n<nLines; ++n )
				{
					if( (cdLine( chLabel[n], wl[n] , &relative , &absolute ))<=0 )
					{
						fprintf(stderr,"did not find %4s",chLabel[n]);
						fprintf(stderr,"%f",wl[n]);
						fprintf(stderr,"\n");
						fprintf(ioDATA,"\ndid not find %4s (",chLabel[n]);
						fprintf(ioDATA,"%f A)",wl[n]);
						fprintf(ioDATA,"\n");
					}
					grid.pred[n] = log10(MAX2(1e-30,relative) ); //prevent zero values which can result in trouble while finding line ratios
				}

#				endif
			}

			/* print the results here if not MPI */
#			ifndef MPI
			/* print exec time and numbers of problems */
			//fprintf(ioDATA,"%i\t%i\t%g\t" , grid.exit_status,	grid.nWarnings , grid.etime);

			/* print grid parameters for this model */
			fprintf(ioDATA,"%.3f\t%.3f\t%.3f\t%.3f\t",grid.hden, temperature,
				grid.eden , grid.cool );

#			if !VERY_QUICK_MODEL
			/* print main set of line intensiies */
			for( n=0; n<nLines; ++n )
			{
				fprintf(ioDATA,"%.4f\t",grid.pred[n] );
			}

#			endif
			fprintf(ioDATA,"\n");
			/* flush this output */
			fflush(ioDATA );
#			endif

#			ifdef MPI
			//MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(&grid,1,message_type,grids,1,message_type,0,MPI_COMM_WORLD);

			if(myrank == 0)
			{
				{
					long int i , nmax ;
					/* on last loop can't expect to get numprocs model results */
					nmax = MIN2( numprocs, nModels-mod);
					for( i=0; i<nmax; ++i )
					{
						/* print exec time and numbers of problems */
						//fprintf(ioDATA,"%i\t%i\t%g\t",grids[i].exit_status, 
						//	grids[i].nWarnings,grids[i].etime );
		
						/* print grid parameters for this model */
						fprintf(ioDATA,
							"%7.3f\t%7.3f\t%7.3f\t%7.3f\t",
							grids[i].hden, grids[i].temperature ,
							grids[i].eden , grids[i].cool  );
#						if !VERY_QUICK_MODEL
						/* print line intensiies */
						for( n=0; n<nLines; ++n )
						{
							fprintf(ioDATA,"%.4f\t",grids[i].pred[n] );
						}

#						endif
						fprintf(ioDATA,"\n");
						/* flush this output */
						fflush(ioDATA );
					}
				}
			}
#			endif
		}
	
		free(xpar );
		free(ypar );
		//free(eden );
		//free(cool );

		/* call MPI_Finalize if MPI is set */
#		ifdef MPI 
		free(grids);
		/* only print message if errors occurred */
		if( nErrorExits )
		{
			fprintf( stderr, "processor %i main exit %i aborts\n" ,myrank, nErrorExits );
			if( DO_PRINT )
			{
				fprintf( ioQQQ, "processor %i main exit %i aborts\n" ,myrank, nErrorExits );
			}
		}
		MPI_Finalize(); 
#		else
		fprintf( stderr, "exit in main with %i aborted models\n" , nErrorExits );
		if( DO_PRINT )
		{
			fprintf( ioQQQ, "exit in main with %i aborted models\n" , nErrorExits );
		}
#		endif

		cdEXIT(exit_type(grid.exit_status));
	}
	catch( bad_alloc )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer "
			 "ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" );
		exit_status = ES_BAD_ALLOC;
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_OUT_OF_RANGE;
	}
	catch( domain_error& e )
	{
		fprintf( ioQQQ, " DISASTER - A vectorized math routine threw a domain_error. Bailing out...\n" );
		fprintf( ioQQQ, " What() = %s", e.what() );
		exit_status = ES_DOMAIN_ERROR;
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		exit_status = ES_BAD_ASSERT;
	}
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT || e.sig() == SIGQUIT )
			{
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
				exit_status = ES_USER_INTERRUPT;
			}
			else if( e.sig() == SIGTERM )
			{
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
				exit_status = ES_TERMINATION_REQUEST;
			}
			else if( e.sig() == SIGILL )
			{
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
				exit_status = ES_ILLEGAL_INSTRUCTION;
			}
			else if( e.sig() == SIGFPE )
			{
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
				exit_status = ES_FP_EXCEPTION;
			}
			else if( e.sig() == SIGSEGV )
			{
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
				exit_status = ES_SEGFAULT;
			}
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
			{
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
				exit_status = ES_BUS_ERROR;
			}
#			endif
			else
			{
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
				exit_status = ES_UNKNOWN_SIGNAL;
			}

		}
	}
#endif
	catch( cloudy_exit& e )
	{
		if( ioQQQ != NULL )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}

	cdPrepareExit(exit_status);

	return exit_status;
}
