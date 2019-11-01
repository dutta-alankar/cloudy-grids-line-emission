# cloudy-grids-line-emission
Generates intensity of Line emission in Astrophysical plasmas as a function of Hydrogen Density and Temperature using CLOUDY

CLOUDY version: 17.01

Using LineList.dat, line intensities of the predicted lines can be found depending on different parameters which are gridded.
MPI is used for running each model on the grid in embarassingly parallel.

Version 17.01 of CLOUDY is available here: 
https://drive.google.com/file/d/1oroKuFMu3Mu1JcKFj-JlNMAxiRhJze6D/view?usp=sharing

Dependency: Any MPI enabled C++ compiler that supports std=c++14 flag (MPI based on c++14 compliant compiler e.g gcc-6.x)

Please be sure that CLOUDY_DATA_PATH environment variable is exported to the Cloudy database location

To compile the code in the repo, do the following:

- Download and untar CLOUDY
- Under the ```sources``` folder, paste the ```sys_mpi_gcc_custom``` folder provided in the repo
- Inside the ```sys_mpi_gcc_custom``` folder type ```make -j <n>``` where ```<n>``` is the number of threads you wish to run concurrently during compilation
- Copy all the files (exclude the folders) from ```sources``` folder and paste it in ```sys_mpi_gcc_custom``` folder
- Also put the ```mpi_grid.cpp``` file provided in the repo in ```sys_mpi_gcc_custom``` folder
- Now inside sys_mpi_gcc_custom folder type 

```mpic++ -O3 -o mycloudy_mpi.exe mpi_grid.cpp -L. -lcloudy -lm -std=c++14 -DMPI_ENABLED```
- If compilation succeeds, ```mycloudy_mpi.exe``` binary will be created
- Copy ```mycloudy_mpi.exe``` and paste it in a folder that also has the ```LineList.dat``` file provided in the repo
- Now from the folder run ```mpiexec -n <n> ./mycloudy_mpi.exe``` where ```<n>``` is the number of parallel processors you wish to deploy

Note that you can also run the code in serial without MPI in any c++14 compliant compiler. 

For example, using g++

```g++ -O3 -o mycloudy_serial.exe mpi_grid.cpp -L. -lcloudy -lm -std=c++14```

You can also look up the comprehensive CLOUDY manual (provided with the package) for any customization of your own.
