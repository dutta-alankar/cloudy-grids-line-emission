This generates an MPI parallelized executable with mpic++ or mpicxx (based on g++)
for running a large grid.

to build enter
make
at the command prompt, to build debug do
make debug

add "-j n" where n is the number of threads to use in the build.

For Mac, seperate Makefile is provided.

To use the MacPorts openmpi mv Makefile_conf.Mac to Makefile.conf

