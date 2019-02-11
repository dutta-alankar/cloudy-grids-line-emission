# cloudy-grids-line-emission
Generates gridded parameters influencing intensity of Line emission in Astrophysical plasmas using CLOUDY

Using LineList.dat, line intensities of the predicted lines can be found depending on different parameters which are gridded.
MPI is used for running each model on the grid in embarassingly parallel
Clouds in Collisional Ionization equilibrium with different densities and temperature is studied.


To compile, paste the .cpp file to where libcloudy.a library is already compiled (using make) and then run,

mpic++ -O3 -o mycloudy_mpi.exe mpi_grid.cpp -L. -lcloudy -lm -std=c++11

CLOUDY version: 17.01
