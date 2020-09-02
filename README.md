# spe10refine

Depends on INMOST (www.inmost.org). Before building this project, build INMOST with COMPILE_EXAMPLES flag enabled. Otherwise it will not provide refinement library.

Download SPE10 data from https://www.spe.org/web/csp/datasets/set02.htm#download and put dat files into the folder with executable.

Run with mpiexec --oversubscribe -np 8 ./spe10grid 0 0 2 20 40 80 120 20 60.
This will take a chunk of SPE10 data, refine it two times and save grid.pvtk file along with 8 vtk files.

Currently the tool is very slow due to general refinement algorithm.

