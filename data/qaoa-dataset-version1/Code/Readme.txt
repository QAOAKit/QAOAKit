The main program is QAOA_BFGS_OMP.f90, the other files have subroutines, parameters, and data to be used by the program.

The code works with gfortran from gcc version 10.2.0 (Homebrew GCC 10.2.0) 

To compile I use "gfortran QAOA_BFGS_OMP.f90 -o q.exe -O3 -fopenmp" then to run "./q.exe"

The program will loop over all connected graphs from a file in the Graphs folder. Each graph will get 50-500 iterations, each iteration does a separate BFGS optimization to optimize <C> from random starting angles.  The program saves the best result from all the optimizations for each graph and outputs results to a file.

The parameters file contains a variety of parameters, most of which don't need to be adjusted.  The most useful parameters are:

save_folder - where to output the save data
graph_file - which graph file to use.  The program will loop over all graphs in the file
n_qubits - number of qubits/vertices in a graph
p_max - the "p" parameter 
loops_param_array - the number of random starts to do with BFGS


The files in the graph folder are generated using the connected graph data from Brendan McKay.  The files list the upper triangle of the adjacency matrix for each non-isomorphic connected graph.
https://users.cecs.anu.edu.au/~bdm/data/graphs.html