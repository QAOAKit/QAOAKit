QAOA_dat_weighted_XXXX

XXXX=graph number

the first column = graph_number
second column = weight_number (which weighted instance out of the ensemble)

proceeding columns are the same as the online dataset - note they are shifted right by one because of the weight_number in the second column. 

The final two columns are the average and standard deviation of the weights for that realization.



weights_XXXX

a list of weights for graph XXXX

the first column is the graph number

The next E columns are the weights of the E edges.

Other columns beyond E+1 are meaningless

These weights were drawn uniformly at random from 1,2