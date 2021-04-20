module parameters

implicit none

character*200, parameter :: save_folder="Test_Graph6/"
character*200, parameter :: graph_file='Graphs/graph6c.txt'
integer, parameter :: n_qubits=6, p_max=3, smaller_p=1
integer, parameter :: loops_param_array(3) = (/ 50, 100, 500 /)
integer, parameter :: min_good_loops =loops_param_array(p_max)

integer, parameter :: dim=2**n_qubits, z_max=2**n_qubits-1,n_edges_max=n_qubits*(n_qubits-1)/2
double precision, parameter :: pi=dacos(-1.d0)
double precision, parameter :: gtol=1.d-10, degen_tol=1.d-7,convergence_tolerance=1.d-7
logical, parameter :: force_angle_ordering=.true.
integer :: graph_num_tot,graph_num
integer :: degeneracy_cz_max

double precision, parameter :: gamma_max = pi, beta_max=pi/4.d0
integer, parameter :: min_it_gamma=0, max_it_gamma=200!for brute force grid search
integer, parameter :: min_it_beta=0, max_it_beta=50

logical :: more_graphs=.true.

integer :: pairs_list(0:n_qubits-1,0:dim/2-1,0:1),basis(0:dim-1,0:n_qubits-1)

integer, allocatable :: n_angles_ma(:)
double precision :: beta(p_max),gamma(p_max),cz_max
double precision, allocatable :: cz_max_save(:), C0(:), C_opt(:), p_C_max(:), beta_opt(:,:), gamma_opt(:,:),p0_C_max(:)
integer, allocatable :: good_loops(:), vertex_degrees(:,:), edges(:,:,:),n_edges(:),cz_vec(:,:)
logical, allocatable :: all_odd_degree(:),all_even_degree(:)


end module parameters