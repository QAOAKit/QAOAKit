
module parameters

implicit none

character*200 :: save_folder='test_weighted/n=5/p=2/'
!character*200 :: save_folder='GHZ/mean/n=8/'
!character*200 :: save_folder='test_GTRI_3_widerintervals/'
!character*200 :: save_folder='Test/Karloff_n=20_p=3_easyangles_noopt/'
!character*200 :: save_folder='even_odd_only/n=7/Full_Generate_Angle_range/n=7_p=3_not_even_odd_2/'
!character*200, parameter ::  save_folder="/Users/5ci/desktop/QAOA/OMP_simpleangles_final/n=2_p=1/"
!character*200, parameter :: save_folder="/Users/5ci/desktop/QAOA/BFGS_TestConvergence_withstdev/graph3c_p=3/"
!character*200, parameter :: save_folder='/Users/5ci/desktop/FinalDataset/graph8c_p=1_2 copy/'
!character*200 :: save_folder='Test/'
character*200, parameter :: graph_file='/Users/5ci/Desktop/Graphs/graph5c.txt'
character*200, parameter :: smaller_p_file='Graphs7_QAOA_p=3_BFGS/QAOA_dat'
character*200, parameter :: old_angles_file='Graphs7_QAOA_p=4_BFGS_test/QAOA_dat'

integer, parameter :: n_qubits=5
integer, parameter :: p_max=2, smaller_p=1

integer, parameter :: min_good_loops =200
integer, parameter :: graph_of_interest=1

integer, parameter :: dim=2**n_qubits, n_edges_max=n_qubits*(n_qubits-1)/2
integer, parameter :: min_it_gamma=0, max_it_gamma=101*40
integer, parameter :: min_it_beta=0, max_it_beta=26*2
double precision, parameter :: pi=dacos(-1.d0)
double precision, parameter :: gtol=1.d-14, degen_tol=1.d-8,convergence_tolerance=1.d-12
logical, parameter :: random_graph=.false.
logical, parameter :: angles_from_smaller_p=.false.
logical, parameter :: old_angles=.false.,directed_angles=.false.
logical, parameter :: force_angle_ordering=.true.
logical, parameter :: even_odd_only=.false., even_odd_angle_dist=.false., not_even_odd_only=.true.
logical, parameter :: no_opt=.true.
logical, parameter :: Optimize_expec_C=.true.,Optimize_p_Cmax=.false.,Optimize_Sq=.false.
logical, parameter :: fix_sign_beta1=.false.,fix_sign_gamma1=.false.,fix_sign_positive=.false.
logical, parameter :: single_graph=.false.
integer, parameter :: single_graph_num = 70
integer :: graph_num_tot,graph_num
integer :: degeneracy_cz_max
logical, parameter :: many_angles=.false.
logical, parameter :: fixed_number_iterations=.true.


double precision, parameter :: gamma_max = pi*1.d0*20.d0, beta_max=pi*1.d0/4.d0*2.d0








double precision, allocatable :: gamma_ma(:), beta_ma(:), gamma_vec_ma(:),gamma_opt_ma(:,:),beta_opt_ma(:,:)
double precision, allocatable :: Sq_opt(:)

double precision :: GHZ_mean_noise=0.d0, GHZ_stdev_noise=0.d0

complex*16 :: psi(0:dim-1) 
complex*16 :: Ub(0:1,0:1)

logical :: more_graphs=.true.





!changed from int!
integer :: cz_vec(0:dim-1)







integer :: pairs_list(0:n_qubits-1,0:dim/2-1,0:1),basis(0:dim-1,0:n_qubits-1)
integer :: n_edges, edges(0:n_edges_max-1,0:1)
integer :: n_graphs_even_odd
!integer :: edges(0:n_edges_max,0:1)
!integer :: n_edges
integer, allocatable :: n_angles_ma(:)
double precision :: beta_a(p_max),gamma_a(p_max),cz_max
double precision, allocatable :: cz_max_save(:), C0(:), C_opt(:), p_C_max(:), beta_opt(:,:), gamma_opt(:,:),p0_C_max(:)
double precision, allocatable :: cz_max_save_weighted(:,:), C0_weighted(:,:), C_opt_weighted(:,:)
double precision, allocatable :: p_C_max_weighted(:,:), beta_opt_weighted(:,:,:), gamma_opt_weighted(:,:,:),p0_C_max_weighted(:,:)
integer :: weight_num

integer, allocatable :: z_max(:)
integer, allocatable :: good_loops(:), vertex_degrees(:,:)
logical, allocatable :: all_odd_degree(:),all_even_degree(:)








logical, parameter :: test_GTRI=.true.








!test graph for comparison with python calculations
logical, parameter :: test=.false.
integer,parameter :: n_edges_test=90
integer :: edges_Test(0:n_edges_test-1,0:1) = &
		 & transpose(reshape((/ 0,7, &	
							  & 0,8, &
							  & 0,9, &      !/),shape(transpose(edges_test))))
							  & 0,13, &
							  & 0,14, &
							  & 0,15, &
							  & 0,16, &
							  & 0,17, &
							  & 0,18, &
							  & 1, 5, &
							  & 1, 6, & 
							  & 1, 9, & 
							  & 1, 11, & 
							  & 1, 12, & 
							  & 1, 15, & 
							  & 1, 16, & 
							  & 1, 17, & 
							  & 1, 19, & 
							  & 2, 4, & 
							  & 2, 6, & 
							  & 2, 8, & 
							  & 2, 10, & 
							  & 2, 12, & 
							  & 2, 14, & 
							  & 2, 16, & 
							  & 2, 18, & 
							  & 2, 19, & 
							  & 3, 4, & 
							  & 3, 5, & 
							  & 3, 7, & 
							  & 3, 10, & 
							  & 3, 11, & 
							  & 3, 13, & 
							  & 3, 17, & 
							  & 3, 18, & 
							  & 3, 19, & 
							  & 4, 9, & 
							  & 4, 11, & 
							  & 4, 12, & 
							  & 4, 13, & 
							  & 4, 14, & 
							  & 4, 18, & 
							  & 4, 19, & 
							  & 5, 8, & 
							  & 5, 10, & 
							  & 5, 12, & 
							  & 5, 13, & 
							  & 5, 15, & 
							  & 5, 17, & 
							  & 5, 19, & 
							  & 6, 7, & 
							  & 6, 10, & 
							  & 6, 11, & 
							  & 6, 14, & 
							  & 6, 15, & 
							  & 6, 16, & 
							  & 6, 19, & 
							  & 7, 10, & 
							  & 7, 11, & 
							  & 7, 14, & 
							  & 7, 15, & 
							  & 7, 17, & 
							  & 7, 18, & 
							  & 8, 10, & 
							  & 8, 12, & 
							  & 8, 13, & 
							  & 8, 15, & 
							  & 8, 16, & 
							  & 8, 18, & 
							  & 9, 11, & 
							  & 9, 12, & 
							  & 9, 13, & 
							  & 9, 14, & 
							  & 9, 16, & 
							  & 9, 17, & 
							  & 10, 15, & 
							  & 10, 18, & 
							  & 10, 19, & 
							  & 11, 14, & 
							  & 11, 17, & 
							  & 11, 19, & 
							  & 12, 13, & 
							  & 12, 16, & 
							  & 12, 19, & 
							  & 13, 17, & 
							  & 13, 18, & 
							  & 14, 16, & 
							  & 14, 18, & 
							  & 15, 16, & 
							  & 15, 17 /),shape(transpose(edges_test))))

!&& transpose(reshape((/ 0,1, 1,2, 2,3, 1,4, 2,4, 1,5, 2,5, 1,6 /),shape(transpose(edges_test))))
end module parameters