include'QAOA_parameters_OMP_mod.f90'
include'QAOA_subroutines_OMP_mod.f90'
include"BFGS_OMP_mod.f90"
Program QAOA
use parameters
use QAOA_subroutines
use BFGS
use omp_lib
implicit none

integer :: i,j,k,x,y,z,n,m,nits
integer :: iter,even,odd,tid
double precision :: fret
double precision :: timef,time0
double precision :: C,p_cz_max,C_test,angles_test(2*p_max)
double precision :: C_max_c, p_max_c, beta_max_c(p_max), gamma_max_c(p_max)
double precision :: prob_cmax
double precision :: angles(2*p_max)
double precision :: tmp,C_graph_avg,C_graph_stdev
double precision :: C_graph_avg_omp, C_graph_stdev_omp

integer :: BFGS_loops=0
integer :: n_threads

double precision :: grad(2*p_max), norm
character*1 :: Mixer='X'

call system('mkdir '//trim(save_folder))
call system('cp QAOA_BFGS_OMP.f90 '//trim(save_folder)//'QAOA_BFGS.f90')
call system('cp QAOA_parameters_OMP_mod.f90 '//trim(save_folder)//'QAOA_parameters_mod.f90')
call system('cp QAOA_subroutines_OMP_mod.f90 '//trim(save_folder)//'QAOA_subroutines_mod.f90')
call system('cp BFGS_OMP_mod.f90 '//trim(save_folder)//'BFGS_mod.f90')

n_threads=omp_get_max_threads()
call OMP_SET_NUM_THREADS(n_threads)
print *, 'n_threads:',n_threads
time0=omp_get_wtime()

!set up the basis for calculations
call enumerate_basis
call calc_pair_list_fast
call init_random_seed()

open(1,file=trim(graph_file),status='old')
call Count_num_graphs
close(1)

allocate(C_opt(graph_num_tot),Beta_opt(graph_num_tot,p_max),Gamma_opt(graph_num_tot,p_max), &
		& C0(graph_num_tot),cz_max_save(graph_num_tot),p_C_max(graph_num_tot), &
		& p0_C_max(graph_num_tot),good_loops(graph_num_tot),vertex_degrees(graph_num_tot,0:n_qubits-1), &
		& all_odd_degree(graph_num_tot),all_even_degree(graph_num_tot),edges(graph_num_tot,0:n_edges_max-1,0:1), &
		& n_edges(graph_num_tot),cz_vec(graph_num_tot,0:dim-1))

edges=0
n_edges=0
cz_vec=0
vertex_degrees=0
open(1,file=trim(graph_file),status='old')
call find_all_edges
close(1)


good_loops=0
vertex_degrees=0
C_opt=0.d0

open(123,file=trim(save_folder)//'C_graph_avg_BFGS',status='unknown')
open(579,file=trim(save_folder)//'degeneracies',status='unknown')

do BFGS_loops=1,min_good_loops

	graph_num=0
	c_graph_avg=0.d0
	C_graph_stdev=0.d0
	
	more_graphs=.true.

!$OMP PARALLEL PRIVATE(angles,fret,C,prob_Cmax,C_graph_avg_omp,C_graph_stdev_omp,tid,graph_num)
tid=omp_get_thread_num()
C_graph_avg_omp=0.d0
C_graph_stdev_omp=0.d0
!$OMP DO SCHEDULE(dynamic)
	do graph_num = 1,graph_num_tot
		if (BFGS_loops .eq. 1) call Calc_cz_vec(BFGS_loops,graph_num)

		call Generate_angles(2*p_max,angles,graph_num)

		call dfpmin(angles,2*p_max,gtol,iter,fret,QAOA_ExpecC,dfunc,graph_num)
		C=-fret*cz_max_save(graph_num)

		call Check_optimum_degeneracies_simple(2*p_max,angles,C,graph_num,tid)
		
		if (cz_max_save(graph_num) .gt. 1.d-8) then
			C_graph_avg_omp= C_graph_avg_omp + C_opt(graph_num)/cz_max_save(graph_num)
			C_graph_stdev_omp = C_graph_stdev_omp + ( C_opt(graph_num) / cz_max_save(graph_num) )**2
		endif
	enddo
!$OMP enddo

!$OMP CRITICAL

	C_graph_avg = C_graph_avg + C_graph_avg_omp
	C_graph_stdev = C_graph_stdev + C_graph_stdev_omp
!$OMP END CRITICAL
!$OMP END PARALLEL

	if (mod(BFGS_loops,10) .eq. 0) then
		timef=omp_get_wtime()
		print *, 'BFGS_loops,time:',BFGS_loops,timef-time0
	endif

	C_graph_avg=C_graph_avg/dble(graph_num_tot)
	C_graph_stdev=C_graph_stdev/dble(graph_num_tot)
	C_graph_stdev = dsqrt( c_graph_stdev - ( (c_graph_avg)**2 ) )

	write(123,*) BFGS_loops, C_graph_avg,C_graph_stdev

enddo

open(2,file=trim(save_folder)//'QAOA_dat',status='unknown')
do graph_num = 1,graph_num_tot
	if (all_even_degree(graph_num)) then
		even = 1
	else
		even = 0
	endif
	if (all_odd_degree(graph_num)) then
		odd = 1
	else
		odd = 0
	endif
	write(2,*) graph_num, cz_max_save(graph_num), C0(graph_num), C_opt(graph_num), p_C_max(graph_num), p_max, &
		& (beta_opt(graph_num,j)/pi,j=1,p_max), (gamma_opt(graph_num,j)/pi,j=1,p_max),even, &
		& odd
enddo

close(2)


print *, 'BFGS_loops:', BFGS_loops
timef=omp_get_wtime()
print*, 'elapsed time:', timef-time0

close(123)
close(579)

end program QAOA

