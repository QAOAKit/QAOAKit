include'QAOA_parameters_mod.f90'
include'QAOA_subroutines_mod.f90'
include"BFGS_mod.f90"
Program QAOA
use parameters
use QAOA_subroutines
use BFGS
implicit none

integer :: i,j,k,x,y,z,n,m,its,q,w
integer :: iter,even,odd
double precision :: fret
double precision :: timef,time0
double precision :: C,p_cz_max,C_test,angles_test(2*p_max),C_tester
double precision :: C_max_c, p_max_c, beta_max_c(p_max), gamma_max_c(p_max)
double precision :: C_max_p, p_max_p, beta_max_p(p_max), gamma_max_p(p_max)
double precision :: prob_cmax
double precision :: angles(2*p_max)
double precision :: tmp,C_graph_avg,C_graph_stdev

character*8 :: fmt= '(I4.4)'
character*5 :: gstring

integer, parameter :: n_weighted_graphs=50
double precision, allocatable :: weights(:,:,:),avg_weights(:,:,:)


logical :: BFGS_not_converged=.true.
integer :: BFGS_loops=0

double precision :: grad(2*p_max), norm
character*1 :: Mixer='X'

call cpu_time(time0)

call system('mkdir '//trim(save_folder))
call system('cp QAOA_weighted_Ruslan.f90 '//trim(save_folder)//'QAOA_weighted_Ruslan.f90')
call system('cp QAOA_parameters_mod.f90 '//trim(save_folder)//'QAOA_parameters_mod.f90')
call system('cp QAOA_subroutines_mod.f90 '//trim(save_folder)//'QAOA_subroutines_mod.f90')
call system('cp BFGS_mod.f90 '//trim(save_folder)//'BFGS_mod.f90')


!set up the basis for calculations
call enumerate_basis
call calc_pair_list_fast

open(1,file=trim(graph_file),status='old')
if (.not. test .and. .not. random_graph) then
	call Count_num_graphs
else
	graph_num_tot=1
endif
close(1)


allocate(C_opt(graph_num_tot),Beta_opt(graph_num_tot,p_max),Gamma_opt(graph_num_tot,p_max), &
		& C0(graph_num_tot),cz_max_save(graph_num_tot),p_C_max(graph_num_tot), &
		& p0_C_max(graph_num_tot),good_loops(graph_num_tot),vertex_degrees(graph_num_tot,0:n_qubits-1), &
		& all_odd_degree(graph_num_tot),all_even_degree(graph_num_tot),Sq_opt(graph_num_tot))

allocate(weights(graph_num_tot,n_weighted_graphs,n_edges_max))

allocate(cz_max_save_weighted(graph_num_tot,n_weighted_graphs), C0_weighted(graph_num_tot,n_weighted_graphs), &
		& C_opt_weighted(graph_num_tot,n_weighted_graphs), p_C_max_weighted(graph_num_tot,n_weighted_graphs), &
		beta_opt_weighted(graph_num_tot,n_weighted_graphs,p_max), gamma_opt_weighted(graph_num_tot,n_weighted_graphs,p_max), &
		p0_C_max_weighted(graph_num_tot,n_weighted_graphs),avg_weights(graph_num_tot,n_weighted_graphs,2))



Sq_opt=log2(2.d0**n_qubits)
P_C_max=0.d0
print *, 'graph_num_tot:', graph_num_tot
good_loops=0
vertex_degrees=0
C_opt=0.d0

C_opt_weighted=0.d0
n_graphs_even_odd=0

call generate_weights(weights,n_weighted_graphs,avg_weights)

open(123,file=trim(save_folder)//'C_graph_avg_BFGS',status='unknown')

do its=1,min_good_loops

	open(1,file=trim(graph_file),status='old')

	graph_num=0
	c_graph_avg=0.d0
	C_graph_stdev=0.d0
	BFGS_loops=BFGS_loops+1

	more_graphs=.true.

	do while (more_graphs)

		graph_num = graph_num+1

		call Read_Adjacency_Matrix			!print *, 'graph_num',graph_num

		if (more_graphs) then

			do weight_num = 1,n_weighted_graphs
					
				if (its .eq. 1) then
					call average_weights(weights,n_weighted_graphs,avg_weights) !(weights,average_weights,n_weighted_graphs)
				endif

				call Calc_cz_vec_weighted(BFGS_loops,weights,n_weighted_graphs)

				call Generate_angles(2*p_max,angles)


				call dfpmin(angles,2*p_max,gtol,iter,fret,QAOA_ExpecC_weighted,dfunc)
				!print *, 'iter:', iter
				C=-fret*cz_max_save_weighted(graph_num,weight_num)

				!print*, 'C,cz_max_save_weighted(graph_num,w):',C,cz_max_save_weighted(graph_num,weight_num)
				do q = 1,p_max
					do while (angles(q) .gt. pi/4.d0) 
						angles(q) = angles(q) - pi/2.d0
					enddo
					do while (angles(q) .lt. -pi/4.d0) 
						angles(q) = angles(q) +pi/2.d0
					enddo 
					do while (angles(q+p_max) .gt. pi) 
						angles(q+p_max) = angles(q+p_max) - 2.d0*pi
					enddo
					do while (angles(q+p_max) .lt. -pi) 
						angles(q+p_max) = angles(q+p_max) +2.d0*pi
					enddo 
					if (angles(1) .gt. 0.d0) angles=-angles
				enddo

				C_tester = -QAOA_ExpecC_weighted(2*p_max,angles)*cz_max_save_weighted(graph_num,weight_num)
				if (dabs(C_tester - C) .gt. 1.d-8) then
					print*, 'error in C_tester! C_tester,C',C_tester,C
					stop
				endif

				if (C .gt. C_opt_weighted(graph_num,weight_num)) then
					C_opt_weighted(graph_num,weight_num) = C
					p_C_max_weighted(graph_num,weight_num) = -QAOA_pmax_weighted(2*p_max,angles)
					Beta_opt_weighted(graph_num,weight_num,1:p_max) = angles(1:p_max)
				    Gamma_opt_weighted(graph_num,weight_num,1:p_max) = angles(p_max+1:2*p_max)
				endif
			enddo
		endif
	enddo

	write(123,*) BFGS_loops, C_graph_avg,C_graph_stdev
	close(1)
	close(11)
enddo

do graph_num = 1,graph_num_tot
	write (gstring,fmt) graph_num
	open(2,file=trim(save_folder)//'QAOA_dat_weighted_'//gstring,status='unknown')
	do w = 1,n_weighted_graphs
		write(2,*) graph_num, w, cz_max_save_weighted(graph_num,w), C0_weighted(graph_num,w), &
			& C_opt_weighted(graph_num,w), p_C_max_weighted(graph_num,w), p_max, &
			& (beta_opt_weighted(graph_num,w,j)/pi,j=1,p_max), (gamma_opt_weighted(graph_num,w,j)/pi,j=1,p_max),&
			& avg_weights(graph_num,w,1),avg_weights(graph_num,w,2)
	enddo
	close(2)
enddo

!open(3,file=trim(save_folder)//'QAOA_dat_opt_p_max',status='unknown')

print *, 'BFGS_loops:', BFGS_loops
call cpu_time(timef)
print*, 'elapsed time:', timef-time0

 close(123)

end program QAOA

