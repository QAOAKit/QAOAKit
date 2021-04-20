include"QAOA_subroutines_UnitaryEv_mod.f90"
include"QAOA_subroutines_graphs_mod.f90"
module QAOA_subroutines
use QAOA_subroutines_UnitaryEv
use QAOA_subroutines_graphs
implicit none

contains


subroutine simplify_angles_even_odd_simple(n_angles,angles,result,graph_number)

	use parameters

	implicit none

	integer :: graph_number
	integer :: n_angles
	double precision :: angles(n_angles),angles_test(n_angles)
	double precision :: result,result_test
	integer :: i,j,k
	logical :: dif_angles,pass

	dif_angles=.false.
	angles_test = angles

	!assumes beta_1 < 0

	!if all even degree, set -pi/2 \leq \gamma \leq pi/2
	if (all_even_degree(graph_number)) then

		do i = 1,n_angles/2
			if (dabs(angles_test(n_angles/2+i)) .gt. pi/2.d0) then
				if (angles_test(n_angles/2+i) .gt. 0.d0) then
					angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
				else
					angles_test(n_angles/2+i) = angles_test(n_angles/2+i) + pi
				endif
				dif_angles=.true.
			endif
		enddo

	!if all odd degree, set -pi/2 \leq \gamma \leq pi/2
	elseif (all_odd_degree(graph_number)) then

		do i = 1,n_angles/2
			!map the gamma angle to have |gamma_1| < \pi/2 while maintaining beta_1 < 0 if (i .eq. 1)
			!this is a conjuction of the odd symemtry and the time reversal symmetry
			if (dabs(angles_test(n_angles/2+i)) .gt. pi/2.d0) then
				if (angles_test(n_angles/2+i) .lt. 0.d0) then
					angles_test(n_angles/2+i) = angles_test(n_angles/2+i) + pi
				else
					angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
				endif
				do j = i,n_angles/2
					angles_test(j) = -angles_test(j)
				enddo
				if (i .eq. 1) angles_test=-angles_test
				dif_angles=.true.
			endif
		enddo

	endif

	if (dif_angles) then
		call Check_new_angle_symmetry(angles_test,result,result_test,pass,graph_number)
		if (pass) then
			angles = angles_test
		else
			print *, 'failed even/odd symmetry'	
			print *, 'result, result_test:',result,result_test
			stop
		endif
	endif

end subroutine simplify_angles_even_odd_simple


subroutine simplify_angles_simple(angles,result,graph_number)

	use parameters
	implicit none

	integer :: graph_number
	integer :: tmp_int
	double precision :: angles(2*p_max),angles_test(2*p_max)
	double precision :: result,result_test
	integer :: i
	logical :: pass

	!put angles into the intervals beta \in [-pi/4,pi/4], gamma \in [-pi,pi]
	angles_test=angles
	do i = 1,p_max
		if ( (angles_test(i) .lt. -0.25d0*pi) .or. (angles_test(i) .gt. 0.25d0*pi) ) then
			if (angles_test(i) .lt. -0.25d0*pi) then
				tmp_int = floor( dabs( (angles_test(i)-0.25*pi)/(pi/2.d0) ) )
				angles_test(i) = angles_test(i) + dble(tmp_int)*(pi/2.d0)
			else
				tmp_int = floor( dabs( (angles_test(i)+0.25*pi)/(pi/2.d0) ) )
				angles_test(i) = angles_test(i) - dble(tmp_int)*(pi/2.d0)
			endif
			if ( (angles_test(i) .lt. -0.25d0*pi) .or. (angles_test(i) .gt. 0.25d0*pi) ) then
				print *, 'error in angles_test!'
				print *, 'angles_test(i)/pi',angles_test(i)/pi
				stop
			endif
		endif
		if ( (angles_test(p_max+i) .lt. -pi) .or. (angles_test(p_max+i) .gt. pi)) then
			if (angles_test(p_max+i) .lt. -pi) then
				tmp_int = floor( dabs( (angles_test(p_max+i)-pi)/(pi*2.d0) ) )
				angles_test(p_max+i) = angles_test(p_max+i) + dble(tmp_int)*(pi*2.d0)
			else
				tmp_int = floor( dabs( (angles_test(p_max+i)+pi)/(pi*2.d0) ) )
				angles_test(p_max+i) = angles_test(p_max+i) - dble(tmp_int)*(pi*2.d0)
			endif
			if ( (angles_test(p_max+i) .lt. -pi) .or. (angles_test(i) .gt. pi) ) then
				print *, 'error in angles_test!'
				print *, 'angles_test(p_max+i)', angles_test(p_max+i)
				stop
			endif
		endif
	enddo

	if (angles_test(1) .gt. 0.d0) angles_test=-angles_test

	call Simplify_angles_even_odd_simple(2*p_max,angles_test,result,graph_number)

	!test these give the same result
	call Check_new_angle_symmetry(angles_test,result,result_test,pass,graph_number)

	if (.not. pass) then
		print *, 'failed to simplify angles'
		print *, 'result,result_test',result,result_test
		print *,'angles, angles_test:',angles,angles_test
		stop 
	else
		angles = angles_test
	endif

end subroutine simplify_angles_simple


subroutine Check_optimum_degeneracies_simple(n_angles,angles,result,graph_number,tid)

	use parameters 

	implicit none

	integer :: graph_number,tid
	integer :: i,j,n_angles
	logical :: degenerate_opt,identical_angles
	double precision :: result, angles(2*p_max)
	integer :: counts_opt,counts_new
	identical_angles=.true.
	degenerate_opt=.false.


	call simplify_angles_simple(angles,result,graph_number)

	if ( (all_even_degree(graph_number) .or. all_odd_degree(graph_number)) .and. &
			& dabs(angles(n_angles/2+1)) .gt. pi/2.d0) then
		print *, 'graph_number, all_even, all_odd:',graph_number, all_even_degree(graph_number), all_odd_degree(graph_number)
		print *, 'angles/pi:',(angles(i)/pi,i=1,n_angles)
		print *, 'angles > pi/2!'
		stop
	endif

	!check for degeneracy
	if (dabs(result - C_opt(graph_number)) .lt. degen_tol) degenerate_opt=.true.
	!if degenerate then try simplifying, check if the angles are the same
	if (degenerate_opt) then
		do i = 1,p_max
			if ( (dabs(angles(i) - beta_opt(graph_number,i)) .gt. 1.d-4) .or. &
					& ( dabs(angles(i+p_max) - gamma_opt(graph_number,i)) .gt. 1.d-4) ) then
				identical_angles=.false.
			endif
		enddo

		if (force_angle_ordering) then
			if ( (.not. identical_angles) ) then
				!choose the angles that best line up with the angle patterns
				!do this by counting how many go into the octant of parameter space 
				!-0.25pi < beta < 0, -pi/2 < gamma < 0
				counts_opt=0
				counts_new=0
				do i = 1,p_max
					if (beta_opt(graph_number,i) .le. 0.d0 .and. gamma_opt(graph_number,i) .le. 0.d0 &
						& .and. gamma_opt(graph_number,i) .ge. -pi/2.d0) then
						counts_opt = counts_opt+1
					endif
					if (angles(i) .le. 0.d0 .and. angles(n_angles/2+i) .le. 0.d0 .and. angles(n_angles/2+i) .ge. -pi/2.d0) then
						counts_new = counts_new+1
					endif
				enddo
				if (counts_new .gt. counts_opt) then
					write(579,*) graph_number,C_opt(graph_number),&
						& (beta_opt(graph_number,j)/pi,j=1,p_max),(gamma_opt(graph_number,j)/pi,j=1,p_max), &
						result, (angles(j)/pi,j=1,2*p_max)
					print *, 'graph_number,C_opt,C',graph_number,C_opt(graph_number),result
					print *, 'degenerate angles/pi:'
					print *, 'counts_old, counts_new:',counts_opt,counts_new
					print *, (beta_opt(graph_number,j)/pi,j=1,p_max),(gamma_opt(graph_number,j)/pi,j=1,p_max)
					print *, (angles(j)/pi,j=1,2*p_max)
					print *, ' '
					beta_opt(graph_number,:) = angles(1:p_max)
					gamma_opt(graph_number,:) = angles(p_max+1:2*p_max)
					C_opt(graph_number)=result
					p_C_max(graph_number) = -QAOA_pmax(2*p_max,angles,graph_number)
				endif
			endif
		endif
	endif

	if (result .gt. (C_opt(graph_number)+convergence_tolerance) ) then

    	C_opt(graph_number) = result
    	Beta_opt(graph_number,1:p_max) = angles(1:p_max)
    	Gamma_opt(graph_number,1:p_max) = angles(p_max+1:2*p_max)
    	p_C_max(graph_number) = -QAOA_pmax(2*p_max,angles,graph_number)
    	good_loops(graph_number) = 0
    else
    	good_loops(graph_number) = good_loops(graph_number)+1
    endif

end subroutine Check_optimum_degeneracies_simple


subroutine Generate_angles(n_angles,angles,graph_number)

	use parameters
	implicit none

	integer :: i
	integer :: n_angles,graph_number
	double precision :: angles(n_angles)
	double precision :: tmp

	do i = 1,p_max
		call random_number(tmp)
		angles(i) = pi/2.d0*(tmp-0.5d0)
		call random_number(tmp)
		angles(p_max+i) = 2.d0*pi*(tmp-0.5d0)
	enddo

end subroutine Generate_angles


subroutine set_angles(x,y,angles)

	!for brute force search of angles:
	!expand x and y in base max_it_gammma and max_it_beta 
	!to map x,y to a p-digit number
	!The p-digit number is then mapped to a point on mesh
	!for the numerical search over angles
	use parameters
	implicit none

	integer :: x,y
	integer :: y0,x0
	integer :: i,j
	double precision :: angles(2*p_max)

	beta=0.d0
	gamma=0.d0
	!write beta(i), gamma(i) on a grid of points
	do i= p_max, 1, -1
		y0=0.d0
		x0=0.d0
		if (i .lt. p_max) then
			do j = p_max, i, -1
				y0 = y0 + beta(j)*(max_it_beta**(j-1))
				x0 = x0 + gamma(j)*(max_it_gamma**(j-1))
			enddo
		endif
		beta(i) = (y-y0)/(max_it_beta**(i-1))
		gamma(i) = ((x-x0)/max_it_gamma**(i-1))
	enddo

	!map the beta(i),gamma(i) grid to the the QAOA angle intervals
	do i = 1,p_max
		beta(i)=pi/2*(dble(beta(i))/dble(max_it_beta-1)-0.5d0)
		gamma(i)=pi*2*(dble(gamma(i))/dble(max_it_gamma-1)-0.5d0)
		angles(i) = beta(i)
		angles(p_max+i) = gamma(i)
	enddo


end subroutine set_angles



SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE init_random_seed


subroutine Calc_cz_vec(BFGS_loops,graph_number)
    !calculate a vector c_z
    !with entries c_z[z] that equal the cost function evaluation for
    !the binary expansion of z
    use parameters
    implicit none
    integer :: z,m
    integer :: BFGS_loops,graph_number
    complex*16 :: state(0:z_max)

    do z =0,2**n_qubits-1
        do m = 0,n_edges(graph_number)-1
            if (basis(z,edges(graph_number,m,0)) .ne. basis(z,edges(graph_number,m,1))) then
            	cz_vec(graph_number,z) = cz_vec(graph_number,z) + 1
            endif
        enddo
    enddo

	if (BFGS_loops .eq. 1) then
		!calculate C0, the <C> for the initial state
		state=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)
		call Calc_expec_C(state,C0(graph_number),p0_C_max(graph_number),graph_number)
		cz_max_save(graph_number) = maxval(cz_vec(graph_number,:))

		call calc_vertex_degrees(graph_number)

		if (all_odd_degree(graph_number)) print *, graph_number,'all odd degree'
		if (all_even_degree(graph_number)) print *, graph_number,'all even degree'

	endif

end subroutine Calc_cz_vec


subroutine Check_new_angle_symmetry(angles,result,result_test,pass,graph_number)

	use parameters
	implicit none

	integer :: graph_number
	double precision :: angles(2*p_max)
	double precision :: result, result_test
	logical :: pass

	result_test = -QAOA_ExpecC(2*p_max,angles,graph_number)*cz_max_save(graph_number)

	if (dabs(result-result_test) .gt. degen_tol) then
		pass = .false.
	else
		pass = .true.
	endif

end subroutine Check_new_angle_symmetry


end module QAOA_subroutines