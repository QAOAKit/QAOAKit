module QAOA_subroutines

implicit none

contains


function Morse(x)

implicit none

double precision :: Morse, x(1)
double precision :: D=1.d0, r0=1.d0, a=1.d0

Morse = D*( exp(-2*a*(x(1)-r0)) - 2*exp(-a*(x(1)-r0)) )

end function Morse


function Morse_test(n,x)

implicit none
integer :: n,m
double precision :: Morse_test, x(2*n)
double precision :: D=1.d0, r0=1.d0, a=1.d0

Morse_test = 0.d0
do m = 1,2*n
	Morse_test = Morse_test + D*( exp(-2*a*(x(m)-r0)) - 2*exp(-a*(x(m)-r0)) )
enddo

end function Morse_test



subroutine QAOA_Algorithm_ma
	use parameters
	implicit none

	integer :: j
	integer :: p

	psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)

    do j = 1,p_max
    	call Calc_gamma_vec_ma
    	call Apply_Uc_ma
    	call Apply_Ub_ma
    enddo

end subroutine QAOA_Algorithm_ma


subroutine QAOA_Algorithm(p,B_angles,C_angles)
	use parameters
	implicit none

	integer :: j
	integer :: p 
	double precision :: b_angles(p),c_angles(p)!beta and gamma

	psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)

    do j = 1,p
    	call construct_Ub(b_angles(j))
    	call Apply_Uc(c_angles(j))
    	call Apply_Ub
    enddo

end subroutine QAOA_Algorithm

function QAOA_ExpecC_ma(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles)
	double precision :: p_cz_max,QAOA_ExpecC_ma


	beta_ma = angles(1:n_qubits*p_max)
	gamma_ma = angles(n_qubits*p_max+1:n_angles)

	Call QAOA_Algorithm_ma
	call Calc_expec_C(QAOA_ExpecC_ma,p_cz_max)

	QAOA_ExpecC_ma = -QAOA_ExpecC_ma/cz_max_save(graph_num)

end function QAOA_ExpecC_ma



function QAOA_ExpecC_weighted(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: p_cz_max,QAOA_ExpecC_weighted
	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)
	call Calc_expec_C(QAOA_ExpecC_weighted,p_cz_max)

	QAOA_ExpecC_weighted = -QAOA_ExpecC_weighted/cz_max_save_weighted(graph_num,weight_num)

	!print*, 'psi(0), psi(dim-1):',psi(0),psi(dim-1)

	!print *, 'angles:', angles
	!print *, 'QAOA_ExpecC:', QAOA_ExpecC
end function QAOA_ExpecC_weighted



function QAOA_ExpecC(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: p_cz_max,QAOA_ExpecC
	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)
	call Calc_expec_C(QAOA_ExpecC,p_cz_max)

	QAOA_ExpecC = -QAOA_ExpecC/cz_max_save(graph_num)

	!print*, 'psi(0), psi(dim-1):',psi(0),psi(dim-1)

	!print *, 'angles:', angles
	!print *, 'QAOA_ExpecC:', QAOA_ExpecC
end function QAOA_ExpecC

function QAOA_Sq(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: p_cz_max,tmp,QAOA_Sq
	integer :: z
	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)
	!call Calc_expec_C(QAOA_ExpecC,p_cz_max)

	QAOA_Sq = 0.d0
	do z = 0,2**n_qubits-1
		tmp = dble(psi(z)*dconjg(psi(z)))
		if (tmp .gt. 1.d-12) QAOA_sq = QAOA_sq - tmp*log2(tmp)
	enddo

	!print *, 'QAOA_Sq:',qaoa_sq

	!QAOA_ExpecC = -QAOA_ExpecC/cz_max_save(graph_num)

end function QAOA_Sq



subroutine generate_weights(weights,n_weighted_graphs,avg_weights)
	use parameters
	implicit none

	integer :: n_weighted_graphs
	double precision :: weights(graph_num_tot,n_weighted_graphs,0:n_edges_max-1)
	double precision :: avg_weights(graph_num_tot,n_weighted_graphs,2)
	logical :: uniform_random=.true.
	integer :: i,j,k
	double precision :: tmp
	logical :: two_weights=.true.
	character*8 :: fmt= '(I4.4)'
	character*5 :: gstring

	avg_weights=0.d0
	weights=0.d0
	do i = 1,graph_num_tot
		do j = 1,n_weighted_graphs
			do k = 0,n_edges_max-1
				call random_number(tmp)
				if (two_weights) weights(i,j,k) = nint(tmp+1.d0)
			enddo
		enddo
	enddo

	do i = 1,graph_num_tot
		write (gstring,fmt) i
		open(9753,file=trim(save_folder)//'weights_'//gstring,status='unknown')
		do j = 1,n_weighted_graphs
			write(9753,*) i,(weights(i,j,k),k=0,n_edges_max-1)
		enddo
		close(9753)
	enddo

end subroutine generate_weights



subroutine average_weights(weights,n_weighted_graphs,avg_weights) !(weights,avg_weights,n_weighted_graphs)
	use parameters
	implicit none
	integer :: k,n_weighted_graphs
	double precision :: weights(graph_num_tot,n_weighted_graphs,0:n_edges_max-1)
	double precision :: avg_weights(graph_num_tot,n_weighted_graphs,2)

	do k = 0,n_edges-1
		avg_weights(graph_num,weight_num,1) = avg_weights(graph_num,weight_num,1) + dble(weights(graph_num,weight_num,k))/dble(n_edges)
		avg_weights(graph_num,weight_num,2) = avg_weights(graph_num,weight_num,2) + dble(weights(graph_num,weight_num,k)**2)/dble(n_edges)
	enddo
	avg_weights(graph_num,weight_num,2) = dsqrt(avg_weights(graph_num,weight_num,2) - avg_weights(graph_num,weight_num,1)**2)

end subroutine average_weights



subroutine Check_optimum_degeneracies_simple(n_angles,angles,result)

	use parameters 

	implicit none

	integer :: i,j,n_angles
	logical :: degenerate_opt,identical_angles,all_negative,ordered_mag
	double precision :: result, angles(2*p_max)
	integer :: counts_opt,counts_new
	identical_angles=.true.
	all_negative=.true.
	degenerate_opt=.false.
	ordered_mag=.false.

	!if (graph_num .eq. 20) then 
	!	print *, 'graph_num, angles/pi, C:',graph_num, angles/pi, result
	!endif
	call simplify_angles_simple(angles,result)

	if (.not. many_angles) then
		!check for degeneracy
		if (dabs(result - C_opt(graph_num)) .lt. degen_tol) degenerate_opt=.true.
		!if degenerate then try simplifying, check if the angles are the same
		if (degenerate_opt) then
			do i = 1,p_max
				if ( (dabs(angles(i) - beta_opt(graph_num,i)) .gt. 1.d-4) .or. &
						& ( dabs(angles(i+p_max) - gamma_opt(graph_num,i)) .gt. 1.d-4) ) then
					identical_angles=.false.
				endif
				if ( (angles(i) .gt. 0.d0) .or. (angles(p_max+i) .gt. 0.d0) ) then
					all_negative=.false.
				endif
			enddo
			if (p_max .gt. 1) then
				if (angles(1) .lt. angles(2) .and. angles(3) .gt. angles(4)) then
					ordered_mag=.true.
				endif
			endif

			!if (graph_num .eq. 20) then 
			!	print *, 'graph_num, angles/pi, C:',graph_num, angles/pi, result
			!	if (all_negative) print *, '**********'
			!	stop
			!endif
			if (force_angle_ordering) then
				!if ( (.not. identical_angles) .and. (ordered_mag) .and. all_negative ) then 
				!	print *, 'graph_num,C_opt,C',graph_num,C_opt(graph_num),result
				!	print *, 'degenerate angles/pi:'
				!	print *, (beta_opt(graph_num,j)/pi,j=1,p_max),(gamma_opt(graph_num,j)/pi,j=1,p_max)
				!	print *, (angles(j)/pi,j=1,2*p_max)
				!	print *, ' '
				!	beta_opt(graph_num,:) = angles(1:p_max)
				!	gamma_opt(graph_num,:) = angles(p_max+1:2*p_max)
				!endif
			!else
				if ( (.not. identical_angles) ) then
					!choose the angles that best line up with the angle patterns
					!do this by counting how many go into the octant of parameter space 
					!-0.25pi < beta < 0, -pi/2 < gamma < 0
					counts_opt=0
					counts_new=0
					do i = 1,p_max
						if (beta_opt(graph_num,i) .le. 0.d0 .and. gamma_opt(graph_num,i) .le. 0.d0 &
							& .and. gamma_opt(graph_num,i) .ge. -pi/2.d0) then
							counts_opt = counts_opt+1
						endif
						if (angles(i) .le. 0.d0 .and. angles(n_angles/2+i) .le. 0.d0 .and. angles(n_angles/2+i) .ge. -pi/2.d0) then
							counts_new = counts_new+1
						endif
					enddo
					if (counts_new .gt. counts_opt) then
						print *, 'graph_num,C_opt,C',graph_num,C_opt(graph_num),result
						print *, 'degenerate angles/pi:'
						print *, 'counts_old, counts_new:',counts_opt,counts_new
						print *, (beta_opt(graph_num,j)/pi,j=1,p_max),(gamma_opt(graph_num,j)/pi,j=1,p_max)
						print *, (angles(j)/pi,j=1,2*p_max)
						print *, ' '
						beta_opt(graph_num,:) = angles(1:p_max)
						gamma_opt(graph_num,:) = angles(p_max+1:2*p_max)
					endif
				endif
			endif
		endif
	endif

	if (result .gt. (C_opt(graph_num)+convergence_tolerance) ) then

    	C_opt(graph_num) = result
    	Beta_opt(graph_num,1:p_max) = angles(1:p_max)
    	Gamma_opt(graph_num,1:p_max) = angles(p_max+1:2*p_max)
    	p_C_max(graph_num) = -QAOA_pmax(2*p_max,angles)
    	good_loops(graph_num) = 0
    	!Sq_opt(graph_num) = QAOA_Sq(2*p_max,angles)
    	!print *, 'Sq_opt:',Sq_opt(graph_num)
    else
    	good_loops(graph_num) = good_loops(graph_num)+1
    endif

end subroutine Check_optimum_degeneracies_simple


subroutine Check_optimum_degeneracies(n_angles,angles,result)

	use parameters 

	implicit none

	integer :: i,j,n_angles
	logical :: degenerate_opt,identical_angles,all_negative,ordered_mag
	double precision :: result, angles(2*p_max)

	identical_angles=.true.
	all_negative=.true.
	degenerate_opt=.false.
	ordered_mag=.false.

	!if (graph_num .eq. 20) then 
	!	print *, 'graph_num, angles/pi, C:',graph_num, angles/pi, result
	!endif

	if (.not. many_angles) then
		!check for degeneracy
		if (dabs(result - C_opt(graph_num)) .lt. 1.d-4) degenerate_opt=.true.
		!if degenerate then try simplifying, check if the angles are the same
		if (degenerate_opt) then
			call simplify_angles(angles,result)
			do i = 1,p_max
				if ( (dabs(angles(i) - beta_opt(graph_num,i)) .gt. 1.d-4) .or. &
						& ( dabs(angles(i+p_max) - gamma_opt(graph_num,i)) .gt. 1.d-4) ) then
					identical_angles=.false.
				endif
				if ( (angles(i) .gt. 0.d0) .or. (angles(p_max+i) .gt. 0.d0) ) then
					all_negative=.false.
				endif
			enddo
			if (p_max .gt. 1) then
				if (angles(1) .lt. angles(2) .and. angles(3) .gt. angles(4)) then
					ordered_mag=.true.
				endif
			endif

			!if (graph_num .eq. 20) then 
			!	print *, 'graph_num, angles/pi, C:',graph_num, angles/pi, result
			!	if (all_negative) print *, '**********'
			!	stop
			!endif
			if (force_angle_ordering) then
				if ( (.not. identical_angles) .and. (ordered_mag) .and. all_negative ) then 
					print *, 'graph_num,C_opt,C',graph_num,C_opt(graph_num),result
					print *, 'degenerate angles/pi:'
					print *, (beta_opt(graph_num,j)/pi,j=1,p_max),(gamma_opt(graph_num,j)/pi,j=1,p_max)
					print *, (angles(j)/pi,j=1,2*p_max)
					print *, ' '
					beta_opt(graph_num,:) = angles(1:p_max)
					gamma_opt(graph_num,:) = angles(p_max+1:2*p_max)
				endif
			else
				if ( (.not. identical_angles) ) then
					print *, 'graph_num,C_opt,C',graph_num,C_opt(graph_num),result
					print *, 'degenerate angles/pi:'
					print *, (beta_opt(graph_num,j)/pi,j=1,p_max),(gamma_opt(graph_num,j)/pi,j=1,p_max)
					print *, (angles(j)/pi,j=1,2*p_max)
					print *, ' '
					beta_opt(graph_num,:) = angles(1:p_max)
					gamma_opt(graph_num,:) = angles(p_max+1:2*p_max)
				endif
			endif
		endif
	endif

	if (result .gt. (C_opt(graph_num)+convergence_tolerance) ) then


    	if (.not. many_angles) call simplify_angles(angles,result)
    	!C_tester = -QAOA_ExpecC(2*p_max,angles)*cz_max
    	!print *, 'dabs(C_tester-C):',C,C_tester,dabs(C_tester-C)
    	C_opt(graph_num) = result
    	if (many_angles) then
    	print*, 'cant do many angles anymore!'
    	stop
    	!	beta_opt_ma(graph_num,1:n_qubits) = angles(1:n_qubits)
    	!	gamma_opt_ma(graph_num,1:n_angles_ma(graph_num)-n_qubits) = angles(n_qubits+1:n_angles_ma(graph_num))
    	else
    		Beta_opt(graph_num,1:p_max) = angles(1:p_max)
    		Gamma_opt(graph_num,1:p_max) = angles(p_max+1:2*p_max)
    	endif
    	p_C_max(graph_num) = -QAOA_pmax(2*p_max,angles)
    	good_loops(graph_num) = 0
    else
    	good_loops(graph_num) = good_loops(graph_num)+1
    	!test reflection symmetry for the final result
    	if ( (good_loops(graph_num) .eq. min_good_loops) .and. (.not. many_angles) ) then
    		do i = 1,p_max
    			angles(i) = beta_opt(graph_num,i)
    			angles(p_max+i) = gamma_opt(graph_num,i)
    		enddo
    		call p1_reflection(angles,C_opt(graph_num))
    		do i = 1,p_max
    			beta_opt(graph_num,i)=angles(i)
    			gamma_opt(graph_num,i)=angles(p_max+i)
    		enddo
    	endif
    endif

end subroutine Check_optimum_degeneracies


subroutine Check_optimum(angles,result)

    use parameters
    
    implicit none

	double precision :: result
	double precision :: angles(2*p_max)
	integer :: i

	if (optimize_expec_C) then

	    if (result .gt. (C_opt(graph_num)+convergence_tolerance)) then
	    	!use angles in simple intervals, with gamma_1 < 0
	    	call simplify_angles(angles,result)
	    	!C_tester = -QAOA_ExpecC(2*p_max,angles)*cz_max
	    	!print *, 'dabs(C_tester-C):',C,C_tester,dabs(C_tester-C)
	    	C_opt(graph_num) = result
	    	Beta_opt(graph_num,1:p_max) = angles(1:p_max)
	    	Gamma_opt(graph_num,1:p_max) = angles(p_max+1:2*p_max)
	    	p_C_max(graph_num) = -QAOA_pmax(2*p_max,angles)
	    	good_loops(graph_num) = 0
	    else
	    	good_loops(graph_num) = good_loops(graph_num)+1
	    	!test reflection symmetry for the final result
	    	if ( (good_loops(graph_num) .eq. min_good_loops) ) then
	    		do i = 1,p_max
	    			angles(i) = beta_opt(graph_num,i)
	    			angles(p_max+i) = gamma_opt(graph_num,i)
	    		enddo
	    		call p1_reflection(angles,C_opt(graph_num))
	    		do i = 1,p_max
	    			beta_opt(graph_num,i)=angles(i)
	    			gamma_opt(graph_num,i)=angles(p_max+i)
	    		enddo
	    	endif
	    endif

	else if (optimize_p_Cmax) then

		if (result .gt. (p_C_max(graph_num)+convergence_tolerance)) then

	    	call simplify_angles(angles,result)
	    	C_opt(graph_num) = -QAOA_ExpecC(2*p_max,angles)*cz_max_save(graph_num)
	    	Beta_opt(graph_num,1:p_max) = angles(1:p_max)
	    	Gamma_opt(graph_num,1:p_max) = angles(p_max+1:2*p_max)
	    	p_C_max(graph_num) = result
	    	good_loops(graph_num) = 0

	    else

	    	good_loops(graph_num) = good_loops(graph_num)+1
	    	!test reflection symmetry for the final result
	    	if ( (good_loops(graph_num) .eq. min_good_loops) ) then
	    		do i = 1,p_max
	    			angles(i) = beta_opt(graph_num,i)
	    			angles(p_max+i) = gamma_opt(graph_num,i)
	    		enddo

	    		call p1_reflection(angles,p_C_max(graph_num))

	    		do i = 1,p_max
	    			beta_opt(graph_num,i)=angles(i)
	    			gamma_opt(graph_num,i)=angles(p_max+i)
	    		enddo
	    	endif

	    endif

	endif

end subroutine Check_optimum

function QAOA_pmax(n_angles,angles)

	use parameters 
	implicit none
	integer :: n_angles,p
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: C,QAOA_pmax

	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)
	call Calc_expec_C(C,QAOA_pmax)

	!print *, 'graph_num,C,QAOA_pmax',graph_num,C,QAOA_pmax
	QAOA_pmax = -QAOA_pmax

end function QAOA_pmax


function QAOA_pmax_weighted(n_angles,angles)

	use parameters 
	implicit none
	integer :: n_angles,p,z
	double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
	double precision :: C,QAOA_pmax_weighted

	p = n_angles/2

	B_angles = angles(1:p)
	C_angles = angles(p+1:2*p)

	Call QAOA_Algorithm(p,B_angles,C_angles)

    QAOA_pmax_weighted=0.d0
    do z = 0,2**n_qubits-1
        if (cz_vec(z) .eq. cz_max_save_weighted(graph_num,weight_num)) then
        	QAOA_pmax_weighted = QAOA_pmax_weighted + realpart(psi(z)*dconjg(psi(z)))
        endif
    enddo

	!print *, 'graph_num,C,QAOA_pmax',graph_num,C,QAOA_pmax
	QAOA_pmax_weighted = -QAOA_pmax_weighted

end function QAOA_pmax_weighted







function QAOA_pmax_ma(n_angles,angles)

	use parameters 
	integer :: n_angles,p
	double precision :: angles(n_angles)
	double precision :: C,QAOA_pmax_ma

	beta_ma = angles(1:n_qubits*p_max)
	gamma_ma = angles(n_qubits*p_max+1:n_angles)

	Call QAOA_Algorithm_ma
	call Calc_expec_C(C,QAOA_pmax_ma)

	QAOA_pmax_ma = -QAOA_pmax_ma

end function QAOA_pmax_ma










subroutine simplify_angles_simple(angles,result)

	use parameters
	implicit none

	double precision :: angles(2*p_max),angles_test(2*p_max)
	double precision :: result,result_test
	integer :: i
	logical :: pass

	!put angles into the intervals beta \in [-pi/4,pi/4], gamma \in [-pi,pi]
	angles_test=angles
	do i = 1,p_max
		do while ( (angles_test(i) .lt. -0.25d0*pi) .or. (angles_test(i) .gt. 0.25d0*pi) )
			if (angles_test(i) .lt. -0.25d0*pi) angles_test(i) = angles_test(i) + 0.5d0*pi
			if (angles_test(i) .gt. 0.25d0*pi) angles_test(i) = angles_test(i) - 0.5d0*pi
		enddo
		do while ( (angles_test(p_max+i) .lt. -pi) .or. (angles_test(p_max+i) .gt. pi))
			if (angles_test(p_max+i) .lt. -pi) angles_test(p_max+i) = angles_test(p_max+i) + 2.d0*pi
			if (angles_test(p_max+i) .gt. pi) angles_test(p_max+i) = angles_test(p_max+i) - 2.d0*pi
		enddo
	enddo

	if (angles_test(1) .gt. 0.d0) angles_test=-angles_test

	call Simplify_angles_even_odd_simple(2*p_max,angles_test,result)

	!test these give the same result
	call Check_new_angle_symmetry(angles_test,result,result_test,pass)

	if (.not. pass) then
		print *, 'failed to simplify angles'
		print *, 'result,result_test',result,result_test
		print *,'angles, angles_test:',angles,angles_test
		stop 
	else
		!print *, 'simplified angles!'
		!print *, 'angles,angles_test', angles,angles_test
		!print *, 'result,result_test:',result,result_test
		angles = angles_test
	endif

end subroutine simplify_angles_simple


subroutine simplify_angles(angles,result)

	use parameters
	implicit none

	double precision :: angles(2*p_max),angles_test(2*p_max)
	double precision :: result,result_test
	integer :: i
	logical :: pass

	!put angles into the intervals beta \in [-pi/4,pi/4], gamma \in [-pi,pi]
	angles_test=angles
	do i = 1,p_max
		do while ( (angles_test(i) .lt. -0.25d0*pi) .or. (angles_test(i) .gt. 0.25d0*pi) )
			if (angles_test(i) .lt. -0.25d0*pi) angles_test(i) = angles_test(i) + 0.5d0*pi
			if (angles_test(i) .gt. 0.25d0*pi) angles_test(i) = angles_test(i) - 0.5d0*pi
		enddo
		do while ( (angles_test(p_max+i) .lt. -pi) .or. (angles_test(p_max+i) .gt. pi))
			if (angles_test(p_max+i) .lt. -pi) angles_test(p_max+i) = angles_test(p_max+i) + 2.d0*pi
			if (angles_test(p_max+i) .gt. pi) angles_test(p_max+i) = angles_test(p_max+i) - 2.d0*pi
		enddo
	enddo

	call Simplify_angles_even_odd(2*p_max,angles_test,result)
	!set the initial gamma angle to be < 0
	!if (angles_test(p_max+1) .gt. 0.d0) angles_test = -angles_test
	if (fix_sign_beta1) then
		!set the initial beta angle to be < 0
		if (fix_sign_positive) then
			if (angles_test(1) .lt. 0.d0) angles_test=-angles_test
		else
			if (angles_test(1) .gt. 0.d0) angles_test=-angles_test
		endif
	endif
	if (fix_sign_gamma1) then
		!set the initial beta angle to be < 0
		if (fix_sign_positive) then
			if (angles_test(p_max+1) .lt. 0.d0) angles_test=-angles_test
		else
			if (angles_test(p_max+1) .gt. 0.d0) angles_test=-angles_test
		endif
	endif
	!test these give the same result
	call Check_new_angle_symmetry(angles_test,result,result_test,pass)

	if (.not. pass) then
		print *, 'failed to simplify angles'
		print *, 'result,result_test',result,result_test
		print *,'angles, angles_test:',angles,angles_test
		stop 
	else
		!print *, 'simplified angles!'
		!print *, 'angles,angles_test', angles,angles_test
		!print *, 'result,result_test:',result,result_test
		angles = angles_test
		

	endif

end subroutine simplify_angles


subroutine Check_new_angle_symmetry(angles,result,result_test,pass)

	use parameters
	implicit none

	double precision :: angles(2*p_max)
	double precision :: result, result_test
	!double precision :: QAOA_ExpecC
	logical :: pass

	!External QAOA_ExpecC

	!test these give the same result
	if (optimize_Expec_C) then
		result_test = -QAOA_ExpecC(2*p_max,angles)*cz_max_save(graph_num)
	else if (optimize_p_Cmax) then
		result_test = -QAOA_pmax(2*p_max,angles)
	endif
	if (dabs(result-result_test) .gt. degen_tol) then
		!print *, 'C,C_test:',C,C_test
		!print *, 'fail. angles/pi:',angles/pi
		!print *, 'fail, angles_test:',angles_test
		pass = .false.
	else
		pass = .true.
	endif
	!print *, 'Delta C, C_test:',C_test - C

end subroutine Check_new_angle_symmetry


subroutine Generate_angles(n_angles,angles)

	use parameters
	implicit none

	integer :: i
	integer :: n_angles
	double precision :: junk(6)
	double precision :: angles(n_angles)
	double precision :: tmp

	if (old_angles) then
		read(12,*) (junk(i),i=1,6),(angles(i),i=1,2*p_max)
		angles=angles*pi
	else if (directed_angles) then
		do i = 1,p_max
			tmp=0.d0
			!call random_number(tmp)
			angles(i) = -pi/8 + tmp*0.1d0*pi
			!call random_number(tmp)
			angles(p_max+i) = -3*pi/8+ tmp*0.1d0*pi
		enddo
	else if (angles_from_smaller_p) then
		read(11,*) (junk(i),i=1,6),(angles(i),i=1,smaller_p),(angles(i),i=p_max+1,p_max+smaller_p)
		angles=angles*pi
		!add a small random perturbation t-o the previous angles
		do i = 1,p_max-1
			call random_number(tmp)
			angles(i) = angles(i) + pi/2.d0*(tmp-0.5d0)/2.d0
			call random_number(tmp)
			angles(p_max+i) = angles(p_max+i) + 2.d0*pi*(tmp-0.5d0)/2.d0
		enddo
			call random_number(tmp)
			angles(p_max) = pi/2.d0*(tmp-0.5d0)
			call random_number(tmp)
			angles(2*p_max) = 2*pi*(tmp-0.5d0)
		!print *, 'angles:', (angles(i),i=1,2*p)
	else if (many_angles) then
		do i = 1,n_qubits
    		call random_number(tmp)
    		angles(i) = beta_max*2.d0*(tmp-0.5d0)
    	enddo
    	beta_ma = angles(1:n_qubits*p_max)
    	do i = 1,n_edges
    		call random_number(tmp)
    		angles(n_qubits*p_max+i) = gamma_max*2.d0*(tmp-0.5d0)
    	enddo
    	gamma_ma(1:n_angles-p_max*n_qubits) = angles(n_qubits*p_max+1:n_angles)
    elseif (even_odd_only .and. even_odd_angle_dist) then
    	do i = 1,p_max
    		call random_number(tmp)
    		if (i .eq. 1) then
    			angles(i) = -pi/4.d0*tmp !-0.25pi < \beta_1 < 0
    		else
    			angles(i) = -pi/2.d0*(tmp-0.5d0)
    		endif
    		call random_number(tmp)
    		angles(p_max+i) = 1.d0*pi*(tmp-0.5d0)
    	enddo
    	!print *, 'angles/pi:',angles/pi
	else
    	do i = 1,p_max
    		call random_number(tmp)
    		angles(i) = pi/2.d0*(tmp-0.5d0)
    		call random_number(tmp)
    		angles(p_max+i) = 2.d0*pi*(tmp-0.5d0)
    	enddo
    	!print *, 'angles/pi:',angles/pi
    endif

end subroutine Generate_angles


subroutine calc_vertex_degrees

	use parameters

	implicit none

	integer :: n,q0,q1

	do n = 0,n_edges-1
		q0 = edges(n,0)
		q1 = edges(n,1)
		vertex_degrees(graph_num,q0) = vertex_degrees(graph_num,q0)+1
		vertex_degrees(graph_num,q1) = vertex_degrees(graph_num,q1)+1
	enddo

	all_odd_degree(graph_num)=.true.
	all_even_degree(graph_num)=.true.
	do n = 0,n_qubits-1
		if (mod(vertex_degrees(graph_num,n),2) .eq. 0) all_odd_degree(graph_num)=.false.
		if (mod(vertex_degrees(graph_num,n),2) .eq. 1) all_even_degree(graph_num)=.false.
	enddo
	!if (all_odd_degree(graph_num)) print *, graph_num,'all odd degree'
	!if (all_even_degree(graph_num)) print *, graph_num,'all even degree'

end subroutine calc_vertex_degrees



subroutine simplify_angles_even_odd_simple(n_angles,angles,result)

	use parameters

	implicit none

	integer :: n_angles
	double precision :: angles(n_angles),angles_test(n_angles)
	double precision :: result,result_test
	integer :: i,j,k
	logical :: dif_angles,pass

	dif_angles=.false.
	angles_test = angles

	!assumes beta_1 < 0

	!if all even degree, set all \gamma_i < 0
	if (all_even_degree(graph_num)) then

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

	!if all odd degree, set all gamma_i < 0 for i > 1
	elseif (all_odd_degree(graph_num)) then

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
		!print *, 'all odd degree symmetry:'
		!print *, 'old angles/pi:', angles/pi
		!print *, 'new angles/pi:', angles_test/pi
		call Check_new_angle_symmetry(angles_test,result,result_test,pass)
		if (pass) then
			!print *, 'satisfied symmetry'
			angles = angles_test
		else
			print *, 'failed even/odd symmetry'	
			print *, 'result, result_test:',result,result_test
			stop
		endif
	endif

end subroutine simplify_angles_even_odd_simple



subroutine Simplify_angles_even_odd(n_angles,angles,result)

	use parameters 

	implicit none

	integer :: n_angles
	double precision :: angles(n_angles),angles_test(n_angles)
	double precision :: result,result_test
	integer :: i,j,k
	logical :: dif_angles,pass

	dif_angles=.false.
	angles_test = angles

	if (all_even_degree(graph_num)) then
		!print *, 'all_even, graph_num', graph_num
		if (fix_sign_gamma1) then
			if (fix_sign_positive) then
				if (angles_test(n_angles/2+1) .gt. pi/2.d0 .and. angles_test(1) .lt. 0.d0) then
					angles_test(n_angles/2+1) = angles_test(n_angles/2+1) - pi
					!angles_test = -angles_test
					dif_angles=.true.
				endif
				if (angles_test(n_angles/2+1) .lt. -pi/2.d0 .and. angles_test(1) .gt. 0.d0) then
					angles_test(n_angles/2+1) = angles_test(n_angles/2+1) + pi
					!angles_test = -angles_test
					dif_angles=.true.
				endif
			else 
				if (angles_test(n_angles/2+1) .lt. -pi/2.d0 .and. angles_test(1) .gt. 0.d0) then
					angles_test(n_angles/2+1) = angles_test(n_angles/2+1) + pi
					!angles_test = -angles_test
					dif_angles=.true.
				endif
			endif
			if (n_angles .gt. 2) then
				do i = 2,n_angles/2
					if (fix_sign_positive) then
						if (angles_test(n_angles/2+i) .gt. 0.d0) then
							angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
							dif_angles=.true.
						endif
					else
						if (angles_test(n_angles/2+i) .lt. 0.d0) then
							angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
							dif_angles=.true.
						endif
					endif
				enddo
			endif



		else if (fix_sign_beta1) then

!!!NOT SURE IF THIS WORKS WITH FIXING THE ANGLES AS POSITIVE!!!
			if (angles_test(1) .gt. 0) then
				angles_test=-angles_test
			endif
			do i = 1,n_angles/2
				!	if (angles_test(n_angles/2+i) .lt. 0.d0) then
				!		angles_test(n_angles/2+i)=angles_test(n_angles/2+i)+pi
				!		dif_angles=.true.
				!	endif
				!else if (angles_test(1) .lt. 0) then
				if (angles_test(n_angles/2+i) .gt. 0.d0) then
					angles_test(n_angles/2+i)=angles_test(n_angles/2+i)-pi
					dif_angles=.true.
				endif
				!endif
			enddo
		endif

		if (dif_angles) then
			!print *, 'all even degree symmetry:'
			!print *, 'old angles:', angles
			!print *, 'new angles:', angles_test
			call Check_new_angle_symmetry(angles_test,result,result_test,pass)
			if (pass) then
				!print *, 'satisfied symmetry'
				angles = angles_test
			else
				print *, 'failed even/odd symmetry'
				print *, 'result, result_test:',result,result_test
				stop
			endif
		endif
	endif

	if (all_odd_degree(graph_num)) then
		!print *, 'all odd, modifying angles'
		if (fix_sign_beta1) then
			do i = 1,n_angles/2
				if (angles_test(1) .gt. 0.d0 .and. angles_test(n_angles/2+i) .lt. 0.d0) then
					angles_test=-angles_test
					dif_angles=.true.
				endif
			enddo
			do i = 1,n_angles/2
				if (angles_test(n_angles/2+i) .gt. 0.d0) then
					angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
					do j = i,n_angles/2
						angles_test(j) = -angles_test(j)
					enddo
					dif_angles=.true.
				endif
			enddo
		elseif (fix_sign_gamma1) then
			if (angles_test(n_angles/2+1) .gt. 0.d0) then
				angles_test=-angles_test
				dif_angles=.true.
			endif
			do i = 1,n_angles/2
				if (angles_test(n_angles/2+i) .gt. 0.d0) then
					angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
					do j = i,n_angles/2
						angles_test(j) = -angles_test(j)
					enddo
					dif_angles=.true.
				endif
			enddo
		else
			if (angles_test(n_angles/2+1) .lt. -pi/2.d0 .and. angles_test(1) .lt. 0.d0) then
				angles_test(n_angles/2+1) = angles_test(n_angles/2+1) + pi
				do i = 1,n_angles/2
					angles_test(i) = -angles_test(i)
				enddo
				angles_test=-angles_test
				dif_angles=.true.
			endif
			if (angles_test(n_angles/2+1) .gt. pi/2.d0 .and. angles_test(1) .gt. 0.d0) then
				angles_test(n_angles/2+1) = angles_test(n_angles/2+1) - pi
				do i = 1,n_angles/2
					angles_test(i) = -angles_test(i)
				enddo
				angles_test=-angles_test
				dif_angles=.true.
			endif
			if (n_angles .gt. 2) then
				do i = 2,n_angles/2
					if (angles_test(n_angles/2+i) .gt. 0.d0) then
						angles_test(n_angles/2+i) = angles_test(n_angles/2+i) - pi
						do j = i,n_angles/2
							angles_test(j) = -angles_test(j)
						enddo
						dif_angles=.true.
					endif
				enddo
			endif
		endif
		if (dif_angles) then
			!print *, 'all odd degree symmetry:'
			print *, 'old angles/pi:', angles/pi
			print *, 'new angles/pi:', angles_test/pi
			call Check_new_angle_symmetry(angles_test,result,result_test,pass)
			if (pass) then
				print *, 'satisfied symmetry'
				angles = angles_test
			else
				print *, 'failed even/odd symmetry'	
				print *, 'result, result_test:',result,result_test
				stop
			endif
		endif
	endif

end subroutine Simplify_angles_even_odd

subroutine p1_reflection(angles,result)

	use parameters
	implicit none

	double precision :: angles(2*p_max),angles_test(2*p_max)
	double precision :: result,result_test
	integer :: i
	logical :: pass

	if (angles(p_max+1) .lt. -(pi/2.d0+1.d-6)) then
		angles_test=angles
		if (angles(1) .gt. 1.d-6) then
			angles_test(p_max+1) = -pi -angles(p_max+1)
			angles_test(1) = -angles(1)
		else
			angles_test(p_max+1) = -pi -angles(p_max+1)
			angles_test(1) = angles(1)
		endif


		call Check_new_angle_symmetry(angles_test,result,result_test,pass)

		if (.not. pass) then
			print *, 'failed p1_reflection'
			print *, 'graph_num,result,result_test,cz_max:',graph_num,result,result_test,cz_max_save(graph_num)
			print *, 'original angles/pi:',angles/pi
			print *, 'new angles/pi:',angles_test/pi
			!if (C_test .gt. C) angles = angles_test
			!stop 
		else
			!print *, 'passed p1_reflection test!'
			!print *, 'graph_num', graph_num
			!print *, 'angles_old', angles
			!print *, 'angles_new', angles_test
			angles = angles_test
		endif

	endif
end subroutine p1_reflection

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

	beta_a=0.d0
	gamma_a=0.d0
	!write beta(i), gamma(i) on a grid of points
	do i= p_max, 1, -1
		y0=0.d0
		x0=0.d0
		if (i .lt. p_max) then
			do j = p_max, i, -1
				y0 = y0 + beta_a(j)*(max_it_beta**(j-1))
				x0 = x0 + gamma_a(j)*(max_it_gamma**(j-1))
			enddo
		endif
		beta_a(i) = (y-y0)/(max_it_beta**(i-1))
		gamma_a(i) = ((x-x0)/max_it_gamma**(i-1))
	enddo

	!map the beta(i),gamma(i) grid to the the QAOA angle intervals
	do i = 1,p_max
		beta_a(i)=beta_max*2*(dble(beta_a(i))/dble(max_it_beta-1)-0.5d0)
		gamma_a(i)=gamma_max*2*(dble(gamma_a(i))/dble(max_it_gamma-1)-0.5d0)
		angles(i) = beta_a(i)
		angles(p_max+i) = gamma_a(i)
	enddo


end subroutine set_angles


subroutine Calc_Z_prod(state)

	use parameters
	implicit none

	double precision :: z0z1,z1z2,z0z2
	integer :: z
	complex*16 :: state(0:2**n_qubits-1)

	z0z1=0.d0
	z0z2=0.d0
	z1z2=0.d0

	do z = 0,2**n_qubits-1
		if (basis(z,0) .eq. basis(z,1)) then
			z0z1 = z0z1 + dble(state(z)*dconjg(state(z)))
		else
			z0z1 = z0z1 - dble(state(z)*dconjg(state(z)))
		endif
		if (basis(z,0) .eq. basis(z,2)) then
			z0z2 = z0z2 + dble(state(z)*dconjg(state(z)))
		else
			z0z2 = z0z2 - dble(state(z)*dconjg(state(z)))
		endif
		if (basis(z,1) .eq. basis(z,2)) then
			z1z2 = z1z2 + dble(state(z)*dconjg(state(z)))
		else
			z1z2 = z1z2 - dble(state(z)*dconjg(state(z)))
		endif
	enddo

	print *, 'z0z1,z1z2,z0z2',z0z1,z1z2,z0z2
end subroutine Calc_Z_prod


subroutine dfunc(n,x,grad,func)

	implicit none

	integer :: n,i
	double precision :: x(n), x_dif_plus(n), x_dif_min(n)
	double precision :: nominal_delta = 1.d-6
	double precision :: h,temp
	double precision :: grad(n)
	double precision :: func

	External func

	!set h exact to numerical precision
	!so h can be represented exactly in base 2
	temp =  x(1) + nominal_delta
 	h = temp-x(1)

 	do i = 1,n
  		x_dif_plus = x
 		x_dif_min = x
 		x_dif_plus(i) = x(i) + h
 		x_dif_min(i) = x(i) - h
 		grad(i) = (func(n,x_dif_plus) - func(n,x_dif_min))/(2.d0*h)
	enddo

end subroutine dfunc

double precision function log2(x)
  implicit none
  double precision, intent(in) :: x

  log2 = log(x) / log(2.d0)
end function

subroutine Count_Num_graphs

	use parameters 
	implicit none

	!the below counts an extra time on the last call
	!to read_adjacency_matrix
	graph_num_tot=-1
	do while (more_graphs)
		call Read_Adjacency_Matrix
		graph_num_tot = graph_num_tot+1
	enddo
	more_graphs=.true.

end subroutine Count_num_graphs


subroutine Read_Adjacency_Matrix

	use parameters
	implicit none

	integer :: n,nn,m
	integer :: iostatus,edge_count
	character (len=n_qubits-1) :: line
	logical :: non_graph_line


	non_graph_line=.true.
	edge_count=-1
	edges=0
	!skip lines that aren't the adjacency matrix entries
	do while (non_graph_line)
		read(1,*,IOSTAT = iostatus) line
		if (iostatus .lt. 0) then
			!print*, 'end of file'
			non_graph_line=.false.
			more_graphs=.false.
		else if (iostatus .gt. 0) then
			print*, 'error reading file:', iostatus
			stop
		endif
		if (line(1:1) .eq. '0' .or. line(1:1) .eq. '1') then
			non_graph_line = .false.
		endif
	enddo
	if (more_graphs) then
		do n = 0,n_qubits-2
			if (n .ne. 0) read(1,*) line
			do nn = n+1,n_qubits-1
				if (line(nn-n:nn-n) .eq. '1') then
					edge_count = edge_count+1
					edges(edge_count,0) = n
					edges(edge_count,1) = nn
				endif
			enddo
		enddo
	endif
	n_edges=edge_count+1

end subroutine Read_Adjacency_Matrix


subroutine Random_edges

	!Generate a graph with random edges
	use parameters
	implicit none

	double precision :: tmp
	integer :: i,j,k,l
	logical :: old_edge=.true.

	call random_number(tmp)

	tmp=0.5d0
	n_edges = nint(tmp*n_edges_max)
	if (n_edges .eq. 0) n_edges =1
	print*,'n_edges,n_edges_max:',n_edges,n_edges_max
	!do i = 1,n_edges-1
	!	edges(i,0)=0
	!	edges(i,1)=i
	!enddo
	if (.true.) then
		do i = 0,n_edges-1
			old_edge=.true.
			do while (old_edge)
				call random_number(tmp)
				j = nint(tmp*n_qubits-0.5d0)
				call random_number(tmp)
				k = nint(tmp*n_qubits-0.5d0)
				old_edge=.false.
				do l = 0,i
					if ( j .eq. k .or. (edges(l,0) .eq. j .and. edges(l,1) .eq. k) .or. &
					   	& (edges(l,0) .eq. k .and. edges(l,1) .eq. j) ) then
					   	old_edge=.true.
					endif
				enddo
				if (.not. old_edge) then
					edges(i,0) = j
					edges(i,1) = k
				endif
			enddo
		enddo
	endif

end subroutine Random_edges


subroutine GHZ_edges

	!Generate a graph with random edges
	use parameters
	implicit none

	double precision :: tmp
	integer :: i,j,k,l

	n_edges = n_qubits
	do i = 0,n_edges-2
		edges(i,0)=i
		edges(i,1)=i+1
	enddo
	edges(n_qubits-1,0) = n_qubits-1
	edges(n_qubits-1,1) = 0

	!print *, 'edges:',edges
end subroutine GHZ_edges

subroutine set_angles_GHZ(angles,stdev_noise,mean_noise)
	use parameters
	implicit none
	double precision :: angles(2*p_max)
	double precision :: stdev_noise,mean_noise
	double precision :: mean_offset_beta,mean_offset_gamma
	double precision :: stdev_error
	integer :: i

	if (mean_noise .gt. 1.d-13) then
		mean_offset_beta=mean_noise
		mean_offset_gamma=mean_noise
		!call Gaussian_variate(0.d0,mean_noise,mean_offset_beta)
		!call Gaussian_variate(0.d0,mean_noise,mean_offset_gamma)
	else
		mean_offset_beta=0.d0
		mean_offset_gamma=0.d0
	endif

    if (p_max .eq. 4) then
        angles(p_max+1)=0.5297*2.d0
        angles(p_max+2)=0.7243*2.d0
        angles(p_max+3)=0.6151*2.d0
        angles(p_max+4)=0.5243*2.d0
    elseif (p_max .eq. 5) then
        angles(p_max+1)=0.5814*2.d0
        angles(p_max+2)=0.6360*2.d0
        angles(p_max+3)=0.5993*2.d0
        angles(p_max+4)=0.7889*2.d0
        angles(p_max+5)=0.5230*2.d0
    elseif (p_max .eq. 6) then
        angles(p_max+1)=0.5466*2.d0
        angles(p_max+2)=0.6902*2.d0
        angles(p_max+3)=0.5946*2.d0
        angles(p_max+4)=0.7276*2.d0
        angles(p_max+5)=0.7212*2.d0
        angles(p_max+6)=0.5452*2.d0
    elseif (p_max .eq. 7) then
        angles(p_max+1)=0.6513*2.d0
        angles(p_max+2)=0.5841*2.d0
        angles(p_max+3)=0.7633*2.d0
        angles(p_max+4)=0.5660*2.d0
        angles(p_max+5)=0.8270*2.d0
        angles(p_max+6)=0.6704*2.d0
        angles(p_max+7)=0.5696*2.d0
    elseif (p_max .eq. 8) then
        angles(p_max+1)=0.5846*2.d0
        angles(p_max+2)=0.6105*2.d0
        angles(p_max+3)=0.7966*2.d0
        angles(p_max+4)=0.6373*2.d0
        angles(p_max+5)=0.7745*2.d0
        angles(p_max+6)=0.6152*2.d0
        angles(p_max+7)=0.7155*2.d0
        angles(p_max+8)=0.5796*2.d0
    elseif (p_max .eq. 9) then
        angles(p_max+1)=0.6064*2.d0
        angles(p_max+2)=0.6632*2.d0
        angles(p_max+3)=0.6660*2.d0
        angles(p_max+4)=0.7773*2.d0
        angles(p_max+5)=0.6904*2.d0
        angles(p_max+6)=0.7133*2.d0
        angles(p_max+7)=0.6302*2.d0
        angles(p_max+8)=0.7780*2.d0
        angles(p_max+9)=0.5232*2.d0
    endif

    do i=1,p_max
    	angles(i) = -angles(2*p_max+1-i)*0.5d0
    enddo
    if (stdev_noise .gt. 1.d-13) then
    	do i = 1,2*p_max
    		call Gaussian_variate(0.d0,stdev_noise,stdev_error)
    		!if (i .le. p_max) stdev_error=stdev_error*0.5d0
    		angles(i) = angles(i)*(1.d0 + stdev_error)
    	enddo
    endif

    if (mean_noise .gt. 1.d-13) then
    	do i = 1,2*p_max
    		if (i .le. p_max) then
    			angles(i) = angles(i)*(1.d0 + mean_offset_beta)!*0.5d0
    		else
    			angles(i) = angles(i)*(1.d0 + mean_offset_gamma)
    		endif
    	enddo
    endif

    !print *, stdev_error
    !print *, 'angles:',angles

end subroutine set_angles_GHZ

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



Subroutine Gaussian_variate(mean,stdev,variate)
	implicit none
	logical :: keepgoing
	double precision, parameter :: x_max=10.d0
	double precision :: mean,stdev
	double precision :: tmp,variate

	keepgoing = .true.
	do while (keepgoing) 
		call random_number(tmp)
		call random_number(variate)
		variate = 2.d0*(variate-0.5d0)*x_max 
		!randomly sample from a Gaussian
		!accept variate based on the probability of variate from the Gaussian
		if ( tmp .le. exp(-variate**2/2)/dsqrt(2.d0*dacos(-1.d0)) )then
			variate = variate*stdev+mean
			keepgoing = .false.
		endif
	enddo

end subroutine Gaussian_variate


subroutine Calc_gamma_vec_ma

	use parameters
	implicit none

	integer :: z,m

	gamma_vec_ma=0.d0
    do z =0,2**n_qubits-1
        do m = 0,n_edges-1
            if (basis(z,edges(m,0)) .ne. basis(z,edges(m,1))) then
           		gamma_vec_ma(z) = gamma_vec_ma(z) + gamma_ma(m+1)
            endif
        enddo
    enddo


end subroutine Calc_gamma_vec_ma


subroutine Calc_cz_vec(BFGS_loops)
    !calculate a vector c_z
    !with entries c_z[z] that equal the cost function evaluation for
    !the binary expansion of z
    use parameters
    implicit none
    integer :: z,m
    integer :: BFGS_loops
    logical :: calc_zmax=.false.

    cz_vec=0.d0

    do z =0,2**n_qubits-1
        do m = 0,n_edges-1
            if (basis(z,edges(m,0)) .ne. basis(z,edges(m,1))) then
            	cz_vec(z) = cz_vec(z) + 1
            	!if (many_angles) then
            	!	gamma_vec_ma(z) = gamma_vec_ma(z) + gamma_ma(m)
            	!endif
            endif
        enddo
    enddo

	if (BFGS_loops .eq. 1) then
		!calculate C0, the <C> for the initial state
		psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)
		call Calc_expec_C(C0(graph_num),p0_C_max(graph_num))
		cz_max_save(graph_num) = maxval(cz_vec)

		call calc_vertex_degrees

		if (all_odd_degree(graph_num)) print *, graph_num,'all odd degree'
		if (all_even_degree(graph_num)) print *, graph_num,'all even degree'

		if (all_even_degree(graph_num) .or. all_odd_degree(graph_num)) n_graphs_even_odd = n_graphs_even_odd+1
		
		if (calc_zmax) then
			P0_C_max(graph_num) = 0.d0
			do z = 0,2**n_qubits-1
				if (cz_vec(z) .eq. maxval(cz_vec)) then
					P0_c_max(graph_num) = P0_C_max(graph_num) +1.d0/dble(dim)
					if (dabs( p0_c_max(graph_num) - 1.d0/dble(dim) ) .lt. 1.d-12) then
						z_max(graph_num) = z
					endif
				endif
			enddo
		endif

	endif

end subroutine Calc_cz_vec


subroutine Calc_cz_vec_weighted(BFGS_loops,weights,n_weighted_graphs)
    !calculate a vector c_z
    !with entries c_z[z] that equal the cost function evaluation for
    !the binary expansion of z
    use parameters
    implicit none
    integer :: n_weighted_graphs
    double precision :: weights(graph_num_tot,n_weighted_graphs,0:n_edges_max-1)
    integer :: z,m
    integer :: BFGS_loops
    logical :: calc_zmax=.false.

    cz_vec=0.d0

    do z =0,2**n_qubits-1
        do m = 0,n_edges-1
            if (basis(z,edges(m,0)) .ne. basis(z,edges(m,1))) then
            	cz_vec(z) = cz_vec(z) + weights(graph_num,weight_num,m)
            	!if (many_angles) then
            	!	gamma_vec_ma(z) = gamma_vec_ma(z) + gamma_ma(m)
            	!endif
            endif
        enddo
    enddo

	if (BFGS_loops .eq. 1) then
		!calculate C0, the <C> for the initial state
		psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)
		call Calc_expec_C(C0_weighted(graph_num,weight_num),p0_C_max_weighted(graph_num,weight_num))
		cz_max_save_weighted(graph_num,weight_num) = maxval(cz_vec)

	endif

end subroutine Calc_cz_vec_weighted




subroutine Calc_cz_vec_weights_GTRI(BFGS_loops)
    !calculate a vector c_z
    !with entries c_z[z] that equal the cost function evaluation for
    !the binary expansion of z
    use parameters
    implicit none
    integer :: z,m
    integer :: BFGS_loops
    logical :: calc_zmax=.false.
    double precision :: weighted_edges(0:14,0:2)
    n_edges=15
    cz_vec=0.d0

    weighted_edges(0,0) = 0
    weighted_edges(0,1) = 1
    weighted_edges(0,2) = 0.71d0
    weighted_edges(1,0) = 0
    weighted_edges(1,1) = 2
    weighted_edges(1,2) = 0.38d0
    weighted_edges(2,0) = 0
    weighted_edges(2,1) = 3
    weighted_edges(2,2) = -0.38d0
    weighted_edges(3,0) = 0
    weighted_edges(3,1) = 4
    weighted_edges(3,2) = -0.71d0
    weighted_edges(4,0) = 0
    weighted_edges(4,1) = 5
    weighted_edges(4,2) = 0.5d0

    weighted_edges(5,0) = 1
    weighted_edges(5,1) = 2
    weighted_edges(5,2) = -0.54d0
    weighted_edges(6,0) = 1
    weighted_edges(6,1) = 3
    weighted_edges(6,2) = 0.54d0
    weighted_edges(7,0) = 1
    weighted_edges(7,1) = 4
    weighted_edges(7,2) = 1.d0
    weighted_edges(8,0) = 1
    weighted_edges(8,1) = 5
    weighted_edges(8,2) = -0.71d0

    weighted_edges(9,0) = 2
    weighted_edges(9,1) = 3
    weighted_edges(9,2) = 0.29d0
    weighted_edges(10,0) = 2
    weighted_edges(10,1) = 4
    weighted_edges(10,2) = 0.54d0
    weighted_edges(11,0) = 2
    weighted_edges(11,1) = 5
    weighted_edges(11,2) = -0.38d0

    weighted_edges(12,0) = 3
    weighted_edges(12,1) = 4
    weighted_edges(12,2) = -0.54d0
    weighted_edges(13,0) = 3
    weighted_edges(13,1) = 5
    weighted_edges(13,2) = 0.38d0

    weighted_edges(14,0) = 4
    weighted_edges(14,1) = 5
    weighted_edges(14,2) = 0.71d0

    do z =0,2**n_qubits-1
        do m = 0,n_edges-1
            if (basis(z,nint(weighted_edges(m,0))) .ne. basis(z,nint(weighted_edges(m,1)))) then
            	cz_vec(z) = cz_vec(z) + weighted_edges(m,2)
            	!if (many_angles) then
            	!	gamma_vec_ma(z) = gamma_vec_ma(z) + gamma_ma(m)
            	!endif
            endif
        enddo
    enddo

	if (BFGS_loops .eq. 1) then
		!calculate C0, the <C> for the initial state
		psi=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)
		call Calc_expec_C(C0(graph_num),p0_C_max(graph_num))
		cz_max_save(graph_num) = maxval(cz_vec)

		call calc_vertex_degrees

		if (all_odd_degree(graph_num)) print *, graph_num,'all odd degree'
		if (all_even_degree(graph_num)) print *, graph_num,'all even degree'

		if (all_even_degree(graph_num) .or. all_odd_degree(graph_num)) n_graphs_even_odd = n_graphs_even_odd+1
		
		if (calc_zmax) then
			P0_C_max(graph_num) = 0.d0
			do z = 0,2**n_qubits-1
				if (cz_vec(z) .eq. maxval(cz_vec)) then
					P0_c_max(graph_num) = P0_C_max(graph_num) +1.d0/dble(dim)
					if (dabs( p0_c_max(graph_num) - 1.d0/dble(dim) ) .lt. 1.d-12) then
						z_max(graph_num) = z
					endif
				endif
			enddo
		endif

	endif

end subroutine Calc_cz_vec_weights_GTRI


subroutine Calc_Pcz_max(Pczmax)

	use parameters
	implicit none

	integer :: z
	double precision :: Pczmax

	pczmax=0.d0
	do z = 0,2**(n_qubits)-1
		if (dabs(dble(cz_vec(z)) -maxval(cz_vec)) .lt. 1.d-8) then
			pczmax=pczmax+1.d0
		endif
		!print *, 'z,cz_vec(z),cz_max,pczmax',z,cz_vec(z),cz_max,pczmax
	enddo
	pczmax=pczmax/dble(2**n_qubits)

end subroutine calc_Pcz_max

subroutine Construct_Ub(b_angle)

	use parameters
	implicit none
	double precision :: b_angle

	Ub(0,0) = dcmplx(dcos(b_angle),0.d0)
	Ub(0,1) = dcmplx(0.d0,-dsin(b_angle))
	Ub(1,0) = dcmplx(0.d0,-dsin(b_angle))
	Ub(1,1) = dcmplx(dcos(b_angle),0.d0)

end subroutine Construct_Ub



subroutine Apply_Uc_ma

	use parameters
	implicit none
	integer :: z

	if (many_angles) then
		do z = 0,2**n_qubits-1
			psi(z) = psi(z) * dcmplx(dcos(gamma_vec_ma(z)),-dsin(gamma_vec_ma(z)))
		enddo
	else
		print *, 'not many angles!'
		stop
	endif

end subroutine Apply_Uc_ma


subroutine Apply_Uc(c_angle)

	use parameters
	implicit none
	integer :: z
	double precision :: c_angle

    do z = 0,2**n_qubits-1
        psi(z) = psi(z)*dcmplx(dcos(cz_vec(z)*c_angle),-dsin(cz_vec(z)*c_angle))
    enddo

end subroutine Apply_Uc


subroutine Calc_expec_C(C,p_cz_max)

	use parameters
	implicit none
	integer :: z
	double precision :: C,p_cz_max

    C=0.d0
    p_cz_max=0.d0
    do z = 0,2**n_qubits-1
        C = C + realpart(psi(z)*dconjg(psi(z)))*dble(cz_vec(z))  
        if (cz_vec(z) .eq. cz_max_save(graph_num)) then
        	p_cz_max = p_cz_max + realpart(psi(z)*dconjg(psi(z)))
        endif
    enddo
end subroutine Calc_expec_C


subroutine Apply_Ub_ma
	use parameters
    implicit none
    integer :: n,z
    complex*16 :: temporary

    !apply Ub to all qubits
    !for each qubit, loop over coupled pairs of states
    !where the basis state of qubit n changes 0->1
    !and all other qubit basis states are constant.
    !The Ub operators are assumed to commute so the
    !order the n are executed in doesn't matter  

    do n = 0,n_qubits-1

    	call construct_Ub(beta_ma(n+1))
        do z = 0, 2**(n_qubits-1)-1
            temporary = Ub(0,0)*psi(pairs_list(n,z,0)) + Ub(0,1)*psi(pairs_list(n,z,1))
            psi(pairs_list(n,z,1)) = Ub(1,1)*psi(pairs_list(n,z,1)) + Ub(1,0)*psi(pairs_list(n,z,0))
            psi(pairs_list(n,z,0)) = temporary
        enddo

    enddo

end subroutine Apply_Ub_ma


subroutine Apply_Ub
	use parameters
    implicit none
    integer :: n,z
    complex*16 :: temporary

    !apply Ub to all qubits
    !for each qubit, loop over coupled pairs of states
    !where the basis state of qubit n changes 0->1
    !and all other qubit basis states are constant.
    !The Ub operators are assumed to commute so the
    !order the n are executed in doesn't matter  

    do n = 0,n_qubits-1

        do z = 0, 2**(n_qubits-1)-1
            temporary = Ub(0,0)*psi(pairs_list(n,z,0)) + Ub(0,1)*psi(pairs_list(n,z,1))
            psi(pairs_list(n,z,1)) = Ub(1,1)*psi(pairs_list(n,z,1)) + Ub(1,0)*psi(pairs_list(n,z,0))
            psi(pairs_list(n,z,0)) = temporary
        enddo

    enddo

end subroutine Apply_Ub


subroutine Print_binary_expansion_string(z,z_str)
	use parameters
	implicit none

	character*10 :: z_str
	integer :: n,z,m

	!print *, 'basis(z):',(basis(z,n),n=0,n_qubits-1)
	z_str = char(48+basis(z,n_qubits-1))
	!print *, '0, z_str',z_str
	do n = 1,n_qubits-1
		m = n_qubits-1-n
		z_str = trim(z_str)//char(48+basis(z,m))
		!print *, 'm,basis(z,m),z_str',m,basis(z,m),z_str
	enddo

	!print *, 'z,z_str',z,z_str

end subroutine Print_binary_expansion_string


subroutine enumerate_basis
	use parameters
	implicit none
	integer :: n,z
	integer :: tmp

	do z = 0,2**n_qubits-1
		tmp=z
		do n = 0,n_qubits-1
			basis(z,n) = mod(tmp,2)
			tmp=tmp/2
		enddo
		!print *, 'z,basis(z,:):',z,basis(z,:)
	enddo

end subroutine enumerate_basis


subroutine calc_pair_list_fast
	!generate a list of pairs of states that are related by switching single
	!qubit states from 0 -> 1
	!resulting list has 2**(n_qubits-1) entries for each qubit n, with pairs
	!of related states where all other qubits are the same and n =0,1
	use parameters
	implicit none
	integer :: z0,n,j
	integer :: tally(0:n_qubits)

	tally=-1
	do z0=0,2**n_qubits-1
    	do n = 0,n_qubits-1
    		if (basis(z0,n) .eq. 0) then
    			tally(n) = tally(n)+1
    			pairs_list(n,tally(n),0) = z0
    			pairs_list(n,tally(n),1) = z0 + 2**n
    		endif
    	enddo
    enddo

end subroutine calc_pair_list_fast


subroutine Calc_pair_list
    !Make a list of basis state numbers that differ in one qubit
    !These will be used later to calculate single-qubit unitaries
    !acting on the quantum state
    use parameters
    implicit none
    integer :: z0,z1,counts,different,j,n
    integer :: tally(0:n_qubits-1)


    tally=-1
    do z0 =0,2**n_qubits-1
        do z1 =z0+1,2**n_qubits-1
            counts=0
            do n = 0,n_qubits-1
                if (basis(z0,n) .ne. basis(z1,n)) then
                    counts=counts+1
                    different=n
                endif
            enddo
            if (counts .eq. 1) then
            	tally(different) = tally(different)+1
                pairs_list(different,tally(different),0) = z0
                pairs_list(different,tally(different),1) = z1
                !print*, 'different,z0,z1,basis(z0,:),basis(z1,:)'
                !print*, different,z0,z1
                !print*, (basis(z0,j),j=1,n_qubits)
               	!print*, (basis(z1,j),j=1,n_qubits)
            endif
        enddo
    enddo

end subroutine Calc_pair_list


end module QAOA_subroutines