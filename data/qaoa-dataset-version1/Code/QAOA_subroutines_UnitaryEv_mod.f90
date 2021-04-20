module QAOA_subroutines_UnitaryEv

contains

function QAOA_pmax(n_angles,angles,graph_number)

    use parameters 
    implicit none
    integer :: n_angles,p,graph_number
    double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
    double precision :: C,QAOA_pmax
    complex*16 :: state(0:z_max)

    p = n_angles/2

    B_angles = angles(1:p)
    C_angles = angles(p+1:2*p)

    Call QAOA_Algorithm(state,p,B_angles,C_angles,graph_number)
    call Calc_expec_C(state,C,QAOA_pmax,graph_number)

    !print *, 'graph_number,C,QAOA_pmax',graph_number,C,QAOA_pmax
    QAOA_pmax = -QAOA_pmax

end function QAOA_pmax

function QAOA_ExpecC(n_angles,angles,graph_number)

    use parameters 
    integer :: n_angles,p,graph_number
    double precision :: angles(n_angles),B_angles(n_angles/2),C_angles(n_angles/2)
    double precision :: p_cz_max,QAOA_ExpecC
    complex*16 :: state(0:z_max)

    p = n_angles/2

    B_angles = angles(1:p)
    C_angles = angles(p+1:2*p)

    Call QAOA_Algorithm(state,p,B_angles,C_angles,graph_number)
    call Calc_expec_C(state,QAOA_ExpecC,p_cz_max,graph_number)

    QAOA_ExpecC = -QAOA_ExpecC/cz_max_save(graph_number)

end function QAOA_ExpecC

subroutine QAOA_Algorithm(state,p,B_angles,C_angles,graph_number)
    use parameters
    implicit none

    integer :: j,graph_number
    integer :: p 
    double precision :: b_angles(p),c_angles(p)!beta and gamma
    complex*16 :: state(0:z_max)

    state=dcmplx(1.d0/dsqrt(dble(dim)),0.d0)

    do j = 1,p
        !call construct_Ub(b_angles(j))
        call Apply_Uc(state,c_angles(j),graph_number)
        call Apply_Ub(state,b_angles(j))
    enddo

end subroutine QAOA_Algorithm

subroutine Calc_expec_C(state,C,p_cz_max,graph_number)

    use parameters
    implicit none
    integer :: z,graph_number
    complex*16 :: state(0:z_max)
    double precision :: C,p_cz_max

    C=0.d0
    p_cz_max=0.d0
    do z = 0,z_max
        C = C + realpart(state(z)*dconjg(state(z)))*dble(cz_vec(graph_number,z))  
        !print *, 'graph_number,z,C',graph_number,C,z
        if (cz_vec(graph_number,z) .eq. cz_max_save(graph_number)) then
            p_cz_max = p_cz_max + realpart(state(z)*dconjg(state(z)))
        endif
    enddo
end subroutine Calc_expec_C

subroutine Apply_Ub(state,angle)
    use parameters
    implicit none
    integer :: n,z
    double precision :: angle
    complex*16 :: temporary
    complex*16 :: state(0:z_max)
    complex*16 :: Ub(0:1,0:1)

    !apply Ub to all qubits
    !for each qubit, loop over coupled pairs of states
    !where the basis state of qubit n changes 0->1
    !and all other qubit basis states are constant.
    !The Ub operators are assumed to commute so the
    !order the n are executed in doesn't matter  

    Ub(0,0) = dcmplx(dcos(angle),0.d0)
    Ub(0,1) = dcmplx(0.d0,-dsin(angle))
    Ub(1,0) = dcmplx(0.d0,-dsin(angle))
    Ub(1,1) = dcmplx(dcos(angle),0.d0)

    do n = 0,n_qubits-1

        do z = 0, 2**(n_qubits-1)-1
            temporary = Ub(0,0)*state(pairs_list(n,z,0)) + Ub(0,1)*state(pairs_list(n,z,1))
            state(pairs_list(n,z,1)) = Ub(1,1)*state(pairs_list(n,z,1)) + Ub(1,0)*state(pairs_list(n,z,0))
            state(pairs_list(n,z,0)) = temporary
        enddo

    enddo

end subroutine Apply_Ub

subroutine Apply_Uc(state,c_angle,graph_number)

    use parameters
    implicit none
    integer :: z,graph_number
    complex*16 :: state(0:z_max)
    double precision :: c_angle

    do z = 0,z_max
        state(z) = state(z)*dcmplx(dcos(dble(cz_vec(graph_number,z))*c_angle),-dsin(dble(cz_vec(graph_number,z))*c_angle))
    enddo

end subroutine Apply_Uc

end module QAOA_subroutines_UnitaryEv
