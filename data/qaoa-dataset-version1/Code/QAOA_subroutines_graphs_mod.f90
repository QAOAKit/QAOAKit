module QAOA_subroutines_graphs

contains 




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

subroutine Count_Num_graphs

	use parameters 
	implicit none

	!the below counts an extra time on the last call
	!to read_adjacency_matrix
	graph_num_tot=-1
	do while (more_graphs)
		call Read_Adjacency_Matrix(.true.)
		graph_num_tot = graph_num_tot+1
	enddo
	more_graphs=.true.

end subroutine Count_num_graphs


subroutine Read_Adjacency_Matrix(first_read)

	use parameters
	implicit none

	integer :: n,nn,m
	integer :: iostatus,edge_count
	character (len=n_qubits-1) :: line
	logical :: non_graph_line,first_read

	non_graph_line=.true.
	edge_count=-1
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
			if (.not. first_read) then
				do nn = n+1,n_qubits-1
					if (line(nn-n:nn-n) .eq. '1') then
						edge_count = edge_count+1
						edges(graph_num,edge_count,0) = n
						edges(graph_num,edge_count,1) = nn
					endif
				enddo
			endif
		enddo
	endif
	if (.not. first_read) n_edges(graph_num)=edge_count+1

end subroutine Read_Adjacency_Matrix

subroutine calc_vertex_degrees(graph_number)

	use parameters

	implicit none

	integer :: n,q0,q1,graph_number

	do n = 0,n_edges(graph_number)-1
		q0 = edges(graph_number,n,0)
		q1 = edges(graph_number,n,1)
		vertex_degrees(graph_number,q0) = vertex_degrees(graph_number,q0)+1
		vertex_degrees(graph_number,q1) = vertex_degrees(graph_number,q1)+1
	enddo

	all_odd_degree(graph_number)=.true.
	all_even_degree(graph_number)=.true.
	do n = 0,n_qubits-1
		if (mod(vertex_degrees(graph_number,n),2) .eq. 0) all_odd_degree(graph_number)=.false.
		if (mod(vertex_degrees(graph_number,n),2) .eq. 1) all_even_degree(graph_number)=.false.
	enddo

end subroutine calc_vertex_degrees

subroutine find_all_edges

	use parameters
	implicit none


	do graph_num = 1, graph_num_tot

		call Read_Adjacency_Matrix(.false.)

	enddo

end subroutine find_all_edges

end module QAOA_subroutines_graphs