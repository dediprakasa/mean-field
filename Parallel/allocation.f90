subroutine allocation

	use global
	implicit none
	
	allocate(wfreq(Nw))
	allocate(id_mat(2*Nsites, 2*Nsites))
	allocate(coord(0:N1D-1, 0:N1D-1, 0:N1D-1))
	allocate(H(2*Nsites, 2*Nsites))
	allocate(w(Nw))
	allocate(sigma(2*Nsites, 2*Nsites, Nw))
	allocate(n_up(Nsites))
	allocate(n_down(Nsites))
	allocate(Gdummy(2*Nsites, 2*Nsites))
	allocate(work(2*2*Nsites))
	allocate(ipiv(2*Nsites))
	allocate(G(2*Nsites, 2*Nsites, Nw))
	allocate(PDOS(2*Nsites, Nw))
	allocate(DOS_up(Nw))
	allocate(DOS_down(Nw))
	allocate(DOS_total(Nw))
	allocate(n_number(2*Nsites))
	allocate(sigma_calc(2*Nsites, 2*Nsites, Nw))
	allocate(del_sigma(2*Nsites, 2*Nsites))
	
	return
	
end subroutine allocation
