	subroutine allocation
	
	use global
	implicit none

	allocate(Gdummy(128,128))
	allocate(G(128,128,Nw))
	allocate(IdMat(128,128))
	allocate(H(128,128))
	allocate(w(Nw))
	allocate(DoS(Nw))
	allocate(PDOS(128,Nw))
	allocate(sigma(128,128))
	allocate(wfreq(Nw))

	return
	
	end subroutine allocation


	subroutine deallocation
	
	use global
	implicit none

	deallocate(Gdummy)
	deallocate(G)
	deallocate(IdMat)
	deallocate(H)
	deallocate(w)
	deallocate(DoS)
	deallocate(PDOS)
	deallocate(sigma)
	deallocate(wfreq)

	return
	
	end subroutine deallocation
