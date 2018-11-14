subroutine deallocation

    use global
    implicit none
    
    deallocate(wfreq)
    deallocate(id_mat)
    deallocate(coord)
    deallocate(H)
    deallocate(sigma)
    deallocate(n_up)
    deallocate(n_down)
    deallocate(Gdummy)
    deallocate(work)
    deallocate(ipiv)
    deallocate(G)
    deallocate(DOS_up)
    deallocate(DOS_down)
    deallocate(DOS_total)
    deallocate(n_number)
    deallocate(sigma_calc)
    deallocate(del_sigma)
    
    return
    
end subroutine deallocation
