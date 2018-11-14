program main

    use global

    implicit none	

! Declaration of variables

    real*8, external :: fermi, rtbis, integrand
    integer :: i, j, counter, l, m, n, ll, mm, nn, iconv, iw
    real*8 :: time_start, time_finish, time_elapsed
    character*50 :: DOS_up_file, DOS_down_file, DOS_Total_file, mu_file, n_file

! Invoke master process
        
    open(unit=10, file="input_main", status="old")
    read(10,*) N1D
    read(10,*) eps
    read(10,*) t1
    read(10,*) t2
    read(10,*) t3
    read(10,*) Nw
    read(10,*) wmin
    read(10,*) wmax
    read(10,*) Nfilling
    read(10,*) U
    read(10,*) T
    read(10,*) eta
    read(10,*) alpha
    read(10,*) Niter
    read(10,*) tol_mu
    read(10,*) tol_sigma
    read(10,*) DOS_up_file
    read(10,*) DOS_down_file
    read(10,*) DOS_total_file
    read(10,*) mu_file
    read(10,*) n_file

    close(10)
    write(*,*) "Finished reading input_main."

    Nsites = N1D**3
    write(*,*) "Nsites", nsites
    Nfilling = Nfilling*real(2*Nsites, 8)
    write(*,*) "Nfilling", nfilling

    call allocation

! Generate weighting factor for Simpson's composite integration 

    dw = (wmax - wmin)/real(Nw-1,8)

    wfreq(1) = dw/3.d0

    do i = 2, Nw-1, 2
        wfreq(i) = 4.d0*dw/3.d0
    end do

    do i = 3, Nw-2, 2
        wfreq(i) = 2.d0*dw/3.d0
    end do

    wfreq(Nw) = dw/3.d0

    do i = 1, Nw
        w(i) = wmin + dw*(i-1)
    end do
    
! Generate identity matrix

    id_mat(:,:) = 0.d0

    do i = 1, 2*Nsites
        id_mat(i,i) = 1.d0
    end do

! Generate Hamiltonian matrix	
    counter = 0
    do l = 0, N1D - 1
    do m = 0, N1D - 1
    do n = 0, N1D - 1
        coord(l, m, n) = counter + 1
        counter = counter + 1
    end do
    end do
    end do

    H(:,:) = 0.d0
    do l = 0, N1D - 1
    do m = 0, N1D - 1
    do n = 0, N1D - 1
        do ll = 0, N1D - 1
        do mm = 0, N1D - 1
        do nn = 0, N1D - 1
            i = coord(l, m, n)
            j = coord(ll, mm, nn)
    
            if (i == j) then
                H(i, j) = eps
            else
                if        ((abs(l - ll) == 1 .and. abs(m - mm) == 0 .and. abs(n - nn) == 0) &
                    & .or. (abs(l - ll) == 0 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0) &
                    & .or. (abs(l - ll) == 0 .and. abs(m - mm) == 0 .and. abs(n - nn) == 1)) then
                    H(i, j) = -t1
                else if   ((abs(l - ll) == 1 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0) &
                    & .or. (abs(l - ll) == 1 .and. abs(m - mm) == 0 .and. abs(n - nn) == 1) &
                    & .or. (abs(l - ll) == 0 .and. abs(m - mm) == 1 .and. abs(n - nn) == 1)) then
                    H(i, j) = -t2
                else if ((abs(l - ll) == 1 .and. abs(m - mm) == 1 .and. abs(n - nn) == 1)) then
                    H(i, j) = -t3
                else
                    H(i, j) = 0.d0
                end if
            end if
        end do
        end do
        end do
    end do
    end do
    end do
                        
    H(Nsites+1:2*Nsites, Nsites+1:2*Nsites) = H(1:Nsites, 1:Nsites)

    ! Initial sigma matrice
    sigma(:, :, :) = 0.d0
    
    do j = 1, Nsites
        n_up(j) = 1.d0
        n_down(j) = 0.d0
    end do
    
    do j = 1, Nsites
        sigma(j, j, :)               = U*n_down(j)
        sigma(j+Nsites, j+Nsites, :) = U*n_up(j) 
    end do
        
!	Start iteration
    if (U == 0.d0) Niter = 1
    
    iconv = 0
    do i = 1, Niter
    
        norm = 0.d0
        write(*,*) 'Iteration number', i
        do iw = 1, Nw ! Iterating over frequency omega	
            Gdummy(:,:) = (w(iw) + im*eta)*id_mat(:,:) - H(:,:) - sigma(:,:,iw)
            
            
            ! Matrix inversion
            call zgetrf(2*Nsites, 2*Nsites, Gdummy, 2*Nsites, ipiv, info)
            call zgetri(2*Nsites, Gdummy, 2*Nsites, ipiv, work, 2*2*Nsites, info)
            
            G(:,:,iw) = Gdummy(:,:)


            traceup = 0.d0
            tracedown = 0.d0

            do j = 1, Nsites
                traceup = traceup + aimag(G(j,j,iw))
            end do
            

            do j = 1+Nsites, 2*Nsites
                tracedown = tracedown + aimag(G(j,j,iw))
            end do

            PDOS(:,:) = 0.d0
            do j=1,2*Nsites
                PDOS(j,iw) = -(1.d0/pi)*aimag(G(j,j,iw))				
            end do

            
            DOS_up(iw)    = -(1.d0/pi)*(traceup)

            DOS_down(iw)  = -(1.d0/pi)*(tracedown)
            DOS_total(iw) = -(1.d0/pi)*(traceup+tracedown)
            norm          = norm + DOS_total(iw)*wfreq(iw)
            
        
        end do
        
        ! Calculate chemical potential, mu
        mu_data = rtbis(integrand, wmin, wmax, tol_mu)
        write(*,*) 'mu =', mu_data
        write(*,*) 'Norm =', norm

        ! Calculate occupation number <n>
        n_number(:) = 0.d0
        do j = 1, 2*Nsites
            do n = 1, Nw
                n_number(j) =  n_number(j) + PDOS(j, n)*fermi(mu, T, w(n))*wfreq(n)
            end do
        end do

    
        sigma_calc(:,:,:) = 0.d0
        error_sigma_init  = 0.d0
        do n = 1, Nw
            do j = 1, Nsites
                sigma_calc(j, j, n)              = n_number(j+Nsites)*U
                sigma_calc(j+Nsites,j+Nsites, n) = n_number(j)*U
            end do
            
            ! Error test
            del_sigma(:,:) 	 = abs(sigma_calc(:,:,n) - sigma(:,:,n))
            error_sigma_init = max(error_sigma, maxval(del_sigma))
        end do
            
            
        write(*,*) 'Error =', error_sigma

        open(unit=11, file=DOS_Total_file, status='unknown')
        open(unit=12, file=DOS_up_file,	   status='unknown')
        open(unit=13, file=DOS_down_file,  status='unknown')
        open(unit=14, file=mu_file, 	   status='unknown')	
                
        do n = 1, Nw
            write(11,*) w(n), DOS_total(n)
            write(12,*) w(n), DOS_up(n)
            write(13,*) w(n), DOS_down(n)
        end do

        close(11)
        close(12)
        close(13)
        close(14)
                
        if (error_sigma .lt. tol_sigma) then
            write(*,*) 'Convergence achieved!'
            iconv = 1
        end if
            
        if (iconv == 0) then
            sigma(:,:,:) = alpha*sigma_calc(:,:,:) + (1 - alpha)*sigma(:,:,:)
        else
            exit
        end if

    end do ! End of self consistency iteration			
    
    call deallocation
    
end program main

!====================================================================

function integrand(x)

    use global
    implicit none
    integer :: i
    real*8 :: integrand, x, total
    real*8, external :: fermi
    
    total = 0.d0
    do i = 1, Nw
        total = total + DOS_total(i)*fermi(x,T,w(i))*wfreq(i)
    end do
    
    integrand = nfilling - (total*2*Nsites/norm)
    
    return
end function integrand













