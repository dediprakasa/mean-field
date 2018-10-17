program main

	use mpi
	use global

	implicit none
	
! For use in MPI environment

	integer, parameter :: nprocmax = 200
	integer :: ierr, myid, nproc, numtasks, numextra, sourceid, destid
	integer :: nstart(0:nprocmax-1), nfinish(0:nprocmax-1), numdata(0:nprocmax-1)
	integer :: status(MPI_STATUS_SIZE)

! Declaration of variables

	real*8, external :: fermi, rtbis, integrand
	integer :: i, j, counter, l, m, n, ll, mm, nn, iconv
	real*8 :: time_start, time_finish, time_elapsed
	character*50 :: DOS_up_file, DOS_down_file, DOS_total_file, mu_file, n_file

	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

! Invoke master process
	
	if (myid == 0) then
		time_start = real(MPI_Wtime(), 4)
	
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
	end if

		!call MPI_BCAST(void *buffer, int count, MPI_Datatype datatype, int root,MPI_Comm comm, ierr)
        call MPI_BCAST(N1D, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
        call MPI_BCAST(eps, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(t1, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(t2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(t3, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(Nw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(wmin, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(wmax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nfilling, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(U, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(T, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(eta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(Niter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(tol_mu, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(tol_sigma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(DOS_up_file, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(DOS_down_file, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(DOS_total_file, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mu_file, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(n_file, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

	Nsites = N1D**3
	Nfilling = Nfilling*real(2*Nsites, 8)

	if (myid == 0) write(*,*) "Finished reading input_main."

	call allocation

! Generate weighting factor for Simpson's composite integration 

	dw = (wmin - wmax)/real(Nw-1,8)

	do i = 1, Nw
		w(i) = wmin + dw*(i-1)
	end do

	wfreq(1) = dw/3.d0

	do i = 2, Nw-1, 2
		wfreq(i) = 4.d0*dw/3.d0
	end do

	do i = 3, Nw-2, 2
		wfreq(1) = 2.d0*dw/3.d0
	end do

	wfreq(Nw) = dw/3.d0
	
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
			end do
		end do
	end do

	H(:,:) = 0
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
				if ((abs(l -ll) == 1 .and. abs(m - mm) == 0 .and. abs(n - nn) == 0) &
					& .or. (abs(l - ll) == 0 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0) &
					& .or. (abs(l - ll) == 0 .and. abs(m - mm) == 0 .and. abs(n - nn) == 1)) then
					H(i, j) = -t1
				else if ((abs(l - ll) == 1 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0) &
					& .or. (abs(l - ll) == 1 .and. abs(m - mm) == 0 .and. abs(n - nn) == 1) &
					& .or. (abs(l - ll) == 0 .and. abs(m - mm) == 1 .and. abs(n - nn) == 1)) then
					H(i, j) = -t2
				else if ((abs(l - ll) == 1 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0)) then
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

	! For U>0, setup initial occupation number with AFM configuration
	n_up(coord(0, 0, 0))   = 1.d0
	n_down(coord(0, 0, 0)) = abs(1.d0 - n_up(coord(1, 1, 1)))
	do n = 1, N1D - 1
		n_up(coord(0, 0 ,n))   = abs(1.d0 - n_up(coord(0, 0, n-1)))
		n_down(coord(0, 0 ,n)) = abs(1.d0 - n_up(coord(0, 0, n)))
	end do
	do m = 1, N1D - 1
	do n = 0 ,N1D - 1
		n_up(coord(0, m, n))   = abs(1.d0 - n_up(coord(0, m-1, n)))
		n_down(coord(0, m, n)) = abs(1.d0 - n_up(coord(0, m, n)))
	end do
	end do
	do l = 1, N1D - 1
	do m = 0, N1D - 1
	do n = 0, N1D - 1
		n_up(coord(l, m, n))   = abs(1.d0 - n_up(coord(l-1, m, n)))
		n_down(coord(l, m, n)) = abs(1.d0 - n_up(coord(l, m, n)))
	end do
	end do
	end do

	
	if (U == 0.d0) then
		n_up = 0.5d0
		n_down = 0.5d0
	end if
	
	do j = 1, Nsites
		sigma(j, j, :) 			     = U*n_down(j)
		sigma(j+Nsites, j+Nsites, :) = U*n_up(j)
	end do
	
! Arrange the job sharing by distributing the real frequency

	numtasks = nw/nproc
	numextra = nw - numtasks*nproc
	
	if (numextra .gt. 0) then
		do i = 0, nproc- 1 - numextra
			nstart(i)	= i*numtasks + 1
			nfinish(i)	= nstart(i) + numtasks - 1
		end do
		
		do i = nproc - numextra, nproc - 1
			nstart(i)	= nfinish(i-1) + 1
			nfinish(i)	= nstart(i) + numtasks
		end do
	
	else
		do i = 0, nproc - 1
			nstart(i)	= i*numtasks + 1
			nfinish(i)	= nstart(i) + numtasks - 1
		end do
	end if

	do i=0,nproc-1
		numdata(i) = nfinish(i) - nstart(i) + 1
	end do 	
	
	
	if (myid == 0) then
		open(unit=55, file='check_jobs_distribution', status='unknown')
			write(55,*) 'Self consistency process with Nw =', Nw
				do i = 0, nproc - 1
					write(55,*) 'id =', 'nstart =', nstart(i), 'nfinish =', nfinish(i)
				end do
		close(55)
		write(*,*) 'Finished writing check_jobs_distribution for the self consistency process!'
	end if
	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	
!	Start iteration
	if (U == 0.d0) Niter = 1
	
	iconv = 0
	
	do i = 1, Niter
		if (myid == 0) then
			write(*,*)
			write(*,*) 'Iteration number', i
		end if
		
		do n = nstart(myid), nfinish(myid) ! Iterating over frequency omega	
			Gdummy(:,:) = (w(n) + im*eta)*id_mat(:,:) - H(:,:) - sigma(:,:,n)
			
			! Matrix inversion
			call zgetrf(2*Nsites, 2*Nsites, Gdummy, 2*Nsites, ipiv, info)
			call zgetri(2*Nsites, Gdummy, 2*Nsites, ipiv, work, 2*2*Nsites, info)

			G(:,:,n) = Gdummy(:,:)
			
			do j = 1, 2*Nsites
				PDOS(j,n) = -(1.d0/pi)*aimag((G(j,j,n)))
			end do
		end do
		
		if (myid /= 0) then
			call MPI_SSEND(PDOS(1,nstart(myid)), 2*Nsites*numdata(myid), &
			&				MPI_DOUBLE_PRECISION, 0, 111, 				 &
			&				MPI_COMM_WORLD, ierr)
		else
			do sourceid = 1, nproc - 1
				! call MPI_RECV(start_address, count, datatype, source rank, tag, communicator, status, ierr)
				 call MPI_RECV(PDOS(1,nstart(sourceid)), 2*Nsites*numdata(sourceid), &
       			 &       		MPI_DOUBLE_PRECISION, sourceid, 111,				&
       			 &				MPI_COMM_WORLD, status, ierr)
       		end do
		 end if
			
		if (myid == 0) then
			norm = 0.d0
			do n = 1, Nw
				DOS_up(n) 	= 0.d0
				DOS_down(n) = 0.d0
				
				do j = 1, Nsites
					DOS_up(n) 	= DOS_up(n) + PDOS(j, n)
					DOS_down(n) = DOS_down(n) + PDOS(j+Nsites, n)			
				end do
				
				DOS_total(n) = DOS_up(n) + DOS_down(n)
				norm 		 = norm + DOS_total(n)*wfreq(n)
			end do
		
			! Calculate chemical potential, mu
			mu_data = rtbis(integrand, wmin, wmax, tol_mu)
			write(*,*) 'mu =', mu

			! Calculate occupation number <n>
			n_number(:) = 0.d0
			do j = 1, 2*Nsites
				do n = 1, Nw
					n_number(j) =  n_number(j) + PDOS(j, n)*fermi(mu, T, w(n))*wfreq(n)
				end do
			end do
		end if
				
		call MPI_BCAST(n_number(1), 2*Nsites,	 &
		&				MPI_DOUBLE_PRECISION, 0, &
		&				MPI_COMM_WORLD, ierr)	
		
		sigma_calc(:,:,:) = 0.d0
		error_sigma_init  = 0.d0
		do n = nstart(myid), nfinish(myid)	
			do j = 1, Nsites
				sigma_calc(j, j, n) 			 = n_number(j+Nsites)*U
				sigma_calc(j+Nsites,j+Nsites, n) = n_number(j)*U
			end do
			
			! Error test
			del_sigma(:,:) 	 = abs(sigma_calc(:,:,n) - sigma(:,:,n))
			error_sigma_init = max(error_sigma, maxval(del_sigma))
		end do			
		
		call MPI_REDUCE(error_sigma_init, error_sigma, 1,	  	  &
		& 				MPI_DOUBLE_PRECISION, MPI_MAX, 0, &
		&				MPI_COMM_WORLD, ierr)
		
		if (myid == 0) then
			write(*,*) 'Error =', error_sigma
			write(*,*) 'Norm =', norm
			open(unit=10, file=DOS_Total_file, status='unknown')
			open(unit=11, file=DOS_up_file,	   status='unknown')
			open(unit=12, file=DOS_down_file,  status='unknown')
			open(unit=13, file=mu_file, 	   status='unknown')	
			
			do n = 1, Nw
				write(10,*) w(n), DOS_total
				write(11,*) w(n), DOS_up
				write(12,*) w(n), DOS_down
			end do
			
			write(13,*) mu

			close(10)
			close(11)
			close(12)
			close(13)
			
			if (error_sigma .lt. tol_sigma) then
				write(*,*) 'Convergence achieved!'
				iconv = 1
			end if
		end if
		
		call MPI_BCAST(iconv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
		
		if (iconv == 0) then
			sigma(:,:,:) = alpha*sigma_calc(:,:,:) + (1 - alpha)*sigma(:,:,:)
		else
			exit
		end if	
	end do ! End of self consistency iteration			
	

! Reconstruct the Green Function using the already converged <n(j)>
	do n = 1, Nw ! Iterating over frequency omega
		sigma_calc(:,:,:) = 0.d0
		do j = 1, Nsites
			sigma_calc(j, j, n) 			 = n_number(j+Nsites)*U
			sigma_calc(j+Nsites,j+Nsites, n) = n_number(j)*U
		end do	
		
		Gdummy(:,:) = (w(n) + im*eta)*id_mat(:,:) - H(:,:) - sigma(:,:,n)
		
		! Matrix inversion
		call zgetrf(2*Nsites, 2*Nsites, Gdummy, 2*Nsites, ipiv, info)
		call zgetri(2*Nsites, Gdummy, 2*Nsites, ipiv, work, 2*2*Nsites, info)

		G(:,:,n) = Gdummy(:,:)
			
		do j = 1, 2*Nsites
			PDOS(j,n) = -(1.d0/pi)*aimag(G(j,j,n))
		end do
	end do

	call Optical_Response(myid)
	
	if (myid == 0) then
        time_finish=real(MPI_Wtime(),4)
        time_elapsed= time_finish-time_start

        write(*,*) 'All calculations finished :) !!'
        write(*,*)
        write(*,*) 'Elapsed time (seconds) =', time_elapsed
        write(*,*) 'Elapsed time (minutes) =', time_elapsed/60.0
        write(*,*) 'Elapsed time (hours) =', time_elapsed/3600.0
		write(*,*)
	end if
	
	call MPI_FINALIZE(ierr)
	
	stop
	
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













