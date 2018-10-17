program main

	use mpi
	use global

	implicit none
	
! For use in MPI environment

	integer, parameter :: nprocmax
	integer :: ierr, myid, numtasks, numextra, sourceid, destid
	integer :: nstart(0:nprocmax-1), nfinish(0:nprocmax-1), numdata(0:nprocmax-1)
	integer :: status(MPI_STATUS_SIZE)

! Declaration of variables

	real*8, external :: fermi, rtbis, integrand
	integer ::
	real*8 :: time_start, time_finish, time_elapsed
	character*50 :: DOS_Up, DOS_down, DOS_total, mu_data, n_data

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
		read(10,*) DOS_Up
		read(10,*) DOS_Down
		read(10,*) DOS_total
		read(10,*) mu_data
		read(10,*) n_data
		
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
        call MPI_BCAST(DOS_Up, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(DOS_Down, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(DOS_total, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mu_data, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(n_data, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

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

	(i,j,k) = ind
	(1,1,1) = 1
	(1,1,2) = 2
	(1,1,3) = 3
	
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
				if ((abs(l - ll) == 1 .and. abs(m - mm) == 0 .and. abs(n - nn) == 0) &
					& .or. (abs(l - ll) == 0 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0) &
					& .or. (abs(l - ll) == 0 .and. abs(m - mm) == 0 .and. abs(n - nn) == 1)) then
					H(i, j) = -t1
				else if ((abs(l - ll) == 1 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0) &
					& .or. (abs(l - ll) == 1 .and. abs(m - mm) == 0 .and. abs(n - nn) == 1) &
					& .or. (abs(l - ll) == 0 .and. abs(m - mm) == 1 .and. abs(n - nn) == 1)) then
					H(i, j) = -t2
				else if ((abs(l - ll) == 1 .and. abs(m - mm) == 1 .and. abs(n - nn) == 0)) then
					H(i, j) = -t3
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
	n_up(ind(0, 0, 0)) = 1.d0
	n_down(ind(0, 0, 0)) = abs(1.d0 - n_up(ind(1, 1, 1)))
	do w = 1, N1D-1
		n_up(ind(0, 0, n)) = abs(1.d0 - n_up(ind(0, 0, n-1)))
		n_down(ind(0, 0, n)) = abs(1.d0 - n_up(ind(0, 0, n)))
	end do
	do v = 1, N1D - 1
	do w = 0, N1D - 1
		n_up(ind(0, m, n)) = abs(1.d0 - n_up(ind(0, m-1, n)))
		n_down(ind(0, m, n)) = abs(1.d0 - n_up(ind(0, m, n)))
	end do
	end do
	do u = 1, N1D - 1
	do v = 0, N1D - 1
	do w = 0, N1D - 1
		n_up(ind(l, m, n))=abs(1.d0 - n_up(ind(l-1, m, n)))
		n_down(ind(l, m, n))=abs(1.d0 - n_up(ind(l, m, n)))
	end do
	end do
	end do

	! For U=0, setup initial occupation number with PM configuration
	if (U == 0.d0) then
		n_up = 0.5d0
		n_down = 0.5d0
	end if

	! Arrange the job sharing by distributing the real frequency jobs to all processors 
	num_tasks = Nw/nproc
	num_extra = Nw - numtasks*nproc

	if (num_extra .gt. 0) then
		do i = 0, nproc - 1 - num_extra
			nstart(i) = i*num_tasks + 1
			nfinish(i) = nstart(i) + num_tasks - 1
		end do
		
		do i = nproc - num_extra, nproc - 1
			nstart(i) = nfinish(i-1) + 1
			nfinish(i) = nstart(i) + num_tasks
		end do
	
	else
		do i = 0, nproc - 1
			nstart(i) = i*num_tasks + 1
			nfinish(i) = nstart(i) + num_tasks - 1
		end do
	end if

	do i = 0, nproc - 1
		num_data(i) = nfinish(i) - nstart(i) + 1
	end do

	! Preparation to start calculation
	if (myid == 0) then
		open(unit=5, file='check_jobs_dist', status='unknown')
		write(5,*) 'For the self-consistency process with Nw =', Nw
		
		do i = 0, nproc - 1
			write(5,*) 'id =', i, 'nstart =', nstart(i), 'nfinish =', nfinish(i)
		end do

		close(5)
		write(*,*) 'Done writing check_jobs_dist for the self-consistency process!'

	
	call MPI_Barrier(MPI_COMM_WORLD, ierr)

	! Start iteration

	if (U == 0.d0) Niter = 1

	iconv = 0

	do i = 1, Niter
		if (myid == 0) then
			write(*,*)
			write(*,*) 'Number of iteration =', i
		end if
	
		do n = nstart(myid), nfinish(myid)

			Gdummy(:,:) = (w(n) + im*eta)*id_mat(:,:) - H(:,:) - sigmac(:,:)

			! Matrix inversion
			call zgetrf(2*Nsites, 2*Nsites, Gdummy, 2*Nsites, ipiv, info)
			call zgetri(2*Nsites, Gdummy, 2*Nsites, ipiv,work, 2*2*Nsites, info)

			G(:,:,n) = Gdummy(:,:)

			do j = 1, 2*Nsites
				PDOS(j,n) = -(1.d0/pi)*aimag(G(j,j,n))
			end do
	end do

	


	















