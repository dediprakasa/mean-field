	subroutine Optical_Response(myid)
	use MPI
	use Global, only: N1D,Nsites,wmin,wmax,Nw,dw,w,G,H,coord,mu,t1,t2,t3,T,pi,im
        implicit none


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Written by: Muhammad Aziz Majidi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! For use in MPI environment

	integer, parameter :: nprocmax = 200
	integer :: ierr, myid, nproc, numtasks, numextra, sourceid, destid
	integer :: nstart(0:nprocmax-1), nfinish(0:nprocmax-1), numdata(0:nprocmax-1)
	integer :: status(MPI_STATUS_SIZE)
    
! Declaration of variables:

	integer :: i,j,k,l,m,n,lp,mp,np
        integer :: nmu,inu,inu1,inu2,numin,numax,it,inumax,isum
        integer :: iw,iwwmin,iwwmax,stepww,iwp,nwpmax,nl1,nl2
	real*8 :: temp,kT,kTmin,kTmax,dkT,delta
        real*8 :: ww,wwmin,wwmax,nu,wwp1,wwp2,wpmax,eps_infty
        real*8 :: deltanu,rsum,coeff
        real*8 :: Ae(2*Nsites,2*Nsites)
	complex*16 :: vxmat(2*Nsites,2*Nsites),vymat(2*Nsites,2*Nsites)
	complex*16 :: vxA1(2*Nsites,2*Nsites),vxA2(2*Nsites,2*Nsites)
	complex*16 :: vyA1(2*Nsites,2*Nsites),vyA2(2*Nsites,2*Nsites)
	complex*16 :: vxAvxA(2*Nsites,2*Nsites),vyAvyA(2*Nsites,2*Nsites)
	complex*16 :: vxAvyA(2*Nsites,2*Nsites)
	complex*16 :: trvxAvxA,trvyAvyA,trvxAvyA
        real*8 :: sigma1_xx(0:2*Nw),sigma1_yy(0:2*Nw),sigma1_xy(0:2*Nw)
	real*8 :: sigma2(0:2*Nw),eps1(0:2*Nw),eps2(0:2*Nw),LF(0:2*Nw)
        real*8 :: dummy_xx,dummy_yy,dummy_xy,re,ima
	real*8 :: dnu(-5*nw:5*nw)
	real*8 :: wp(0:20*nw),wwp(0:20*nw),sigma1int(0:20*nw),dwp,wk,wj,prod
        real*8, external :: fermi
	real*8 :: e,hbar,a,opt_cond_prefactor,e0

	complex*16 :: cmat(2*Nsites,2*Nsites)
	character*50 :: s1xxfile,s1yyfile,s1xyfile,s2file,eps1file,eps2file,LFfile

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


! Read parameters from the inputfile
	if (myid==0) then	
	open(unit=10,file="input_opt_response",status="old")
	read(10,*) wwmin ! in eV
	read(10,*) wwmax ! in eV
	read(10,*) stepww
	read(10,*) eps_infty
	read(10,*) s1xxfile
 	read(10,*) s1yyfile
	read(10,*) s1xyfile
	read(10,*) s2file
	read(10,*) eps1file
	read(10,*) eps2file
	read(10,*) LFfile
	close(10)
	end if
        
        call MPI_BCAST(wwmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(wwmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(stepww,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	e = 1.60217662d-19 ! in Coulomb
	hbar = 1.0545718d-34 ! in m^2 kg / s
	e0 = 8.854187817d-12 ! vacuum permittivity in farad/m
	a = 5.d-10 ! lattice constant in m
	opt_cond_prefactor=pi*e**2/(hbar*a)/real(Nsites,8) ! in Ohm^-1 m^-1


	Temp=T
        kT=8.617332478d-5*temp ! in eV

	! Defining real-frequency variable and their Simpson integration weights

	! Defining v-matrices
	vxmat=0.d0
	vymat=0.d0
	do l=0,N1D-1
	do m=0,N1D-1
	do n=0,N1D-1
		do lp=0,N1D-1
		do mp=0,N1D-1
		do np=0,N1D-1		
			if ((abs(l-lp)==1).and.(abs(m-mp)==0).and.(abs(n-np)==0)) then
				vxmat(coord(l,m,n),coord(lp,mp,np))=sign(1,l-lp)*im*H(coord(l,m,n),coord(lp,mp,np))	
			end if
			if ((abs(l-lp)==0).and.(abs(m-mp)==1).and.(abs(n-np)==0)) then
				vymat(coord(l,m,n),coord(lp,mp,np))=sign(1,m-mp)*im*H(coord(l,m,n),coord(lp,mp,np))	
			end if
		end do
		end do
		end do
	end do
	end do
	end do

	vxmat(Nsites+1:2*Nsites,Nsites+1:2*Nsites)=vxmat(1:Nsites,1:Nsites)
	vymat(Nsites+1:2*Nsites,Nsites+1:2*Nsites)=vymat(1:Nsites,1:Nsites)

	if (myid==0) then
	write(*,*)
        write(*,*) 'Calculating the Optical Conductivity ...'
        write(*,*) 'temp(K)=',T
        write(*,*) 'mu(eV)=',mu
	write(*,*)
	end if

        iwwmin=int(wwmin/dw)/stepww
        iwwmin=iwwmin*stepww

        iwwmax=int(wwmax/dw)/stepww
        iwwmax=(iwwmax+1)*stepww


! Arrange the job sharing by distributing the real frequency jobs to all processors 

	numtasks=(iwwmax-iwwmin)/stepww + 1
	rsum=0.d0
	do i=0,nproc-1
		rsum=rsum+1.d0/real(1+i,8)
	end do
	coeff=real(numtasks,8)/rsum
	do i=0,nproc-1
		numdata(i)=int(coeff/(real(1+i,8)))
	end do
	isum=0
	do i=0,nproc-1
		isum=isum+numdata(i)
	end do
	numextra=numtasks-isum	

	if (numextra.gt.0) then
		do i=0,numextra-1
			numdata(i)=numdata(i)+1
		end do
	end if

	nstart(0)=iwwmin
	nfinish(0)=iwwmin+(numdata(0)-1)*stepww
	do i=1,nproc-2
		nstart(i)=nfinish(i-1)+stepww
		nfinish(i)=nstart(i)+(numdata(i)-1)*stepww
	end do
	nstart(nproc-1)=nfinish(nproc-2)+stepww
	nfinish(nproc-1)=iwwmax

        if (myid==0) then

      	open(unit=777, file='check_jobs_dist', status='old', position='append')
	write(777,*)
      	write(777,*) 'For OptCond calculation with number of freqs =',numtasks
       	do i=0,nproc-1 
       		write(777,*) 'id=',i,'nstrt=',nstart(i),'nfin=',nfinish(i)
      	end do
       	close(777)
      	write(*,*) 'DONE writing check_jobs_dist for OptCond calculation!'
      
        write(*,*)
        write(*,*) 'From w =',iwwmin*dw,' eV'
        write(*,*) 'To   w =',iwwmax*dw,' eV'
        write(*,*) 'With increment dw =',stepww*dw,' eV'

        write(*,*)
        write(*,*) '          w(eV)         ','        opt_cond_xx      ',&
	&	'        opt_cond_yy      ','        opt_cond_xy      '
	
!	open(unit=25,file=s1xxfile,status="replace") 
!	open(unit=26,file=s1yyfile,status="replace") 
!	open(unit=27,file=s1xyfile,status="replace") 
!	! This is simply to clear the content of any existing OptCondfile of the same name. 
!	close(25)
!	close(26)
!	close(27)

	end if ! myid=0

        do n=1,nw
           if (w(n)>mu) then
              nmu=n
              exit
           end if
        end do

        if (myid==0) write(*,*)


!        delta=1.0d-3*dw !<---- probably too small
        delta=dw         !<---- try change it with this.
        
        if (kT.gt.0.0d0) then
           it=16
        else
           it=0
        end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        do iw=nstart(myid),nfinish(myid),stepww    ! Start varying photon frequency

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ww=real(iw,8)*dw 

	sigma1_xx(iw)=0.d0   
	sigma1_yy(iw)=0.d0   
	sigma1_xy(iw)=0.d0   


!        if (mu+ww+8.0d0*kT .ge. w(nw)) then
!        write (*,*) "There is not enough self-energy data to compute for" 
!        write (*,*) "w >=",ww,'eV .'
!        exit
!        end if


!        if (mu+ww+8.0d0*kT .ge. w(nw)) goto 100


!>>>>>>>>>>>>>>>>>>>>>>>  This part is for T=0:  >>>>>>>>>>>>>>>>>>>>>>>>>

        if (it .eq. 0) then


        if (iw .eq. 0) then

           dnu(0)=1.0d0

        else

           dnu(:)=dw          ! Trapezoidal integration weights
           dnu(-iw)=0.5d0*dw  ! less accurate but more flexible
           dnu(0)=0.5d0*dw    ! since the number of integration
                                       ! points vary.

!           Simpson integration weights might work better,
!           but requires odd number of points (>=5) to proceed.
!           dnu(-iw)=dw/3.0d0
!           do inu=-iw+1,-1,2
!              dnu(inu)=dw*4.0d0/3.0d0
!           end do
!           do inu=-iw+2,-2,2
!               dnu(inu)=dw*2.0d0/3.0d0
!           end do
!           dnu(0)=dw/3.0d0  

        end if


        do inu=-iw,0

        nu=mu+real(inu,8)*dw

        inu1=int( (nu-w(1))/dw ) + 1
        inu2=int( (nu+ww-w(1))/dw ) + 1

        if ((inu1 .lt. 1).or.(inu2 .ge. nw)) cycle
          
        cmat(:,:)=G(:,:,inu1)+ &
      &        (G(:,:,inu1+1)-G(:,:,inu1))/dw&
      &        *(nu-w(inu1))

	Ae=-dimag(cmat)/pi

	vxA1=matmul(vxmat,Ae)
	vyA1=matmul(vymat,Ae)

        cmat(:,:)=G(:,:,inu2)+ &
      &        (G(:,:,inu2+1)-G(:,:,inu2))/dw&
      &        *(nu+ww-w(inu2))

	Ae=-dimag(cmat)/pi

	vxA2=matmul(vxmat,Ae)
	vyA2=matmul(vymat,Ae)
        
	vxAvxA=matmul(vxA1,vxA2)
	vyAvyA=matmul(vyA1,vyA2)
	vxAvyA=matmul(vxA1,vyA2)

	trvxAvxA=0.d0
	trvyAvyA=0.d0
	trvxAvyA=0.d0

	do i=1,2*Nsites
		trvxAvxA=trvxAvxA+vxAvxA(i,i)
		trvyAvyA=trvyAvyA+vyAvyA(i,i)
		trvxAvyA=trvxAvyA+vxAvyA(i,i)
	end do

        if (iw .eq. 0) then
           dummy_xx=trvxAvxA
           dummy_yy=trvyAvyA
           dummy_xy=trvxAvyA
	else
           dummy_xx=trvxAvxA/ww
           dummy_yy=trvyAvyA/ww
           dummy_xy=trvxAvyA/ww
        end if

        sigma1_xx(iw)=sigma1_xx(iw)+dummy_xx*dnu(inu)
        sigma1_yy(iw)=sigma1_yy(iw)+dummy_yy*dnu(inu)
        sigma1_xy(iw)=sigma1_xy(iw)+dummy_xy*dnu(inu)

        end do     ! end do inu (only inu=0 done if iw=0)


!<<<<<<<<<<<<<<<<<<<<<<<<<<<   End T=0 part  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        else

!>>>>>>>>>>>>>>>>>>>>>>>  This part is for T > 0:  >>>>>>>>>>>>>>>>>>>>>>>>>


!       Case 1:  if w <= 16*kT
!       ----------------------


        if (ww .le. 16.0d0*kT) then

        dkT=(16.0d0*kT+ww)/real(2*it,8)

        dnu(0)=dkT/3.0d0
        do inu=1,2*it-1,2
           dnu(inu)=dkT*4.0d0/3.0d0
        end do
        do inu=2,2*it-2,2
            dnu(inu)=dkT*2.0d0/3.0d0
        end do
        dnu(2*it)=dkT/3.0d0


        do inu=0,2*it

        nu=mu-ww-8.0d0*kT+real(inu,8)*dkT

        inu1=int( (nu-w(1))/dw ) + 1
        inu2=int( (nu+ww-w(1))/dw ) + 1

        if ((inu1 .lt. 1).or.(inu2 .ge. nw)) cycle
          
        cmat(:,:)=G(:,:,inu1)+ &
      &        (G(:,:,inu1+1)-G(:,:,inu1))/dw&
      &        *(nu-w(inu1))

	Ae=-dimag(cmat)/pi

	vxA1=matmul(vxmat,Ae)
	vyA1=matmul(vymat,Ae)

        cmat(:,:)=G(:,:,inu2)+ &
      &        (G(:,:,inu2+1)-G(:,:,inu2))/dw&
      &        *(nu+ww-w(inu2))

	Ae=-dimag(cmat)/pi

	vxA2=matmul(vxmat,Ae)
	vyA2=matmul(vymat,Ae)

	vxAvxA=matmul(vxA1,vxA2)
	vyAvyA=matmul(vyA1,vyA2)
	vxAvyA=matmul(vxA1,vyA2)

	trvxAvxA=0.d0
	trvyAvyA=0.d0
	trvxAvyA=0.d0

	do i=1,2*Nsites
		trvxAvxA=trvxAvxA+vxAvxA(i,i)
		trvyAvyA=trvyAvyA+vyAvyA(i,i)
		trvxAvyA=trvxAvyA+vxAvyA(i,i)
	end do
 

        if (iw .eq. 0) then
        dummy_xx=trvxAvxA*fermi(mu,temp,nu)*(1.0d0-fermi(mu,temp,nu))/kT
        dummy_yy=trvyAvyA*fermi(mu,temp,nu)*(1.0d0-fermi(mu,temp,nu))/kT
        dummy_xy=trvxAvyA*fermi(mu,temp,nu)*(1.0d0-fermi(mu,temp,nu))/kT
        else
        dummy_xx=trvxAvxA*&
      &     (fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        dummy_yy=trvyAvyA*&
      &     (fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        dummy_xy=trvxAvyA*&
      &     (fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        end if

        sigma1_xx(iw)=sigma1_xx(iw)+dummy_xx*dnu(inu)
        sigma1_yy(iw)=sigma1_yy(iw)+dummy_yy*dnu(inu)
        sigma1_xy(iw)=sigma1_xy(iw)+dummy_xy*dnu(inu)


        end do     ! end do inu 


        else

!       Case 2: if w > 16*kT
!       --------------------

!       First, take care of the "tails" region 


        dkT=kT

        dnu(0)=dkT/3.0d0
        do inu=1,it-1,2
           dnu(inu)=dkT*4.0d0/3.0d0
        end do
        do inu=2,it-2,2
            dnu(inu)=dkT*2.0d0/3.0d0
        end do
        dnu(it)=dkT/3.0d0


        do m=1,2

           ! m=1 for region: mu-w-8*kT < nu < mu-w+8*kT  ("front tail")
           ! m=2 for region: mu-8*kT < nu < mu+8*kT      ("back tail")

        if (m.eq.1) then
           kTmin=mu-ww-8.0d0*kT
           kTmax=mu-ww+8.0d0*kT
        else
           kTmin=mu-8.0d0*kT
           kTmax=mu+8.0d0*kT
        end if


        do inu=0,it

        nu=kTmin+real(inu,8)*(kTmax-kTmin)/real(it,8)

        inu1=int( (nu-w(1))/dw ) + 1
        inu2=int( (nu+ww-w(1))/dw ) + 1

        if ((inu1 .lt. 1).or.(inu2 .ge. nw)) cycle
          
        cmat(:,:)=G(:,:,inu1)+ &
      &        (G(:,:,inu1+1)-G(:,:,inu1))/dw&
      &        *(nu-w(inu1))

	Ae=-dimag(cmat)/pi

	vxA1=matmul(vxmat,Ae)
	vyA1=matmul(vymat,Ae)

        cmat(:,:)=G(:,:,inu2)+ &
      &        (G(:,:,inu2+1)-G(:,:,inu2))/dw&
      &        *(nu+ww-w(inu2))

	Ae=-dimag(cmat)/pi

	vxA2=matmul(vxmat,Ae)
	vyA2=matmul(vymat,Ae)
        
	vxAvxA=matmul(vxA1,vxA2)
	vyAvyA=matmul(vyA1,vyA2)
	vxAvyA=matmul(vxA1,vyA2)

	trvxAvxA=0.d0
	trvyAvyA=0.d0
	trvxAvyA=0.d0

	do i=1,2*Nsites
		trvxAvxA=trvxAvxA+vxAvxA(i,i)
		trvyAvyA=trvyAvyA+vyAvyA(i,i)
		trvxAvyA=trvxAvyA+vxAvyA(i,i)
	end do
      
        
        dummy_xx=trvxAvxA*(fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        dummy_yy=trvyAvyA*(fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        dummy_xy=trvxAvyA*(fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww

        sigma1_xx(iw)=sigma1_xx(iw)+dummy_xx*dnu(inu)
        sigma1_yy(iw)=sigma1_yy(iw)+dummy_yy*dnu(inu)
        sigma1_xy(iw)=sigma1_xy(iw)+dummy_xy*dnu(inu)

        end do     ! end do inu 

        end do     ! end do m


!       Now take care of the integration region due to finite w, 
!       i.e. the region : mu-w+8*kT < nu <  mu-8*kT

        if ( (ww-16.0d0*kT) .le. max(dw,kT) ) then 
           inumax=2
        else
           inumax=2*(1+int( (ww-16.0d0*kT)/(max(dw,kT)) ))
        end if
       
        deltanu=(ww-16.0d0*kT)/real(inumax,8)

        dnu(0)=deltanu/3.0d0
        do inu=1,inumax-1,2
           dnu(inu)=deltanu*4.0d0/3.0d0
        end do
        do inu=2,inumax-2,2
           dnu(inu)=deltanu*2.0d0/3.0d0
        end do
        dnu(inumax)=deltanu/3.0d0  


        do inu=0,inumax
        
        nu= mu-ww+8.0d0*kT+real(inu,8)*deltanu

        inu1=int( (nu-w(1))/dw ) + 1
        inu2=int( (nu+ww-w(1))/dw ) + 1

        if ((inu1 .lt. 1).or.(inu2 .ge. nw)) cycle
          
        cmat(:,:)=G(:,:,inu1)+ &
      &        (G(:,:,inu1+1)-G(:,:,inu1))/dw&
      &        *(nu-w(inu1))

	Ae=-dimag(cmat)/pi

	vxA1=matmul(vxmat,Ae)
	vyA1=matmul(vymat,Ae)

        cmat(:,:)=G(:,:,inu2)+ &
      &        (G(:,:,inu2+1)-G(:,:,inu2))/dw&
      &        *(nu+ww-w(inu2))

	Ae=-dimag(cmat)/pi

	vxA2=matmul(vxmat,Ae)
	vyA2=matmul(vymat,Ae)
        
	vxAvxA=matmul(vxA1,vxA2)
	vyAvyA=matmul(vyA1,vyA2)
	vxAvyA=matmul(vxA1,vyA2)

	trvxAvxA=0.d0
	trvyAvyA=0.d0
	trvxAvyA=0.d0
	do i=1,2*Nsites
		trvxAvxA=trvxAvxA+vxAvxA(i,i)
		trvyAvyA=trvyAvyA+vyAvyA(i,i)
		trvxAvyA=trvxAvyA+vxAvyA(i,i)
	end do


        dummy_xx=trvxAvxA*&
      &     (fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        dummy_yy=trvyAvyA*&
      &     (fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww
        dummy_xy=trvxAvyA*&
      &     (fermi(mu,temp,nu)-fermi(mu,temp,nu+ww))/ww

        sigma1_xx(iw)=sigma1_xx(iw)+dummy_xx*dnu(inu)
        sigma1_yy(iw)=sigma1_yy(iw)+dummy_yy*dnu(inu)
        sigma1_xy(iw)=sigma1_xy(iw)+dummy_xy*dnu(inu)


        end do         ! end do inu


        end if	! (ww .le. 16.0d0*kT or not)


        end if	! (T=0 or T>0)


!<<<<<<<<<<<<<<<<<<<<<<<<<<<   End T > 0 part  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	if (myid==0) then

        	write(*,*) ww, 1.d-5*opt_cond_prefactor*sigma1_xx(iw),&
		& 	1.d-5*opt_cond_prefactor*sigma1_yy(iw),&
		&	1.d-5*opt_cond_prefactor*sigma1_xy(iw)
		! Here, sigma1 values are stored in unit of 10^3 Ohm cm^-1

	end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        end do ! iw    ! End variation of photon frequency
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       	call MPI_Barrier(MPI_COMM_WORLD,ierr)

	if (myid.ne.0) then
           call MPI_SSEND(sigma1_xx(nstart(myid)),numdata(myid),&
       &          MPI_DOUBLE_PRECISION,0,110,MPI_COMM_WORLD,ierr)
           else
        do sourceid=1, nproc-1
          call MPI_RECV(sigma1_xx(nstart(sourceid)),numdata(sourceid),&
       &       MPI_DOUBLE_PRECISION,sourceid,110,MPI_COMM_WORLD,status,ierr)
        end do
	end if

	if (myid.ne.0) then
           call MPI_SSEND(sigma1_yy(nstart(myid)),numdata(myid),&
       &          MPI_DOUBLE_PRECISION,0,120,MPI_COMM_WORLD,ierr)
           else
        do sourceid=1, nproc-1
          call MPI_RECV(sigma1_yy(nstart(sourceid)),numdata(sourceid),&
       &       MPI_DOUBLE_PRECISION,sourceid,120,MPI_COMM_WORLD,status,ierr)
        end do
	end if

	if (myid.ne.0) then
           call MPI_SSEND(sigma1_xy(nstart(myid)),numdata(myid),&
       &          MPI_DOUBLE_PRECISION,0,130,MPI_COMM_WORLD,ierr)
           else
        do sourceid=1, nproc-1
          call MPI_RECV(sigma1_xy(nstart(sourceid)),numdata(sourceid),&
       &       MPI_DOUBLE_PRECISION,sourceid,130,MPI_COMM_WORLD,status,ierr)
        end do
	end if

!---------------------------------------------------------------------------------------
!	The rest of the calculation process below is done by master only 
!---------------------------------------------------------------------------------------

	if (myid==0) then

	do iw=nstart(1),iwwmax,stepww
       	ww=real(iw,8)*dw 
        	write(*,*) ww, 1.d-5*opt_cond_prefactor*sigma1_xx(iw),&
		& 	1.d-5*opt_cond_prefactor*sigma1_yy(iw),&
		&	1.d-5*opt_cond_prefactor*sigma1_xy(iw)
		! Here, sigma1 values are stored in unit of 10^3 Ohm cm^-1

	end do

       	open(unit=30,file=s1xxfile,status='replace')
        open(unit=31,file=s1yyfile,status='replace')
        open(unit=32,file=s1xyfile,status='replace')
	do iw=iwwmin,iwwmax,stepww
        ww=real(iw,8)*dw 
		! Here, sigma1 values are stored in unit of 10^3 Ohm cm^-1
        	write(30,*) ww, 1.d-5*opt_cond_prefactor*sigma1_xx(iw)
        	write(31,*) ww, 1.d-5*opt_cond_prefactor*sigma1_yy(iw)
        	write(32,*) ww, 1.d-5*opt_cond_prefactor*sigma1_xy(iw)
	end do
	close(30)
	close(31)
	close(32)


! Performing Lagrange interpolation of the optical conductivity data

	write(*,*)
	write(*,*) 'Interpolating sigma1_xx(w) ...'
	open(unit=33,file='sigma1int.dat',status='unknown')

	wpmax=iwwmax*dw
	nwpmax=8*iwwmax
	dwp=wpmax/real(nwpmax,8)
	do n=0,nwpmax
		wp(n)=n*dwp
	end do

	! Do linear interpolation
	do n=0,nwpmax
	sigma1int(n)=sigma1_xx(0)+(sigma1_xx(1)-sigma1_xx(0))/dw*wp(n)
	nl1=n
	if (wp(n).gt.dw) exit
	end do
	do n=nwpmax,0,-1
	sigma1int(n)=sigma1_xx(iwwmax)-(sigma1_xx(iwwmax)-sigma1_xx(iwwmax-1))/&
	&	dw*(iwwmax*dw-wp(n))
	nl2=n
	if (wp(n).lt.((iwwmax-1)*dw)) exit
	end do

	! Do cubic Lagrange interpolation

	do n=nl1+1,nl2-1
		sigma1int(n)=0.d0
		i=int(wp(n)/dw)
		if (i.lt.1) cycle
		if (i.gt.iwwmax-1) exit
		do j=i-1,i+2
			wj=j*dw
			prod=1.d0
			do k=i-1,i+2
				if (k==j) cycle
				wk=k*dw
				prod=prod*(wp(n)-wk)/(wj-wk)
			end do
			sigma1int(n)=sigma1int(n)+prod*sigma1_xx(j)
		end do
	end do

	do n=0,nwpmax
		write(33,*) wp(n),1.d-5*opt_cond_prefactor*sigma1int(n) ! sigma1int values are stored in unit of 10^3 Ohm cm^-1
	end do

	close(33)

	write(*,*) 'Interpolation of sigma1_xx(w) DONE.'
	write(*,*)


! Defining Simpson integration weights for computing sigma2 through KK transformation

	write(*,*) 'Calculating sigma2(w) ...'
	open(unit=35,file=s2file,status='unknown')

	wwp(0)=dwp/3.d0
	do n=1,nwpmax-1,2
		wwp(n)=4.d0*dwp/3.d0
	end do
	do n=2,nwpmax-2,2
		wwp(n)=2.d0*dwp/3.d0
	end do
	wwp(nwpmax)=dwp/3.d0


! Now performing Kramers-Kronig transformation to calculate sigma2(w)
! Here, we apply the KK formula:
! sigma2(w)=-(2*w/pi) integral(from 0 to infnty) sigma(w')/(w'^2-w^2) dw'

	sigma2(0)=0.d0

	do iw=1,iwwmax

	ww=real(iw,8)*dw
	sigma2(iw)=0.d0

	do n=0,nwpmax
	if (wp(n)==ww) cycle
	sigma2(iw)=sigma2(iw)-(2.d0*ww/pi)*sigma1int(n)/(wp(n)**2 - ww**2)*wwp(n) ! sigma2 values are stored in unit of 10^3 Ohm cm^-1
	end do

	end do

	do iw=0,iwwmax
	ww=real(iw,8)*dw
	write(35,*) ww,1.d-5*opt_cond_prefactor*sigma2(iw)
	end do

	close(35)

	write(*,*) 'Calculation of sigma2(w) DONE.'
	write(*,*)

! Now calculating epsilon1(w): epsilon1(w) = 1.d0-sigma2(w)/(epsilon0*ww)

	write(*,*) 'Calculating epsilon1(w) ...'
	open(unit=40,file=eps1file,status='unknown')

	do iw=iwwmax,1,-1
	ww=real(iw,8)*dw
	eps1(iw)=1.d0-opt_cond_prefactor*sigma2(iw)/(e0*ww*1.5193d15) ! 1 eV = 1.5193d15 rad/s
	end do

	eps1(0)=2.d0*eps1(1)-eps1(2) ! eps1(0) is obtained by extrapolation

	do iw=0,iwwmax
	ww=real(iw,8)*dw
	write(40,*) ww,eps1(iw)
	end do

	close(40)

	write(*,*) 'Calculation of epsilon1(w) DONE.'
	write(*,*)

! Now calculating epsilon2(w): epsilon2(w) = sigma1_xx(w)/(epsilon0*ww)

	write(*,*) 'Calculating epsilon2(w) ...'
	open(unit=45,file=eps2file,status='unknown')

	eps2(0)=0.d0
	do iw=1,iwwmax
	ww=real(iw,8)*dw
	eps2(iw)=opt_cond_prefactor*sigma1_xx(iw)/(e0*ww*1.5193d15) ! 1 eV = 1.5193d-15 rad/s
	end do

	do iw=0,iwwmax
	ww=real(iw,8)*dw
	write(45,*) ww,eps2(iw)
	end do

	close(45)

	write(*,*) 'Calculation epsilon2(w) DONE.'
	write(*,*)

! Now calculating Loss Function(w): LF(w) = epsilon2(w)/(epsilon1^2(w)+epsilon2^2(w))

	write(*,*) 'Calculating Loss Function ...'
	open(unit=50,file=LFfile,status='unknown')

	do iw=0,iwwmax

	ww=real(iw,8)*dw
	LF(iw) = eps2(iw)/(eps1(iw)**2 + eps2(iw)**2)

	write(50,*) ww,LF(iw)
	end do

	close(50)

	write(*,*) 'Calculation of Loss Function DONE.'
	write(*,*)

	end if ! myid

!---------------------------------------------------------------------------------------
!	END of the block of the last processes done by master only 
!---------------------------------------------------------------------------------------

       	call MPI_Barrier(MPI_COMM_WORLD,ierr)

	RETURN
	
	END SUBROUTINE Optical_Response


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






