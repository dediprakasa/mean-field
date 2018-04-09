	program Percobaan
	
	use global
	implicit none
	
	REAL*8,EXTERNAL::fermi,rtbis,integrand
	integer, allocatable::ind(:,:,:)
	
!	Arrays and variables for LAPACK
	complex*16, dimension(2*128)::work
	integer, dimension(128)::ipiv
	integer::info
	
	eps=0.D0
	t1=1.D0
	t2=0.D0
	Nw=501
	wmin=-5.D0
	wmax=5.D0
	nfilling=64.D0
		
	call allocation
	allocate(ind(0:3,0:3,0:3))	
	!	Real frequency Simpson weights
	dw=(wmax-wmin)/REAL(Nw-1,8)
	wfreq(1)=dw/3.D0
	DO i=2,Nw-1,2
		wfreq(i)=4.D0*dw/3.D0
	END DO
	DO i=3,Nw-2,2
		wfreq(i)=2.D0*dw/3.D0
	END DO
	wfreq(Nw)=dw/3.D0


	! Matriks sigma awal
	sigma(:,:)=(0.D0,0.D0)


	
	!	Membuat Matriks Identitas
	IdMat(:,:)=0.D0
	do i=1,128
		IdMat(i,i)=1.D0
	end do
	
	!	Generate Hamiltonian matrix
	counter=0
	DO l=0,3
	DO m=0,3
	DO n=0,3
		ind(l,m,n)=counter+1
		counter=counter+1
	END DO
	END DO
	END DO
	
	H(:,:)=0.D0
	DO l=0,3
	DO m=0,3
	DO n=0,3
		DO lp=0,3
		DO mp=0,3
		DO np=0,3
			i=ind(l,m,n)
			j=ind(lp,mp,np)
			
			IF(i==j) THEN
				H(i,j)=eps
			ELSE
				IF((ABS(l-lp)==1 .AND. ABS(m-mp)==0 .AND. ABS(n-np)==0) .OR. &
					&(ABS(m-mp)==1 .AND. ABS(l-lp)==0 .AND. ABS(n-np)==0) .OR. &
					&(ABS(n-np)==1 .AND. ABS(l-lp)==0 .AND. ABS(m-mp)==0)) THEN
					H(i,j)=-t1
				ELSE IF((ABS(l-lp)==2 .AND. ABS(m-mp)==0 .AND. ABS(n-np)==0) .OR. &
					&(ABS(m-mp)==2 .AND. ABS(l-lp)==0 .AND. ABS(n-np)==0) .OR. &
					&(ABS(n-np)==2 .AND. ABS(l-lp)==0 .AND. ABS(m-mp)==0)) THEN
					H(i,j)=-t2
				ELSE
					H(i,j)=0.D0
				END IF
			END IF
		END DO
		END DO
		END DO
	END DO
	END DO
	END DO
	
	H(65:128,65:128)=H(1:64,1:64)
	
	open(unit=11,file='DoS_vs_w.dat',status='unknown')
    open(unit=12,file='mu.dat',status='unknown')
		
	!	Hitung [G(w)]
	dw=(wmax-wmin)/(Nw-1) !Rentang antar titik
	norm=0.D0
	do iw=1,Nw
	w(iw)=wmin+dw*(iw-1)
	
	Gdummy(:,:)=(w(iw)+ii*0.1D0)*IdMat(:,:)-H(:,:)-sigma(:,:)
	
	!	   Matrix Inversion
	call zgetrf(128,128,Gdummy,128,ipiv,info)
	call zgetri(128,Gdummy,128,ipiv,work,2*128,info)
	
	G(:,:,iw)=Gdummy(:,:)
	
	trace=0.D0
	do j=1,128
	trace=trace+aimag(G(j,j,iw))
	PDOS(j,iw)=-(1.D0/pi)*aimag(G(j,j,iw))
	
	end do
	
	DOS(iw)=-(1.D0/pi)*trace

!	Calculate chemical potential
	
	norm=norm+DOS(iw)*wfreq(iw)
	
	
	write(11,*)w(iw),DoS(iw)
	end do
	
	mu=rtbis(integrand,wmin,wmax,tolerance)
	print*, mu
    write(12,*)mu,0.D0
    write(12,*)mu,300.D0
	close(11)

	deallocate(ind)
	call deallocation

	
	end program
	

	FUNCTION integrand(x)
	
	USE global
	IMPLICIT NONE
	REAL*8::integrand,x,tot
	REAL*8,EXTERNAL::fermi

	tot=0.D0
	DO i=1,Nw
		tot=tot+DOS(i)*fermi(x,T,w(i))*wfreq(i) ! Pemberat integrasi
	END DO
	integrand=nfilling-(tot*128.D0/norm) ! Normalisasi, dari 128 ke 127,...
	
	RETURN

	
	END FUNCTION integrand
