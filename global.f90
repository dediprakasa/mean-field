	MODULE global
	
!	Mathematical and physical constants
	COMPLEX*16,PARAMETER::ii=CMPLX(0.D0,1.D0)
	REAL*8,PARAMETER::PI=2.D0*ASIN(1.D0)
	REAL*8,PARAMETER::tolerance=1.D-5, T=0.D0!, x=0.D0,

	
!	Variables
	REAL*8::eps,t1,t2,wmin,wmax,dw,trace,dws,nfilling,norm,mu
	INTEGER::i,j,l,m,n,lp,mp,np,counter,iw,Nw
		
	
!	Matrices
	COMPLEX*16,ALLOCATABLE::G(:,:,:), Gdummy(:,:)
	REAL*8,ALLOCATABLE::IdMat(:,:),H(:,:),w(:),DoS(:),PDOS(:,:), sigma(:,:), wfreq(:)
!	integer, allocatable::ind(:,:,:)

		
	END MODULE global

