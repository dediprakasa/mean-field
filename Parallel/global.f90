module global

! Mathematical and physical constant
	complex*16, parameter :: im = cmplx(0.d0, 1.d0)
	real*8, 	parameter :: pi = 2.d0*asin(1.d0)
	
! Variables
	real*8 	:: eps, t1, t2, t3, T, wmin, wmax, dw, nfilling, U, norm, mu
	real*8 	:: error_sigma_init, error_sigma, alpha, tol_mu, tol_sigma
	integer :: Nw, N1D, Nsites, nproc
	
! Matrices
	real*8, 	dimension(:),	  allocatable :: wfreq, n_up, n_down, DOS_up, DOS_down, DOS_total, n_number
	real*8, 	dimension(:,:),	  allocatable :: id_mat, H, w, PDOS, del_sigma
	real*8,		dimension(:,:,:), allocatable :: sigma, sigma_calc
	complex*16, dimension(:,:),	  allocatable :: Gdummy
	complex*16, dimension(:,:,:), allocatable :: G
	integer, 	dimension(:,:,:), allocatable :: ind
	
! Arrays and variables for LAPACK
	complex*16, dimension(:),	  allocatable :: work ! dimension(2*2*Nsites)
	integer, 	dimension(:),	  allocatable :: ipiv ! dimension(2*Nsites)
	
end module global
