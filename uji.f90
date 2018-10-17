program uji

	implicit none

	integer :: numtasks, numextra, nw, nproc, i
	integer, allocatable :: nfinish(:), nstart(:), numdata(:)

	allocate(nfinish(20))
	allocate(nstart(20))
	allocate(numdata(20))

	nw = 20
	nproc = 6

	  numtasks=nw/nproc
      numextra=nw-numtasks*nproc
            
        if (numextra.gt.0) then
        do i=0,nproc-1-numextra
           nstart(i)=i*numtasks+1
           nfinish(i)=nstart(i)+numtasks-1
        end do
        do i=nproc-numextra,nproc-1
           nstart(i)=nfinish(i-1)+1
           nfinish(i)=nstart(i)+numtasks
        end do
        else
        do i=0,nproc-1
           nstart(i)=i*numtasks+1
           nfinish(i)=nstart(i)+numtasks-1
        end do
        end if

        do i=0,nproc-1
        numdata(i)=nfinish(i)-nstart(i)+1
        end do 
			
		do i = 0, nproc-1
		write(*,*)
		write(*,*) 'i = ', i, 'nstart, nfinish, numdata = ', nstart(i), ' ', nfinish(i), ' ', numdata(i)
		end do

	deallocate(nfinish)
	deallocate(nstart)
	deallocate(numdata)
	
end program
