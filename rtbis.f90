	FUNCTION rtbis(func,x1,x2,xacc)
	IMPLICIT NONE
	REAL*8,INTENT(IN)::x1,x2,xacc
	REAL*8::rtbis
	INTERFACE
		FUNCTION func(x)
		IMPLICIT NONE
		REAL*8,INTENT(IN)::x
		REAL*8::func
		END FUNCTION func
	END INTERFACE
	INTEGER,PARAMETER::MAXIT=5000
	INTEGER::j
	REAL*8::dx,f,fmid,xmid
	fmid=func(x2)
	f=func(x1)
	IF(f*fmid>=0.D0) CALL nerror('rtbis: root must be bracketed')
	IF(f<0.D0) THEN
		rtbis=x1
		dx=x2-x1
	ELSE
		rtbis=x2
		dx=x1-x2
	END IF
	DO j=1,MAXIT
		dx=dx*0.5D0
		xmid=rtbis+dx
		fmid=func(xmid)
		IF(fmid<=0.D0) rtbis=xmid
		IF((ABS(dx)<xacc).OR.(fmid==0.D0)) RETURN
	END DO
	CALL nerror('rtbis: too many bisections')
	END FUNCTION rtbis
	
	
	SUBROUTINE nerror(string)
	CHARACTER(LEN=*),INTENT(IN)::string
	WRITE(*,*)'nerror: ',string
	STOP 'program terminated by nerror'
	END SUBROUTINE nerror