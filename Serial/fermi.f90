	FUNCTION fermi(mu,T,nu)
	
	IMPLICIT NONE
	REAL*8::eng,T,beta,mu,nu,fermi,kT

	kT=8.617332478d-5*T ! kT is in eV, while T is in Kelvin
	
	IF(kT==0.D0) THEN
		IF(nu.LT.mu)fermi=1.D0
		IF(nu==mu)fermi=0.5D0
		IF(nu.GT.mu)fermi=0.D0
	ELSE
		beta=1.D0/kT
		eng=beta*(nu-mu)
		IF(eng.GT.0.D0) THEN
			fermi=EXP(-eng)/(1.D0+EXP(-eng))
		ELSE
			fermi=(1.D0/(1.D0+EXP(eng)))
		END IF
	END IF
	
	RETURN
	
	END FUNCTION fermi
