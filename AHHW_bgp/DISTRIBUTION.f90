	
SUBROUTINE DISTRIBUTION(app,inf0,outf)

    USE MOD_PARAMETERS
    USE MPI
    IMPLICIT NONE

    ! ==================================================================== !
    ! Input Variables
    ! ==================================================================== !

    REAL (8), DIMENSION(m,n,p), INTENT(IN) :: app
    REAL (8), DIMENSION(m,n), INTENT(IN)   :: inf0
    
    ! ==================================================================== !
    ! Output Variables
    ! ==================================================================== !

    REAL (8), INTENT(OUT), DIMENSION(m,n) :: outf 

    ! ==================================================================== !
    ! Defining local variables
    ! ==================================================================== !

    REAL (8) :: kk0, &
                kk1, &
                n0

    REAL (8), DIMENSION(m,n) :: f, &
                                f0

    INTEGER :: z, z_prev, found_z
	
	INTEGER, DIMENSION(m,n,p) :: zz

    ! Loop Variables

    REAL (8) :: d2, &       ! Initial value for the loop's convergence measure.
                tol2        ! Tolerance.
    INTEGER ::  ite , &     ! Iterations
                maxite
    INTEGER ::  i, j, q, ii, jj ! Loop indexes

    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!

    DO i = 1,n
		DO q = 1,p
			DO j = 1,m
				kk1 = app(j,i,q)
				IF (kk1 .eq. amin) THEN
					zz(j,i,q) = 1					
					z_prev    = 1
				ELSEIF (kk1 .ge. amax) THEN
					zz(j,i,q) = m
					z_prev    = m
				ELSE				
					IF (j .EQ. 1) THEN
						z_prev = 1
					ENDIF					    
					found_z = 0
					jj      = z_prev-1
					DO WHILE (found_z .EQ. 0)
						jj = jj + 1
						IF (kk1 .LE. ap(MIN(jj+1,m))) THEN
							IF ((kk1-ap(jj)) .LE. (ap(MIN(jj+1,m))-kk1) ) THEN
								z = jj
							ELSE
								z = jj+1
							ENDIF
							found_z = 1		
						ENDIF	
					ENDDO					
					z_prev = z					
					IF (ap(z) .le. kk1) THEN
						z = z + 1
					ENDIF
					zz(j,i,q) = z					
				ENDIF
			ENDDO			
		ENDDO
	ENDDO
	
	f0 = inf0

    ite    = 0
    d2     = 1.0D0
	
	IF (MORE_PRECISE .EQ. 0) THEN
		tol2 = 1e-6
	ELSEIF (MORE_PRECISE .EQ. 1) THEN
		tol2 = 1e-8
	ELSEIF (MORE_PRECISE .EQ. 2) THEN
		tol2 = 5e-9
	ENDIF

    maxite = 1e5
    
    DO WHILE (d2 .ge. tol2 .AND. ite .lt. maxite)
        ite = ite + 1
        IF (ite .gt. (maxite-1)) THEN
             ER      = 1
             ER_TYPE = 1
        ENDIF
        f  = 0.0D0
        DO i = 1,n
            DO q = 1,p
                DO j = 1,m
                    kk0 = ap(j)
                    kk1 = app(j,i,q)
                    IF (kk1 .eq. amin) THEN
                        DO ii = 1,n
                            f(1,ii) = f(1,ii) + f0(j,i)*Ps(i,ii)*Pe(q)
                        ENDDO
                    ELSEIF (kk1 .ge. amax) THEN
                        DO ii = 1,n
                            f(m,ii) = f(m,ii) + f0(j,i)*Ps(i,ii)*Pe(q)
                        ENDDO
                    ELSE
						z  = zz(j,i,q)
                        n0 = (kk1 - ap(z-1))/(ap(z)-ap(z-1))
                        DO ii = 1,n
                            f(z,ii) = f(z,ii) + n0*f0(j,i)*Ps(i,ii)*Pe(q)
                            f(z-1,ii) = f(z-1,ii) + (1.0D0-n0)*f0(j,i)*Ps(i,ii)*Pe(q)
                        ENDDO
                    ENDIF
                ENDDO
				
            ENDDO
        ENDDO
        DO i = 1,n
            f(:,i) = SSPs(1,i)*f(:,i)/(SUM(f(:,i)))
        ENDDO

        d2 = MAXVAL(ABS(f-f0)/(1.0D0+ABS(f0)))
        f0 = f
        
		time_end_bound = MPI_WTIME()

		IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
			ER = 1
			ER_TYPE = 7
			d2  = 0.0D0
		ENDIF
		
    ENDDO	
	
	outf = f

END SUBROUTINE DISTRIBUTION	

! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################
! ###############################################################################################

SUBROUTINE DISTRIBUTION0(app,inf0,outf)

    USE MOD_PARAMETERS
    USE MPI
    IMPLICIT NONE

    ! ==================================================================== !
    ! Input Variables
    ! ==================================================================== !

    REAL (8), DIMENSION(m,n,p), INTENT(IN) :: app
    REAL (8), DIMENSION(m,n), INTENT(IN)   :: inf0

    ! ==================================================================== !
    ! Output Variables
    ! ==================================================================== !

    REAL (8), INTENT(OUT), DIMENSION(m,n) :: outf

    ! ==================================================================== !
    ! Defining local variables
    ! ==================================================================== !

    REAL (8) :: kk0, &
                kk1, &
                n0

    REAL (8), DIMENSION(m,n) :: f, &
                                f0

    INTEGER :: z, z_prev, found_z

    ! Loop Variables

    REAL (8) :: d2, &       ! Initial value for the loop's convergence measure.
                tol2        ! Tolerance.
    INTEGER ::  ite , &     ! Iterations
                maxite
    INTEGER ::  i, j, ii, jj    ! Loop indexes

    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!

    f0 = inf0

    ite    = 0
    d2     = 1.0D0
    tol2   = 1e-5
    maxite = 1e5

    DO WHILE (d2 .ge. tol2 .AND. ite .lt. maxite)
        ite = ite + 1
        IF (ite .gt. (maxite-1)) THEN
             ER      = 1
             ER_TYPE = 1
        ENDIF
        f  = 0.0D0
        DO i = 1,n
            DO j = 1,m
                kk0 = ap(j)
				kk1 = app(j,i,median_e_ind)
                IF (kk1 .le. amin) THEN
                    DO ii = 1,n
                        f(1,ii) = f(1,ii) + f0(j,i)*Ps(i,ii)
                    ENDDO
					z_prev = 1
                ELSEIF (kk1 .ge. amax) THEN
                    DO ii = 1,n
                        f(m,ii) = f(m,ii) + f0(j,i)*Ps(i,ii)
                    ENDDO
					z_prev = m
                ELSE
                    IF (j .EQ. 1) THEN
						z_prev = 1
					ENDIF					    
					found_z = 0
					jj      = z_prev-1
					DO WHILE (found_z .EQ. 0)
						jj = jj + 1
						IF (kk1 .LE. ap(MIN(jj+1,m))) THEN
							IF ((kk1-ap(jj)) .LE. (ap(MIN(jj+1,m))-kk1) ) THEN
								z = jj
							ELSE
								z = jj+1
							ENDIF
							found_z = 1		
						ENDIF	
					ENDDO
					z_prev = z
                    IF (ap(z) .le. kk1) THEN
                        z = z + 1
                    ENDIF
                    n0 = (kk1 - ap(z-1))/(ap(z)-ap(z-1))
                    DO ii = 1,n
                        f(z,ii) = f(z,ii) + n0*f0(j,i)*Ps(i,ii)
                        f(z-1,ii) = f(z-1,ii) + (1.0D0-n0)*f0(j,i)*Ps(i,ii)
                    ENDDO
                ENDIF

            ENDDO
        ENDDO
        DO i = 1,n
            f(:,i) = SSPs(1,i)*f(:,i)/(SUM(f(:,i)))
        ENDDO

        d2 = MAXVAL(ABS(f-f0)/(1.0D0+ABS(f0)))
        f0 = f

		time_end_bound = MPI_WTIME()

		IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
			ER = 1
			ER_TYPE = 7
			d2  = 0.0D0
		ENDIF

    ENDDO

    outf = f

END SUBROUTINE DISTRIBUTION0