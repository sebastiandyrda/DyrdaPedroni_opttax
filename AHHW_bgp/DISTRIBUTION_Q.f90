SUBROUTINE DISTRIBUTION_Q(app,f,aa_const,bb_const,out_fq)

    USE MOD_PARAMETERS
	USE MPI
    IMPLICIT NONE

    ! ==================================================================== !
    ! Input Variables
    ! ==================================================================== !

    REAL (8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL (8), DIMENSION(m,n,p), INTENT(IN) :: app, aa_const, bb_const
    
    ! ==================================================================== !
    ! Output Variables
    ! ==================================================================== !

    REAL (8), INTENT(OUT), DIMENSION(m,n,n_q) :: out_fq 

    ! ==================================================================== !
    ! Defining local variables
    ! ==================================================================== !

    REAL (8) :: kk0, qq0, &
                kk1, qq1, &
                n0,  nq0, &
				phi

    REAL (8), DIMENSION(:,:,:), ALLOCATABLE :: fq, fq0
	
	REAL (8), DIMENSION(:,:,:,:), ALLOCATABLE :: q_pol

    INTEGER :: z, zq, z_prev, found_z, z_test
	
	INTEGER, DIMENSION(m,n,p)   :: zz	
	INTEGER, DIMENSION(m,n,p,n_q) :: zzq
	
	REAL (8), DIMENSION(2,2) :: pp
    REAL (8), DIMENSION(2)   :: ppp

    ! Loop Variables

    REAL (8) :: d2, &       ! Initial value for the loop's convergence measure.
                tol2        ! Tolerance.
    INTEGER ::  ite , &     ! Iterations
                maxite
    INTEGER ::  i, j, ii, jj, q, qq, i_q  ! Loop indexes

    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	ALLOCATE(fq(m,n,n_q),fq0(m,n,n_q),q_pol(m,n,p,n_q))
	
    DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				DO i_q = 1,n_q
					q_pol(i,j,q,i_q) =  aa_const(i,j,q)+bb_const(i,j,q)*q_grid(i_q)
				ENDDO
			ENDDO
		ENDDO
	ENDDO	
	
	! zz(j,i) is the index of ap that is closest to app(j,i)
	
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
	
	! zzq(j,i,q,i_q) is the index of q_grid that is closest to q_pol(j,i,q,i_q)
	
	DO i = 1,n
		DO j = 1,m
		    DO q = 1,p
				DO i_q = 1,n_q
					qq1 = q_pol(j,i,q,i_q)
					IF (qq1 .le. q_grid(1)) THEN
						zzq(j,i,q,i_q) = 1					
						z_prev         = 1
					ELSEIF (qq1 .ge. q_grid(n_q)) THEN
						zzq(j,i,q,i_q) = m
						z_prev     	   = m
					ELSE                    
						! IF (bb_const(j,i,q) .GE. 0.0D0) THEN					
							! IF (i_q .EQ. 1) THEN
								! z_prev = 1
							! ENDIF					    
							! found_z = 0
							! qq      = z_prev-1
							! DO WHILE (found_z .EQ. 0)
								! qq = qq + 1
								! IF (qq1 .LT. q_grid(qq+1)) THEN
									! z = qq
									! found_z = 1		
								! ENDIF	
							! ENDDO	
							! z_prev = z					
						! ELSE	
							! IF (i_q .EQ. 1) THEN
								! z_prev = p
							! ENDIF					    
							! found_z = 0
							! qq      = z_prev+1
							! DO WHILE (found_z .EQ. 0)
								! qq = qq - 1
								! IF (qq1 .LT. q_grid(MAX(qq-1,1))) THEN								
									! z = qq								
									! found_z = 1		
								! ENDIF		
							! ENDDO					
							! z_prev = z					
						! ENDIF	

						! IF (z .NE. MINLOC(ABS(q_grid-qq1),1)) THEN
							! WRITE(*,'(I5,I5,I5,I5,I5,I5,F20.10,F20.10,F20.10)') j,i,q,i_q,z, MINLOC(ABS(q_grid-qq1),1), qq1, q_grid(z), q_grid(MINLOC(ABS(q_grid-qq1),1))
							! READ(*,*)
						! ENDIF
						
						! IF (q_grid(z) .le. qq1) THEN
							! z = z + 1
						! ENDIF
						! zzq(j,i,q,i_q) = z
						
						z = MINLOC(ABS(q_grid-qq1),1)
						IF (q_grid(z) .le. qq1) THEN
							z = z + 1
						ENDIF
						zzq(j,i,q,i_q) = z
						
					ENDIF
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	
	DO i = 1,n
		DO j = 1,m
		    DO i_q = 1,n_q
				fq0(j,i,i_q) = f(j,i)/DBLE(n_q)
			ENDDO
		ENDDO
	ENDDO
	
    ite    = 0
    d2     = 1.0D0
	
	IF (MORE_PRECISE .EQ. 0) THEN
		tol2 = 1e-6
	ELSEIF (MORE_PRECISE .EQ. 1) THEN
		tol2 = 1e-8
	ENDIF
	
    maxite = 1e5
	
	DO WHILE (d2 .ge. tol2 .AND. ite .lt. maxite)
        
        ite = ite + 1
		
		IF (ite .gt. (maxite-1)) THEN
             ER      = 1
             ER_TYPE = 6
        ENDIF
        
        fq  = 0.0D0

        DO i_q = 1,n_q
            DO i = 1,n
                DO j = 1,m
					DO q = 1,p
						kk0 = ap(j)
						kk1 = app(j,i,q)
						qq0 = q_grid(i_q)
						qq1 = q_pol(j,i,q,i_q)
						IF ((kk1 .eq. amin) .and. (qq1 .le. q_grid(1))) THEN
							DO ii = 1,n
								fq(1,ii,1) = fq(1,ii,1) + fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF ((kk1 .eq. amin) .and. (qq1 .ge. q_grid(n_q))) THEN
							DO ii = 1,n
								fq(1,ii,n_q) = fq(1,ii,n_q) + fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF ((kk1 .ge. amax) .and. (qq1 .le. q_grid(1))) THEN
							DO ii = 1,n
								fq(m,ii,1) = fq(m,ii,1) + fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF ((kk1 .ge. amax) .and. (qq1 .ge. q_grid(n_q))) THEN
							DO ii = 1,n
								fq(m,ii,n_q) = fq(m,ii,n_q) + fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF (qq1 .le. q_grid(1)) THEN
							z = zz(j,i,q)
							ppp(1) = ((ap(z)  -kk1)/(ap(z)-ap(z-1)))
							ppp(2) = ((kk1-ap(z-1))/(ap(z)-ap(z-1)))
							DO ii = 1,n
								fq(z-1,ii,1) = fq(z-1,ii,1) + ppp(1)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(z  ,ii,1) = fq(z  ,ii,1) + ppp(2)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF (qq1 .ge. q_grid(n_q)) THEN
							z  = zz(j,i,q)
							ppp(1) = ((ap(z)  -kk1)/(ap(z)-ap(z-1)))
							ppp(2) = ((kk1-ap(z-1))/(ap(z)-ap(z-1)))
							DO ii = 1,n
								fq(z-1,ii,n_q) = fq(z-1,ii,n_q) + ppp(1)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(z  ,ii,n_q) = fq(z  ,ii,n_q) + ppp(2)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF (kk1 .eq. amin) THEN
							zq = zzq(j,i,q,i_q)
							ppp(1) = ((q_grid(zq)  -qq1)/(q_grid(zq)-q_grid(zq-1)))
							ppp(2) = ((qq1-q_grid(zq-1))/(q_grid(zq)-q_grid(zq-1)))
							DO ii = 1,n
								fq(1,ii,zq-1) = fq(1,ii,zq-1) + ppp(1)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(1,ii,zq  ) = fq(1,ii,zq  ) + ppp(2)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSEIF (kk1 .ge. amax) THEN
							zq = zzq(j,i,q,i_q)
							ppp(1) = ((q_grid(zq)  -qq1)/(q_grid(zq)-q_grid(zq-1)))
							ppp(2) = ((qq1-q_grid(zq-1))/(q_grid(zq)-q_grid(zq-1)))
							DO ii = 1,n
								fq(m,ii,zq-1) = fq(m,ii,zq-1) + ppp(1)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(m,ii,zq  ) = fq(m,ii,zq  ) + ppp(2)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ELSE
							z  = zz(j,i,q)
							zq = zzq(j,i,q,i_q)
							pp(1,1) = ((ap(z)  -kk1)/(ap(z)-ap(z-1)))*((q_grid(zq)  -qq1)/(q_grid(zq)-q_grid(zq-1)))
							pp(1,2) = ((ap(z)  -kk1)/(ap(z)-ap(z-1)))*((qq1-q_grid(zq-1))/(q_grid(zq)-q_grid(zq-1)))
							pp(2,1) = ((kk1-ap(z-1))/(ap(z)-ap(z-1)))*((q_grid(zq)  -qq1)/(q_grid(zq)-q_grid(zq-1)))
							pp(2,2) = ((kk1-ap(z-1))/(ap(z)-ap(z-1)))*((qq1-q_grid(zq-1))/(q_grid(zq)-q_grid(zq-1)))
							DO ii = 1,n
								fq(z-1,ii,zq-1) = fq(z-1,ii,zq-1) + pp(1,1)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(z-1,ii,zq  ) = fq(z-1,ii,zq  ) + pp(1,2)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(z  ,ii,zq-1) = fq(z  ,ii,zq-1) + pp(2,1)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
								fq(z  ,ii,zq  ) = fq(z  ,ii,zq  ) + pp(2,2)*fq0(j,i,i_q)*Ps(i,ii)*Pe(q)
							ENDDO
						ENDIF
					ENDDO
                ENDDO
            ENDDO
        ENDDO
		
        DO i = 1,n
            fq(:,i,:) = SSPs(1,i)*fq(:,i,:)/(SUM(fq(:,i,:)))
        ENDDO
		
        phi = 0.0D0
        d2 = MAXVAL(ABS(fq-fq0)/(1.0D0+ABS(fq0)))
		
        fq0 = phi*fq0 + (1.0D0-phi)*fq	
		
		time_end_bound = MPI_WTIME()
		IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
			ER = 1
			ER_TYPE = 7
			d2  = 0.0D0
		ENDIF
		
		!PRINT *, 'CONV DIST = ', d2

    ENDDO
	
    out_fq = fq
	
	DEALLOCATE(fq,fq0,q_pol)

END SUBROUTINE DISTRIBUTION_Q