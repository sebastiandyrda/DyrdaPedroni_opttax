SUBROUTINE DECISION_RULES(incp0,tauk,tauh,tauc,LS,r,w,outC,outA)

USE MOD_PARAMETERS
USE MPI
IMPLICIT NONE

REAL(8), DIMENSION(m,n,p), INTENT(IN)  :: incp0
REAL(8), INTENT(IN)                    :: r,w,tauk,tauh,tauc,LS

REAL(8), DIMENSION(m,n,p), INTENT(OUT) :: outC,outA

REAL(8), DIMENSION(m,n) :: EMU

REAL(8), DIMENSION(m,n,p) :: cc,cc1,aa,app,cp,cp0

REAL(8) :: d1 , &        ! Initial value for the loop's convergence measure.
           tol1,&        ! Tolerance.
           phi, &
           d_prev
INTEGER :: ite, maxite
INTEGER :: i, j, jj, q, qq  ! Loop indexes

LOGICAL r8vec_ascends_strictly


!-----------------------------------------------------------------------------------!
! EXECUTABLE
!-----------------------------------------------------------------------------------!

d_prev  = 0.0D0
d1      = 1.0D0         ! Initial value for the loop's convergence measure.
maxite  = 500000
ite     = 0

IF (MORE_PRECISE .EQ. 0) THEN
	tol1 = 1e-6
ELSEIF (MORE_PRECISE .EQ. 1) THEN
	tol1 = 5e-7
ELSEIF (MORE_PRECISE .EQ. 2) THEN
	tol1 = 8e-9
ENDIF

! ==================================================================== !
! Initial guess for the policy function for consumption next period
! c'(a',s')
! ==================================================================== !

cp0 = incp0

DO WHILE ((d1 .gt. tol1) .AND. (ite .le. maxite))

	ite = ite + 1

	IF (ite .gt. (maxite-1)) THEN
		ER = 1
		ER_TYPE = 2
	ENDIF

! ==================================================================== !
! Expected marginal utility tomorrow: 
!       
!  EMU(a',s) = beta*(1+(1-tauk)*r)*E[u'(c'(a',s'))|s]
!
! ==================================================================== !

    EMU = 0.0D0
	
	DO i = 1,m
        DO j = 1,n
			DO jj = 1,n
				DO qq = 1,p
					EMU(i,j) = EMU(i,j) + beta*(1.0D0+(1.0D0-tauk)*r)*Ps(j,jj)*Pe(qq)*&
							   MIN((((1.0D0-gam))/gam)*(((1.0D0+tauc))/((1.0D0-tauh)*w*(s(jj)+e(qq)*mult_e(jj))))*(cp0(i,jj,qq)**(-(sigma/((1.0D0-gam)*(1.0D0-sigma))))),&
								   cp0(i,jj,qq)**((((1.0D0-sigma)*gam-1.0D0)/((1.0D0-gam)*(1.0D0-sigma)))))**((1.0D0-gam)*(1.0D0-sigma))
				ENDDO                    
			ENDDO
        ENDDO
    ENDDO	
	 

! ==================================================================== !
! Using the Euler equation to solve for cc(a',s) [observe that this is
! not the policy function yet]:
!       
!       u'(cc(a',s)) = EMU(a',s)
!
! ==================================================================== !
	
	DO j = 1,n
        DO q = 1,p
            cc1(:,j,q) = ((((gam/((1.0D0-gam)))*(((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))/((1.0D0+tauc))))**&
			             ((1.0D0-gam)*(1.0D0-sigma)))*EMU(:,j))**(-1.0D0/sigma)
        ENDDO
	ENDDO
	
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				IF (cc1(i,j,q) .LE. (gam/((1.0D0-gam)))*(((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))/((1.0D0+tauc)))) THEN
					cc(i,j,q) = cc1(i,j,q)
				ELSE
					cc(i,j,q) = EMU(i,j)**(1.0D0/((1.0D0-sigma)*gam-1.0D0))
				ENDIF
			ENDDO
		ENDDO
	ENDDO	

 ! ==================================================================== !
 ! Using the budget constraint, we solve for a(a',s) [the inverse of 
 ! the policy function for assets]:
 !
 !       cc(a',s) + a' = (1+(1-tauk)*r)a(a',s) + (1-tauh)*w*s + LS
 !
 ! ==================================================================== !
 
    DO i = 1,m
       DO j = 1,n
           DO q = 1,p
			   aa(i,j,q) = MIN((((1.0D0+tauc)*cc(i,j,q)/gam) + ap(i) - (1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)) - LS)/(1.0D0+(1.0D0-tauk)*r),&
							   (((1.0D0+tauc)*cc(i,j,q)+ ap(i)- LS)/(1.0D0+(1.0D0-tauk)*r)))
           ENDDO
       ENDDO
    ENDDO	

     DO j =1,n
         DO i = 1,m
             DO q = 1,p
                 IF (ISNAN(aa(i,j,q))) THEN
                     ER = 1
					 ER_TYPE = 11
                     d1 = 0.0D0
                 ENDIF
             ENDDO
         ENDDO
     ENDDO

! ==================================================================== !
! Finding a''(a',s') [ignoring the borrowing constraint]
! ==================================================================== !
  
    DO i = 1,n
        DO q = 1,p
            IF ( .not. r8vec_ascends_strictly (m,aa(:,i,q)) ) THEN
                CALL sub_myinterp1(aa(:,i,q),ap,ap,m,app(:,i,q))
            ELSE
                CALL interp_linear(1,m,aa(:,i,q),ap,m,ap,app(:,i,q))
            ENDIF
        ENDDO
    ENDDO	
	
! ==================================================================== !
! Taking care of the borrowing constraint for a''
! ==================================================================== !

    app = MAX(app,amin)
	
! ==================================================================== !
! Solving for the policy function c'(a',s') 
! ==================================================================== !

    DO i = 1,n
        DO q = 1,p            
			DO j = 1,m
				IF (app(j,i,q) .GE. (1.0D0+(1.0D0-tauk)*r)*ap(j)+LS-(gam/((1.0D0-gam)))*(1.0D0-tauh)*w*(s(i)+e(q)*mult_e(i))) THEN
					cp(j,i,q) = gam*((1.0D0+(1.0D0-tauk)*r)*ap(j) + (1.0D0-tauh)*w*(s(i)+e(q)*mult_e(i)) + LS - app(j,i,q))/(1.0D0+tauc)
				ELSE
					cp(j,i,q) = ((1.0D0+(1.0D0-tauk)*r)*ap(j) + LS - app(j,i,q))/(1.0D0+tauc)
				ENDIF
			ENDDO
        ENDDO
    ENDDO
    
! ==================================================================== !
! Updating convergence measure and policy function for consumption
! ==================================================================== !    
    
    d1  = MAXVAL(ABS(cp0-cp)/(1.0D0+ABS(cp)))
            
    IF (d1 .lt. d_prev) THEN
       phi = 0.9D0 - 0.5D0*MAX(MIN((DBLE(ite-5)/20.0D0),1.0D0),0.0D0)
    ELSE
       phi = 0.9D0
	ENDIF
	    
    d_prev = d1
    
    cp0 = phi*cp0+(1.0D0-phi)*cp
	
	time_end_bound = MPI_WTIME()

    IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
		ER = 1
		ER_TYPE = 7
		d1  = 0.0D0
	ENDIF
                       
ENDDO  ! end of while for policy function for consumption

outA = app
outC = cp

END SUBROUTINE DECISION_RULES

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

SUBROUTINE DECISION_RULES0(tauk,tauh,tauc,LS,r,w,outC)


    USE MOD_PARAMETERS
    USE MPI
    IMPLICIT NONE

    REAL(8), INTENT(IN)                    :: r,w,tauk,tauh,tauc,LS

    REAL(8), DIMENSION(m,n,p), INTENT(OUT) :: outC

    REAL(8), DIMENSION(m,n) :: EMU

    REAL(8), DIMENSION(m,n) :: cc,cc1,aa,app,cp,cp0
	
	REAL(8) :: M_e

    REAL(8) :: d1 , &        ! Initial value for the loop's convergence measure.
               tol1,&        ! Tolerance.
               phi
    INTEGER :: ite, maxite
    INTEGER :: i, j, jj, q  ! Loop indexes

    LOGICAL r8vec_ascends_strictly

    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!

    d1      = 1.0D0         ! Initial value for the loop's convergence measure.
    maxite  = 200000
    ite     = 0
    tol1    = 1e-4
    phi     = 0.9D0
	
	M_e = 0.0D0
	DO q = 1,p
		M_e = Pe(q)*e(q)
	ENDDO	
	
    ! ==================================================================== !
    ! Initial guess for the policy function for consumption next period
    ! c'(a',s')
    ! ==================================================================== !

    DO i = 1,m
        DO j = 1,n
			cp0(i,j) = gam*(((1.0D0-tauk)*r)*ap(i) + (1.0D0-tauh)*w*(s(j)+M_e) + LS)/(1.0D0+tauc)
        ENDDO
    ENDDO
		
   DO WHILE ((d1 .gt. tol1) .AND. (ite .le. maxite))
 
   ite = ite + 1
 
   IF (ite .gt. (maxite-1)) THEN
       ER = 1
       ER_TYPE = 2
   ENDIF
 
   ! ==================================================================== !
   ! Expected marginal utility tomorrow:
   !
   !  EMU(a',s) = beta*(1+(1-tauk)*r)*E[u'(c'(a',s'))|s]
   !
   ! ==================================================================== !
 
   EMU = 0.0D0
	
	DO i = 1,m
       DO j = 1,n
		   DO jj = 1,n
				EMU(i,j) = EMU(i,j) + beta*(1.0D0+(1.0D0-tauk)*r)*Ps(j,jj)*&
						   MIN((((1.0D0-gam))/gam)*(((1.0D0+tauc))/((1.0D0-tauh)*w*(s(jj)+M_e)))*(cp0(i,jj)**(-(sigma/((1.0D0-gam)*(1.0D0-sigma))))),&
							   cp0(i,jj)**((((1.0D0-sigma)*gam-1.0D0)/((1.0D0-gam)*(1.0D0-sigma)))))**((1.0D0-gam)*(1.0D0-sigma))                    
		   ENDDO           
       ENDDO
   ENDDO	
 
   ! ==================================================================== !
   ! Using the Euler equation to solve for cc(a',s) [observe that this is
   ! not the policy function yet]:
   !
   !       u'(cc(a',s)) = EMU(a',s)
   !
   ! ==================================================================== !
   
	DO j = 1,n
		cc1(:,j) = ((((gam/((1.0D0-gam)))*(((1.0D0-tauh)*w*(s(j)+M_e))/((1.0D0+tauc))))**&
					 ((1.0D0-gam)*(1.0D0-sigma)))*EMU(:,j))**(-1.0D0/sigma)
	ENDDO
	
	DO i = 1,m
		DO j = 1,n
			IF (cc1(i,j) .LE. (gam/((1.0D0-gam)))*(((1.0D0-tauh)*w*(s(j)+M_e))/((1.0D0+tauc)))) THEN
				cc(i,j) = cc1(i,j)
			ELSE
				cc(i,j) = EMU(i,j)**(1.0D0/((1.0D0-sigma)*gam-1.0D0))
			ENDIF
		ENDDO
	ENDDO	
 
    ! ==================================================================== !
    ! Using the budget constraint, we solve for a(a',s) [the inverse of
    ! the policy function for assets]:
    !
    !       cc(a',s) + a' = (1+(1-tauk)*r)a(a',s) + (1-tauh)*w*s + LS
    !
    ! ==================================================================== !
 
	DO i = 1,m
       DO j = 1,n
			aa(i,j) = MIN((((1.0D0+tauc)*cc(i,j)/gam) + ap(i) - (1.0D0-tauh)*w*(s(j)+M_e) - LS)/(1.0D0+(1.0D0-tauk)*r),&
						   (((1.0D0+tauc)*cc(i,j)+ ap(i)- LS)/(1.0D0+(1.0D0-tauk)*r)))
       ENDDO
   ENDDO	
 
   DO j =1,n
       DO i = 1,m
           IF (ISNAN(aa(i,j))) THEN
               ER = 1
				ER_TYPE = 12
               d1 = 0.0D0
           ENDIF
       ENDDO
   ENDDO
 
   ! ==================================================================== !
   ! Finding a''(a',s') [ignoring the borrowing constraint]
   ! ==================================================================== !
 
    DO i = 1,n
       IF ( .not. r8vec_ascends_strictly (m,aa(:,i)) ) THEN
           CALL sub_myinterp1(aa(:,i),ap,ap,m,app(:,i))
       ELSE
           CALL interp_linear(1,m,aa(:,i),ap,m,ap,app(:,i))
       ENDIF
    ENDDO
 
   
    ! ==================================================================== !
    ! Taking care of the borrowing constraint for a''
    ! ==================================================================== !
 
    app = MAX(app,amin)
	
	! ==================================================================== !
	! Solving for the policy function c'(a',s') 
	! ==================================================================== !
		
	DO i = 1,n         
		DO j = 1,m
			IF (app(j,i) .GE. (1.0D0+(1.0D0-tauk)*r)*ap(j)+LS-(gam/((1.0D0-gam)))*(1.0D0-tauh)*w*(s(i)+M_e)) THEN
				cp(j,i) = gam*((1.0D0+(1.0D0-tauk)*r)*ap(j) + (1.0D0-tauh)*w*(s(i)+M_e) + LS - app(j,i))/(1.0D0+tauc)
			ELSE
				cp(j,i) = ((1.0D0+(1.0D0-tauk)*r)*ap(j) + LS - app(j,i))/(1.0D0+tauc)
			ENDIF
		ENDDO
    ENDDO
 
    ! ==================================================================== !
    ! Updating convergence measure and policy function for consumption
    ! ==================================================================== !
 
	d1  = MAXVAL(ABS(cp0-cp)/(1.0D0+ABS(cp)))
	cp0 = phi*cp0+(1.0D0-phi)*cp
	
	
	time_end_bound = MPI_WTIME()

    IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
		ER = 1
		ER_TYPE = 7
		d1  = 0.0D0
	ENDIF
 
	ENDDO  ! end of while for policy function for consumption

	DO j =1,n
		DO i = 1,m
			DO q = 1,p
				outC(i,j,q) = MAX(cp0(i,j) + 0.1D0*e(q),1E-2)
			ENDDO
		ENDDO
	ENDDO

END SUBROUTINE DECISION_RULES0