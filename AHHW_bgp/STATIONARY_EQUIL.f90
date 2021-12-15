	SUBROUTINE STATIONARY_EQUIL(in_tauk,in_tauh,out_BY,out_K,out_H,out_C,in_r,out_w,out_f,out_a_pol,out_c_pol,out_h_pol)

    USE MOD_PARAMETERS
    USE MPI
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: in_tauk, in_tauh, in_r
    
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: out_BY, out_K, out_H, out_C, out_w
	
	REAL(8), DIMENSION(m,n), INTENT(OUT) :: out_f
	
	REAL(8), DIMENSION(m,n,p), INTENT(OUT) :: out_a_pol, out_c_pol, out_h_pol
	
	! ==================================================================== !
    ! LOCAL
	! ==================================================================== !

    REAL(8) :: H0, &
               w, &     
               Y, &     
               amean, &
               K, &
               C, &
		       H, &
			   BY0, &
			   BY, &
			   tauk, &
			   tauh, &
			   tauc, &
			   r, &
			   LS, &
			   MGR, &
			   K_H
   
    REAL(8), DIMENSION(m,n)   :: f0, &
                                 f

    REAL(8), DIMENSION(m,n,p) :: cp0, &   ! Initial guess for the policy function for consumption next period c'(a',s')
                                 app, &
                                 cp, &
                                 hp, &
                                 hp0

    ! Outer loop
    REAL(8) ::  d, d_prev, &  ! Initial value for the loop's convergence measure.
                tol     , &                ! Tolerance.
                phi
    INTEGER ::  ite, maxite, iii ! Iterations
    INTEGER ::  i, j, q          ! Loop indexes
	
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
    
	ER = 0
	
	tauk = in_tauk    
    tauh = in_tauh
    r    = in_r
    
	
    tauc = tauc0
	
    IF (MORE_PRECISE .EQ. 0) THEN
		tol = 1e-6
	ELSEIF (MORE_PRECISE .EQ. 1) THEN
        tol = 1e-8
        !tol = 5e-6
	ENDIF
    maxite = 100000
	
	! ==================================================================== !
    ! Implementing the Modified Golden Rule
    ! ==================================================================== !
    	
    !MGR =  (1.0D0-beta)/beta
    
    !MGR =  (1.0D0-beta)/beta + DBLE(ITE_MGR-1)*MGR_ADJ
	
    !r = MGR
    !r = 0.04828830051711787D0

    ! DP interest rate:
    !r  = 0.0481723988D0
		
	K_H = ((r+delta)/alpha)**(1.0D0/(alpha-1.0D0))

    w = (1.0D0-alpha)*(K_H**(alpha))

    ! BY  = (amean - K)
	
	!!!!!!!!!!########################
	! BY0 = 0.615D0  !-2.473368907D0*tauh**2 -3.201018332D0*tauk**2 -9.532265906D0*tauh -14.31612215D0*tauk +16.64734741D0
    ! H0  = 0.445D0  !-0.079388487D0*tauh**2 +0.021183669D0*tauk**2 -0.116451752D0*tauh -0.063180378D0*tauk +0.402074125D0	 
    ! BY0 = 2.5223866730D0
    ! H0  = 0.3999420622D0
    !!!!!!!!!!########################
    
    IF (ite_bob .LE. 1) THEN

        BY0 = 1.6363813446D0	 
        H0  = 0.40193D0
    ELSE
		BY0 = BY0_PREV 
		H0  = H0_PREV 
    ENDIF    
	
    ! ==================================================================== !
    ! ==================================================================== !
    ! Start of outer loop to find K
    ! ==================================================================== !
    ! ==================================================================== !

    ite = 0
    iii = 0
    d = 1.0D0
    d_prev = 0.0D0

    DO WHILE ((d .gt. tol) .and. (ite .le. maxite))

        ite = ite + 1

        IF (ite .gt. (maxite-1)) THEN
             ER      = 1
             ER_TYPE = 3
        ENDIF        
        
        ! ==================================================================== !
        ! Calculating Government Consumption (G) and Output (Y)
        ! ==================================================================== !
        
		K = K_H*H0	
		
        Y = (K**alpha)*(H0**(1.0D0-alpha))
		
		amean = K + (BY0*Y)
		
		C = w*H0 + r*(amean - BY0*Y)  - G
		
		LS = tauk*r*amean - r*BY0*Y + tauh*w*H0 + tauc*C - G
		

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Start of loop to find the policy function c(a,s) that solves the
        ! Household's problem
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! ==================================================================== !
        ! Initial guess for the policy function for consumption next period
        ! c'(a',s')
        ! ==================================================================== !

        IF (ite .eq. 1) THEN
            CALL DECISION_RULES0(tauk,tauh,tauc,LS,r,w,cp0)
        ELSE
            cp0 = cp
        ENDIF

        CALL DECISION_RULES(cp0,tauk,tauh,tauc,LS,r,w,cp,app)
		
		! ==================================================================== !
        ! Labor hour levels
        ! ==================================================================== !
		
        DO i = 1,m
			DO j = 1,n
				DO q = 1,p
					hp(i,j,q) = MAX(1.0D0-((1.0D0-gam)/gam)*(((1.0D0+tauc)*cp(i,j,q))/((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))),0.0D0)
				ENDDO
			ENDDO
        ENDDO   
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! End of loop to find the policy function
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Start of loop to find the ergodic distribution of (a,s)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! ==================================================================== !
        ! Initial guess for stationary distribution for (a,s)
        ! ==================================================================== !
         
        IF (ite .eq. 1) THEN            
			DO i = 1,m
				f0(i,:) = SSPs(1,:)/DBLE(m)
			ENDDO
			CALL DISTRIBUTION0(app,f0,f0)
        ELSE
            f0 = f
        ENDIF

        DO j =1,n
            DO i = 1,m
                IF (ISNAN(f0(i,j))) THEN
                     ER      = 1
                     ER_TYPE = 3
                     d       = 0.0D0
                ENDIF
            ENDDO
        ENDDO

        CALL DISTRIBUTION(app,f0,f)
        
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! End of loop to find the ergodic distribution of (a,s)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Calculating supply of assets
        amean = SUM(MATMUL(TRANSPOSE(f),ap))

        ! Market Clearing for Assets
        BY = (amean-K)/Y

        ! Market Clearing for Labor
        H = 0.0D0
        DO i = 1,m
			DO j = 1,n
				DO q = 1,p
					H = H + f(i,j)*Pe(q)*((s(j)+e(q)*mult_e(j))*hp(i,j,q)) 
				ENDDO
			ENDDO
		ENDDO

        d = MAX(ABS(BY-BY0)/(1.0D0+ABS(BY0)),ABS(H-H0)/(1.0D0+ABS(H0)))

        IF (d .lt. d_prev) THEN
            phi = 0.99D0 - 0.09D0*MAX(MIN((REAL(iii-10)/25.0D0),1.0D0),0.0D0)
            iii = iii + 1
        ELSE
            phi = 0.99D0
            iii = 0
        ENDIF
        
		!WRITE(*,'(A,F15.10,I5,I5,F10.5,F10.5,F10.5,F10.5,F10.5,F10.5,F10.5,F10.5,F10.5)') ' INITIAL ',d, ER, ER_TYPE, BY, BY0, K, H, H0, w, r, beta*(1.0D0+r), phi
		
        d_prev = d
        
		BY0 = phi*BY0+(1.0D0-phi)*BY
		H0  = phi*H0 +(1.0D0-phi)*H
		
		time_end_bound = MPI_WTIME()	
		IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
			ER = 1
			ER_TYPE = 7
		ENDIF

        IF (ISNAN(BY0)) THEN
            ER = 1
			ER_TYPE = 3
        ENDIF

        IF (ER .EQ. 1) THEN
            d = 0.0D0
        ENDIF

    ENDDO

    ! ==================================================================== !
    ! ==================================================================== !
    ! End of outer loop to find K
    ! ==================================================================== !
    ! ==================================================================== !
		
	! PRINT *, 'G/Y  = ',G/Y
	! PRINT *, 'LS/Y = ',LS/Y
	! PRINT *, 'B/Y  = ',BY
	! PRINT *, 'K0   = ',K
	! PRINT *, 'H0   = ',H0
	
    !WRITE(*,'(A,F10.5,F10.5,F20.10,F10.5,F10.5,F10.5,F20.10,F20.10,I5,I5)') 'INITIAL ', tauk, tauh, LS/Y, BY, K, H, w, (amean-K), r, ER, ER_TYPE
    WRITE(*,'(A,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)') 'INITIAL ', tauk, tauh, LS/Y, BY, K/Y, H, w , r, beta*(1+r)
    
    !fig_temp(ite_tauk_mgr,:) = [tauk, tauh, LS/Y, BY, K/Y, H, w , r, beta*(1+r)]

    DO j = 1,n
        DO q = 1,p
            hp_thr(j,q)  =  minloc(hp(:,j,q),1,MASK=(hp(:,j,q) .gt. 0.0D0))
            !if (hp_thr(j,q) .eq. 1) THEN
                !print *, j,q,tauk, tauh, ER
                !print *, hp(:,j,q)
            !end if
            ap_thr(j,q)  =  minloc(app(:,j,q),1,MASK=(app(:,j,q) .gt. amin))
        END DO
    END DO        
    
    BY0_PREV = BY
	H0_PREV  = H0    

	
	out_BY = BY
    out_K = K
	out_H = H
	out_C = C
	out_w = w
	out_f = f
	out_a_pol = app
	out_c_pol = cp
	out_h_pol = hp    
        
END SUBROUTINE STATIONARY_EQUIL
