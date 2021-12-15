SUBROUTINE POLICY_FUNCTION_ITERATION(in_tauk,in_tauh,in_BY,in_K,in_H,in_C,in_r,in_w,in_f,in_a_pol,in_c_pol,in_h_pol,out_aa_const,out_bb_const)

    USE MOD_PARAMETERS
    USE MPI
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: in_tauk, in_tauh, in_BY, in_K, in_H, in_C, in_r, in_w
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: in_f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: in_a_pol, in_c_pol, in_h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), DIMENSION(m,n,p), INTENT(OUT) :: out_aa_const, out_bb_const
	
	! ==================================================================== !
    ! LOCAL
	! ==================================================================== !

    REAL(8) :: tauk, tauh, BY, K, H, C, r, w, part_c_ap0, part_n_ap0, part_cp_ap0, part_np_ap0
	
	REAL(8), DIMENSION(1) :: aa_const_p, bb_const_p, cp, hp
	
	REAL(8), DIMENSION(m,n,p,n,p) :: part_cp_ap, part_np_ap, uccp, ucnp
	
    REAL(8), DIMENSION(m,n) :: f
	
	REAL(8), DIMENSION(m,n,p) :: a_pol, c_pol, h_pol, ucc, ucn, aa_const, bb_const, aa_const_0, bb_const_0, EE1, EE2, EE3, EE4, part_c_ap, part_n_ap, &
	                             der_a_pol
	
    ! Outer loop
    REAL(8) ::  d, tol
    INTEGER ::  ite, maxite, i, j, jj, q, qq  

    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
    
	tauk = in_tauk    
    tauh = in_tauh
    BY   = in_BY
	K    = in_K
	H    = in_H
	C    = in_C
	r    = in_r
	w    = in_w
	f    = in_f
	a_pol = in_a_pol
	c_pol = in_c_pol
	h_pol = in_h_pol
	
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p				
				CALL PARTIAL_DERIVATIVE_C_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_ap(i,j,q))
				CALL PARTIAL_DERIVATIVE_N_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_ap(i,j,q))
				ucc(i,j,q) = -gam*(gam*sigma+1.0D0-gam)*((((c_pol(i,j,q)**gam)*((1.0D0-h_pol(i,j,q))**(1.0D0-gam)))**(1.0D0-sigma))/(c_pol(i,j,q)**2.0D0))
				ucn(i,j,q) = (sigma-1.0D0)*gam*(1.0D0-gam)*((((c_pol(i,j,q)**gam)*((1.0D0-h_pol(i,j,q))**(1.0D0-gam)))**(1.0D0-sigma))/(c_pol(i,j,q)*(1.0D0-h_pol(i,j,q))))
			ENDDO
		ENDDO
	ENDDO
	
	EE3 = 0.0D0
	EE4 = 0.0D0
	
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				DO jj = 1,n
					DO qq = 1,p
						CALL interp_linear(1,m,ap,c_pol(:,jj,qq),1,a_pol(i,j,q),cp(1))
						!CALL interp_linear(1,m,ap,h_pol(:,jj,qq),1,a_pol(i,j,q),hp(1))
						if (hp_thr(jj,qq) .gt. 1) then
						CALL interp_linear(1,hp_thr(jj,qq),ap(1:hp_thr(jj,qq)),h_pol(1:hp_thr(jj,qq),jj,qq),1,a_pol(i,j,q),hp(1))
						else
							hp(1)  = 0.0D0
						end if							
						hp(1) = MAX(hp(1),0.0D0)						
						CALL PARTIAL_DERIVATIVE_CP_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,jj,qq,part_cp_ap(i,j,q,jj,qq))
						CALL PARTIAL_DERIVATIVE_NP_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,jj,qq,part_np_ap(i,j,q,jj,qq))
						uccp(i,j,q,jj,qq) = -gam*(gam*sigma+1.0D0-gam)*((((cp(1)**gam)*((1.0D0-hp(1))**(1.0D0-gam)))**(1.0D0-sigma))/(cp(1)**2.0D0))
						ucnp(i,j,q,jj,qq) = (sigma-1.0D0)*gam*(1.0D0-gam)*((((cp(1)**gam)*((1.0D0-hp(1))**(1.0D0-gam)))**(1.0D0-sigma))/(cp(1)*(1.0D0-hp(1))))
						EE3(i,j,q) = EE3(i,j,q) + (s(jj)+e(qq)*mult_e(jj))*part_np_ap(i,j,q,jj,qq)*Ps(j,jj)*Pe(qq)
						EE4(i,j,q) = EE4(i,j,q) + part_cp_ap(i,j,q,jj,qq)*Ps(j,jj)*Pe(qq)
					ENDDO  
				ENDDO				
			ENDDO
		ENDDO
	ENDDO
		
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				IF (a_pol(i,j,q) .GT. amin) THEN
					CALL PARTIAL_DERIVATIVE_C_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_ap0)
					CALL PARTIAL_DERIVATIVE_N_AP_0(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_ap0)
				    CALL PARTIAL_DERIVATIVE_CP_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,j,q,part_cp_ap0)
					CALL PARTIAL_DERIVATIVE_NP_AP_0(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,j,q,part_np_ap0)                    
					aa_const_0(i,j,q) = (beta*tauk*r + (tauh*w)*(part_n_ap0+beta*part_np_ap0) + tauc0*(part_c_ap0+beta*part_cp_ap0))/&
					                    (ucn(i,j,q)*(part_n_ap0+beta*part_np_ap0) + ucc(i,j,q)*(part_c_ap0+beta*part_cp_ap0))
					!bb_const_0(i,j,q) = 1.0D0+(1.0D0-tauk)*r
				ELSE
					aa_const_0(i,j,q) = 0.0D0
					bb_const_0(i,j,q) = 0.0D0			
				ENDIF
			ENDDO
		ENDDO
	ENDDO	


	DO i = 2,m-1
		DO j = 1,n
			DO q = 1,p
				bb_const_0(i,j,q) = 0.5D0*(a_pol(i+1,j,q)-a_pol(i,j,q))/(ap(i+1)-ap(i)) + 0.5D0*(a_pol(i,j,q)-a_pol(i-1,j,q))/(ap(i)-ap(i-1))
			ENDDO
		ENDDO		
	ENDDO
	bb_const_0(1,:,:) = bb_const_0(2,:,:)
	bb_const_0(m,:,:) = bb_const_0(m-1,:,:)	


	!bb_const = bb_const_0
	
	
    IF (MORE_PRECISE .EQ. 0) THEN
		tol = 1e-6
	ELSEIF (MORE_PRECISE .EQ. 1) THEN
		!tol = 1e-8
		tol = 1e-12
	ENDIF

    maxite  = 100000
    ite     = 0
	d       = 1.0D0
	

	DO WHILE ((d .gt. tol) .and. (ite .le. maxite))

        ite = ite + 1

        IF (ite .gt. (maxite-1)) THEN
             ER      = 1
             ER_TYPE = 5
        ENDIF
		
		!EE1 = 0.0D0
		EE2 = 0.0D0
		
		DO i = 1,m
			DO j = 1,n
				DO q = 1,p
					DO jj = 1,n
						DO qq = 1,p
							!CALL interp_linear(1,m,ap,aa_const_0(:,jj,qq),1,a_pol(i,j,q),aa_const_p(1))
							CALL interp_linear(1,m,ap,bb_const_0(:,jj,qq),1,a_pol(i,j,q),bb_const_p(1))
							!EE1(i,j,q) = EE1(i,j,q) + beta*aa_const_p(1)*(uccp(i,j,q,jj,qq)*part_cp_ap(i,j,q,jj,qq)+ucnp(i,j,q,jj,qq)*part_np_ap(i,j,q,jj,qq))*Ps(j,jj)*Pe(qq)
							EE2(i,j,q) = EE2(i,j,q) + beta*(1.0D0+(1.0D0-tauk)*r-bb_const_p(1))*(uccp(i,j,q,jj,qq)*part_cp_ap(i,j,q,jj,qq)+ucnp(i,j,q,jj,qq)*part_np_ap(i,j,q,jj,qq))*Ps(j,jj)*Pe(qq)
						ENDDO  
					ENDDO				
				ENDDO
			ENDDO
		ENDDO

        DO i = 1,m
			DO j = 1,n
				DO q = 1,p
					IF (a_pol(i,j,q) .GT. amin) THEN
						!aa_const(i,j,q) = (EE1(i,j,q)-beta*tauk*r-tauh*w*((s(j)+e(q)*mult_e(j))*part_n_ap(i,j,q)+beta*EE3(i,j,q))-tauc0*(part_c_ap(i,j,q)+beta*EE4(i,j,q)))/(-(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q))+EE2(i,j,q))
						bb_const(i,j,q) = (-(1.0D0+(1.0D0-tauk)*r)*(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q)))/(-(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q))+EE2(i,j,q))	    
					ELSE
						!aa_const(i,j,q) = 0.0D0
						bb_const(i,j,q) = 0.0D0	
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		
		!d  = MAX(MAXVAL(ABS(aa_const-aa_const_0)/(1.0D0+ABS(aa_const_0))),MAXVAL(ABS(bb_const-bb_const_0)/(1.0D0+ABS(bb_const_0))))
		d  = MAXVAL(ABS(bb_const-bb_const_0)/(1.0D0+ABS(bb_const_0)))
		
		!PRINT *, d, ite
		
		time_end_bound = MPI_WTIME()	
		IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
			ER = 1
			ER_TYPE = 7
			d  = 0.0D0
		ENDIF
		
		!aa_const_0 = aa_const
		bb_const_0 = bb_const

	ENDDO	


	maxite  = 100000
    ite     = 0
	d       = 1.0D0
	

	DO WHILE ((d .gt. tol) .and. (ite .le. maxite))

        ite = ite + 1

        IF (ite .gt. (maxite-1)) THEN
             ER      = 1
             ER_TYPE = 5
        ENDIF
		
		EE1 = 0.0D0
		EE2 = 0.0D0
		
		DO i = 1,m
			DO j = 1,n
				DO q = 1,p
					DO jj = 1,n
						DO qq = 1,p
							CALL interp_linear(1,m,ap,aa_const_0(:,jj,qq),1,a_pol(i,j,q),aa_const_p(1))
							CALL interp_linear(1,m,ap,bb_const_0(:,jj,qq),1,a_pol(i,j,q),bb_const_p(1))
							EE1(i,j,q) = EE1(i,j,q) + beta*aa_const_p(1)*(uccp(i,j,q,jj,qq)*part_cp_ap(i,j,q,jj,qq)+ucnp(i,j,q,jj,qq)*part_np_ap(i,j,q,jj,qq))*Ps(j,jj)*Pe(qq)
							EE2(i,j,q) = EE2(i,j,q) + beta*(1.0D0+(1.0D0-tauk)*r-bb_const_p(1))*(uccp(i,j,q,jj,qq)*part_cp_ap(i,j,q,jj,qq)+ucnp(i,j,q,jj,qq)*part_np_ap(i,j,q,jj,qq))*Ps(j,jj)*Pe(qq)
						ENDDO  
					ENDDO				
				ENDDO
			ENDDO
		ENDDO

        DO i = 1,m
			DO j = 1,n
				DO q = 1,p
					IF (a_pol(i,j,q) .GT. amin) THEN
						aa_const(i,j,q) = (EE1(i,j,q)-beta*tauk*r-tauh*w*((s(j)+e(q)*mult_e(j))*part_n_ap(i,j,q)+beta*EE3(i,j,q))-tauc0*(part_c_ap(i,j,q)+beta*EE4(i,j,q)))/(-(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q))+EE2(i,j,q))
						!bb_const(i,j,q) = (-(1.0D0+(1.0D0-tauk)*r)*(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q)))/(-(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q))+EE2(i,j,q))	    
					ELSE
						aa_const(i,j,q) = 0.0D0
						!bb_const(i,j,q) = 0.0D0	
					ENDIF
				ENDDO
			ENDDO
		ENDDO
		
		d  = MAXVAL(ABS(aa_const-aa_const_0)/(1.0D0+ABS(aa_const_0)))
		!d  = MAXVAL(ABS(bb_const-bb_const_0)/(1.0D0+ABS(bb_const_0)))
		
		!PRINT *, d, ite
		
		time_end_bound = MPI_WTIME()	
		IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
			ER = 1
			ER_TYPE = 7
			d  = 0.0D0
		ENDIF
		
		aa_const_0 = aa_const
		!bb_const_0 = bb_const

	ENDDO		


    ! DO WHILE ((d .gt. tol) .and. (ite .le. maxite))

    !     ite = ite + 1

    !     IF (ite .gt. (maxite-1)) THEN
    !          ER      = 1
    !          ER_TYPE = 5
    !     ENDIF
		
	! 	EE1 = 0.0D0
	! 	EE2 = 0.0D0
		
	! 	DO i = 1,m
	! 		DO j = 1,n
	! 			DO q = 1,p
	! 				DO jj = 1,n
	! 					DO qq = 1,p
	! 						CALL interp_linear(1,m,ap,aa_const_0(:,jj,qq),1,a_pol(i,j,q),aa_const_p(1))
	! 						CALL interp_linear(1,m,ap,bb_const_0(:,jj,qq),1,a_pol(i,j,q),bb_const_p(1))
	! 						EE1(i,j,q) = EE1(i,j,q) + beta*aa_const_p(1)*(uccp(i,j,q,jj,qq)*part_cp_ap(i,j,q,jj,qq)+ucnp(i,j,q,jj,qq)*part_np_ap(i,j,q,jj,qq))*Ps(j,jj)*Pe(qq)
	! 						EE2(i,j,q) = EE2(i,j,q) + beta*(1.0D0+(1.0D0-tauk)*r-bb_const_p(1))*(uccp(i,j,q,jj,qq)*part_cp_ap(i,j,q,jj,qq)+ucnp(i,j,q,jj,qq)*part_np_ap(i,j,q,jj,qq))*Ps(j,jj)*Pe(qq)
	! 					ENDDO  
	! 				ENDDO				
	! 			ENDDO
	! 		ENDDO
	! 	ENDDO

    !     DO i = 1,m
	! 		DO j = 1,n
	! 			DO q = 1,p
	! 				IF (a_pol(i,j,q) .GT. amin) THEN
	! 					aa_const(i,j,q) = (EE1(i,j,q)-beta*tauk*r-tauh*w*((s(j)+e(q)*mult_e(j))*part_n_ap(i,j,q)+beta*EE3(i,j,q))-tauc0*(part_c_ap(i,j,q)+beta*EE4(i,j,q)))/(-(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q))+EE2(i,j,q))
	! 					bb_const(i,j,q) = (-(1.0D0+(1.0D0-tauk)*r)*(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q)))/(-(ucc(i,j,q)*part_c_ap(i,j,q)+ucn(i,j,q)*part_n_ap(i,j,q))+EE2(i,j,q))	    
	! 				ELSE
	! 					aa_const(i,j,q) = 0.0D0
	! 					bb_const(i,j,q) = 0.0D0	
	! 				ENDIF
	! 			ENDDO
	! 		ENDDO
	! 	ENDDO
		
	! 	d  = MAX(MAXVAL(ABS(aa_const-aa_const_0)/(1.0D0+ABS(aa_const_0))),MAXVAL(ABS(bb_const-bb_const_0)/(1.0D0+ABS(bb_const_0))))
		
	! 	!PRINT *, d, ite
		
	! 	time_end_bound = MPI_WTIME()	
	! 	IF ((time_end_bound - time_start_bound) .GE. TIME_BOUND) THEN
	! 		ER = 1
	! 		ER_TYPE = 7
	! 		d  = 0.0D0
	! 	ENDIF
		
	! 	aa_const_0 = aa_const
	! 	bb_const_0 = bb_const

	! ENDDO
	
	! DO i = 2,m-1
	! 	DO j = 1,n
	! 		DO q = 1,p
    !         der_a_pol(i,j,q) = 0.5D0*(a_pol(i+1,j,q)-a_pol(i,j,q))/(ap(i+1)-ap(i)) + 0.5D0*(a_pol(i,j,q)-a_pol(i-1,j,q))/(ap(i)-ap(i-1))
	! 		WRITE(*,'(A,F20.10,F20.10,F20.10)') 'TEST DERIV', der_a_pol(i,j,q), bb_const_0(i,j,q), ABS(der_a_pol(i,j,q)-bb_const_0(i,j,q))
	! 		ENDDO
	! 	ENDDO		
	! ENDDO	

	! DO j = 1,n
	! DO q = 1,p
	! DO i = 1,m        
	! 	WRITE(*,'(A,I5,I5,I5,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)') 'PolFun', i,j,q, aa_const(i,j,q), bb_const(i,j,q), h_pol(i,j,q), ap(i), a_pol(i,j,q), c_pol(i,j,q) 
	! END DO
	! END DO
	! END DO
	
	out_aa_const = aa_const_0
	out_bb_const = bb_const_0
	
        
END SUBROUTINE POLICY_FUNCTION_ITERATION
    
    
    

