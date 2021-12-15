SUBROUTINE RESIDUAL(in_tauk,in_tauh,in_BY,in_K,in_H,in_C,in_r,in_w,in_f,in_a_pol,in_c_pol,in_h_pol,in_aa_const,in_bb_const,in_fq,in_kappa,out_RES1,out_RES2)

    USE MOD_PARAMETERS
    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: in_tauk, in_tauh, in_BY, in_K, in_H, in_C, in_r, in_w, in_kappa
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: in_f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: in_a_pol, in_c_pol, in_h_pol, in_aa_const, in_bb_const
	
	REAL(8), DIMENSION(m,n,n_q), INTENT(IN) :: in_fq
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: out_RES1, out_RES2
	
	! ==================================================================== !
    ! LOCAL
	! ==================================================================== !

    REAL(8) :: tauk, tauh, BY, K, H, C, r, w, kappa, Y, B, part_c_r, part_n_r, part_c_w, part_n_w, uc, ucc, ucn, SS1, SS2, SS3, SS4
   
    REAL(8), DIMENSION(m,n) :: f
	
	REAL(8), DIMENSION(m,n,p) :: a_pol, c_pol, h_pol, aa_const, bb_const
	
	REAL(8), DIMENSION(m,n,n_q) :: fq
	
	REAL (8), DIMENSION(:,:,:,:), ALLOCATABLE :: q_pol

    ! Outer loop
    REAL(8) ::  d, tol
    INTEGER ::  ite, maxite, i, j, jj, q, qq, i_q

    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
    
	ALLOCATE(q_pol(m,n,p,n_q))
	
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
	aa_const = in_aa_const
	bb_const = in_bb_const		
	fq       = in_fq
	kappa    = in_kappa   
	
	Y = (K**alpha)*(H**(1.0D0-alpha))
	B = BY*Y
    
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				DO i_q = 1,n_q
					q_pol(i,j,q,i_q) =  aa_const(i,j,q)+bb_const(i,j,q)*q_grid(i_q)
				ENDDO
			ENDDO
		ENDDO
	ENDDO		
	
	SS1 = 0.0D0
	SS2 = 0.0D0
	SS3 = 0.0D0
	SS4 = 0.0D0
	
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				uc  = (((c_pol(i,j,q)**gam)*((1.0D0-h_pol(i,j,q))**(1.0D0-gam)))**(1.0D0-sigma))*(gam/c_pol(i,j,q))
				ucc = -gam*(gam*sigma+1.0D0-gam)*((((c_pol(i,j,q)**gam)*((1.0D0-h_pol(i,j,q))**(1.0D0-gam)))**(1.0D0-sigma))/(c_pol(i,j,q)**2.0D0))
				ucn = (sigma-1.0D0)*gam*(1.0D0-gam)*((((c_pol(i,j,q)**gam)*((1.0D0-h_pol(i,j,q))**(1.0D0-gam)))**(1.0D0-sigma))/(c_pol(i,j,q)*(1.0D0-h_pol(i,j,q))))
				CALL PARTIAL_DERIVATIVE_C_R(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_r)
				CALL PARTIAL_DERIVATIVE_C_W(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_w)
				CALL PARTIAL_DERIVATIVE_N_R(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_r)
				CALL PARTIAL_DERIVATIVE_N_W(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_w)
				SS3 = SS3 + (tauh*w*part_n_r*(s(j)+e(q)*mult_e(j))+tauc0*part_c_r)*f(i,j)*Pe(q)
				SS4 = SS4 + (tauh*w*part_n_w*(s(j)+e(q)*mult_e(j))+tauc0*part_c_w)*f(i,j)*Pe(q)
				DO i_q = 1,n_q
					SS1 = SS1 + ((ap(i)*uc/(kappa*(1.0D0+tauc0)))+q_grid(i_q)*uc+(q_grid(i_q)*(1.0D0+(1.0D0-tauk)*r)-q_pol(i,j,q,i_q))*(ucc*part_c_r+ucn*part_n_r))*fq(i,j,i_q)*Pe(q)
					SS2 = SS2 + (((s(j)+e(q)*mult_e(j))*h_pol(i,j,q)*uc/(kappa*(1.0D0+tauc0)))+(q_grid(i_q)*(1.0D0+(1.0D0-tauk)*r)-q_pol(i,j,q,i_q))*(ucc*part_c_w+ucn*part_n_w))*fq(i,j,i_q)*Pe(q)
				ENDDO
			ENDDO
		ENDDO
	ENDDO
	
	out_RES1 = 1.0D0 - ((SS1+SS3)/(K+B))
	out_RES2 = 1.0D0 - ((SS2+SS4)/H)	
	
	DEALLOCATE(q_pol)
        
END SUBROUTINE RESIDUAL
    
    
    
