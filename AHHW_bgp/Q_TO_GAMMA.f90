SUBROUTINE Q_TO_GAMMA(in_tauk,in_tauh,in_BY,in_K,in_H,in_C,in_r,in_w,in_f,in_a_pol,in_c_pol,in_h_pol,in_aa_const,in_bb_const,in_fq,out_gamma)

    USE MOD_PARAMETERS
    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: in_tauk, in_tauh, in_BY, in_K, in_H, in_C, in_r, in_w
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: in_f, in_a_pol, in_c_pol, in_h_pol, in_aa_const, in_bb_const
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: in_fq
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: out_gamma
	
	! ==================================================================== !
    ! LOCAL
	! ==================================================================== !

    REAL(8) :: tauk, tauh, BY, K, H, C, r, w, part_c_T, part_n_T, uc, ucc, SS1, Ss2, SS3
   
    REAL(8), DIMENSION(m,n) :: f, a_pol, c_pol, h_pol, aa_const, bb_const
	
	REAL(8), DIMENSION(m,n,p) :: fq, q_pol

    ! Outer loop
    REAL(8) ::  d, tol
    INTEGER ::  ite, maxite, i, j, jj, q     

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
	aa_const = in_aa_const
	bb_const = in_bb_const		
	fq       = in_fq
    
	DO i = 1,m
		DO j = 1,n
			DO q = 1,p
				q_pol(i,j,q) =  aa_const(i,j)+bb_const(i,j)*q_grid(q)
			ENDDO
		ENDDO
	ENDDO		
	
	SS1 = 0.0D0
	SS2 = 0.0D0
	SS3 = 0.0D0
	
	DO i = 1,m
		DO j = 1,n
		    uc  = c_pol(i,j)**(-sigma)
			ucc = -sigma*(c_pol(i,j)**(-sigma-1.0D0))
			CALL PARTIAL_DERIVATIVE_C_T(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,part_c_T)
			CALL PARTIAL_DERIVATIVE_N_T(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,part_n_T)
			SS1 = SS1 + uc*f(i,j)
			SS3 = SS3 + part_n_T*s(j)*f(i,j)
			DO q = 1,p
				SS2 = SS2 + (q_grid(q)*(1.0D0+(1.0D0-tauk)*r)-q_pol(i,j,q))*ucc*part_c_T*fq(i,j,q)
			ENDDO
		ENDDO
	ENDDO
	
	out_gamma = SS1/(1.0D0-SS2- tauh*w*SS3)	
        
END SUBROUTINE Q_TO_GAMMA
    
    
    
