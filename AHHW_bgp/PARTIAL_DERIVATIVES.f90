SUBROUTINE PARTIAL_DERIVATIVE_CP_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,jj,qq,part_cp_ap)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q,jj,qq
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_cp_ap	
	
	! ==================================================================== !
	! LOCAL
	! ==================================================================== !
	
	REAL(8), DIMENSION(1) :: cp, hp
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	!CALL interp_linear(1,m,ap,h_pol(:,jj,qq),1,a_pol(i,j,q),hp(1))
	if (hp_thr(jj,qq) .gt. 1) then	
		CALL interp_linear(1,hp_thr(jj,qq),ap(1:hp_thr(jj,qq)),h_pol(1:hp_thr(jj,qq),jj,qq),1,a_pol(i,j,q),hp(1))
	else
		hp(1)  = 0.0D0
	end if		
	hp(1) = MAX(hp(1),0.0D0)	
	
	IF (hp(1) .GT. 0.0D0) THEN	
		part_cp_ap = (1.0D0+(1.0D0-tauk)*r)*gam/(1.0D0+tauc0)
	ELSE
		part_cp_ap = (1.0D0+(1.0D0-tauk)*r)/(1.0D0+tauc0)
	ENDIF
	
        
END SUBROUTINE PARTIAL_DERIVATIVE_CP_AP

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_C_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_ap)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_c_ap	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_c_ap = -gam/(1.0D0+tauc0)
	ELSE
		part_c_ap = -1.0D0/(1.0D0+tauc0)
	ENDIF
		
	
END SUBROUTINE PARTIAL_DERIVATIVE_C_AP

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_C_T(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_T)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_c_T	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_c_T = gam/(1.0D0+tauc0)
	ELSE
		part_c_T = 1.0D0/(1.0D0+tauc0)
	ENDIF	
	
	
END SUBROUTINE PARTIAL_DERIVATIVE_C_T

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_N_T(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_T)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_n_T	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!	
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_n_T = -(1.0D0-gam)/((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))
	ELSE
		part_n_T = 0.0D0
	ENDIF			
		
        
END SUBROUTINE PARTIAL_DERIVATIVE_N_T

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_C_R(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_r)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_c_r	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!	
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_c_r = gam*ap(i)/(1.0D0+tauc0)
	ELSE
		part_c_r = ap(i)/(1.0D0+tauc0)
	ENDIF		
		
        
END SUBROUTINE PARTIAL_DERIVATIVE_C_R

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_N_R(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_r)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_n_r	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!	
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_n_r = -(1.0D0-gam)*ap(i)/((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))
	ELSE
		part_n_r = 0.0D0
	ENDIF	
		
        
END SUBROUTINE PARTIAL_DERIVATIVE_N_R

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_C_W(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_c_w)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_c_w	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!	
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_c_w = gam*(s(j)+e(q)*mult_e(j))/(1.0D0+tauc0)
	ELSE
		part_c_w = 0.0D0
	ENDIF		
		
        
END SUBROUTINE PARTIAL_DERIVATIVE_C_W

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_N_W(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_w)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_n_w	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!	
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_n_w = (gam-h_pol(i,j,q))/((1.0D0-tauh)*w)
	ELSE
		part_n_w = 0.0D0
	ENDIF
		
        
END SUBROUTINE PARTIAL_DERIVATIVE_N_W

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_NP_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,jj,qq,part_np_ap)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q,jj,qq
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_np_ap	
	
	! ==================================================================== !
	! LOCAL
	! ==================================================================== !
	
	REAL(8), DIMENSION(1) :: cp, hp
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	!CALL interp_linear(1,m,ap,h_pol(:,jj,qq),1,a_pol(i,j,q),hp(1))
	if (hp_thr(jj,qq) .gt. 1) then	
		CALL interp_linear(1,hp_thr(jj,qq),ap(1:hp_thr(jj,qq)),h_pol(1:hp_thr(jj,qq),jj,qq),1,a_pol(i,j,q),hp(1))
	else
		hp(1)  = 0.0D0
	end if	
	hp(1) = MAX(hp(1),0.0D0)	
	
	IF (hp(1) .GT. 0.0D0) THEN	
		!part_np_ap = -(1.0D0-gam)*(1.0D0+(1.0D0-tauk)*r)/((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))
		part_np_ap = -(1.0D0-gam)*(1.0D0+(1.0D0-tauk)*r)/((1.0D0-tauh)*w*(s(jj)+e(qq)*mult_e(jj)))
	ELSE
		part_np_ap = 0.0D0
	ENDIF	
	
        
END SUBROUTINE PARTIAL_DERIVATIVE_NP_AP

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_N_AP(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_ap)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_n_ap	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_n_ap = (1.0D0-gam)/((1.0D0-tauh)*w*(s(j)+e(q)*mult_e(j)))
	ELSE
		part_n_ap = 0.0D0
	ENDIF
		
	
END SUBROUTINE PARTIAL_DERIVATIVE_N_AP

! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_N_AP_0(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,part_n_ap)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_n_ap	
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	IF (h_pol(i,j,q) .GT. 0.0D0) THEN	
		part_n_ap = (1.0D0-gam)/((1.0D0-tauh)*w)
	ELSE
		part_n_ap = 0.0D0
	ENDIF
	
END SUBROUTINE PARTIAL_DERIVATIVE_N_AP_0
	
! ==================================================================== !
! ==================================================================== !
! ==================================================================== !

SUBROUTINE PARTIAL_DERIVATIVE_NP_AP_0(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,i,j,q,jj,qq,part_np_ap)

    USE MOD_PARAMETERS    
    IMPLICIT NONE
    
	! ==================================================================== !
    ! INPUTS
	! ==================================================================== !

    REAL(8), INTENT(IN) :: tauk, tauh, BY, K, H, C, r, w
	
	INTEGER, INTENT(IN) :: i,j,q,jj,qq
    
	REAL(8), DIMENSION(m,n), INTENT(IN) :: f
	
	REAL(8), DIMENSION(m,n,p), INTENT(IN) :: a_pol, c_pol, h_pol
	
	! ==================================================================== !
	! OUTPUTS
	! ==================================================================== !
	
	REAL(8), INTENT(OUT) :: part_np_ap	
	
	! ==================================================================== !
	! LOCAL
	! ==================================================================== !
	
	REAL(8), DIMENSION(1) :: cp, hp
	       
    !-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
	
	!CALL interp_linear(1,m,ap,h_pol(:,jj,qq),1,a_pol(i,j,q),hp(1))
	if (hp_thr(jj,qq) .gt. 1) then	
		CALL interp_linear(1,hp_thr(jj,qq),ap(1:hp_thr(jj,qq)),h_pol(1:hp_thr(jj,qq),jj,qq),1,a_pol(i,j,q),hp(1))
	else
		hp(1)  = 0.0D0
	end if		
	hp(1) = MAX(hp(1),0.0D0)	
	
	IF (hp(1) .GT. 0.0D0) THEN	
		part_np_ap = -(1.0D0-gam)*(1.0D0+(1.0D0-tauk)*r)/((1.0D0-tauh)*w)
	ELSE
		part_np_ap = 0.0D0
	ENDIF	
	
        
END SUBROUTINE PARTIAL_DERIVATIVE_NP_AP_0
		
	




