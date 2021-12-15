SUBROUTINE MAIN_FUNCTION(in_tauk,in_tauh,in_r,out_RES1,out_RES2,out_RES3)

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
	
	REAL(8), INTENT(OUT) :: out_RES1, out_RES2, out_RES3
	
	! ==================================================================== !
    ! LOCAL
	! ==================================================================== !
	
	REAL(8) :: tauk, tauh, BY, K, H, C, r, w, kappa, q_L, q_H, RES1, RES2, AA
	
	REAL(8) :: start_se, finish_se, start_pfi, finish_pfi, start_dist, finish_dist
	
	REAL(8), DIMENSION(m,n)   :: f
	
	REAL(8), DIMENSION(m,n,p) :: a_pol, c_pol, h_pol, aa_const, bb_const
	
	!REAL(8), DIMENSION(m,n,n_q) :: fq
	REAL(8), ALLOCATABLE :: fq(:,:,:)	
	
	REAL(8), DIMENSION(m*n*p) :: q_grid_0, q_grid_1, int_grid_0

	real(8), allocatable :: q_aux(:), q_left(:), q_right(:), q_vec(:), q_nodups(:)
	integer, allocatable :: key(:), key_fig(:)	
	
	!REAL(8), DIMENSION(n_q) :: int_grid
	REAL(8) :: PercDev_MGR, PercDev_RES1, PercDev_RES2, MGR
	REAL(8), DIMENSIOn(3) :: tauk_vec
	
	INTEGER :: i, j, q, info, ii, ind0, status

	INTEGER :: size_dups, N_ITE_MGR, ITE_TAU
	
	INTEGER :: ierror, rank
	
	!-----------------------------------------------------------------------------------!
    ! EXECUTABLE
    !-----------------------------------------------------------------------------------!
    
	! ==================================================================== !
    ! Stationary distribution for s
    ! ==================================================================== !
    
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)	! obtain rank for each node
	
	time_start_bound = MPI_WTIME()

	if (ct .eq. 0) THEN
		ite_bob = ite_bob + 1
	end if	
	
    SSPs = MATMUL(Ps,Ps)
    DO i = 1,500
        SSPs = MATMUL(SSPs,Ps)
    ENDDO
	median_e_ind = 0
    DO q = 1,p
        IF (SUM(Pe(1:q)) .LE. 0.5D0) THEN
            median_e_ind = median_e_ind+1
        ENDIF
    ENDDO
	median_e_ind = MAX(median_e_ind,1)
	
	
	tauk = in_tauk
	tauh = in_tauh
	r    = in_r
	
	!tauk = tauk0
	!tauh = tauh0

	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! Best fit DP Method
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! tauk   = 0.266091653024420D0 
	! tauh   = 0.390978064124877D0
	! BY     = 1.541D0	
	
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! Best fit AHHW Method - r_DP
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! tauk   = 0.245397676087836D0 
	! tauh   = 0.382408847700694D0
	! tauk_vec = [0.245397676D0, 0.245397676D0, 0.265397676D0]

	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	! Best fit AHHW Method - r_MGR
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!   tauk   = 0.211755829D0
    !   tauh   = 0.389209838D0
	! tauk   = 0.179432728D0
	! tauh   = 0.396340371D0
	! tauk_vec = [0.201755829D0, 0.211755829D0, 0.225397676D0]	

	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! Best fit AHHW Method - 3 var
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tauk   = 0.2456580920D0 
	tauh   = 0.3825493020D0
	r      = 0.048189108D0

	! tauk_vec = [0.245397676D0, 0.245397676D0, 0.265397676D0]	

	! ==================================================================== !
    ! Stationary Equilibrium
	! ==================================================================== !

	! N_ITE_MGR = 22

	! allocate(fig_inputs(size(tauk_vec)*(2*N_ITE_MGR),9))
	! DO i  = 1,size(tauk_vec) 

	! ite_tauk_mgr = 0	
	! allocate(fig_temp(2*N_ITE_MGR,9),key_fig(2*N_ITE_MGR))	
	! MGR_ADJ = 0.00025D0	
	! DO ITE_MGR = 1,N_ITE_MGR      
	
	! 	ite_tauk_mgr = ite_tauk_mgr + 1
	! 	tauk     = tauk_vec(i)
	! 	start_se = MPI_WTIME()	
	! 	CALL STATIONARY_EQUIL(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol)
	! 	finish_se = MPI_WTIME()	
	
	! ENDDO                 

	! MGR_ADJ = -0.00025D0
	! DO ITE_MGR = 1,N_ITE_MGR      
		
	! 	ite_tauk_mgr = ite_tauk_mgr + 1
	! 	tauk     = tauk_vec(i)
	! 	start_se = MPI_WTIME()	
	! 	CALL STATIONARY_EQUIL(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol)
	! 	finish_se = MPI_WTIME()	
	! 	! !PRINT '("SEQ  Time = ",f12.1," seconds.")',finish_se-start_se
		
	! ENDDO                

	! key_fig = (/ (j,j=1,2*N_ITE_MGR) /)
	! call DLASRT2('I', 2*N_ITE_MGR, fig_temp(:,9), key_fig, info)
	! do j=1,8
	! 	fig_temp(:,j) = fig_temp(key_fig,j)
	! end do

	! do j=1,2*N_ITE_MGR
	! 	WRITE(*,'(A,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)') 'fig_temp check ', fig_temp(j,:)
	! end do	

	! fig_inputs(((2*N_ITE_MGR)*(i-1)+1):(2*N_ITE_MGR)*i,:) = fig_temp(:,:)

	! deallocate(fig_temp,key_fig)
	
	
	! END DO	

	! do j=1,(size(tauk_vec)*(2*N_ITE_MGR))
	! 	WRITE(*,'(A,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)') 'fig_inputs ', fig_inputs(j,:)
	! end do	



	DO ITE_TAU = 1,31

	!tauk  = 0.266091653024420D0 + (DBLE(ITE_TAU - 6)/5.0D0)*0.05D0
	!tauk  = 0.2456580920D0 + (DBLE(ITE_TAU - 16)/15.0D0)*0.05D0
	tauh  = 0.3825493020D0 + (DBLE(ITE_TAU - 16)/15.0D0)*0.05D0

	start_se = MPI_WTIME()	
		CALL STATIONARY_EQUIL(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol)
		!CALL STATIONARY_EQUIL_DEBT(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol)
	finish_se = MPI_WTIME()		



	!STOP                 !!!!!!!!!!!!!!!!!!!!!
	
	! ==================================================================== !
    ! Policy Function Iteration to find q'(a,e,q)
    ! ==================================================================== !	
	
	IF (ER .EQ. 0) THEN
		start_pfi = MPI_WTIME()
		CALL POLICY_FUNCTION_ITERATION(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,aa_const,bb_const)
		finish_pfi = MPI_WTIME()
		!PRINT '("PFI  Time = ",f12.1," seconds.")',finish_pfi-start_pfi
	ENDIF	
	
	q_L = MINVAL(aa_const/(1.0D0-bb_const))
	q_H = MAXVAL(aa_const/(1.0D0-bb_const))
	
	IF (GRID_TYPE .EQ. 0) THEN
    !    CALL GRID(0,n_q,q_L,q_H,q_grid)
	! ELSEIF (GRID_TYPE .EQ. 1) THEN   
	!    CALL GRID(1,n_q,q_L,q_H,q_grid)
	! ELSEIF (GRID_TYPE .EQ. 2) THEN   
	!    CALL GRID(2,n_q,q_L,q_H,q_grid)
	! ELSEIF (GRID_TYPE .EQ. 3) THEN 
    !     ii = 0
	! 	DO i = 1,m
	! 		DO j = 1,n
	! 			DO q = 1,p
	! 				ii = ii + 1
	! 				q_grid_0(ii) = aa_const(i,j,q)/(1.0D0-bb_const(i,j,q))
	! 			ENDDO
	! 		ENDDO
	! 	ENDDO		
	!     CALL DLASRT('I',m*n*p,q_grid_0,info)
	! 	ii = 0
	! 	DO i = 1,m*n*p-1
	! 		IF (q_grid_0(i) .NE. q_grid_0(i+1)) THEN
	! 		    ii = ii+1
	! 			q_grid_1(ii) = q_grid_0(i) 
	! 		ENDIF
	! 	ENDDO
	! 	DO i = 1,ii
	! 		int_grid_0(i) = DBLE(i-1)
	! 	ENDDO
	! 	DO i = 1,n_q
	! 		int_grid(i) = (DBLE(ii-1)/DBLE(n_q-1))*DBLE(i-1)
	! 	ENDDO
	! 	CALL interp_linear(1,ii,int_grid_0(1:ii),q_grid_1(1:ii),n_q,int_grid,q_grid)
	! 	ind0 = MINLOC(ABS(q_grid),1)
	! 	q_grid(ind0) = 0.0D0

	ELSE IF (GRID_TYPE .EQ. 4) THEN

		allocate(q_aux(n_aux),q_left(n_aux),q_right(n_aux),stat = status)
		if (status /= 0) then
			print*, 'Allocation not successfull'       
			stop
		end if 		
		allocate(q_vec(2*n_aux),q_nodups(2*n_aux),stat = status)
		if (status /= 0) then
			print*, 'Allocation not successfull'       
			stop
		end if 

		! Create a symmetric grid

		CALL grid(1,n_aux,0.0D0,abs(q_L),q_aux)
		q_left    = -q_aux
		! DO i=1,(n_aux)
		! 	PRINT *, 'q_left 1', i, q_left(i)
		! END DO	

		call reverse_array(n_aux,q_left,q_left)

		! DO i=1,(n_aux)
		! 	PRINT *, 'q_left 2', i, q_left(i)
		! END DO			

		CALL grid(1,n_aux,0.0D0,q_H,q_aux)
		q_right   =  q_aux
		q_vec     =  [q_left, q_right]

		! DO i=1,(2*n_aux)
		! 	PRINT *, 'q_vec', i, q_vec(i)
		! END DO			
		
		! Remove duplicates
		call remove_dups_real(q_vec,2*n_aux,q_nodups,size_dups)
		n_q    = size_dups

		allocate(q_grid(n_q),key(n_q), stat = status)
		if (status /= 0) then
			print*, 'Allocation q_grid not successfull'       
			stop
		end if

		q_grid    = q_nodups(1:n_q)

		! Sort in increasing order
		!call DLASRT2('I', p, q_grid, key, info)
		
	ENDIF

	allocate(fq(m,n,n_q),stat = status)
	if (status /= 0) then
		print*, 'Allocation fq not successfull'       
		stop
	end if	

	! ! ==================================================================== !
    ! ! Computing the stationary distribution of (a,e,q) associated with 
	! ! policy function a'(a,e) and q'(a,e,q)
    ! ! ==================================================================== !    
    
	IF (ER .EQ. 0) THEN
		start_dist = MPI_WTIME()
		CALL DISTRIBUTION_Q(a_pol,f,aa_const,bb_const,fq)  
		finish_dist = MPI_WTIME()	
		!PRINT '("DIST  Time = ",f12.1," seconds.")',finish_dist-start_dist
	ENDIF   
	
    ! ==================================================================== !
    ! Computing gamma
    ! ==================================================================== !	
	
	IF (ER .EQ. 0) THEN
		CALL Q_TO_KAPPA(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,aa_const,bb_const,fq,kappa)
	ENDIF
	
	! ==================================================================== !
    ! Computing the residual of the 3 remaining equations
    ! ==================================================================== !
    
	IF (ER .EQ. 0) THEN
		CALL RESIDUAL(tauk,tauh,BY,K,H,C,r,w,f,a_pol,c_pol,h_pol,aa_const,bb_const,fq,kappa,RES1,RES2)
	ELSE 
		RES1 = 10.0D0
		RES2 = 10.0D0
	ENDIF
	
	out_RES1 = RES1
	out_RES2 = RES2
	out_RES3 =  (1.0D0-beta*(1.0D0+r))

    MGR =  (1.0D0-beta)/beta
	!PercDev_MGR  = ((r - MGR)/MGR)*100
	PercDev_MGR  = (1.0D0-beta*(1.0D0+r))*100
	PercDev_RES1 = (RES1)*100
	PercDev_RES2 = (RES2)*100

	!WRITE(*,'(A,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10,F20.10)') 'ResFig', tauk,tauh,BY,PercDev_MGR, PercDev_RES1, PercDev_RES2, sqrt((PercDev_MGR**2)+(PercDev_RES1**2)+(PercDev_RES2**2))
    
    IF (ER .EQ. 0) THEN
		WRITE(*,'(F15.10,F15.10,F15.10,F15.10,F15.10,F15.10,F15.10,F15.10,A10,E15.6,E15.6,E15.6,E15.6,A10,F10.1,A10,I6,A,F10.1,A,F10.1,A,F10.1,I10)') &
             		tauk, tauh, BY, K, H, C, r, w,&
		           ' F = ', out_RES1,out_RES2,out_RES3,(out_RES1**2.0D0) + (out_RES2**2.0D0)+ (out_RES3**2.0D0), &
                   ' TIME = ', time_end_bound - time_start_bound, ' RANK = ', rank, &
                   ' TIME SE = ', finish_se-start_se, &
				   ' TIME PFI = ', finish_pfi-start_pfi, &
				   ' TIME DIST = ', finish_dist-start_dist, ct
	ELSE
		IF (ER_TYPE .eq. 1) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: DISTRIBUTION", ' RANK = ', rank
		ELSEIF (ER_TYPE .eq. 2) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: DECISION RULES", ' RANK = ', rank
		ELSEIF (ER_TYPE .eq. 3) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: STATIONARY EQUILIBRIUM", ' RANK = ', rank
		ELSEIF (ER_TYPE .eq. 4) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: STATIONARY EQUILIBRIUM", ' RANK = ', rank
		ELSEIF (ER_TYPE .eq. 5) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: POLICY FUNCTION ITERATION", ' RANK = ', rank
		ELSEIF (ER_TYPE .eq. 6) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: DISTRIBUTION Q", ' RANK = ', rank
		ELSEIF (ER_TYPE .eq. 7) THEN
			WRITE(*,"(A35,A,I8)") " ERROR: TIME BOUND", ' RANK = ', rank
		ELSE
			WRITE(*,"(A31,I4,A,I8)") " ERROR TYPE = ", ER_TYPE,' RANK = ', rank
		ENDIF
	ENDIF

	deallocate(q_aux,q_left,q_right,q_vec,q_nodups,q_grid,key,fq,stat = status)

	END DO

	



	STOP
    
END SUBROUTINE MAIN_FUNCTION
    
    