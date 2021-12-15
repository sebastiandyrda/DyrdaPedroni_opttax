PROGRAM MAIN

!-----------------------------------------------------------------------------
! PROGRAM     : MAIN.f90
! PROGRAMMER  : Sebastian Dyrda, Marcelo Pedroni, University of Minnesota
! DESCRIPTION : Global solver
! DATE        : 08/08/13
!-----------------------------------------------------------------------------

 USE MOD_PARAMETERS
 USE BOBYQA_MODULE
 USE MPI
 IMPLICIT NONE

   
!---------------------------------------------------------------------------------------------------------------------------!
! VARIABLES RELATED TO MPI
!---------------------------------------------------------------------------------------------------------------------------!
	INTEGER 		   :: ierror			     ! returns error messages from the mpi subroutines 
	INTEGER 		   :: rank					 ! identification number of each processor
	INTEGER 		   :: nproc        			 ! total number of processors
	INTEGER, PARAMETER :: id_master=0            ! the identification of a master thread	
	INTEGER            :: b_each                 ! number of iterations in bobyqua parts done by each thread
	INTEGER            :: d_start, d_end,&
						  b_start, b_end         ! indeces used to denote the start and end point for each thread's iterations
	REAL(8)            :: time_start,time_end, time_start_clean, time_end_clean   ! variables to measure time
	INTEGER 		   :: id_now, tag            ! index for the send/receive loop
	INTEGER, DIMENSION(MPI_STATUS_SIZE)          :: status
	
!---------------------------------------------------------------------------------------------------------------------------!	
! SOBOL SEQUENCE RELATED VARIABLES
!---------------------------------------------------------------------------------------------------------------------------!

	REAL(8), DIMENSION(NVAR)        :: draw, draw_B            
	INTEGER                         :: d, dd  ! indeces in the Sobol loops
	REAL(8), DIMENSION(NVAR)        :: XL,XU  ! lower and upper bounds
    REAL(8), DIMENSION(NVAR,d_each) :: X	  ! stores the Sobol draws	
    	
!---------------------------------------------------------------------------------------------------------------------------!	
! GLOBAL LOOP VARIABLES
!---------------------------------------------------------------------------------------------------------------------------!	

	REAL(8) 							 :: EW        ! expected number of local minima                             
	INTEGER    							 :: W, W_bound, W_prev
	INTEGER    							 :: glob_ite
	REAL(8), DIMENSION(d_each)           :: F
	REAL(8), DIMENSION(NVAR,2*N_r)       :: X_master
	REAL(8), DIMENSION(2*N_r)            :: F_master
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: MF_recv
	REAL(8), ALLOCATABLE, DIMENSION(:)   :: F_recv
	REAL(8), ALLOCATABLE, DIMENSION(:,:) :: X_LM_recv
	REAL(8), DIMENSION(N_r)              :: F_r
    REAL(8), DIMENSION(NVAR,N_r)         :: X_r	
    INTEGER, DIMENSION(d_each) 	 	     :: KEY_each 
	INTEGER, DIMENSION(2*N_r)            :: KEY
	INTEGER, DIMENSION(N_r) 	 	     :: KEY1
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: KEY2
	INTEGER 				             :: INFO, i, j, ii
	REAL(8), DIMENSION(NVAR)             :: X_glob_opt, XX_LM_min, XX_LM_max, XXL_old, XXU_old
	REAL(8)								 :: F_glob_opt
	REAL(8)					           	 :: omega_n, meas_B, crit_dist, sigma_alg, alpha_r, FF, DIFF, X1_TEMP
	REAL(8)                              :: X_BISEC, X_BISEC_L, X_BISEC_U, DIFF_BISEC
	REAL(8)                              :: X_LM_PCHIP(1)
	
!---------------------------------------------------------------------------------------------------------------------------!	
! LOCAL SOLVER VARIABLES
!---------------------------------------------------------------------------------------------------------------------------!		
	
	INTEGER 				                                           :: iter, N_CLUSTER, NBR, NF
	INTEGER, DIMENSION(N_r)                                            :: MF, NLM   ! NLM - number of local minima
	REAL(8), DIMENSION(N_r)                                            :: F_LM_0 
	REAL(8), DIMENSION((2*NVAR+1+5)*(2*NVAR+1+NVAR)+3*NVAR*(NVAR+5)/2) :: WSPACE
    REAL(8), DIMENSION(NVAR, N_r)               					   :: X_LM_0, X_LM_AC0
	REAL(8), DIMENSION(NVAR, max_glob_ite*N_r)  					   :: X_LM, X_LM_AC
	REAL(8), DIMENSION(max_glob_ite*N_r) 	    					   :: F_LM
	REAL(8), DIMENSION(NVAR,2*max_glob_ite*N_r) 					   :: CLUSTER
	REAL(8), DIMENSION(max_glob_ite*N_r,N_r)  						   :: DIST_LM	
	
!---------------------------------------------------------------------------------------------------------------------------!
! INITIALIZE THE MPI ENVIRONMENT
!---------------------------------------------------------------------------------------------------------------------------!    
	CALL MPI_INIT(ierror)                         	! initialize MPI 
    time_start = MPI_WTIME()	
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)	! obtain rank for each node
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)	! obtain the number of nodes

    
	ERRORS = 0
	
	time_start_bound = MPI_WTIME()
			
	IF (N_r .NE. np) THEN 
	   PRINT *, 'N_r should be equal to np'
	ENDIF
		
!---------------------------------------------------------------------------------------------------------------------------!	  
! SET THE BOUNDS FOR THE VARIABLES
!---------------------------------------------------------------------------------------------------------------------------!	 

	ERRORS = 0
    
	! Normalize the variables to be between 0 and 1  
	XU = 1.0D0
	XL = 0.0D0
	
!---------------------------------------------------------------------------------------------------------------------------!		
	
    ! XXU = (/  1.0D0,  1.0D0 /)
	! XXL = (/ -1.0D0, -1.0D0 /)
	! tauk   = 0.245397676087836D0 
	! tauh   = 0.382408847700694D0
    XXU = (/ 0.345397676087836D0, 0.482408847700694D0, 0.0581723988D0 /)
	XXL = (/ 0.145397676087836D0, 0.282408847700694D0, 0.0381723988D0 /)	

	

	
	W_bound = 8   ! how many local minima to use to adjust the bounds XXL and XXU	
	
!----------------------------------------------------------------------------------------------------------------------------!
! SET THE VARIABLES FOR THE GLOBAL SOLVER
!----------------------------------------------------------------------------------------------------------------------------!	
    
    omega_n   = (3.141592653589793238462643383279D0**(DBLE(NVAR)/2.0D0))/GAMMA(1.0D0+(DBLE(NVAR)/2.0D0))
    meas_B    = product(XU-XL)
    sigma_alg = 2.0D0 
    alpha_r   = DBLE(N_r)/DBLE(IMAX)

!----------------------------------------------------------------------------------------------------------------------------!
! DRAW A SOBOL SEQUENCE (ONLY MASTER THREAD EXECUTES THIS PART)
!----------------------------------------------------------------------------------------------------------------------------!
	
    IF (rank .EQ. id_master) THEN     
		CALL INSOBL(NVAR,2**29)
    ENDIF

!---------------------------------------------------------------------------------------------------------------------------!
! Set W = 0 (number of local minima), glob_ite = 0 (global iteration) and start the GLOBAL loop
!---------------------------------------------------------------------------------------------------------------------------!

	W 		   = 0
	glob_ite   = 0  ! This is k from the paper
	gi         = 0
    EW         = 0.0D0
	
	! Compute the start and the end point for iterations for each thread
	d_start = rank*d_each+1 
	d_end   = min((rank+1)*d_each,IMAX)
	IF (d_end .lt. d_start) THEN
	    PRINT *, 'd_end is less than d_start' 
	    STOP
	ENDIF
	

    ! Determine the number of BOBYQA evaluations done by each thread in each global iteration
	DO b_each=1,N_r
		IF (b_each*nproc .GE. N_r) EXIT
	END DO

    ! Compute the start and the end point for iterations for each thread
	b_start = rank*b_each+1 
	b_end   = min((rank+1)*b_each,N_r)
	IF (b_end .lt. b_start) THEN
	    PRINT *, 'b_end is less than b_start' 
	    STOP
	END IF
	
	DO WHILE (glob_ite .lt. max_glob_ite)
	    
        IF (rank .EQ. id_master) THEN 
			glob_ite  = glob_ite + 1
			gi        = gi + 1
		END IF 
		CALL MPI_BCAST(glob_ite,1,mpi_integer,id_master,MPI_COMM_WORLD,ierror)
		CALL MPI_BCAST(gi,1,mpi_integer,id_master,MPI_COMM_WORLD,ierror)	
        
		tag = 1111
        IF (rank .EQ. id_master) THEN     
			
			!Draw and send to each other thread
			DO i = 1,nproc-1
				DO d = 1,d_each
					CALL I4_SOBOL(NVAR,draw)
					X(:,d) = draw
					! X(1,d) = XXL(1) + (XXU(1) - XXL(1))*DBLE(d-1)/39.0D0
					! X(2,d) = XXL(2) + (XXU(2) - XXL(2))*DBLE(i-1)/39.0D0
				ENDDO
				CALL MPI_SEND (X,NVAR*d_each,MPI_REAL8,i,tag,MPI_COMM_WORLD,ierror)
        	ENDDO
			
			!Draw for master thread
			DO d = 1,d_each
			    CALL I4_SOBOL(NVAR,draw)
				X(:,d) = draw
				! X(1,d) = XXL(1) + (XXU(1) - XXL(1))*DBLE(d-1)/39.0D0
				! X(2,d) = XXL(2) + (XXU(2) - XXL(2))*DBLE(40-1)/39.0D0				
			ENDDO		
			
		ELSE
		    CALL MPI_RECV (X,NVAR*d_each,MPI_REAL8,0,tag,MPI_COMM_WORLD,status,ierror)
        ENDIF 	

		crit_dist = ((meas_B*sigma_alg*(log(DBLE(glob_ite)*DBLE(N_r))/(omega_n*DBLE(glob_ite)*DBLE(N_r)))))**(1.0D0/DBLE(NVAR))

        !---------------------------------------------------------------------------------------------------------------------------!	
		! COMPUTING FUNCTION FOR EACH DRAW
		!---------------------------------------------------------------------------------------------------------------------------!
				
		ct = 0
		ite_bob = 0	
        DO d = 1,d_each		
            ct = ct + 1
            CALL CALFUN(NVAR,X(:,d),F(d))     
		ENDDO		
				
        !---------------------------------------------------------------------------------------------------------------------------!	
		! EACH THREAD ORDERS ITS RESULTS	
		!---------------------------------------------------------------------------------------------------------------------------!
		
		! Order the vector F in the increasing order:
		KEY_each = (/ (d,d=1,d_each) /)
		CALL DLASRT2('I',d_each,F,KEY_each,INFO)

		! Order the sample points vector in a corresponding way:
		X = X(:,KEY_each)
		
		!---------------------------------------------------------------------------------------------------------------------------!	
		! GATHER AND ORDER THE RESULTS FROM THE SOBOL EVALUATIONS	
		!---------------------------------------------------------------------------------------------------------------------------!
  		
		tag = 2222
		
		IF (rank .EQ. id_master) THEN
			F_master(1)   = F(1)
			X_master(:,1) = X(:,1)
			DO i = 1,nproc-1
				CALL MPI_RECV(F_master(i+1),1,MPI_REAL8,i,tag,MPI_COMM_WORLD,status,ierror)
				CALL MPI_RECV(X_master(:,i+1),NVAR,MPI_REAL8,i,tag,MPI_COMM_WORLD,status,ierror)
			ENDDO            
		ELSE
			CALL MPI_SEND(F(1),1,MPI_REAL8,id_master,tag,MPI_COMM_WORLD,ierror)
			CALL MPI_SEND(X(:,1),NVAR,MPI_REAL8,id_master,tag,MPI_COMM_WORLD,ierror)
		ENDIF
		
		DO j = 2,d_each
		
			IF (rank .EQ. id_master) THEN
				F_master(N_r+1)   = F(j)
				X_master(:,N_r+1) = X(:,j)
				DO i = 1,nproc-1
					CALL MPI_RECV(F_master(N_r+i+1),1,MPI_REAL8,i,tag,MPI_COMM_WORLD,status,ierror)
					CALL MPI_RECV(X_master(:,N_r+i+1),NVAR,MPI_REAL8,i,tag,MPI_COMM_WORLD,status,ierror)
				ENDDO            
			ELSE
				CALL MPI_SEND(F(j),1,MPI_REAL8,id_master,tag,MPI_COMM_WORLD,ierror)
				CALL MPI_SEND(X(:,j),NVAR,MPI_REAL8,id_master,tag,MPI_COMM_WORLD,ierror)
			ENDIF
			
			IF (rank .EQ. id_master) THEN
			
				! Order the vector F in the increasing order:
				KEY = (/ (d,d=1,2*N_r) /)
				CALL DLASRT2('I',2*N_r,F_master,KEY,INFO)
		
				! Order the sample points vector in a corresponding way:
				X_master = X_master(:,KEY)
				
			ENDIF
		
        ENDDO
		  		
	    IF (rank .EQ. id_master) THEN 
	
	        ! Define the reduced sample set
	        X_r = X_master(:,1:N_r)
            F_r = F_master(1:N_r)	
				
	    ENDIF		
	
        CALL MPI_BCAST(X_r,NVAR*N_r,mpi_double_precision,id_master,MPI_COMM_WORLD,ierror)        
		
	    !---------------------------------------------------------------------------------------------------------------------------!	
	    ! START THE LOCAL SOLVER PART	
	    !---------------------------------------------------------------------------------------------------------------------------!	
        
        MF = 0

	    ! RUN BOBYQA FOR THE REDUCED SET
		ct = 0
		ite_bob = 0
	    DO iter = b_start,b_end	
		    X_LM_0(:,iter) = X_r(:,iter)
			CALL BOBYQA(NVAR,NPT,X_LM_0(:,iter),XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,NF,FF,CALFUN)
			F_LM_0(iter) = FF
			IF (NF .EQ. MAXFUN) THEN
				MF(iter) = 1
			ELSE
			    MF(iter) = 0
			ENDIF
		ENDDO			
		
	    !---------------------------------------------------------------------------------------------------------------------------!	
	    ! GATHER AND ORGANIZE THE RESULTS FROM BOBYQUA	
	    !---------------------------------------------------------------------------------------------------------------------------!		
	    IF (rank .EQ. id_master) THEN
		    ALLOCATE(F_recv(b_each),X_LM_recv(NVAR,b_each),MF_recv(b_each),stat=ierror)
		    IF (ierror/=0) THEN
			    PRINT *, 'Allocation error in F_id'
			    STOP
		    ENDIF	
		    DO id_now =1,nproc-1
				tag = id_now + 100
				CALL MPI_RECV(F_recv,b_each,mpi_double_precision,id_now,tag,MPI_COMM_WORLD,status,ierror)				
				F_LM_0((id_now*b_each+1):(id_now+1)*b_each) = F_recv
				CALL MPI_RECV(X_LM_recv,b_each*NVAR,mpi_double_precision,id_now,tag,MPI_COMM_WORLD,status,ierror)
				X_LM_0(:,(id_now*b_each+1):(id_now+1)*b_each) = X_LM_recv
				CALL MPI_RECV(MF_recv,b_each,MPI_INTEGER,id_now,tag,MPI_COMM_WORLD,status,ierror)				
				MF((id_now*b_each+1):(id_now+1)*b_each) = MF_recv
		    ENDDO
		    DEALLOCATE(F_recv,X_LM_recv,MF_recv,stat=ierror)
		    IF (ierror .NE. 0) THEN
			    PRINT *, 'Allocation error in F_id'
		    STOP
		    ENDIF	   
	    ELSE
		    tag = rank + 100
		    CALL MPI_SEND(F_LM_0(b_start:b_end),b_each,mpi_double_precision,id_master,tag,MPI_COMM_WORLD,ierror)
		    CALL MPI_SEND(X_LM_0(:,b_start:b_end),b_each*NVAR,mpi_double_precision,id_master,tag,MPI_COMM_WORLD,ierror)					  
            CALL MPI_SEND(MF(b_start:b_end),b_each,MPI_INTEGER,id_master,tag,MPI_COMM_WORLD,ierror)					  
	    ENDIF
	    
		!---------------------------------------------------------------------------------------------------------------------------!	
	    ! ORDER RESULTS FROM BOBYQUA	
	    !---------------------------------------------------------------------------------------------------------------------------!
		
	    IF (rank .EQ. id_master) THEN 
	
		    ! Order the vector F_LM_0 in the increasing order.
		    KEY1 = (/ (d,d=1,N_r) /)
		    CALL DLASRT2('I',N_r,F_LM_0,KEY1,INFO)

		    ! Order the sample points vector in a corresponding way
		    X_LM_0 = X_LM_0(:,KEY1)
            MF     = MF(KEY1)
			
			! Convert into the actual parameter values
			DO j = 1,N_r 
				FORALL (i=1:NVAR) X_LM_AC0(i,j) = XXL(i)+X_LM_0(i,j)*(XXU(i)-XXL(i))
			ENDDO
		
		    ! Check which points generate new local minima. Check on the domain and check
		    ! whether the BOBYQA hit maxite:
		    
            NLM = 1 !initially all are considered new local minima

            IF (glob_ite .eq. 1) THEN            

		        DO j = 1,(N_r-1)
			        IF (NLM(j) .eq. 1) THEN
			        DO i = (j+1),N_r
                        IF (sqrt(sum((X_LM_AC0(:,i)-X_LM_AC0(:,j))**2.0D0)) .le. criterium) THEN
					        NLM(i) = 0
                        ENDIF
                        IF (MF(i) .eq. 1) THEN
					        NLM(i) = 0
                        ENDIF
			        ENDDO
			        ENDIF
		        ENDDO    
		        DO i = 1,N_r
			        IF (NLM(i) .eq. 1) THEN
				        X_LM(:,SUM(NLM(1:i))) = X_LM_0(:,i)
				        F_LM(SUM(NLM(1:i)))   = F_LM_0(i)
			        ENDIF			
		        ENDDO	
               W = SUM(NLM(:)) 
                
                ! Convert into the actual parameter values
                DO j = 1,W
                    FORALL (i=1:NVAR) X_LM_AC(i,j) = XXL(i)+X_LM(i,j)*(XXU(i)-XXL(i))
                ENDDO            

            ELSE
            
                DO j = 1,(N_r-1)
		        IF (NLM(j) .eq. 1) THEN
			        DO i = (j+1),N_r
                        IF (sqrt(sum((X_LM_AC0(:,i)-X_LM_AC0(:,j))**2.0D0)) .le. criterium) THEN
					        NLM(i) = 0
                        ENDIF
                        IF (MF(i) .eq. 1) THEN
					        NLM(i) = 0
                        ENDIF
			        ENDDO
			        ENDIF
		        ENDDO

                DO i = 1,N_r
				    DO j = 1,W
					    DIST_LM(j,i) = sqrt(sum((X_LM_AC0(:,i)-X_LM_AC(:,j))**2.0D0))
				    ENDDO
			    ENDDO

			    DO i = 1,N_r
                    IF  (MINVAL(DIST_LM(1:W,i)) .le. criterium) THEN
                        IF  (F_LM_0(i) .lt. F_LM(MINLOC(DIST_LM(1:W,i),1))) THEN
                        IF  (NLM(i) .eq. 1) THEN
						    X_LM_AC(:,MINLOC(DIST_LM(1:W,i),1)) = X_LM_AC0(:,i)                            
					        F_LM(MINLOC(DIST_LM(1:W,i),1))   = F_LM_0(i)
					    ENDIF 
                        ENDIF 
                        NLM(i) = 0                                   
                    ENDIF
			    ENDDO

			    DO i = 1,N_r
                   IF (NLM(i) .eq. 1) THEN
				        X_LM_AC(:,W+SUM(NLM(1:i))) = X_LM_AC0(:,i)					    
					    F_LM(W+SUM(NLM(1:i)))   = F_LM_0(i)
                    ENDIF
			    ENDDO
				
                W = W+SUM(NLM(:))
                
            ENDIF
	
        ENDIF	
	
	    IF (rank .EQ. id_master) THEN
		
            ! Order the vector F_LM in the increasing order:
		    ALLOCATE(KEY2(W),stat=ierror)
            KEY2 = (/ (d,d=1,W) /)
		    CALL DLASRT2('I',W,F_LM(1:W),KEY2,INFO)

            ! Order the X_LM points vector in a corresponding way:		    
            X_LM_AC(:,1:W) = X_LM_AC(:,KEY2)
			
            DEALLOCATE(KEY2,stat=ierror)      

		    !EW = ((DBLE(W)*((DBLE(glob_ite)*DBLE(N_r))-1.0D0))/((DBLE(glob_ite)*DBLE(N_r))-DBLE(W)-2.0D0))
			
			EW = ((DBLE(SUM(NLM(:)))*((DBLE(N_r))-1.0D0))/((DBLE(N_r))-DBLE(SUM(NLM(:)))-2.0D0))
		
		    PRINT *, " "
		    WRITE( *, '(a,F16.4,a,i8)' ) 'EW =', EW, '  W =', W
		    PRINT *, " "	
                        
            DO i = 1,NVAR
                XX_LM_max(i) = MAXVAL(X_LM_AC(i,1:MIN(W_bound,W))) 
                XX_LM_min(i) = MINVAL(X_LM_AC(i,1:MIN(W_bound,W)))
                PRINT *, 'XX_LM: ',XX_LM_min(i),XX_LM_max(i)
            ENDDO

            XXL_old = XXL
            XXU_old = XXU

		DO i = 1,NVAR
            IF ((XX_LM_max(i)-XX_LM_min(i)) .gt. 0.5D0*(XXU_old(i)-XXL_old(i))) THEN
                XXU(i) = XX_LM_max(i) + 0.30D0*(XX_LM_max(i)-XX_LM_min(i))
                XXL(i) = XX_LM_min(i) - 0.30D0*(XX_LM_max(i)-XX_LM_min(i))
            ELSEIF (XX_LM_min(i) .gt. (XXU_old(i)-0.1D0*(XXU_old(i)-XXL_old(i)))) THEN
                XXU(i) = XXU_old(i) + 0.8D0*(XXU_old(i)-XXL_old(i))
                XXL(i) = XXU_old(i) - 0.2D0*(XXU_old(i)-XXL_old(i))
            ELSEIF (XX_LM_max(i) .lt. (XXL_old(i)+0.1D0*(XXU_old(i)-XXL_old(i)))) THEN
                XXU(i) = XXL_old(i) + 0.2D0*(XXU_old(i)-XXL_old(i))
                XXL(i) = XXL_old(i) - 0.8D0*(XXU_old(i)-XXL_old(i))
            ELSE
                XXU(i) = XX_LM_max(i) + 0.25D0*(XXU_old(i)-XXL_old(i))
                XXL(i) = XX_LM_min(i) - 0.25D0*(XXU_old(i)-XXL_old(i))
			ENDIF
        ENDDO

        DO i = 1,NVAR
            PRINT *, 'BOUNDS: ',XXL(i),XXU(i)
        ENDDO			
			
		
		PRINT *, '================================'
		PRINT *, '================================'
		PRINT *, 'W =',W, 'EW =',EW	
		PRINT *, '================================'
		PRINT *, '================================'
		
	    ENDIF
	    !Send the an glob_ite value to the other threads
	    CALL MPI_BCAST(W,1,mpi_integer,id_master,MPI_COMM_WORLD,ierror) 
	    CALL MPI_BCAST(EW,1,mpi_double_precision,id_master,MPI_COMM_WORLD,ierror)
        CALL MPI_BCAST(XXU,NVAR,mpi_double_precision,id_master,MPI_COMM_WORLD,ierror)	
        CALL MPI_BCAST(XXL,NVAR,mpi_double_precision,id_master,MPI_COMM_WORLD,ierror)	
    	
    ENDDO

	
!---------------------------------------------------------------------------------------------------------------------------!	
! GET THE RESULTS AND SAVE THEM IN THE .DAT FILES
!---------------------------------------------------------------------------------------------------------------------------!		
	
	IF (rank .EQ. id_master) THEN	    
		
		F_glob_opt    = MINVAL(F_LM(1:W))
		X_glob_opt(:) = X_LM_AC(:,MINLOC(F_LM(1:W),1))
		
		PRINT *, 'F_glob_opt is', F_glob_opt, 'send by thread', rank
		PRINT *, 'X_glob_opt    is', X_glob_opt   , 'send by thread', rank	
		PRINT *, 'Number of global iterations', glob_ite
		
		time_end = MPI_WTIME()
		PRINT *, 'MPI elapse time =',time_end-time_start

	END IF 
	
    CALL CALFUN(NVAR,X_glob_opt,F_glob_opt)

	IF (rank .EQ. id_master) THEN
	PRINT *, '================================'
	PRINT *, '================================'
	PRINT *, 'THE CODE IS COMPLETED'	
	PRINT *, '================================'
	PRINT *, '================================'	
	END IF	
   
    CALL MPI_FINALIZE(ierror)
	
	CONTAINS
	
	SUBROUTINE CALFUN (NV,in_X,out_F)

		USE MOD_PARAMETERS

		IMPLICIT NONE
		REAL(8),INTENT(OUT) :: out_F
		INTEGER,INTENT(IN) :: NV
		REAL(8),DIMENSION(:),INTENT(IN)  :: in_X

		REAL(8) :: OBJECTIVE
		REAL(8), DIMENSION(NV) :: X, RES
		INTEGER  :: i
				
        X = XXL + in_X*(XXU-XXL)
		
		CALL MAIN_FUNCTION(X(1),X(2),X(3),RES(1),RES(2),RES(3)) 
		
		IF (ER .EQ. 0) THEN
			out_F = (RES(1)**2.0D0) + (RES(2)**2.0D0)+ (RES(3)**2.0D0)
		ELSE
		    out_F = 1000000.0D0
		ENDIF		
		  
	END SUBROUTINE	
    	
END PROGRAM MAIN