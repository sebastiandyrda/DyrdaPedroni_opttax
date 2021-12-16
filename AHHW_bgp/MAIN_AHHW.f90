PROGRAM MAIN_AHHW

    !--------------------------------------------------------------------------------------!
    ! PROGRAM     : MAIN_AHHW.f90
    ! PROGRAMMER  : Sebastian Dyrda (U Toronto), Marcelo Pedroni (U Amsterdam)
    ! DESCRIPTION : This codes evaluates residuals from the FOCs in the Ramsey
    !               steady-state derived based on the Ackigoz et. al (2018) method
    !               for BGP preferences. See Online Appendix for Dyrda, Pedroni(2021)
    !               for derivations.
    !               http://dyrda.info/files/Dyrda_Pedroni_2021_Appendix.pdf
    ! DATE        : 15/12/2021
    !--------------------------------------------------------------------------------------!
    
     USE MOD_PARAMETERS
     USE MPI
     IMPLICIT NONE

     REAL(8) :: tauk, tauh, r, res1_dp, res2_dp, res3_dp, &
                res1_ahhw, res2_ahhw, res3_ahhw
	INTEGER  :: ierror			     ! returns error messages from the mpi subroutines 
	INTEGER  :: rank					 ! identification number of each processor
	INTEGER  :: nproc        			 ! total number of processors                
                
    !--------------------------------------------------------------------------------------!
    ! INITIALIZE THE MPI ENVIRONMENT
    !--------------------------------------------------------------------------------------!    
    CALL MPI_INIT(ierror)                         	! initialize MPI 
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)	! obtain rank for each node
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)	! obtain the number of nodes
            
    ERRORS = 0
	ER     = 0
    ct     = 0
    time_start_bound = MPI_WTIME()                

    !--------------------------------------------------------------------------------------!
	! DP Method optimal policy
    !--------------------------------------------------------------------------------------!
    tauk   = 0.266091653024420D0 
    tauh   = 0.390978064124877D0
    r      = 0.0481723988D0
    
    CALL MAIN_FUNCTION(tauk,tauh,r,res1_dp,res2_dp,res3_dp)

    !--------------------------------------------------------------------------------------!
	! AHHW Method, DP calibration
    !--------------------------------------------------------------------------------------!
	tauk   = 0.2456580920D0 
	tauh   = 0.3825493020D0
	r      = 0.048189108D0  

    CALL MAIN_FUNCTION(tauk,tauh,r,res1_ahhw,res2_ahhw,res3_ahhw) 

    !--------------------------------------------------------------------------------------!
    write(*,'(A5,A50)')  ' ','==========================================================='
    write(*,'(A5,A50)')  ' ', 'Residuals Comparison for DP calibration     '
    write(*,'(A5,A50)')  ' ','==========================================================='	
    write(*,'(A5,A15)')  ' ','DP Method'
    write(*,'(A5,A50)')  ' ','-----------------------------------------------------------'	    
    write(*,'(A5,3(A5,F10.6))')  ' ', 'R1 ', res1_dp,   'R2 ', res2_dp,   'R3 ', res3_dp
    write(*,'(A5,A50)')  ' ','-----------------------------------------------------------'	        
    write(*,'(A5,A15)')  ' ','AHHW Method'
    write(*,'(A5,A50)')  ' ','-----------------------------------------------------------'	        
    write(*,'(A5,3(A5,F10.6))')  ' ', 'R1 ', res1_ahhw, 'R2 ', res2_ahhw, 'R3 ', res3_ahhw
    write(*,'(A5,A50)')  ' ','==========================================================='     
  

       
    CALL MPI_FINALIZE(ierror)

END PROGRAM MAIN_AHHW


