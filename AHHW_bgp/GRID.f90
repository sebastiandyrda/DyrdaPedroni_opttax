SUBROUTINE GRID(n_case,len_grid,min,max,grid_out)
    ! Returns various types of grid
    ! Inputs:
    !           n_case: selects grid type
    !           len_grid : grid length
    !           min, max : smallest and largest grid node
    ! Output:
    !           grid_out
	
    IMPLICIT NONE

    INTEGER, INTENT(in)   :: n_case, len_grid
    REAL(8), INTENT(in)   :: min, max
    REAL(8), INTENT(out), DIMENSION(len_grid) :: grid_out

    INTEGER  :: i1
    real(8)  :: den



    IF (n_case<0 .OR. n_case>3) THEN
       PRINT *, "Error, the first argument of grid has to be between 0 and 3"
       STOP
    ENDIF
    den = DBLE(len_grid) - 1.d0

    !PRINT *, 'len_grid', len_grid

    SELECT CASE(n_case)

       CASE(0)
          ! Uniform grid
          grid_out(len_grid) = max
          grid_out(1) = min
          DO i1 = 2, len_grid
             grid_out(i1) = grid_out(1) + DBLE(i1-1)/den*&
                  (grid_out(len_grid)-grid_out(1))
          END DO

       CASE(1)
            ! Exponential grid
		    IF (min .GT. -1.0D0) THEN
			    grid_out(len_grid) = LOG(max+1.0D0)
			    grid_out(1) = LOG(min+1.0D0)
			    DO i1 = 2, len_grid
				    grid_out(i1) = grid_out(1) + DBLE(i1-1)/den*&
					  (grid_out(len_grid)-grid_out(1))
			    END DO
			    grid_out = EXP(grid_out)-1.0D0
			ELSE
				grid_out(len_grid) = LOG(max-min+1.0D0)
			    grid_out(1) = LOG(1.0D0)
			    DO i1 = 2, len_grid
				    grid_out(i1) = grid_out(1) + DBLE(i1-1)/den*&
					  (grid_out(len_grid)-grid_out(1))
			    END DO
			    grid_out = EXP(grid_out)-1.0D0
				grid_out = grid_out+min
			ENDIF			

       CASE(2)
            ! Double exponential grid
		    IF (min .GT. -1.0D0) THEN
			    grid_out(len_grid) = LOG(LOG(max+1.0D0)+1.0D0)
				grid_out(1) = LOG(LOG(min+1.0D0)+1.0D0)
				DO i1 = 2, len_grid
					grid_out(i1) = grid_out(1) + DBLE(i1-1)/den*&
					  (grid_out(len_grid)-grid_out(1))
				END DO
				grid_out = EXP(EXP(grid_out)-1.0D0)-1.0D0
			ELSE
				grid_out(len_grid) = LOG(LOG(max-min+1.0D0)+1.0D0)
				grid_out(1) = LOG(LOG(1.0D0)+1.0D0)
				DO i1 = 2, len_grid
					grid_out(i1) = grid_out(1) + DBLE(i1-1)/den*&
					  (grid_out(len_grid)-grid_out(1))
				END DO
				grid_out = EXP(EXP(grid_out)-1.0D0)-1.0D0
				grid_out = grid_out+min
			ENDIF		
          

       CASE(3)
          ! Triple exponential grid
          grid_out(len_grid) = LOG(LOG(LOG(max+1)+1)+1)
          grid_out(1) = LOG(LOG(LOG(min+1)+1)+1)
          DO i1 = 2, len_grid
             grid_out(i1) = grid_out(1) + DBLE(i1-1)/den*&
                  (grid_out(len_grid)-grid_out(1))
          END DO
          grid_out = EXP(EXP(EXP(grid_out)-1)-1)-1

       END SELECT

END SUBROUTINE GRID