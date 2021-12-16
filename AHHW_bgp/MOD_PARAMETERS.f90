MODULE MOD_PARAMETERS

IMPLICIT NONE
SAVE

! ==================================================================== !		

INTEGER, PARAMETER  :: MORE_PRECISE = 1 

INTEGER, PARAMETER  :: GRID_TYPE = 4 ! 0 FOR LINEAR, 1 FOR LOG, 2 FOR DOUBLE LOG, 3 FOR STATIONARY POINTS METHOD

REAL(8), PARAMETER  :: TIME_BOUND = 108000.0D0

INTEGER :: ITE_MGR

REAL(8)  :: MGR_ADJ

! ==================================================================== !
! Parameters
! ==================================================================== !

REAL (8), PARAMETER ::  alpha	=	0.3780000000000000D0,	&
						beta	=	0.9537640920000000D0,	&
						sigma	=	2.0694102003728100D0,	&
						gam		=	0.5103869964000000D0,	&
						BY00	=	0.6150000000000000D0,	&
						delta	=	0.1042332268370610D0,	&
						s1		=	0.5803074850000000D0,	&
						s2		=	1.1525636286000000D0,	&
						s3		=	1.9261521752000000D0,	&
						s4		=	27.222932490355200D0,	&
						Ps11	=	0.9940266521000000D0,	&
						Ps12	=	0.0024244852000000D0,	&
						Ps13	=	0.0035220882000000D0,	&
						Ps21	=	0.0191389086000000D0,	&
						Ps22	=	0.9794480388000000D0,	&
						Ps23	=	0.0013244613000000D0,	&
						Ps31	=	0.0226370719000000D0,	&
						Ps32	=	0.0001224279000000D0,	&
						Ps33	=	0.9771916274000000D0,	&
						Ps41	=	0.0001611225000000D0,	&
						Ps42	=	0.0004117773000000D0,	&
						Ps43	=	0.0120005967000000D0,	&
						e1		=	-0.574036347700000D0,	&
						e2		=	-0.232049432900000D0,	&
						e3		=	0.1138448392000000D0,	&
						e4		=	0.1329957689000000D0,	&
						e5		=	0.8173100202000000D0,	&
						e6		=	1.2450563548000000D0,	&
						Pe1		=	0.2627937879000000D0,	&
						Pe2		=	0.0032936749000000D0,	&
						Pe3		=	0.5564683969000000D0,	&
						Pe4		=	0.0005142282000000D0,	&
						Pe5		=	0.0006054404000000D0,	&
						amin_	=	-0.078125514700000D0,	&
						amax_	=	500.00000000000000D0,	&
						LS0		=	0.0882761188067851D0,	&
						G       =   0.06866947216161348D0,&
						eta_transitory	=	1.1533470554000000D0

! ==================================================================== !		
! For More Elaborate Income Process
! ==================================================================== !		

 REAL (8), PARAMETER :: Ps14 = 1.0D0-Ps11-Ps12-Ps13, &
                        Ps24 = 1.0D0-Ps21-Ps22-Ps23, &
                        Ps34 = 1.0D0-Ps31-Ps32-Ps33, &
                        Ps44 = 1.0D0-Ps41-Ps42-Ps43, &
                        Pe6  = 1.0D0-Pe1 -Pe2 -Pe3 - Pe4 - Pe5
 
 INTEGER, PARAMETER  :: n = 4, &
                        p = 6
 
 REAL (8), PARAMETER, DIMENSION (n) ::  s = (/s1 , s2, s3, s4/)
 
 REAL (8), PARAMETER, DIMENSION (n) :: mult_e = s**eta_transitory
 
 
 REAL (8), PARAMETER, DIMENSION (n,n) :: Ps = reshape( (/ Ps11 , Ps21 , Ps31 , Ps41,  &    ! This representation is
                                                          Ps12 , Ps22 , Ps32 , Ps42,  &    ! the transpose of the real one
                                                          Ps13 , Ps23 , Ps33 , Ps43,  &
                                                          Ps14 , Ps24 , Ps34 , Ps44 /) , shape(Ps))
 
 REAL (8), PARAMETER, DIMENSION (p) ::  e  = (/e1  , e2 , e3 , e4 , e5 , e6 /)
 REAL (8), PARAMETER, DIMENSION (p) ::  Pe = (/Pe1 , Pe2, Pe3, Pe4, Pe5, Pe6 /)

! ==================================================================== !
! Stationary distribution for s
! ==================================================================== !

REAL (8), DIMENSION (n,n) :: SSPs

! ==================================================================== !
! Grid for a'
! ==================================================================== !

INTEGER, PARAMETER  :: m  = 200  ! Points in the grid for a'
!INTEGER, PARAMETER  :: m  = 400  ! Points in the grid for a'

INTEGER :: i_grid

REAL (8), PARAMETER, DIMENSION (m) :: ap  = (/ ( exp(((DBLE(i_grid)-1.0D0)/DBLE(m-1))*&
                                             log(exp(((DBLE(i_grid)-1.0D0)/DBLE(m-1))*&
                                             log(1.0D0+amax_-amin_))))-1.0D0+amin_ , i_grid = 1,m ) /)

REAL (8), PARAMETER :: amin   = ap(1), &
					   amax   = ap(m)
					   
! ==================================================================== !
! Grid for q
! ==================================================================== !

! INTEGER, PARAMETER  :: n_q = 1000 !m*n*p

! REAL(8), DIMENSION(n_q) :: q_grid  

INTEGER, PARAMETER  :: n_aux = 2000
!INTEGER, PARAMETER  :: n_aux = 4000
INTEGER             :: n_q

!REAL(8), DIMENSION(p) :: q_grid
REAL(8), ALLOCATABLE :: q_grid(:)  


! ==================================================================== !
! Initial tax system and Initial guess for interest rates
! ==================================================================== !

REAL (8), PARAMETER :: tauh0 = 0.225D0, &
                       tauk0 = 0.415D0, &
                       tauc0 = 0.047D0

! ==================================================================== !
! Count of errors and transitions computed
! ==================================================================== !

INTEGER :: ERRORS     ,&     ! shared variable
           ER         ,&     ! private variable
           ER_TYPE           ! private variable

! ==================================================================== !
! Global optimization
! ==================================================================== !

INTEGER, PARAMETER :: NVAR = 3

REAL(8), DIMENSION(NVAR) :: XXU, XXL


INTEGER, PARAMETER :: IMAX         = 400*20,        &           ! number of (Sobol) draws for each global iteration    [EXP07]
                      N_r          = 1,           &           ! number of BOBYQAs run in each global iteration
                      max_glob_ite = 1                          ! maximum number of iterations
					  
INTEGER, PARAMETER  :: np          = N_r
INTEGER, PARAMETER  :: d_each      = IMAX/np				  

REAL (8), PARAMETER :: criterium = sqrt(NVAR*(0.01D0**2.0D0)) ! distance to decide if new local minimum

! ==================================================================== !
! Local Optimization (BOBYQA)
! ==================================================================== !

REAL (8), PARAMETER :: RHOBEG   = 5.0D-2, &
                       RHOEND   = 5.0D-8 
INTEGER, PARAMETER ::  IPRINT   = 2, &           
		               NPT      = 2*NVAR+1, &
		               MAXFUN   = 1000, &    ! maximum number of iterations in BOBYQA
                       USE_PREV = 2*NVAR+2   ! when BOBYQA starts to use previous KK and AC as initial guess


REAL(8) :: time_start_bound, time_end_bound
REAL(8), ALLOCATABLE :: fig_inputs(:,:), fig_temp(:,:)

REAL(8) ::     BY0_PREV, H0_PREV

INTEGER :: ct, gi, median_e_ind, ite_bob, ite_tauk_mgr

INTEGER, DIMENSION(n,p)  :: hp_thr
INTEGER, DIMENSION(n,p)  :: ap_thr


END MODULE MOD_PARAMETERS
