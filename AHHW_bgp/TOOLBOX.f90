! module toolbox

! implicit none
!
      SUBROUTINE DLASRT2(ID, N, D, KEY, INFO)
!
!  -- ScaLAPACK routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      INTEGER            KEY(*)
      DOUBLE PRECISION   D(*)
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2!!32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER!1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  KEY     (input/output) INTEGER array, dimension (N)
!          On entry, KEY contains a key to each of the entries in D()
!          Typically, KEY(I) = I for all I
!          On exit, KEY is permuted in exactly the same manner as
!          D() was permuted from input to output
!          Therefore, if KEY(I) = I for all I upon input, then
!          D_out(I) = D_in(KEY(I))
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT, TMPKEY
      DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
!
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRT2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.1 ) RETURN
!
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.GT.0 ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
         IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                     TMPKEY = KEY( J )
                     KEY( J ) = KEY( J-1 )
                     KEY( J-1 ) = TMPKEY
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
!
         ELSE
!
!           Sort into increasing order
!
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                     TMPKEY = KEY( J )
                     KEY( J ) = KEY( J-1 )
                     KEY( J-1 ) = TMPKEY
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
!
         END IF
!
      ELSE IF( ENDD-START.GT.SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1.LT.D2 ) THEN
            IF( D3.LT.D1 ) THEN
               DMNMX = D1
            ELSE IF( D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3.LT.D2 ) THEN
               DMNMX = D2
            ELSE IF( D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
!
         IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J ).LT.DMNMX ) GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I ).GT.DMNMX ) GO TO 80
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               TMPKEY = KEY( J )
               KEY( J ) = KEY( I )
               KEY( I ) = TMPKEY
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
!
!           Sort into increasing order
!
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J ).GT.DMNMX ) GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I ).LT.DMNMX ) GO TO 110
            IF( I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               TMPKEY = KEY( J )
               KEY( J ) = KEY( I )
               KEY( I ) = TMPKEY
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 ) GO TO 10
!
!
      RETURN
!
!     End of DLASRT2
!
      END

function enorm ( n, x )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) enorm

  enorm = sqrt ( sum ( x(1:n)**2 ))

  return
end

function enorm2 ( n, x )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!
!! ENORM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This routine was named ENORM.  It has been renamed "ENORM2",
!    and a simplified routine has been substituted.
!
!    The Euclidean norm is computed by accumulating the sum of
!    squares in three different sums.  The sums of squares for the
!    small and large components are scaled so that no overflows
!    occur.  Non-destructive underflows are permitted.  Underflows
!    and overflows do not occur in the computation of the unscaled
!    sum of squares for the intermediate components.
!
!    The definitions of small, intermediate and large components
!    depend on two constants, RDWARF and RGIANT.  The main
!    restrictions on these constants are that RDWARF!!2 not
!    underflow and RGIANT!!2 not overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM2, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) agiant
  real    ( kind = 8 ) enorm2
  integer ( kind = 4 ) i
  real    ( kind = 8 ) rdwarf
  real    ( kind = 8 ) rgiant
  real    ( kind = 8 ) s1
  real    ( kind = 8 ) s2
  real    ( kind = 8 ) s3
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xabs
  real    ( kind = 8 ) x1max
  real    ( kind = 8 ) x3max

  rdwarf = sqrt ( tiny ( rdwarf ) )
  rgiant = sqrt ( huge ( rgiant ) )

  s1 = 0.0D+00
  s2 = 0.0D+00
  s3 = 0.0D+00
  x1max = 0.0D+00
  x3max = 0.0D+00
  agiant = rgiant / real ( n, kind = 8 )

  do i = 1, n

    xabs = abs ( x(i) )

    if ( xabs <= rdwarf ) then

      if ( x3max < xabs ) then
        s3 = 1.0D+00 + s3 * ( x3max / xabs )**2
        x3max = xabs
      else if ( xabs /= 0.0D+00 ) then
        s3 = s3 + ( xabs / x3max )**2
      end if

    else if ( agiant <= xabs ) then

      if ( x1max < xabs ) then
        s1 = 1.0D+00 + s1 * ( x1max / xabs )**2
        x1max = xabs
      else
        s1 = s1 + ( xabs / x1max )**2
      end if

    else

      s2 = s2 + xabs!!2

    end if

  end do
!
!  Calculation of norm.
!
  if ( s1 /= 0.0D+00 ) then

    enorm2 = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )

  else if ( s2 /= 0.0D+00 ) then

    if ( x3max <= s2 ) then
      enorm2 = sqrt ( s2 * ( 1.0D+00 + ( x3max / s2 ) * ( x3max * s3 ) ) )
    else
      enorm2 = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
    end if

  else

    enorm2 = x3max * sqrt ( s3 )

  end if

  return
end


      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER!1
!  CB      (input) CHARACTER!1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. INTA.GE.145 .AND. INTA.LE.153 .OR. INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. INTB.GE.145 .AND. INTB.LE.153 .OR. INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END

            SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER!(!)
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
!     ..
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
!
      STOP
!
 9999 FORMAT( ' !! On entry to ', A, ' parameter number ', I2, ' had ', 'an illegal value' )
!
!     End of XERBLA
!
    END
    
subroutine sub_myinterp1(x,f_x,xp,length_x,interp_value)

! From Villaverde
					 
	implicit none

	integer, intent(in) :: length_x

	real (8), dimension(length_x), intent(in) :: x, f_x, xp
	real (8), dimension(length_x), intent(out) :: interp_value
	
	integer, dimension(1) :: x_min
	integer :: index_k, x_index
	real (8) :: t
	
	
	do index_k = 1,length_x
	
		x_min = minloc(abs(x-xp(index_k)))
		x_index = x_min(1)
	
		if (xp(index_k)<x(x_index)) x_index = x_index-1
		x_index = max(x_index,1)
		x_index = min(x_index,length_x-1)
		t = (f_x(x_index+1)-f_x(x_index))/(x(x_index+1)-x(x_index))
		interp_value(index_k) = t*(xp(index_k)-x(x_index))+f_x(x_index)
	
	enddo

    end subroutine sub_myinterp1
    
subroutine sub_myinterp2(x,f_x,xp,length_x,interp_value)

	implicit none

	integer, intent(in) :: length_x

	real (8), dimension(length_x), intent(in) :: x, f_x
	real (8), intent(in) :: xp
	real (8), intent(out) :: interp_value
	
	integer, dimension(1) :: x_min
	integer :: x_index

	real (8) :: t
	
	x_min   = minloc(abs(x-xp))
	x_index = x_min(1)
	if (xp<x(x_index)) x_index = x_index-1
	x_index = max(x_index,1)
	x_index = min(x_index,length_x-1)
	 
	t = (f_x(x_index+1)-f_x(x_index))/(x(x_index+1)-x(x_index))
  
  interp_value = t*(xp-x(x_index))+f_x(x_index)

end subroutine sub_myinterp2    


subroutine reverse_array(dim,inArray,outArray)

      implicit none
      
      integer, intent(in) :: dim
      real(8), intent(in),  dimension(dim)  :: inArray
      real(8), intent(out), dimension(dim)  :: outArray
      
      real(8), dimension(dim) :: a
      integer                 :: Head, Tail
      real(8) :: Temp
      
            a    = inArray
            Head = 1                             ! start with the beginning
            Tail = dim                           ! start with the end
            do                                   
            if (Head >= Tail)  exit              !    if Head crosses Tail, exit
            Temp    = a(Head)                 ! otherwise, swap them
            a(Head) = a(Tail)
            a(Tail) = Temp
            Head    = Head + 1                ! move forward
            Tail    = Tail - 1                ! move backward
            end do                               
      
            outArray = a      
      end subroutine
      
      
      subroutine remove_dups_real(a,n,b,k)
      
      implicit none
      
      integer, intent(in)                :: n
      real(8), intent(in),  dimension(n) :: a
      real(8), intent(out), dimension(n) :: b
      integer, intent(out)               :: k
      
      integer :: i,j
      
      k     = 1
      b(1)  = a(1)
      
      outer: do i=2,n
            do j=1,k
                  if (abs(b(j) - a(i)) .le. 1D-12) then
                     cycle outer    
                  end if
            end do 
            k  = k + 1
            b(k)  = a(i)
      
      end do outer

end subroutine  

!end module toolbox

    

