MODULE Cubic_Spline_GCV

!     ALGORITHM 642 COLLECTED ALGORITHMS FROM ACM.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-06-04  Time: 15:44:36
 
!     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.12, NO. 2,
!     JUN., 1986, P. 150.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


!   SUBROUTINE NAME     - CUBGCV

!--------------------------------------------------------------------------

! COMPUTER            - VAX/DOUBLE

! AUTHOR              - M.F.HUTCHINSON
!                       CSIRO DIVISION OF MATHEMATICS AND STATISTICS
!                       P.O. BOX 1965
!                       CANBERRA, ACT 2601
!                       AUSTRALIA

! LATEST REVISION     - 15 AUGUST 1985

! PURPOSE             - CUBIC SPLINE DATA SMOOTHER

! USAGE               - CALL CUBGCV (X, F, DF, N, Y, C, IC, VAR, JOB, SE, WK, IER)

! ARGUMENTS    X      - VECTOR OF LENGTH N CONTAINING THE ABSCISSAE OF THE N
!                         DATA POINTS (X(I),F(I)) I=1..N. (INPUT) X
!                         MUST BE ORDERED SO THAT X(I) < X(I+1).
!              F      - VECTOR OF LENGTH N CONTAINING THE ORDINATES
!                         (OR FUNCTION VALUES) OF THE N DATA POINTS (INPUT).
!              DF     - VECTOR OF LENGTH N. (INPUT/OUTPUT)
!                         DF(I) IS THE RELATIVE STANDARD DEVIATION OF THE
!                         ERROR ASSOCIATED WITH DATA POINT I.
!                         EACH DF(I) MUST BE POSITIVE.  THE VALUES IN DF ARE
!                         SCALED BY THE SUBROUTINE SO THAT THEIR MEAN SQUARE
!                         VALUE IS 1, AND UNSCALED AGAIN ON NORMAL EXIT.
!                         THE MEAN SQUARE VALUE OF THE DF(I) IS RETURNED
!                         IN WK(7) ON NORMAL EXIT.
!                         IF THE ABSOLUTE STANDARD DEVIATIONS ARE KNOWN,
!                         THESE SHOULD BE PROVIDED IN DF AND THE ERROR VARIANCE
!                         PARAMETER VAR (SEE BELOW) SHOULD THEN BE SET TO 1.
!                         IF THE RELATIVE STANDARD DEVIATIONS ARE UNKNOWN,
!                         SET EACH DF(I)=1.
!              N      - NUMBER OF DATA POINTS (INPUT).
!                         N MUST BE >= 3.
!              Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y IS A VECTOR OF LENGTH N.
!                         C IS AN N-1 BY 3 MATRIX.  THE VALUE OF THE SPLINE
!                         APPROXIMATION AT T IS
!                         S(T) = ((C(I,3)*D + C(I,2))*D + C(I,1))*D + Y(I)
!                         WHERE X(I) <= T < X(I+1) AND D = T - X(I).
!              IC     - ROW DIMENSION OF MATRIX C EXACTLY AS SPECIFIED IN THE
!                         DIMENSION STATEMENT IN THE CALLING PROGRAM. (INPUT)
!              VAR    - ERROR VARIANCE. (INPUT/OUTPUT)
!                         IF VAR IS NEGATIVE (I.E. UNKNOWN) THEN THE SMOOTHING
!                         PARAMETER IS DETERMINED BY MINIMIZING THE GENERALIZED
!                         CROSS VALIDATION AND AN ESTIMATE OF THE ERROR
!                         VARIANCE IS RETURNED IN VAR.
!                         IF VAR IS NON-NEGATIVE (I.E. KNOWN) THEN THE
!                         SMOOTHING PARAMETER IS DETERMINED TO MINIMIZE AN
!                         ESTIMATE, WHICH DEPENDS ON VAR, OF THE TRUE MEAN
!                         SQUARE ERROR, AND VAR IS UNCHANGED.
!                         IN PARTICULAR, IF VAR IS ZERO, THEN AN INTERPOLATING
!                         NATURAL CUBIC SPLINE IS CALCULATED.
!                         VAR SHOULD BE SET TO 1 IF ABSOLUTE STANDARD
!                         DEVIATIONS HAVE BEEN PROVIDED IN DF (SEE ABOVE).
!              JOB    - JOB SELECTION PARAMETER. (INPUT)
!                       JOB = 0 SHOULD BE SELECTED IF POINT STANDARD ERROR
!                         ESTIMATES ARE NOT REQUIRED IN SE.
!                       JOB = 1 SHOULD BE SELECTED IF POINT STANDARD ERROR
!                         ESTIMATES ARE REQUIRED IN SE.
!              SE     - VECTOR OF LENGTH N CONTAINING BAYESIAN STANDARD
!                         ERROR ESTIMATES OF THE FITTED SPLINE VALUES IN Y.
!                         SE IS NOT REFERENCED IF JOB=0. (OUTPUT)
!              WK     - WORK VECTOR OF LENGTH 7*(N + 2).  ON NORMAL EXIT THE
!                         FIRST 7 VALUES OF WK ARE ASSIGNED AS FOLLOWS:-

!                         WK(1) = SMOOTHING PARAMETER (= RHO/(RHO + 1))
!                         WK(2) = ESTIMATE OF THE NUMBER OF DEGREES OF
!                                 FREEDOM OF THE RESIDUAL SUM OF SQUARES
!                         WK(3) = GENERALIZED CROSS VALIDATION
!                         WK(4) = MEAN SQUARE RESIDUAL
!                         WK(5) = ESTIMATE OF THE TRUE MEAN SQUARE ERROR
!                                 AT THE DATA POINTS
!                         WK(6) = ESTIMATE OF THE ERROR VARIANCE
!                         WK(7) = MEAN SQUARE VALUE OF THE DF(I)

!                         IF WK(1)=0 (RHO=0) AN INTERPOLATING NATURAL CUBIC
!                         SPLINE HAS BEEN CALCULATED.
!                         IF WK(1)=1 (RHO=INFINITE) A LEAST SQUARES REGRESSION
!                         LINE HAS BEEN CALCULATED.
!                         WK(2) IS AN ESTIMATE OF THE NUMBER OF DEGREES OF
!                         FREEDOM OF THE RESIDUAL WHICH REDUCES TO THE
!                         USUAL VALUE OF N-2 WHEN A LEAST SQUARES REGRESSION
!                         LINE IS CALCULATED.
!                         WK(3),WK(4),WK(5) ARE CALCULATED WITH THE DF(I)
!                         SCALED TO HAVE MEAN SQUARE VALUE 1.  THE UNSCALED
!                         VALUES OF WK(3),WK(4),WK(5) MAY BE CALCULATED BY
!                         DIVIDING BY WK(7).
!                         WK(6) COINCIDES WITH THE OUTPUT VALUE OF VAR IF
!                         VAR IS NEGATIVE ON INPUT.  IT IS CALCULATED WITH
!                         THE UNSCALED VALUES OF THE DF(I) TO FACILITATE
!                         COMPARISONS WITH A PRIORI VARIANCE ESTIMATES.

!              IER    - ERROR PARAMETER. (OUTPUT)
!                       TERMINAL ERROR
!                         IER = 129, IC IS LESS THAN N-1.
!                         IER = 130, N IS LESS THAN 3.
!                         IER = 131, INPUT ABSCISSAE ARE NOT
!                           ORDERED SO THAT X(I) < X(I+1).
!                         IER = 132, DF(I) IS NOT POSITIVE FOR SOME I.
!                         IER = 133, JOB IS NOT 0 OR 1.

! PRECISION/HARDWARE  - DOUBLE

! REQUIRED ROUTINES   - SPINT1, SPFIT1, SPCOF1, SPERR1

! REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY THE SUBROUTINE
!              IS PROPORTIONAL TO N.  THE SUBROUTINE USES AN ALGORITHM
!              DEVELOPED BY M.F. HUTCHINSON AND F.R. DE HOOG,
!              'SMOOTHING NOISY DATA WITH SPLINE FUNCTIONS',
!              NUMER. MATH. (IN PRESS)

!-----------------------------------------------------------------------

SUBROUTINE cubgcv(x, f, df, n, y, c, ic, var, job, se, wk, ier)

!---SPECIFICATIONS FOR ARGUMENTS---

REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN OUT)  :: f(:)
REAL (dp), INTENT(OUT)     :: df(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: y(:)
INTEGER, INTENT(IN)        :: ic
REAL (dp), INTENT(IN OUT)  :: c(ic,3)
REAL (dp), INTENT(IN OUT)  :: var
INTEGER, INTENT(IN)        :: job
REAL (dp), INTENT(IN OUT)  :: se(:)
REAL (dp), INTENT(OUT)     :: wk(0:n+1,7)
INTEGER, INTENT(OUT)       :: ier

!---SPECIFICATIONS FOR LOCAL VARIABLES---
REAL (dp)  :: delta, ERR, gf1, gf2, gf3, gf4, r1, r2, r3, r4,  &
              avh, avdf, avar, stat(6), p, q

REAL (dp), PARAMETER  :: ratio = 2.0_dp, tau = 1.618033989_dp,   &
                         zero = 0.0_dp, one = 1.0_dp
INTEGER    :: i

!---INITIALIZE---
ier = 133
IF (job >= 0 .AND. job <= 1) THEN
  CALL spint1(x, avh, f, df, avdf, n, y, c, ic, wk, wk(0,4), ier)
  IF (ier == 0) THEN
    avar = var
    IF (var > zero) avar = var * avdf * avdf
    
!---CHECK FOR ZERO VARIANCE---
    IF (var == zero) THEN
      r1 = zero
    ELSE
      
!---FIND LOCAL MINIMUM OF GCV OR THE EXPECTED MEAN SQUARE ERROR---
      r1 = one
      r2 = ratio * r1
      CALL spfit1(x, avh, df, n, r2, p, q, gf2, avar, stat, y, c, ic, wk,  &
                  wk(0,4), wk(0,6), wk(0,7))
      10 CALL spfit1(x, avh, df, n, r1, p, q, gf1, avar, stat, y, c, ic, wk,  &
                     wk(0,4), wk(0,6), wk(0,7))
      IF (gf1 <= gf2) THEN
        
!---EXIT IF P ZERO---
        IF (p <= zero) GO TO 40
        r2 = r1
        gf2 = gf1
        r1 = r1 / ratio
        GO TO 10
      END IF
      
      r3 = ratio * r2
      20 CALL spfit1(x, avh, df, n, r3, p, q, gf3, avar, stat, y, c, ic, wk,  &
                     wk(0,4), wk(0,6), wk(0,7))
      IF (gf3 <= gf2) THEN
        
!---EXIT IF Q ZERO---
        IF (q <= zero) GO TO 40
        r2 = r3
        gf2 = gf3
        r3 = ratio * r3
        GO TO 20
      END IF
      
      r2 = r3
      gf2 = gf3
      delta = (r2-r1) / tau
      r4 = r1 + delta
      r3 = r2 - delta
      CALL spfit1(x, avh, df, n, r3, p, q, gf3, avar, stat, y, c, ic, wk,  &
                  wk(0,4), wk(0,6), wk(0,7))
      CALL spfit1(x, avh, df, n, r4, p, q, gf4, avar, stat, y, c, ic, wk,  &
                  wk(0,4), wk(0,6), wk(0,7))
      
!---GOLDEN SECTION SEARCH FOR LOCAL MINIMUM---
      30 IF (gf3 <= gf4) THEN
        r2 = r4
        gf2 = gf4
        r4 = r3
        gf4 = gf3
        delta = delta / tau
        r3 = r2 - delta
        CALL spfit1(x, avh, df, n, r3, p, q, gf3, avar, stat, y, c, ic, wk,  &
                    wk(0,4), wk(0,6), wk(0,7))
      ELSE
        
        r1 = r3
        gf1 = gf3
        r3 = r4
        gf3 = gf4
        delta = delta / tau
        r4 = r1 + delta
        CALL spfit1(x, avh, df, n, r4, p, q, gf4, avar, stat, y, c, ic, wk,  &
                    wk(0,4), wk(0,6), wk(0,7))
      END IF
      ERR = (r2-r1) / (r1+r2)
      IF (ERR*ERR+one > one .AND. ERR > 1.0D-6) GO TO 30
      r1 = (r1+r2) * 0.5D0
    END IF
    
!---CALCULATE SPLINE COEFFICIENTS---
    CALL spfit1(x, avh, df, n, r1, p, q, gf1, avar, stat, y, c, ic, wk,  &
                wk(0,4), wk(0,6), wk(0,7))
    40 CALL spcof1(x, avh, f, df, n, p, q, y, c, ic, wk(0,6), wk(0,7))
    
!---OPTIONALLY CALCULATE STANDARD ERROR ESTIMATES---
    IF (var < zero) THEN
      avar = stat(6)
      var = avar / (avdf*avdf)
    END IF
    IF (job == 1) CALL sperr1(x, avh, df, n, wk, p, avar, se)
    
!---UNSCALE DF---
    df(1:n) = df(1:n) * avdf
    
!--PUT STATISTICS IN WK---
    DO  i = 0, 5
      wk(i,1) = stat(i+1)
    END DO
    wk(5,1) = stat(6) / (avdf*avdf)
    wk(6,1) = avdf * avdf
    GO TO 70
  END IF
END IF

!---CHECK FOR ERROR CONDITION---
!     IF (IER .NE. 0) CONTINUE
70 RETURN
END SUBROUTINE cubgcv



SUBROUTINE spint1(x, avh, y, dy, avdy, n, a, c, ic, r, t, ier)

! INITIALIZES THE ARRAYS C, R AND T FOR ONE DIMENSIONAL CUBIC
! SMOOTHING SPLINE FITTING BY SUBROUTINE SPFIT1.  THE VALUES
! DF(I) ARE SCALED SO THAT THE SUM OF THEIR SQUARES IS N
! AND THE AVERAGE OF THE DIFFERENCES X(I+1) - X(I) IS CALCULATED
! IN AVH IN ORDER TO AVOID UNDERFLOW AND OVERFLOW PROBLEMS IN SPFIT1.

! SUBROUTINE SETS IER IF ELEMENTS OF X ARE NON-INCREASING,
! IF N IS LESS THAN 3, IF IC IS LESS THAN N-1 OR IF DY(I) IS
! NOT POSITIVE FOR SOME I.

!---SPECIFICATIONS FOR ARGUMENTS---

REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(OUT)     :: avh
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN OUT)  :: dy(:)
REAL (dp), INTENT(OUT)     :: avdy
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: a(:)
INTEGER, INTENT(IN)        :: ic
REAL (dp), INTENT(OUT)     :: c(ic,3)
REAL (dp), INTENT(OUT)     :: r(0:n+1,3)
REAL (dp), INTENT(OUT)     :: t(0:n+1,2)
INTEGER, INTENT(OUT)       :: ier


!---SPECIFICATIONS FOR LOCAL VARIABLES---
INTEGER    :: i
REAL (dp)  :: e, f, g, h
REAL (dp), PARAMETER  :: zero = 0.0_dp

!---INITIALIZATION AND INPUT CHECKING---
ier = 0
IF (n >= 3) THEN
  IF (ic < n-1) GO TO 60
  
!---GET AVERAGE X SPACING IN AVH---
  g = zero
  DO  i = 1, n - 1
    h = x(i+1) - x(i)
    IF (h <= zero) GO TO 70
    g = g + h
  END DO
  avh = g / (n-1)
  
!---SCALE RELATIVE WEIGHTS---
  g = zero
  DO  i = 1, n
    IF (dy(i) <= zero) GO TO 80
    g = g + dy(i) * dy(i)
  END DO
  avdy = SQRT(g/n)
  
  DO  i = 1, n
    dy(i) = dy(i) / avdy
  END DO
  
!---INITIALIZE H,F---
  h = (x(2)-x(1)) / avh
  f = (y(2)-y(1)) / h
  
!---CALCULATE A,T,R---
  DO  i = 2, n - 1
    g = h
    h = (x(i+1)-x(i)) / avh
    e = f
    f = (y(i+1)-y(i)) / h
    a(i) = f - e
    t(i,1) = 2.0D0 * (g+h) / 3.0D0
    t(i,2) = h / 3.0D0
    r(i,3) = dy(i-1) / g
    r(i,1) = dy(i+1) / h
    r(i,2) = -dy(i) / g - dy(i) / h
  END DO
  
!---CALCULATE C = R'*R---
  r(n,2) = zero
  r(n,3) = zero
  r(n+1,3) = zero
  DO  i = 2, n - 1
    c(i,1) = r(i,1) * r(i,1) + r(i,2) * r(i,2) + r(i,3) * r(i,3)
    c(i,2) = r(i,1) * r(i+1,2) + r(i,2) * r(i+1,3)
    c(i,3) = r(i,1) * r(i+2,3)
  END DO
  RETURN
END IF

!---ERROR CONDITIONS---
ier = 130
RETURN

60 ier = 129
RETURN

70 ier = 131
RETURN

80 ier = 132
RETURN
END SUBROUTINE spint1



SUBROUTINE spfit1(x, avh, dy, n, rho, p, q, fun, var, stat, a, c, ic, r, t, &
                  u, v)

! FITS A CUBIC SMOOTHING SPLINE TO DATA WITH RELATIVE WEIGHTING DY FOR A
! GIVEN VALUE OF THE SMOOTHING PARAMETER RHO USING AN ALGORITHM BASED ON THAT
! OF C.H. REINSCH (1967), NUMER. MATH. 10, 177-183.

! THE TRACE OF THE INFLUENCE MATRIX IS CALCULATED USING AN ALGORITHM DEVELOPED
! BY M.F.HUTCHINSON AND F.R.DE HOOG (NUMER. MATH., IN PRESS), ENABLING THE
! GENERALIZED CROSS VALIDATION AND RELATED STATISTICS TO BE CALCULATED
! IN ORDER N OPERATIONS.

! THE ARRAYS A, C, R AND T ARE ASSUMED TO HAVE BEEN INITIALIZED
! BY THE SUBROUTINE SPINT1.  OVERFLOW AND UNDERFLOW PROBLEMS ARE
! AVOIDED BY USING P=RHO/(1 + RHO) AND Q=1/(1 + RHO) INSTEAD OF
! RHO AND BY SCALING THE DIFFERENCES X(I+1) - X(I) BY AVH.

! THE VALUES IN DF ARE ASSUMED TO HAVE BEEN SCALED SO THAT THE
! SUM OF THEIR SQUARED VALUES IS N.  THE VALUE IN VAR, WHEN IT IS
! NON-NEGATIVE, IS ASSUMED TO HAVE BEEN SCALED TO COMPENSATE FOR
! THE SCALING OF THE VALUES IN DF.

! THE VALUE RETURNED IN FUN IS AN ESTIMATE OF THE TRUE MEAN SQUARE WHEN VAR IS
! NON-NEGATIVE, AND IS THE GENERALIZED CROSS VALIDATION WHEN VAR IS NEGATIVE.

!---SPECIFICATIONS FOR ARGUMENTS---

REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN OUT)  :: avh
REAL (dp), INTENT(IN)      :: dy(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: rho
REAL (dp), INTENT(OUT)     :: p
REAL (dp), INTENT(OUT)     :: q
REAL (dp), INTENT(OUT)     :: fun
REAL (dp), INTENT(IN)      :: var
REAL (dp), INTENT(OUT)     :: stat(6)
REAL (dp), INTENT(IN)      :: a(:)
INTEGER, INTENT(IN)        :: ic
REAL (dp), INTENT(IN)      :: c(ic,3)
REAL (dp), INTENT(OUT)     :: r(0:n+1,3)
REAL (dp), INTENT(IN)      :: t(0:n+1,2)
REAL (dp), INTENT(OUT)     :: u(0:n+1)
REAL (dp), INTENT(OUT)     :: v(0:n+1)


!---LOCAL VARIABLES---
INTEGER    :: i
REAL (dp)  :: e, f, g, h, zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, rho1

!---USE P AND Q INSTEAD OF RHO TO PREVENT OVERFLOW OR UNDERFLOW---
rho1 = one + rho
p = rho / rho1
q = one / rho1
IF (rho1 == one) p = zero
IF (rho1 == rho) q = zero

!---RATIONAL CHOLESKY DECOMPOSITION OF P*C + Q*T---
f = zero
g = zero
h = zero
DO  i = 0, 1
  r(i,1) = zero
END DO
DO  i = 2, n - 1
  r(i-2,3) = g * r(i-2,1)
  r(i-1,2) = f * r(i-1,1)
  r(i,1) = one / (p*c(i,1) + q*t(i,1) - f*r(i-1,2) - g*r(i-2,3))
  f = p * c(i,2) + q * t(i,2) - h * r(i-1,2)
  g = h
  h = p * c(i,3)
END DO

!---SOLVE FOR U---
u(0) = zero
u(1) = zero
DO  i = 2, n - 1
  u(i) = a(i) - r(i-1,2) * u(i-1) - r(i-2,3) * u(i-2)
END DO
u(n) = zero
u(n+1) = zero
DO  i = n - 1, 2, -1
  u(i) = r(i,1) * u(i) - r(i,2) * u(i+1) - r(i,3) * u(i+2)
END DO

!---CALCULATE RESIDUAL VECTOR V---
e = zero
h = zero
DO  i = 1, n - 1
  g = h
  h = (u(i+1)-u(i)) / ((x(i+1)-x(i))/avh)
  v(i) = dy(i) * (h-g)
  e = e + v(i) * v(i)
END DO
v(n) = dy(n) * (-h)
e = e + v(n) * v(n)

!---CALCULATE UPPER THREE BANDS OF INVERSE MATRIX---
r(n,1) = zero
r(n,2) = zero
r(n+1,1) = zero
DO  i = n - 1, 2, -1
  g = r(i,2)
  h = r(i,3)
  r(i,2) = -g * r(i+1,1) - h * r(i+1,2)
  r(i,3) = -g * r(i+1,2) - h * r(i+2,1)
  r(i,1) = r(i,1) - g * r(i,2) - h * r(i,3)
END DO

!---CALCULATE TRACE---
f = zero
g = zero
h = zero
DO  i = 2, n - 1
  f = f + r(i,1) * c(i,1)
  g = g + r(i,2) * c(i,2)
  h = h + r(i,3) * c(i,3)
END DO
f = f + two * (g+h)

!---CALCULATE STATISTICS---
stat(1) = p
stat(2) = f * p
stat(3) = n * e / (f*f)
stat(4) = e * p * p / n
stat(6) = e * p / f
IF (var < zero) THEN
  stat(5) = stat(6) - stat(4)
  fun = stat(3)
ELSE
  
  stat(5) = MAX(stat(4) - two*var*stat(2)/n + var, zero)
  fun = stat(5)
END IF
RETURN
END SUBROUTINE spfit1



SUBROUTINE sperr1(x, avh, dy, n, r, p, var, se)

! CALCULATES BAYESIAN ESTIMATES OF THE STANDARD ERRORS OF THE FITTED
! VALUES OF A CUBIC SMOOTHING SPLINE BY CALCULATING THE DIAGONAL ELEMENTS
! OF THE INFLUENCE MATRIX.

!---SPECIFICATIONS FOR ARGUMENTS---

REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: avh
REAL (dp), INTENT(IN)      :: dy(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: r(0:n+1,3)
REAL (dp), INTENT(IN)      :: p
REAL (dp), INTENT(IN)      :: var
REAL (dp), INTENT(OUT)     :: se(:)


!---SPECIFICATIONS FOR LOCAL VARIABLES---
INTEGER    :: i
REAL (dp)  :: f, g, h, f1, g1, h1
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

!---INITIALIZE---
h = avh / (x(2)-x(1))
se(1) = one - p * dy(1) * dy(1) * h * h * r(2,1)
r(1,1) = zero
r(1,2) = zero
r(1,3) = zero

!---CALCULATE DIAGONAL ELEMENTS---
DO  i = 2, n - 1
  f = h
  h = avh / (x(i+1)-x(i))
  g = -f - h
  f1 = f * r(i-1,1) + g * r(i-1,2) + h * r(i-1,3)
  g1 = f * r(i-1,2) + g * r(i,1) + h * r(i,2)
  h1 = f * r(i-1,3) + g * r(i,2) + h * r(i+1,1)
  se(i) = one - p * dy(i) * dy(i) * (f*f1+g*g1+h*h1)
END DO
se(n) = one - p * dy(n) * dy(n) * h * h * r(n-1,1)

!---CALCULATE STANDARD ERROR ESTIMATES---
DO  i = 1, n
  se(i) = SQRT(MAX(se(i)*var, zero)) * dy(i)
END DO
RETURN
END SUBROUTINE sperr1



SUBROUTINE spcof1(x, avh, y, dy, n, p, q, a, c, ic, u, v)

! CALCULATES COEFFICIENTS OF A CUBIC SMOOTHING SPLINE FROM
! PARAMETERS CALCULATED BY SUBROUTINE SPFIT1.

!---SPECIFICATIONS FOR ARGUMENTS---

REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN OUT)  :: avh
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN)      :: dy(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: p
REAL (dp), INTENT(IN)      :: q
REAL (dp), INTENT(OUT)     :: a(:)
INTEGER, INTENT(IN)        :: ic
REAL (dp), INTENT(OUT)     :: c(ic,3)
REAL (dp), INTENT(OUT)     :: u(0:n+1)
REAL (dp), INTENT(IN)      :: v(0:n+1)

!---SPECIFICATIONS FOR LOCAL VARIABLES---
INTEGER    :: i
REAL (dp)  :: h, qh

!---CALCULATE A---
qh = q / (avh*avh)
DO  i = 1, n
  a(i) = y(i) - p * dy(i) * v(i)
  u(i) = qh * u(i)
END DO

!---CALCULATE C---
DO  i = 1, n - 1
  h = x(i+1) - x(i)
  c(i,3) = (u(i+1) - u(i)) / (3.0_dp*h)
  c(i,1) = (a(i+1) - a(i)) / h - (h*c(i,3) + u(i)) * h
  c(i,2) = u(i)
END DO
RETURN
END SUBROUTINE spcof1

END MODULE Cubic_Spline_GCV

