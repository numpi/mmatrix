      SUBROUTINE DMMTRI( N, A, LDA, WORK, LWORK, INFO )

C     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, LWORK, N
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
C     ..

C     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     ..
C     .. LOCAL SCALARS ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
C     ..
C     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMM, DGEMV, DTRSM, DTRTRI,
     $                   XERBLA
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
C     ..
C     .. EXECUTABLE STATEMENTS ..
C
C     TEST THE INPUT PARAMETERS.
C
      INFO = 0
C     If NB from ILAENV is not optimal, we can load it to input 
      NB = ILAENV( 1, 'DMMTRI', ' ', N, -1, -1, -1 )
      LWKOPT = MAX( 1, N*NB )
      WORK( 1 ) = LWKOPT
C
      LQUERY = ( LWORK .EQ. -1 )
      IF ( N .LT. 0 ) THEN
         INFO = -1
      ELSE IF ( LDA .LT. MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF ( LWORK .LT. MAX( 1, N ) .AND. .NOT. LQUERY ) THEN
         INFO = -5
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DMMTRI', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
         RETURN
      END IF
C
C     QUICK RETURN IF POSSIBLE
C
      IF ( N .EQ. 0 )
     $   RETURN
C
C     FORM INV(U). IF INFO > 0 FROM DTRTRI, THEN U IS SINGULAR
C
      CALL DTRTRI( 'UPPER', 'NON-UNIT', N, A, LDA, INFO )
      IF ( INFO .GT. 0 )
     $   RETURN
C
      NBMIN = 2
      LDWORK = N
      IF ( NB .GT. 1 .AND. NB .LT. N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF ( LWORK .LT. IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DMMTRI', ' ',
     $                   N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
C
C     SOLVE INV(A)*L = INV(U) FOR INV(A)
C
      IF ( NB .LT. NBMIN .OR. NB .GE. N ) THEN
C
C        UNBLOCKED CODE
C
         DO 20 J = N, 1, -1
C
C           COPY COLUMN OF L TO WORK AND ZERO IT
C
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
C
C           COMPUTE COLUMN OF INV(A)
C
            IF ( J .LT. N )
     $         CALL DGEMV( 'N', N, N-J, -ONE,
     $                     A( 1, J+1 ), LDA,
     $                     WORK( J+1 ), 1,
     $                     ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
C
C        BLOCKED CODE
C
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
C
C           COPY BLOCK COLUMN OF L TO WORK AND ZERO IT
C
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I + ( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
C
C           COMPUTE BLOCK COLUMN OF INV(A)
C
            IF ( J+JB .LE. N )
     $         CALL DGEMM( 'N', 'N', N, JB, N-J-JB+1,
     $                     -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK,
     $                     ONE, A( 1, J ), LDA )
C
            CALL DTRSM( 'R', 'L', 'N', 'U', N, JB,
     $                  ONE, WORK( J ), LDWORK,
     $                  A( 1, J ), LDA )
   50    CONTINUE
      END IF
C
      WORK( 1 ) = IWS
      RETURN
C
C     END OF DMMTRI
C
      END
