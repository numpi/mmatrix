C *********************************************************************
C *  DMLTRFF (LEFT BLOCKED GTH using DMLTRF2 inside)                  *
C *  DMLTRFF computes a recursive blocked GTH factorization for a     *
C *  given M-matrix A and vectors V and W satisfying V'*A = W'.       *
C *                                                                   *
C *  Recursive blocked factorization of an M-by-N matrix:             *
C *                                                                   *
C *        A = L * U                                                  *
C *                                                                   *
C *  where:                                                           *
C *     - L is lower triangular with unit diagonal (M x K)            *
C *     - U is upper triangular (K x N)                               *
C *     - K = min(M,N)                                                *
C *********************************************************************
      SUBROUTINE DMLTRFF( M, N, A, LDA, V, W, WORK, NB, INFO )

C     .. Scalar Arguments ..
      INTEGER            M, N, LDA, NB, INFO

C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), V( * ), W( * ), WORK( * )

C     .. Local Scalars ..
      INTEGER            J, JB, K, M2, N2, IINFO

C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )

C     .. External Subroutines ..
      EXTERNAL           DMLTRF2, DGEMM, DGEMV, DTRSM, DTRSV, DCOPY
      EXTERNAL           XERBLA
      
C     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
C ------------------------------------------------------------
C  Test input parameters
      INFO = 0
      IF ( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF ( LDA .LT. MAX(1,M) ) THEN
         INFO = -4
      ELSE IF ( NB .LT. 1 ) THEN
         INFO = -7
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DMLTRFF', -INFO )
         RETURN
      END IF

      K = MIN( M, N )
      IF ( K .EQ. 0 ) RETURN

      IF ( NB .EQ. 1 .OR. NB .GE. K ) THEN
         CALL DMLTRF2( M, N, A, LDA, V, W, WORK, INFO )
         RETURN
      END IF
C ------------------------------------------------------------
C Blocked loop
      DO 10 J = 1, K, NB

         JB = MIN( NB, K - J + 1 )
         M2 = M - (J + JB) + 1
         N2 = N - (J + JB) + 1
C ------------------------------------------------------------
C Compute w_hat = W(J:J+JB-1) - V(J+JB:M)' * A(J+JB:M, J:J+JB-1)
         CALL DCOPY( JB, W(J), 1, WORK, 1 )
         IF ( M2 .GT. 0 ) THEN
            CALL DGEMV( 'T', M2, JB, -ONE,
     $                  A(J+JB, J), LDA,
     $                  V(J+JB), 1,
     $                  ONE, WORK, 1 )
         END IF
C ------------------------------------------------------------
C Factor diagonal block (LEFT recursive)
         CALL DMLTRF2( JB, JB, A(J,J), LDA,
     $                 V(J), WORK, WORK(JB+1), IINFO )

         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1
C ------------------------------------------------------------
C U12 = L11 \ A12  (L11 is unit lower triangular)
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L', 'L', 'N', 'U',
     $                  JB, N2, ONE,
     $                  A(J,J), LDA,
     $                  A(J,J+JB), LDA )
         END IF
C ------------------------------------------------------------
C L21 = A21 / U11  (U11 is upper triangular)
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R', 'U', 'N', 'N',
     $                  M2, JB, ONE,
     $                  A(J,J), LDA,
     $                  A(J+JB,J), LDA )
         END IF
C ------------------------------------------------------------
C Schur complement: A22 = A22 - L21 * U12
         IF ( M2 .GT. 0 .AND. N2 .GT. 0 ) THEN
            CALL DGEMM( 'N', 'N', M2, N2, JB, -ONE,
     $                  A(J+JB,J), LDA,
     $                  A(J,J+JB), LDA,
     $                  ONE, A(J+JB,J+JB), LDA )
         END IF
C ------------------------------------------------------------
C Update W(J+JB:N) = W(J+JB:N) - (W(J:J+JB-1) / U11) * U12
         IF ( N2 .GT. 0 ) THEN
            CALL DCOPY( JB, W(J), 1, WORK, 1 )
            CALL DTRSV( 'U', 'T', 'N',
     $                  JB, A(J,J), LDA, WORK, 1 )
            CALL DGEMV( 'T', JB, N2, -ONE,
     $                  A(J,J+JB), LDA,
     $                  WORK, 1,
     $                  ONE, W(J+JB), 1 )
         END IF

   10 CONTINUE

      RETURN
      END

C *********************************************************************
C  DMLTRF2 - Recursive GTH factorization for left case
C *********************************************************************
      RECURSIVE SUBROUTINE DMLTRF2(M, N, A, LDA, V, W, WORK, INFO)
      INTEGER            M, N, LDA, INFO
      DOUBLE PRECISION   A(LDA,*), V(*), W(*), WORK(*)

      INTEGER            I, IINFO, K, M1, M2, N1, N2
      DOUBLE PRECISION   ONE, ZERO, PIVOT
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )

      EXTERNAL           DGEMV, DGEMM, DTRSM, DTRSV, DCOPY, XERBLA
      INTRINSIC          MIN, MAX
C ------------------------------------------------------------
C  Test input parameters
      INFO = 0
      IF (M .LT. 0) THEN
         INFO = -1
      ELSE IF (N .LT. 0) THEN
         INFO = -2
      ELSE IF (LDA .LT. MAX(1,M)) THEN
         INFO = -4
      END IF
      IF (INFO .NE. 0) THEN
         CALL XERBLA('DMLTRF2', -INFO)
         RETURN
      END IF

C  Quick return if possible
      K = MIN(M,N)
      IF (K .EQ. 0) RETURN

C Base case
      IF (K .EQ. 1) THEN
         PIVOT = W(1) / V(1)
         A(1,1) = PIVOT
         RETURN
      END IF
C ------------------------------------------------------------
C Recursive case 
      N1 = K/2
      N2 = N - N1
      M1 = N1
      M2 = M - N1
C ------------------------------------------------------------
C Compute w_hat for top-left block
      CALL DCOPY(N1, W, 1, WORK, 1)
      IF (M2 .GT. 0) THEN
         CALL DGEMV('T', M2, N1, -ONE,
     $              A(N1+1,1), LDA,
     $              V(N1+1), 1,
     $              ONE, WORK, 1)
      END IF
C ------------------------------------------------------------
C Recursively factor A11
      CALL DMLTRF2(N1, N1, A, LDA, V, WORK, WORK(N1+1), IINFO)
      IF (INFO .EQ. 0 .AND. IINFO .GT. 0) INFO = IINFO
C ------------------------------------------------------------
C U12 = L11 \ A12
      IF (N2 .GT. 0) THEN
         CALL DTRSM('L', 'L', 'N', 'U',
     $              N1, N2, ONE,
     $              A, LDA,
     $              A(1, N1+1), LDA)
      END IF
C ------------------------------------------------------------
C L21 = A21 / U11
      IF (M2 .GT. 0) THEN
         CALL DTRSM('R', 'U', 'N', 'N',
     $              M2, N1, ONE,
     $              A, LDA,
     $              A(N1+1, 1), LDA)
      END IF
C ------------------------------------------------------------
C Schur complement: A22 = A22 - L21 * U12
      IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
         CALL DGEMM('N', 'N',
     $              M2, N2, N1, -ONE,
     $              A(N1+1, 1), LDA,
     $              A(1, N1+1), LDA,
     $              ONE, A(N1+1, N1+1), LDA)
      END IF
C ------------------------------------------------------------
C Update W 
      CALL DCOPY(N1, W, 1, WORK, 1)
      CALL DTRSV('U', 'T', 'N',
     $           N1, A, LDA, WORK, 1)
      IF (N2 .GT. 0) THEN
         CALL DGEMV('T', N1, N2, -ONE,
     $              A(1, N1+1), LDA,
     $              WORK, 1,
     $              ONE, W(N1+1), 1)
      END IF
C ------------------------------------------------------------
C Recursively factor A22
      IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
         CALL DMLTRF2(M2, N2, A(N1+1, N1+1), LDA,
     $                V(N1+1), W(N1+1), WORK, IINFO)
         IF (INFO .EQ. 0 .AND. IINFO .GT. 0)
     $        INFO = IINFO + N1
      END IF

      RETURN
      END