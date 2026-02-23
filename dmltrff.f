C *********************************************************************
C *  DMLTRFF (LEFT BLOCKED GTH using DMLTRF2 inside)                  *
C *  DMLTRFF computes a recursive blocked GTH factorization for a     *
C *  given M-matrix A and vectors V and W satisfying W^T = V^T * A.   *
C *                                                                   *
C *  Recursive blocked factorization of an M-by-N matrix:             *
C *                                                                   *
C *        A = L * U                                                  *
C *                                                                   *
C *  where:                                                           *
C *     - L is lower triangular (M x K)                               *
C *     - U is upper triangular with unit diagonal (K x N)            *
C *     - K = min(M,N)                                                *
C *  In this case, U has unit diagonal.                               *
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
      EXTERNAL           DMLTRF2, DGEMM, DGEMV, DTRSM, DTRSV
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
         M2 = M - J - JB + 1
         N2 = N - J - JB + 1
C ------------------------------------------------------------
C WORK = W(J:J+JB-1)
         CALL DCOPY( JB, W(J), 1, WORK, 1 )
C ------------------------------------------------------------
C WORK -= V(J+JB:N)^T * A(J+JB:N, J:J+JB-1)
         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'T', N2, JB, -ONE,
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
C U12 = L11^{-1} * A12
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L', 'L', 'N', 'N',
     $                  JB, N2, ONE,
     $                  A(J,J), LDA,
     $                  A(J,J+JB), LDA )
         END IF
C ------------------------------------------------------------
C L21 = A21 * U11^{-1}  (U unit upper)
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R', 'U', 'N', 'U',
     $                  M2, JB, ONE,
     $                  A(J,J), LDA,
     $                  A(J+JB,J), LDA )
         END IF
C ------------------------------------------------------------
C Schur complement
         IF ( M2 .GT. 0 .AND. N2 .GT. 0 ) THEN
            CALL DGEMM( 'N', 'N', M2, N2, JB, -ONE,
     $                  A(J+JB,J), LDA,
     $                  A(J,J+JB), LDA,
     $                  ONE, A(J+JB,J+JB), LDA )
         END IF
C ------------------------------------------------------------
C W(J:J+JB-1) = W(J:J+JB-1) * U11^{-1}
         CALL DTRSV( 'U', 'T', 'U',
     $               JB, A(J,J), LDA, W(J), 1 )
C ------------------------------------------------------------
C W(J+JB:N) -= W(J:J+JB-1) * U12
         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'T', JB, N2, -ONE,
     $                  A(J,J+JB), LDA,
     $                  W(J), 1,
     $                  ONE, W(J+JB), 1 )
         END IF

   10 CONTINUE

      RETURN
      END

C *********************************************************************
C  DMLTRF2  
C *********************************************************************
      RECURSIVE SUBROUTINE DMLTRF2(M, N, A, LDA, V, W, WORK, INFO)
      INTEGER            M, N, LDA, INFO
      DOUBLE PRECISION   A(LDA,*), V(*), W(*), WORK(*)

      INTEGER            I, IINFO, K, M1, M2, N1, N2
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )

      EXTERNAL           DGEMV, DGEMM, DTRSM, DTRSV, XERBLA
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
         A(1,1) = W(1)/V(1)
         RETURN
      END IF
C ------------------------------------------------------------
C Recursive case
      N1 = K/2
      N2 = N - N1
      M1 = N1
      M2 = M - N1
C ------------------------------------------------------------
      CALL DCOPY(N1, W, 1, WORK, 1)
C ------------------------------------------------------------
C W(1:N1) -= V(N1+1:N)^T * A(N1+1:N,1:N1)
      IF (N2 .GT. 0) THEN
         CALL DGEMV('T', N2, N1, -ONE,
     $              A(N1+1,1), LDA,
     $              V(N1+1),1,
     $              ONE, W,1)
      END IF
C ------------------------------------------------------------
C Recursively factor A11
      CALL DMLTRF2(N1,N1,A,LDA,V,W,WORK(N1+1),IINFO)
      IF (INFO .EQ. 0 .AND. IINFO .GT. 0) INFO = IINFO
C ------------------------------------------------------------
      CALL DCOPY(N1, WORK,1,W,1)
C ------------------------------------------------------------
C U12 = L11^{-1} * A12
      IF (N2 .GT. 0) THEN
         CALL DTRSM('L','L','N','N',
     $              N1,N2,ONE,
     $              A,LDA,
     $              A(1,N1+1),LDA)
      END IF
C ------------------------------------------------------------
C L21 = A21 * U11^{-1}
      IF (M2 .GT. 0) THEN
         CALL DTRSM('R','U','N','U',
     $              M2,N1,ONE,
     $              A,LDA,
     $              A(N1+1,1),LDA)
      END IF
C ------------------------------------------------------------
C Schur complement
      IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
         CALL DGEMM('N','N',
     $              M2,N2,N1,-ONE,
     $              A(N1+1,1),LDA,
     $              A(1,N1+1),LDA,
     $              ONE,A(N1+1,N1+1),LDA)
      END IF
C ------------------------------------------------------------
C W(1:N1) = W(1:N1) * U11^{-1}
      CALL DTRSV('U','T','U',
     $           N1,A,LDA,W,1)
C ------------------------------------------------------------
C W(N1+1:N) -= W(1:N1) * U12
      IF (N2 .GT. 0) THEN
         CALL DGEMV('T',N1,N2,-ONE,
     $              A(1,N1+1),LDA,
     $              W,1,
     $              ONE,W(N1+1),1)
      END IF
C ------------------------------------------------------------
C Recursively factor A22
      IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
         CALL DMLTRF2(M2,N2,A(N1+1,N1+1),LDA,
     $                V(N1+1),W(N1+1),WORK,IINFO)
         IF (INFO .EQ. 0 .AND. IINFO .GT. 0)
     $        INFO = IINFO + N1
      END IF

      RETURN
      END