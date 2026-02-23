C *********************************************************************
C *  DMLTRF (LEFT BLOCKED GTH using DMLTF2 inside)                    *
C *  DMLTRF computes a blocked GTH factorization for a given          *
C *  M-matrix A and vectors V and W such that W^T = V^T * A.          *
C *                                                                   *
C *  Blocked factorization of an M-by-N matrix:                       *
C *                                                                   *
C *        A = L * U                                                  *
C *                                                                   *
C *  where:                                                           *
C *     - L is lower triangular (M x K)                               *
C *     - U is upper triangular with unit diagonal (K x N)            *
C *     - K = min(M,N)                                                *
C *                                                                   *
C *  In this case, U has unit diagonal.                               *
C *********************************************************************
      SUBROUTINE DMLTRF( M, N, A, LDA, V, W, WORK, NB, INFO )

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
      EXTERNAL           DMLTF2
      EXTERNAL           DGEMM, DGEMV, DTRSM, DTRSV
      
C     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX

C ------------------------------------------------------------
C  Test the input parameters.
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
      IF ( INFO .NE. 0 ) RETURN

      K = MIN( M, N )
      IF ( K .EQ. 0 ) RETURN

C ------------------------------------------------------------
C Use unblocked if trivial
      IF ( NB .EQ. 1 .OR. NB .GE. K ) THEN
         CALL DMLTF2( M, N, A, LDA, V, W, INFO )
         RETURN
      END IF
C ------------------------------------------------------------
C Main blocked loop
      DO 10 J = 1, K, NB

         JB = MIN( NB, K - J + 1 )

         M2 = M - (J + JB) + 1
         N2 = N - (J + JB) + 1
C ------------------------------------------------------------
C WORK = W(J:J+JB-1)
         WORK(1:JB) = W(J:J+JB-1)
C ------------------------------------------------------------
C WORK -= V(J+JB:N)^T * A(J+JB:N,J:J+JB-1)
C (row times matrix)
         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'T',
     $           N2, JB,
     $           -ONE,
     $           A( J+JB, J ), LDA,
     $           V( J+JB ), 1,
     $           ONE,
     $           WORK, 1 )
         END IF
C ------------------------------------------------------------
C Factor diagonal block 
         CALL DMLTF2( JB, JB,
     $        A( J, J ), LDA,
     $        V( J ), WORK, IINFO )

         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1
C ------------------------------------------------------------
C U12 = L11^{-1} * A12
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L','L','N','N',
     $           JB, N2,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J, J+JB ), LDA )
         END IF
C ------------------------------------------------------------
C L21 = A21 * U11^{-1} (U11 is unit upper)
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R','U','N','U',
     $           M2, JB,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J+JB, J ), LDA )
         END IF
C ------------------------------------------------------------
C Schur complement
         IF ( M2 .GT. 0 .AND. N2 .GT. 0 ) THEN
            CALL DGEMM( 'N','N',
     $           M2, N2, JB,
     $           -ONE,
     $           A( J+JB, J ), LDA,
     $           A( J, J+JB ), LDA,
     $           ONE,
     $           A( J+JB, J+JB ), LDA )
         END IF
C ------------------------------------------------------------
C W(J:J+JB-1) = W(J:J+JB-1) * U11^{-1}
C (row solve with upper unit)
        CALL DTRSV( 'U','T','U',
     $        JB,
     $        A( J, J ), LDA,
     $        W( J ), 1 )
C ------------------------------------------------------------
C W(J+JB:N) -= W(J:J+JB-1) * U12
         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'T',
     $           JB, N2,
     $           -ONE,
     $           A( J, J+JB ), LDA,
     $           W( J ), 1,
     $           ONE,
     $           W( J+JB ), 1 )
         END IF

 10   CONTINUE

      RETURN
      END