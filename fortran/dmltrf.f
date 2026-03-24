C ------------------------------------------------------------
C DMLTRF: Blocked GTH factorization using DMLTF2.
C The routine computes a factorization of an M-by-N matrix A
C  A = L * U
C  where L  is unit lower triangular, U  is upper triangular and 
C  V and W are row vectors such that V * A = W
C  --------------------------------------------------------------------
C  The matrix and vectors are partitioned as
C         A = [ A11  A12 ]
C             [ A21  A22 ]
C         V = [ v1 v2 ],  W = [ w1  w2 ]

C The process corresponds to the following steps:
C          w_hat = w1 - v2 * A21
C          [L11, U11] = DMLTF2(A11, v1, w_hat)
C          U12 = L11^{-1} * A12
C          L21 = A21 * U11^{-1}
C          A22 = A22 - L21 * U12
C          w2 = w2 - (w1 / U11) * U12
C  --------------------------------------------------------------------
C   BLAS Level Usage
C       DGEMV   : compute w_hat
C       DTRSM   : triangular solves for U12 and L21
C       DGEMM   : Schur complement update
C       DTRSV   : compute w1 / U11
C  --------------------------------------------------------------------

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
      EXTERNAL           DGEMM, DGEMV, DTRSM, DTRSV, DCOPY

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
C ------------------------------------------------------------
C  Quick return if possible
      K = MIN( M, N )
      IF ( K .EQ. 0 ) RETURN
C ------------------------------------------------------------
C  Use unblocked algorithm if block size is trivial
      IF ( NB .EQ. 1 .OR. NB .GE. K ) THEN
         CALL DMLTF2( M, N, A, LDA, V, W, INFO )
         RETURN
      END IF
C ------------------------------------------------------------
C  Main loop
      DO 10 J = 1, K, NB
         JB = MIN( NB, K - J + 1 )

         M2 = M - (J + JB) + 1
         N2 = N - (J + JB) + 1
C ------------------------------------------------------------
C   Compute w_hat = w1 - v2' * A21
C ------------------------------------------------------------
         CALL DCOPY( JB, W( J ), 1, WORK, 1 )
         IF ( M2 .GT. 0 ) THEN
            CALL DGEMV( 'T',
     $           M2, JB,
     $           -ONE,
     $           A( J+JB, J ), LDA,
     $           V( J+JB ), 1,
     $           ONE,
     $           WORK, 1 )
         END IF
C ------------------------------------------------------------
C  Factor diagonal block using unblocked GTH
C ------------------------------------------------------------
         CALL DMLTF2( JB, JB,
     $        A( J, J ), LDA,
     $        V( J ), WORK, IINFO )

         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1
C ------------------------------------------------------------
C  Compute U12 = L11^{-1} * A12
C ------------------------------------------------------------
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L','L','N','U',
     $           JB, N2,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J, J+JB ), LDA )
         END IF
C ------------------------------------------------------------
C  Compute L21 = A21 * U11^{-1}
C ------------------------------------------------------------
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R','U','N','N',
     $           M2, JB,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J+JB, J ), LDA )
         END IF
C ------------------------------------------------------------
C  Update  A22 = A22 - L21 * U12
C ------------------------------------------------------------
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
C   Update W: w2 = w2 - (w1 / U11) * U12
C ------------------------------------------------------------
         IF ( N2 .GT. 0 ) THEN
            CALL DCOPY( JB, W( J ), 1, WORK, 1 )
            CALL DTRSV( 'U', 'T', 'N',
     $                  JB, A( J, J ), LDA, WORK, 1 )
            CALL DGEMV( 'T',
     $                  JB, N2,
     $                  -ONE,
     $                  A( J, J+JB ), LDA,
     $                  WORK, 1,
     $                  ONE,
     $                  W( J+JB ), 1 )
         END IF
 10   CONTINUE

      RETURN
      END