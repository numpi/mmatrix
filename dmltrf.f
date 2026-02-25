C *********************************************************************
C *  DMLTRF (LEFT BLOCKED GTH using DMLTF2 inside)                    *
C *  DMLTRF computes a blocked GTH factorization for a given          *
C *  M-matrix A and vectors V and W such that V'*A = W'.              *
C *                                                                   *
C *  Blocked factorization of an M-by-N matrix:                       *
C *                                                                   *
C *        A = L * U                                                  *
C *                                                                   *
C *  where:                                                           *
C *     - L is lower triangular with unit diagonal (M x K)            *
C *     - U is upper triangular (K x N)                               *
C *     - K = min(M,N)                                                *
C *                                                                   *
C *  In this case, L has unit diagonal.                               *
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
      EXTERNAL           DGEMM, DGEMV, DTRSM, DCOPY
      
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
C Compute w_hat = W(J:J+JB-1) - V(J+JB:M)' * A(J+JB:M, J:J+JB-1)
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
C Factor diagonal block using left GTH
         CALL DMLTF2( JB, JB,
     $        A( J, J ), LDA,
     $        V( J ), WORK, IINFO )

         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1
C ------------------------------------------------------------
C U12 = L11 \ A12  (L11 is unit lower triangular)
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L','L','N','U',
     $           JB, N2,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J, J+JB ), LDA )
         END IF
C ------------------------------------------------------------
C L21 = A21 / U11  (U11 is upper triangular)
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R','U','N','N',
     $           M2, JB,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J+JB, J ), LDA )
         END IF
C ------------------------------------------------------------
C Schur complement: A22 = A22 - L21 * U12
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
C Update W(J+JB:N) = W(J+JB:N) - (W(J:J+JB-1) / U_II) * U12
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