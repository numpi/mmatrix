      SUBROUTINE DMLTRFF( M, N, A, LDA, V, W, NB, INFO )
C ------------------------------------------------------------
C Blocked GTH factorization using DMLTRF2.
C
C v*A = w  (v,w row vectors)
C
C Partition:
C
C     A = [ A11 A12 ]
C         [ A21 A22 ]
C
C     v = [ v1 v2 ],  w = [ w1 w2 ]
C
C Steps:
C 1) w1 = w1 - v2*A21 (It seems this computation is not needed for the first block, so lines 64-69 can be removed.)
C 2) [L11, U11, w_out] = DMLTRF2(A11, v1, w1)
C 3) U12 = L11 \ A12
C 4) L21 = A21 / U11
C 5) A22 = A22 - L21*U12
C 6) w2 = w2 - w_out * U12
C ------------------------------------------------------------
      INTEGER            M, N, LDA, NB, INFO
      DOUBLE PRECISION   A(LDA,*), V(*), W(*)

      INTEGER            J, JB, K, M2, N2, IINFO
      DOUBLE PRECISION   ONE
      PARAMETER          (ONE = 1.0D0)

      EXTERNAL           DMLTRF2, DGEMM, DGEMV, DTRSM
      INTRINSIC          MIN, MAX
C ------------------------------------------------------------
C Test input parameters
      INFO = 0
      IF ( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF ( LDA .LT. MAX(1,M) ) THEN
         INFO = -4
      ELSE IF ( NB .LT. 1 ) THEN
         INFO = -6
      END IF

      IF ( INFO .NE. 0 ) RETURN
C ------------------------------------------------------------
C Factorization size
      K = MIN( M, N )
      IF ( K .EQ. 0 ) RETURN
C ------------------------------------------------------------
C Use unblocked code if NB small
      IF ( NB .EQ. 1 .OR. NB .GE. K ) THEN
         CALL DMLTRF2( M, N, A, LDA, V, W, INFO )
         RETURN
      END IF
C ------------------------------------------------------------
C Main blocked loop
      DO 10 J = 1, K, NB
         JB = MIN( NB, K-J+1 )
         M2 = M - J - JB + 1
         N2 = N - J - JB + 1
C ------------------------------------------------------------
C Compute w1 = w1 - v2* A21
C ------------------------------------------------------------
C         IF ( M2 .GT. 0 ) THEN
C            CALL DGEMV( 'T', M2, JB, -ONE,
C     $                  A(J+JB,J), LDA,
C     $                  V(J+JB), 1,
C     $                  ONE, W(J), 1 )
C         END IF
C ------------------------------------------------------------
C Factor diagonal block
C ------------------------------------------------------------
         CALL DMLTRF2( JB, JB, A(J,J), LDA,
     $                 V(J), W(J), IINFO )

         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1
C ------------------------------------------------------------
C Compute U12 = L11^{-1} * A12
C ------------------------------------------------------------
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L','L','N','U',
     $                  JB, N2, ONE,
     $                  A(J,J), LDA,
     $                  A(J,J+JB), LDA )
         END IF
C ------------------------------------------------------------
C Compute L21 = A21 * U11^{-1}
C ------------------------------------------------------------
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R','U','N','N',
     $                  M2, JB, ONE,
     $                  A(J,J), LDA,
     $                  A(J+JB,J), LDA )
         END IF
C ------------------------------------------------------------
C Update Schur complement
C ------------------------------------------------------------
         IF ( M2 .GT. 0 .AND. N2 .GT. 0 ) THEN
            CALL DGEMM( 'N','N',
     $                  M2, N2, JB, -ONE,
     $                  A(J+JB,J), LDA,
     $                  A(J,J+JB), LDA,
     $                  ONE,
     $                  A(J+JB,J+JB), LDA )
         END IF
C ------------------------------------------------------------
C Update w2 using w_out returned by DMLTRF2
C ------------------------------------------------------------
         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'T',
     $                  JB, N2, -ONE,
     $                  A(J,J+JB), LDA,
     $                  W(J), 1,
     $                  ONE,
     $                  W(J+JB), 1 )
         END IF

 10   CONTINUE

      RETURN
      END