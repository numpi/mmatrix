      SUBROUTINE DMLTRFF( M, N, A, LDA, V, W, WORK, NB, INFO )
C ------------------------------------------------------------
C Blocked GTH factorization using DMLTRF2 
C The routine computes a componentwise accurate LU factorization
C A = L * U, where L is unit lower triangular and U is upper triangular.
C The algorithm assumes the relation v' * A = w' where v > 0, w >= 0
C Block partition of matrix A and vectors v, w: 
C     A = [ A11  A12 ]
C         [ A21  A22 ]
C     v = [ v1  v2 ], w = [ w1  w2 ]
C Steps:
C 1) Save original w1
C 2) w1_hat = w1 - v2*A21
C 3) [L11,U11] = DMLTRF2(A11, v1, w1_hat)
C 4) U12 = L11 \ A12
C 5) L21 = A21 / U11
C 6) A22 = A22 - L21*U12
C 7) w2 = w2 - (w1/U11) * U12
C ------------------------------------------------------------
      INTEGER            M, N, LDA, NB, INFO
      DOUBLE PRECISION   A(LDA,*), V(*), W(*), WORK(*)

      INTEGER            J, JB, K, M2, N2, IINFO
      DOUBLE PRECISION   ONE
      PARAMETER          (ONE = 1.0D0)

      EXTERNAL           DMLTRF2, DGEMM, DGEMV, DTRSM, DCOPY
      INTRINSIC          MIN, MAX
C ------------------------------------------------------------
C Parameter check
C ------------------------------------------------------------
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
C Factorization size
C ------------------------------------------------------------
      K = MIN(M,N)
      IF (K .EQ. 0) RETURN
C ------------------------------------------------------------
C Use recursive routine when blocking is not beneficial
C ------------------------------------------------------------
      IF (NB .EQ. 1 .OR. NB .GE. K) THEN
         CALL DMLTRF2(M,N,A,LDA,V,W,INFO)
         RETURN
      END IF
C ------------------------------------------------------------
C Main loop
C ------------------------------------------------------------
      DO 10 J = 1, K, NB

         JB = MIN(NB, K-J+1)
         M2 = M - J - JB + 1
         N2 = N - J - JB + 1
C ------------------------------------------------------------
C Step 1: Save original w1: WORK = w(I)
C ------------------------------------------------------------
         CALL DCOPY(JB, W(J), 1, WORK, 1)
C ------------------------------------------------------------
C Step 2: w1_hat = w1 - v2*A21
C ------------------------------------------------------------
         IF (M2 .GT. 0) THEN
            CALL DGEMV('T',M2,JB,-ONE,
     $                 A(J+JB,J),LDA,
     $                 V(J+JB),1,
     $                 ONE,
     $                 W(J),1)
         END IF
C ------------------------------------------------------------
C Step 3: Factor diagonal block [L11,U11] = DMLTRF2(A11,v1,w1_hat)
C ------------------------------------------------------------
         CALL DMLTRF2(JB,JB,
     $                A(J,J),LDA,
     $                V(J),
     $                W(J),
     $                IINFO)

         IF (INFO .EQ. 0 .AND. IINFO .GT. 0)
     $        INFO = IINFO + J - 1
C ------------------------------------------------------------
C Step 4: U12 = L11^{-1} * A12
C ------------------------------------------------------------
         IF (N2 .GT. 0) THEN
            CALL DTRSM('L','L','N','U',
     $                 JB,N2,ONE,
     $                 A(J,J),LDA,
     $                 A(J,J+JB),LDA)
         END IF
C ------------------------------------------------------------
C Step 5: L21 = A21 * U11^{-1}
C ------------------------------------------------------------
         IF (M2 .GT. 0) THEN
            CALL DTRSM('R','U','N','N',
     $                 M2,JB,ONE,
     $                 A(J,J),LDA,
     $                 A(J+JB,J),LDA)
         END IF
C ------------------------------------------------------------
C Step 6: Schur complement A22 = A22 - L21*U12
C ------------------------------------------------------------
         IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
            CALL DGEMM('N','N',
     $                 M2,N2,JB,-ONE,
     $                 A(J+JB,J),LDA,
     $                 A(J,J+JB),LDA,
     $                 ONE,
     $                 A(J+JB,J+JB),LDA)
         END IF
C ------------------------------------------------------------
C Step 7: w2 = w2 - (w1/U11)*U12
C ------------------------------------------------------------
C Compute w1/U11 
         CALL DTRSM('R','U','N','N',
     $              1,JB,ONE,
     $              A(J,J),LDA,
     $              WORK,1)

C Update w2
         IF (N2 .GT. 0) THEN
            CALL DGEMV('T',
     $                 JB,N2,-ONE,
     $                 A(J,J+JB),LDA,
     $                 WORK,1,
     $                 ONE,
     $                 W(J+JB),1)
         END IF

 10   CONTINUE

      RETURN
      END
