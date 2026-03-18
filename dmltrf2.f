      RECURSIVE SUBROUTINE DMLTRF2(M, N, A, LDA, V, W, INFO)
C ------------------------------------------------------------
C  DMLTRF2 - Recursive GTH (LEFT CASE) 
C
C  Computes a componentwise accurate LU factorization of the
C  M-by-N block of matrix A using the GTH trick.
C  The algorithm assumes the relation v' * A = w' where v > 0, w >= 0
C  The factorization produced is
C        A = L * U
C  where L is unit lower triangular and U is upper triangular.
C ------------------------------------------------------------
C  Recursive block partition
C      A = [ A11  A12 ]
C          [ A21  A22 ]
C      v = [ v1  v2 ]
C      w = [ w1  w2 ]
C  Steps
C  1) [L11,U11] = DMLTRF2(A(:,1:n1), v, w1)
C  2) Compute  U12 = L11^{-1} * A12
C  3) Update w:  w2 = w2 - w1 * U12
C  4) Update Schur complement A22 = A22 - L21 * U12
C  5) Recursively factor the submatrix A22: DMLTRF2(A22, v2, w2)
C ------------------------------------------------------------

      INTEGER M, N, LDA, INFO
      DOUBLE PRECISION A(LDA,*), V(*), W(*)

      INTEGER N1, N2, M2, IINFO, I

      DOUBLE PRECISION ONE, ZERO, PIV, TEMP, SFMIN
      PARAMETER (ONE = 1.0D0, ZERO = 0.0D0)

      DOUBLE PRECISION DDOT, DLAMCH
      EXTERNAL DDOT, DLAMCH
      EXTERNAL DSCAL, DGEMM, DGEMV, DTRSM

      INFO = 0
C ------------------------------------------------------------
C Quick return if matrix has zero columns
C ------------------------------------------------------------
      IF (N .EQ. 0) RETURN
C ------------------------------------------------------------
C Base case: N=1
C ------------------------------------------------------------
      IF (N .EQ. 1) THEN
         IF (M .GT. 1) THEN
            TEMP = DDOT(M-1, A(2,1), 1, V(2), 1)
         ELSE
            TEMP = ZERO
         END IF
         PIV = (W(1) - TEMP) / V(1)
         A(1,1) = PIV
         IF (M .GT. 1) THEN
            SFMIN = DLAMCH('S')
            IF (ABS(PIV) .GE. SFMIN) THEN
               CALL DSCAL(M-1, ONE/PIV, A(2,1), 1)
            ELSE
               DO I = 2, M
                  A(I,1) = A(I,1) / PIV
               END DO
            END IF
         END IF
         W(1) = W(1) / PIV
         RETURN
      END IF
C ------------------------------------------------------------
C Split matrix
C ------------------------------------------------------------
      N1 = N/2
      N2 = N - N1
      M2 = M - N1
C ------------------------------------------------------------
C Step 1:  [L11,U11] = DMLTRF2(A(:,1:n1), v, w1)
C ------------------------------------------------------------
      CALL DMLTRF2(M, N1, A, LDA, V, W, IINFO)

      IF (IINFO .NE. 0) THEN
         INFO = IINFO
         RETURN
      END IF
C ------------------------------------------------------------
C Step 2: Compute U12 = L11^{-1} * A12
C ------------------------------------------------------------
      CALL DTRSM('L','L','N','U',
     $           N1, N2, ONE,
     $           A, LDA,
     $           A(1,N1+1), LDA)
C ------------------------------------------------------------
C Step 3: Update  w2 = w2 - w1 * U12
C ------------------------------------------------------------
      CALL DGEMV('T',
     $           N1, N2,
     $           -ONE,
     $           A(1,N1+1), LDA,
     $           W, 1,
     $           ONE,
     $           W(N1+1), 1)
C ------------------------------------------------------------
C Step 4: Update Schur complement A22 = A22 - L21 * U12
C ------------------------------------------------------------
      IF (M2 .GT. 0) THEN
         CALL DGEMM('N','N',
     $              M2, N2, N1,
     $              -ONE,
     $              A(N1+1,1), LDA,
     $              A(1,N1+1), LDA,
     $              ONE,
     $              A(N1+1,N1+1), LDA)
      END IF
C -----------------------------------------------------------
C Step 5: Factor trailing block recursively
C ------------------------------------------------------------
      CALL DMLTRF2(M2, N2,
     $             A(N1+1,N1+1), LDA,
     $             V(N1+1),
     $             W(N1+1),
     $             IINFO)

      IF (IINFO .NE. 0) INFO = IINFO + N1

      RETURN
      END
