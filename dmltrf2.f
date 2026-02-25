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
C Recursive case - split at first half
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
C Update W for bottom-right block
C First, compute temp = W(1:N1) / U11
      CALL DCOPY(N1, W, 1, WORK, 1)
      CALL DTRSV('U', 'T', 'N',
     $           N1, A, LDA, WORK, 1)
C Then update W(N1+1:N) -= temp * U12
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