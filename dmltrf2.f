C  DMLTRF2 - Recursive GTH factorization (LEFT CASE)
C ------------------------------------------------------------
      RECURSIVE SUBROUTINE DMLTRF2(M, N, A, LDA, V, W, INFO)
      INTEGER            M, N, LDA, INFO
      DOUBLE PRECISION   A(LDA,*), V(*), W(*)

      INTEGER            I, IINFO, N2, M2
      DOUBLE PRECISION   ONE
      PARAMETER          (ONE = 1.0D0)

      EXTERNAL           DGEMM, DGEMV, DTRSM

      INFO = 0
C ------------------------------------------------------------
C Quick return
      IF (N .EQ. 0) RETURN
C Base case
      IF (N .EQ. 1) THEN
         W(1) = W(1) / A(1,1)
C L = M / pivot
         IF (M .GT. 1) THEN
            DO I = 2, M
               A(I,1) = A(I,1) / A(1,1)
            END DO
         END IF
         RETURN
      END IF
C ------------------------------------------------------------
C Split
      N2 = N/2
      M2 = M - N2
C ------------------------------------------------------------
C First block
      CALL DMLTRF2(M, N2, A, LDA, V, W, IINFO)

      IF (IINFO .NE. 0) THEN
         INFO = IINFO
         RETURN
      END IF
C ------------------------------------------------------------
C U12 = L11 \ A12
      CALL DTRSM('L','L','N','U',
     $           N2, N-N2, ONE,
     $           A, LDA,
     $           A(1,N2+1), LDA)
C ------------------------------------------------------------
C Update w
      CALL DGEMV('T',
     $           N2, N-N2,
     $           -ONE,
     $           A(1,N2+1), LDA,
     $           W, 1,
     $           ONE,
     $           W(N2+1), 1)
C ------------------------------------------------------------
C Schur complement A22 = A22 - L21 * U12
      IF (M2 .GT. 0) THEN
         CALL DGEMM('N','N',
     $              M2, N-N2, N2,
     $              -ONE,
     $              A(N2+1,1), LDA,
     $              A(1,N2+1), LDA,
     $              ONE,
     $              A(N2+1,N2+1), LDA)
      END IF
C ------------------------------------------------------------
C Second block
      CALL DMLTRF2(M2,
     $             N-N2,
     $             A(N2+1,N2+1),
     $             LDA,
     $             V(N2+1),
     $             W(N2+1),
     $             IINFO)
      IF (IINFO .NE. 0) INFO = IINFO + N2
      RETURN
      END