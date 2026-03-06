      RECURSIVE SUBROUTINE DMLTRF2(M, N, A, LDA, V, W, INFO)
C  DMLTRF2 - Recursive GTH factorization (LEFT CASE)
C ------------------------------------------------------------    
      INTEGER M, N, LDA, INFO
      DOUBLE PRECISION A(LDA,*), V(*), W(*)

      INTEGER N1, N2, M2, IINFO, I
      DOUBLE PRECISION ONE,PIV
      PARAMETER (ONE=1.0D0)

      EXTERNAL DGEMM, DGEMV, DTRSM

      INFO = 0
C ------------------------------------------------------------
C Quick return
      IF (N .EQ. 0) RETURN
C Base case
      IF (N .EQ. 1) THEN
         PIV = W(1) / V(1)
         A(1,1) = PIV

         IF (M .GT. 1) THEN
            DO I = 2,M
               A(I,1) = A(I,1) / PIV
            END DO
         END IF

         RETURN
      END IF
C -------------------------
C Split
      N1 = N/2
      N2 = N-N1
      M2 = M-N1
C -------------------------
C First block
      CALL DMLTRF2(M, N1, A, LDA, V, W, IINFO)
      IF(IINFO.NE.0) THEN
         INFO=IINFO
         RETURN
      END IF
C -------------------------
C U12 = L11 \ A12
      CALL DTRSM('L','L','N','U',
     $           N1,N2,ONE,
     $           A,LDA,
     $           A(1,N1+1),LDA)
C -------------------------
C Update W
      CALL DGEMV('T',
     $           N1,N2,
     $           -ONE,
     $           A(1,N1+1),LDA,
     $           W,1,
     $           ONE,
     $           W(N1+1),1)
C -------------------------
C Schur complement
      IF(M2.GT.0) THEN
         CALL DGEMM('N','N',
     $              M2,N2,N1,
     $              -ONE,
     $              A(N1+1,1),LDA,
     $              A(1,N1+1),LDA,
     $              ONE,
     $              A(N1+1,N1+1),LDA)
      END IF
C -------------------------
C Second block
      CALL DMLTRF2(M2,N2,
     $             A(N1+1,N1+1),LDA,
     $             V(N1+1),
     $             W(N1+1),
     $             IINFO)

      IF(IINFO.NE.0) INFO = IINFO + N1

      RETURN
      END