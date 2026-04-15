
      SUBROUTINE DMMLTF(N, A, LDA, S, PIV, INFO)
C-----------------------------------------------------------------------
C  DMMLTF: Double-precision M-Matrix LDU Factorization
C
C  This subroutine implements LDU factorization from:
C
C  J. Demmel and P. Koev,
C  "Accurate SVDs of weakly diagonally dominant M-matrices,"
C  Numerische Mathematik, 98:99–104, 2004.
C
C  The algorithm computes an LDU factorization of a diagonally dominant
C  M-matrix using only its off-diagonal entries and row sums.
C
C  Input:
C     N     - matrix size
C     A     - matrix storing off-diagonal entries 
C     LDA   - leading dimension of A
C     S     - row sums vector
C
C  Output:
C     A     - contains L (unit lower) and U (upper) factors
C     S     - updated row sums
C     PIV   - pivot indices
C     INFO  - = 0 successful
C-----------------------------------------------------------------------
      INTEGER N, LDA, INFO
      INTEGER PIV(*)
      DOUBLE PRECISION A(LDA,*), S(*)

      INTEGER I, J, K, KP
      DOUBLE PRECISION MAXVAL, TEMP, DIAG

      EXTERNAL DSWAP, DGER

      INFO = 0
C-----------------------------------------------------------------------
C Main loop
      DO K = 1, N

C Construct diagonal entries from row sums
         DO I = K, N
            TEMP = S(I)
            DO J = K, N
               IF (J .NE. I) TEMP = TEMP - A(I,J)
            END DO
            A(I,I) = TEMP
         END DO

C Pivot (largest diagonal entry)
         KP = K
         MAXVAL = ABS(A(K,K))

         DO I = K+1, N
            IF (ABS(A(I,I)) .GT. MAXVAL) THEN
               MAXVAL = ABS(A(I,I))
               KP = I
            END IF
         END DO

         PIV(K) = KP

         IF (KP .NE. K) THEN
            CALL DSWAP(N, A(K,1), LDA, A(KP,1), LDA)
            TEMP = S(K)
            S(K) = S(KP)
            S(KP) = TEMP
         END IF

         DIAG = A(K,K)
         IF (DIAG .EQ. 0.0D0) THEN
            INFO = K
            RETURN
         END IF

         DO I = K+1, N
            A(I,K) = A(I,K) / DIAG
         END DO

         CALL DAXPY(N-K, -S(K), A(K+1,K), 1, S(K+1), 1)


         IF (K .LT. N) THEN
            CALL DGER(N-K, N-K,
     &                -1.0D0,
     &                A(K+1,K), 1,
     &                A(K,K+1), LDA,
     &                A(K+1,K+1), LDA)
         END IF

      END DO

      RETURN
      END