       SUBROUTINE DMMLTF(N, A, LDA, U, V, PIV, INFO)
C-----------------------------------------------------------------------
C  DMMLTF: Double-precision M-Matrix LDU Factorization 
C
C  Input:
C     N     - matrix size
C     A     - matrix storing off-diagonal entries 
C     LDA   - leading dimension of A
C     U     -  vector u
C     V     - right-hand side vector v
C
C  Output:
C     A     - contains L (unit lower) and U (upper) factors
C     U - permuted vector u
C     V - updated vector v
C     PIV   - pivot indices
C     INFO  - = 0 successful
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N, LDA, INFO
      INTEGER PIV(*)
      DOUBLE PRECISION A(LDA,*), U(*), V(*)

      INTEGER I, J, K, KP
      DOUBLE PRECISION MAXVAL, TEMP, DIAG, S

      EXTERNAL DSWAP, DGER, DAXPY

      INFO = 0
C-----------------------------------------------------------------------
C Main loop
      DO K = 1, N

C Construct diagonal entries from u and v
         DO I = K, N
            S = 0.0D0
            DO J = K, N
               IF (J .NE. I) S = S + A(I,J) * U(J)
            END DO
            A(I,I) = (V(I) - S) / U(I)
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
            CALL DSWAP(N, A(1,K), 1, A(1,KP), 1)
C           Swap vectors
            TEMP = U(K)
            U(K) = U(KP)
            U(KP) = TEMP
            
            TEMP = V(K)
            V(K) = V(KP)
            V(KP) = TEMP
         END IF

         DIAG = A(K,K)
         IF (DIAG .EQ. 0.0D0) THEN
            INFO = K
            RETURN
         END IF

         DO I = K+1, N
            A(I,K) = A(I,K) / DIAG
         END DO

C        Update V
         CALL DAXPY(N-K, -V(K), A(K+1,K), 1, V(K+1), 1)

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
