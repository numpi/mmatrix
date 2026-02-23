C *********************************************************************
C *  DMLTF2 computes an unblocked GTH factorization for a given       *
C *  M-matrix A and vectors V and W such that W^T = V^T * A.          *
C *                                                                   *
C *  Blocked factorization of an M-by-N matrix:                       *
C *                                                                   *
C *        A = L * U                                                  *
C *                                                                   *
C *  where:                                                           *
C *     - L is lower triangular (M x K)                               *
C *     - U is upper triangular with unit diagonal (K x N)            *
C *     - K = min(M,N)                                                *
C *                                                                   *
C *  In this case, U has unit diagonal.                               *
C *********************************************************************
      SUBROUTINE DMLTF2( M, N, A, LDA, V, W, INFO )

C     .. Scalar Arguments ..
      INTEGER            M, N, LDA, INFO

C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), V( * ), W( * )

C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )

C     .. Local Scalars ..
      INTEGER            J, I, K, MJ, NJ
      DOUBLE PRECISION   PIVOT, TEMP

C     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT

C     .. External Subroutines ..
      EXTERNAL           DGER, DAXPY, XERBLA

C     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
C  --------------------------------------------------------------------
C     Test the input parameters.
      INFO = 0
      IF ( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF ( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF ( LDA .LT. MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DMLTF2', -INFO )
         RETURN
      END IF
C ------------------------------------------------------------
C     Quick return if possible
      IF ( M .EQ. 0 .OR. N .EQ. 0 ) RETURN

      K = MIN( M, N )
C ------------------------------------------------------------
C     Main GTH loop  
      DO 10 J = 1, K
C        Compute pivot:
C        L(j,j) = ( w(j) - V(j+1:M)^T * A(j+1:M,j) ) / V(j)

         MJ = M - J
         IF ( MJ .GT. 0 ) THEN
            TEMP = DDOT( MJ, V( J+1 ), 1, A( J+1, J ), 1 )
            PIVOT = ( W(J) - TEMP ) / V(J)
         ELSE
            PIVOT = W(J) / V(J)
         END IF

         A( J, J ) = PIVOT

         IF ( PIVOT .EQ. ZERO .AND. INFO .EQ. 0 ) THEN
            INFO = J
         END IF
C ------------------------------------------------------------
C        Form U(j,j+1:N) = A(j,j+1:N) / pivot

         NJ = N - J
         IF ( NJ .GT. 0 ) THEN
            DO 20 I = 1, NJ
               A( J, J+I ) = A( J, J+I ) / PIVOT
   20       CONTINUE
         END IF
C ------------------------------------------------------------
C        Update W(j+1:N) = W(j+1:N) - W(j) * U(j,j+1:N)

         IF ( NJ .GT. 0 ) THEN
            CALL DAXPY( NJ, -W(J), A( J, J+1 ), LDA, W( J+1 ), 1 )
         END IF
C ------------------------------------------------------------
C Update trailing submatrix: A(j+1:M, j+1:N) -= A(j+1:M,j) * U(j,j+1:N)

         IF ( MJ .GT. 0 .AND. NJ .GT. 0 ) THEN
            CALL DGER( MJ, NJ, -ONE,
     $                 A( J+1, J ), 1,
     $                 A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF

   10 CONTINUE

      RETURN

C     End of DMLTF2

      END