      SUBROUTINE DMLTF2( M, N, A, LDA, V, W, INFO )
*
*  DMLTF2 computes an unblocked GTH factorization of a general
*  M-by-N M-matrix A using vectors V and W where V'*A = W'.
*
*  The factorization has the form
*
*     A = L * U
*
*  where:
*     L is M-by-K  
*     U is K-by-N  
*     K = MIN(M,N)
*
*  --------------------------------------------------------------------
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA, INFO
*
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), V( * ), W( * )
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
*
*     .. Local Scalars ..
      INTEGER            J, I, K, MJ, NJ
      DOUBLE PRECISION   PIVOT, TEMP
*
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*
*     .. External Subroutines ..
      EXTERNAL           DGER, DAXPY, XERBLA
*
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*
*  --------------------------------------------------------------------
*
*     Test the input parameters.
*
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
*
*     Quick return if possible
*
      IF ( M .EQ. 0 .OR. N .EQ. 0 ) RETURN
*
*     Compute number of steps
*
      K = MIN( M, N )
*
*     Main GTH loop 
*
      DO 10 J = 1, K
*
*        Compute pivot 
*
         MJ = M - J
         IF ( MJ .GT. 0 ) THEN
            TEMP = DDOT( MJ, V( J+1 ), 1, A( J+1, J ), 1 )
            PIVOT = ( W(J) - TEMP ) / V(J)
         ELSE
            PIVOT = W(J) / V(J)
         END IF
*
         A( J, J ) = PIVOT
*
*        Test for singularity
*
         IF ( PIVOT .EQ. ZERO .AND. INFO .EQ. 0 ) THEN
            INFO = J
         END IF
*
*        Compute L(j+1:M, j)
*      
         MJ = M - J
         IF ( MJ .GT. 0 ) THEN
            DO 20 I = 1, MJ
               A( J+I, J ) = A( J+I, J ) / PIVOT
   20       CONTINUE
         END IF
*
*        Update W(j+1:N) = W(j+1:N) - (W(j)/PIVOT) * U(j,j+1:N)
*
         NJ = N - J
         IF ( NJ .GT. 0 ) THEN
            CALL DAXPY( NJ, -W(J)/PIVOT, A( J, J+1 ), LDA,
     $                  W( J+1 ), 1 )
         END IF
*
*        Update  A(j+1:M, j+1:N) -= L21 * U12
*
         IF ( MJ .GT. 0 .AND. NJ .GT. 0 ) THEN
            CALL DGER( MJ, NJ, -ONE,
     $                 A( J+1, J ), 1,
     $                 A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
*
   10 CONTINUE
*
      DO 30 J = 1, K
          W(J) = W(J) / A(J,J)
   30     CONTINUE
      RETURN
*
*     End of DMLTF2
*
      END
