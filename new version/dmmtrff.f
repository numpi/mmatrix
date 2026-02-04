C======================================================================
C     DMMTRFF: Blocked version of the GTH algorithm with DMMTRF2 inside
C======================================================================
      SUBROUTINE DMMTRFF( M, N, A, LDA, V, W, WORK, NB, INFO )
      INTEGER            M, N, LDA, NB, INFO
      DOUBLE PRECISION   A( LDA, * ), V( * ), W( * ), WORK( * )
      
      INTEGER            J, JB, K, M2, N2, IINFO
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      
      EXTERNAL           DMMTRF2, DGEMM, DGEMV, DTRSM, DTRSV
      EXTERNAL           XERBLA
      INTRINSIC          MIN, MAX
      
C     Test input parameters
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
      IF ( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DMMTRFF', -INFO )
         RETURN
      END IF
      
      K = MIN( M, N )
      IF ( K .EQ. 0 ) RETURN
      
      IF ( NB .EQ. 1 .OR. NB .GE. K ) THEN
         CALL DMMTRF2( M, N, A, LDA, V, W, WORK, INFO )
         RETURN
      END IF
      
C     Blocked algorithm
      DO 10 J = 1, K, NB
         JB = MIN( NB, K - J + 1 )
         M2 = M - J - JB + 1
         N2 = N - J - JB + 1
         
C        Compute WORK(1:JB) = w(J:J+JB-1) - A(J:J+JB-1,J+JB:N)*v(J+JB:N)
         DO 20 I = 1, JB
            WORK(I) = W(J+I-1)
   20    CONTINUE
         
         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'N', JB, N2, -ONE,
     $                  A(J, J+JB), LDA,
     $                  V(J+JB), 1,
     $                  ONE, WORK, 1 )
         END IF
         
C        Factor the diagonal block
         CALL DMMTRF2( JB, JB, A(J, J), LDA,
     $                 V(J), WORK, WORK(JB+1), IINFO )
         
         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1
         
C        Update trailing submatrix
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L', 'L', 'N', 'U', JB, N2, ONE,
     $                  A(J, J), LDA,
     $                  A(J, J+JB), LDA )
         END IF
         
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R', 'U', 'N', 'N', M2, JB, ONE,
     $                  A(J, J), LDA,
     $                  A(J+JB, J), LDA )
         END IF
         
         IF ( M2 .GT. 0 .AND. N2 .GT. 0 ) THEN
            CALL DGEMM( 'N', 'N', M2, N2, JB, -ONE,
     $                  A(J+JB, J), LDA,
     $                  A(J, J+JB), LDA,
     $                  ONE, A(J+JB, J+JB), LDA )
         END IF
         
C        Update W  
         CALL DTRSV( 'L', 'N', 'U', JB,
     $               A(J, J), LDA, W(J), 1 )
         
         IF ( M2 .GT. 0 ) THEN
            CALL DGEMV( 'N', M2, JB, -ONE,
     $                  A(J+JB, J), LDA,
     $                  W(J), 1,
     $                  ONE, W(J+JB), 1 )
         END IF
         
   10 CONTINUE
      
      RETURN
      END

C======================================================================
C     DMMTRF2  
C======================================================================
      RECURSIVE SUBROUTINE DMMTRF2(M, N, A, LDA, V, W, WORK, INFO)
      INTEGER            M, N, LDA, INFO
      DOUBLE PRECISION   A(LDA, *), V(*), W(*), WORK(*)
      
      INTEGER            I, IINFO, J, K, M1, M2, N1, N2
      DOUBLE PRECISION   PIVOT, TEMP
      
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          (ONE = 1.0D0, ZERO = 0.0D0)
      
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT, DAXPY, DGER, DGEMM, DGEMV, DTRSM, DTRSV
      EXTERNAL           XERBLA
      
      INTRINSIC          MAX, MIN
      
C     Test input parameters
      INFO = 0
      IF (M .LT. 0) THEN
         INFO = -1
      ELSE IF (N .LT. 0) THEN
         INFO = -2
      ELSE IF (LDA .LT. MAX(1, M)) THEN
         INFO = -4
      END IF
      IF (INFO .NE. 0) THEN
         CALL XERBLA('DMMTRF2', -INFO)
         RETURN
      END IF
      
C     Quick return if possible
      K = MIN(M, N)
      IF (K .EQ. 0) RETURN
      
      IF (M .EQ. 0 .OR. N .EQ. 0) RETURN
      
      IF (M .EQ. 1) THEN
         A(1,1) = W(1) / V(1)
         IF (A(1,1) .EQ. ZERO) INFO = 1
         RETURN
      ELSE IF (N .EQ. 1) THEN
         A(1,1) = W(1) / V(1)
         SFMIN = DLAMCH('S')
         IF (ABS(A(1,1)) .GE. SFMIN) THEN
            CALL DSCAL(M-1, ONE/A(1,1), A(2,1), 1)
         ELSE
            DO I = 1, M-1
               A(1+I, 1) = A(1+I, 1) / A(1, 1)
            END DO
         END IF
         RETURN
      END IF
      
C     Recursive case
      N1 = K / 2
      N2 = N - N1
      M1 = N1
      M2 = M - N1
      
      CALL DCOPY( N1, W, 1, WORK, 1 )

C     Update W(1:N1) = W(1:N1) - A(1:N1, N1+1:N) * V(N1+1:N)
      IF (N2 .GT. 0) THEN
         CALL DGEMV('N', N1, N2, -ONE,
     $              A(1, N1+1), LDA,
     $              V(N1+1), 1,
     $              ONE, W, 1)
      END IF
      
C     Recursively factor A11
      CALL DMMTRF2(N1, N1, A, LDA, V, W, WORK(N1+1), IINFO)
      IF (INFO .EQ. 0 .AND. IINFO .GT. 0) INFO = IINFO
      
      CALL DCOPY( N1, WORK, 1, W, 1 )

C     Solve U12 = L11^{-1} * A12
      IF (N2 .GT. 0) THEN
         CALL DTRSM('L', 'L', 'N', 'U',
     $              N1, N2, ONE,
     $              A, LDA,
     $              A(1, N1+1), LDA)
      END IF
      
C     Solve L21 = A21 * U11^{-1}
      IF (M2 .GT. 0) THEN
         CALL DTRSM('R', 'U', 'N', 'N',
     $              M2, N1, ONE,
     $              A, LDA,
     $              A(N1+1, 1), LDA)
      END IF
      
C     Update Schur complement: A22 = A22 - L21 * U12
      IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
         CALL DGEMM('N', 'N',
     $              M2, N2, N1, -ONE,
     $              A(N1+1, 1), LDA,
     $              A(1, N1+1), LDA,
     $              ONE, A(N1+1, N1+1), LDA)
      END IF
      
C     Update W(1:N1) = L11^{-1} * W(1:N1)
      CALL DTRSV('L', 'N', 'U',
     $           N1, A, LDA, W, 1)
      
C     Update W(N1+1:M) = W(N1+1:M) - L21 * W(1:N1)
      IF (M2 .GT. 0) THEN
         CALL DGEMV('N', M2, N1, -ONE,
     $              A(N1+1, 1), LDA,
     $              W, 1, ONE,
     $              W(N1+1), 1)
      END IF
      
C     Recursively factor A22
      IF (M2 .GT. 0 .AND. N2 .GT. 0) THEN
         CALL DMMTRF2(M2, N2, A(N1+1, N1+1), LDA,
     $                V(N1+1), W(N1+1), WORK, IINFO)
         IF (INFO .EQ. 0 .AND. IINFO .GT. 0) THEN
            INFO = IINFO + N1
         END IF
      END IF
      
      RETURN
      END

