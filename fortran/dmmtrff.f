
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
      
      IF ( NB .GE. K ) THEN
         CALL DMMTRF2( M, N, A, LDA, V, W, WORK, INFO )
         RETURN
      END IF
      
C     Blocked algorithm
      DO 10 J = 1, K, NB
         JB = MIN( NB, K - J + 1 )
         M2 = M - J - JB + 1
         N2 = N - J - JB + 1
         
C        Compute WORK(1:JB) = w(J:J+JB-1) - A(J:J+JB-1,J+JB:N)*v(J+JB:N)
C         DO 20 I = 1, JB
C            WORK(I) = W(J+I-1)
C   20    CONTINUE
          CALL DCOPY(JB, W(J), 1, WORK, 1)
         
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
      
