      SUBROUTINE DMMTRF( M, N, A, LDA, V, W, WORK, NB, INFO )
C
C  Blocked GTH factorization for a general M-by-N M-matrix.
C
C     A = L * U
C
C  where:
C     L is M-by-K unit lower trapezoidal
C     U is K-by-N upper trapezoidal
C     K = min(M,N)
C
      INTEGER            M, N, LDA, NB, INFO
      DOUBLE PRECISION   A( LDA, * ), V( * ), W( * ), WORK( * )

      INTEGER            J, JB, K, M2, N2, IINFO
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )

      EXTERNAL           DMMTF2
      EXTERNAL           DGEMM, DGEMV, DTRSM, DTRSV
      INTRINSIC          MIN, MAX

C------------------------------------------------------------
C Argument check
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
      IF ( INFO .NE. 0 ) RETURN

      K = MIN( M, N )
      IF ( K .EQ. 0 ) RETURN

C------------------------------------------------------------
C Use unblocked code if NB is small
      IF ( NB .EQ. 1 .OR. NB .GE. K ) THEN
         CALL DMMTF2( M, N, A, LDA, V, W, INFO )
         RETURN
      END IF

C------------------------------------------------------------
C Main blocked loop (over pivot blocks)
      DO 10 J = 1, K, NB

         JB = MIN( NB, K - J + 1 )

C        Remaining rows and columns
         M2 = M - (J + JB) + 1
         N2 = N - (J + JB) + 1

C------------------------------------------------------------
C Compute WORK = w(J:J+JB-1) - A(J:J+JB-1,J+JB:N)*v(J+JB:N)
         WORK(1:JB) = W(J:J+JB-1)

         IF ( N2 .GT. 0 ) THEN
            CALL DGEMV( 'N',
     $           JB, N2,
     $           -ONE,
     $           A( J, J+JB ), LDA,
     $           V( J+JB ), 1,
     $           ONE,
     $           WORK, 1 )
         END IF

C------------------------------------------------------------
C Factor diagonal block using unblocked GTH
         CALL DMMTF2( JB, JB,
     $        A( J, J ), LDA,
     $        V( J ), WORK, IINFO )

         IF ( INFO .EQ. 0 .AND. IINFO .GT. 0 )
     $        INFO = IINFO + J - 1

C------------------------------------------------------------
C U12 = L11^{-1} * A12
         IF ( N2 .GT. 0 ) THEN
            CALL DTRSM( 'L','L','N','U',
     $           JB, N2,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J, J+JB ), LDA )
         END IF

C------------------------------------------------------------
C L21 = A21 * U11^{-1}
         IF ( M2 .GT. 0 ) THEN
            CALL DTRSM( 'R','U','N','N',
     $           M2, JB,
     $           ONE,
     $           A( J, J ), LDA,
     $           A( J+JB, J ), LDA )
         END IF

C------------------------------------------------------------
C Schur complement: A22 -= L21 * U12
         IF ( M2 .GT. 0 .AND. N2 .GT. 0 ) THEN
            CALL DGEMM( 'N','N',
     $           M2, N2, JB,
     $           -ONE,
     $           A( J+JB, J ), LDA,
     $           A( J, J+JB ), LDA,
     $           ONE,
     $           A( J+JB, J+JB ), LDA )
         END IF

C------------------------------------------------------------
C Update W(J:J+JB-1) = L11^{-1} * W
         CALL DTRSV( 'L','N','U',
     $        JB,
     $        A( J, J ), LDA,
     $        W( J ), 1 )

C------------------------------------------------------------
C Update W(J+JB:M) -= L21 * W(J:J+JB-1)
         IF ( M2 .GT. 0 ) THEN
            CALL DGEMV( 'N',
     $           M2, JB,
     $           -ONE,
     $           A( J+JB, J ), LDA,
     $           W( J ), 1,
     $           ONE,
     $           W( J+JB ), 1 )
         END IF

 10   CONTINUE

      RETURN
      END
