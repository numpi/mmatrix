      SUBROUTINE DMMSVD( N, L, LDL, D, U, LDU,
     &                   S, UMAT, LDUF, VMAT, LDVF,
     &                   WORK, LWORK, IWORK,
     &                   PIV,
     &                   INFO )
C-----------------------------------------------------------------------
C  DMMSVD: SVD for a double-precision matrix 
C
C  This subroutine computes the Singular Value Decomposition (SVD)
C  of a matrix G given in the form:
C
C        G = L * D * U
C
C  where:
C     L - unit lower triangular matrix
C     D - diagonal matrix
C     U - upper triangular matrix
C
C  This routine implements Algorithm 3.1 from:
C
C  J. Demmel, M. Gu, S. Eisenstat, I. Slapnicar, K. Veselic, Z. Drmac,
C  "Computing the singular value decomposition with high relative accuracy,"
C  Linear Algebra and its Applications, 299 (1999), 21–80. 
C
C  Algorithm outline:
C
C     1. Compute QR factorization with pivoting of (L*D): XD = Q * R * P
C     2. Compute W = R * (P^T * U)
C     3. Compute SVD of W using one-sided Jacobi (DGESVJ): W = U_w * Σ * V^T
C     4. Compute U = Q * U_w
C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER            N, LDL, LDU, LDUF, LDVF, LWORK, INFO
      INTEGER            IWORK(*)
      INTEGER            PIV(*)

      DOUBLE PRECISION   L(LDL,*), D(*), U(LDU,*)
      DOUBLE PRECISION   UMAT(LDUF,*), VMAT(LDVF,*), S(*)
      DOUBLE PRECISION   WORK(*)

      INTEGER            I, J, K, INFO2
      INTEGER            JPIV(N)

      DOUBLE PRECISION   XD(N,N), Q(N,N)
      DOUBLE PRECISION   R(N,N), W(N,N), UP(N,N)
      DOUBLE PRECISION   TAU(N)

      EXTERNAL           DGEQP3, DORGQR, DGEMM, DGESVJ, DSWAP

      INFO = 0
C-----------------------------------------------------------------------
C Build permutation matrix from  PIV
C-----------------------------------------------------------------------
      DO J = 1, N
         DO I = 1, N
            Q(I,J) = 0.0D0
         END DO
         Q(J,J) = 1.0D0
      END DO
      DO K = N, 1, -1
         IF (PIV(K) .NE. K) THEN
            CALL DSWAP(N, Q(1,K), 1, Q(1,PIV(K)), 1)
         END IF
      END DO
C-----------------------------------------------------------------------
C Compute XD = (P^T L) * D   
C-----------------------------------------------------------------------
      DO J = 1, N
         DO I = 1, N
            XD(I,J) = 0.0D0
         END DO
      END DO
C---------compute Q^T * L----------------------------------------------- 
      DO J = 1, N
         DO K = 1, N
            DO I = 1, N
               XD(I,J) = XD(I,J) + Q(K,I) * L(K,J)
            END DO
         END DO
      END DO
      DO J = 1, N
         DO I = 1, N
            XD(I,J) = XD(I,J) * D(J)
         END DO
      END DO
C-----------------------------------------------------------------------
C QR with pivoting
C-----------------------------------------------------------------------
      CALL DGEQP3(N, N, XD, N, JPIV, TAU,
     &            WORK, LWORK, INFO2)
      IF (INFO2 .NE. 0) THEN
         INFO = 1
         RETURN
      END IF
C--------- Form Q explicitly--------------------------------------------
      Q = XD
      CALL DORGQR(N, N, N, Q, N, TAU,
     &            WORK, LWORK, INFO2)

      IF (INFO2 .NE. 0) THEN
         INFO = 2
         RETURN
      END IF
C--------- Extract R----------------------------------------------------
      DO I = 1, N
         DO J = 1, N
            R(I,J) = 0.0D0
         END DO
      END DO
      DO I = 1, N
         DO J = I, N
            R(I,J) = XD(I,J)
         END DO
      END DO
C-----------------------------------------------------------------------
C Apply QR permutation
C-----------------------------------------------------------------------
      DO J = 1, N
         DO I = 1, N
            UP(I,J) = U(I, JPIV(J))
         END DO
      END DO
      CALL DGEMM('N','N',
     &           N, N, N,
     &           1.0D0,
     &           R, N,
     &           UP, N,
     &           0.0D0,
     &           W, N)
C-----------------------------------------------------------------------
C  SVD OF W
C-----------------------------------------------------------------------
      CALL DGESVJ('G','U','V',
     &            N, N, W, N,
     &            S,
     &            N,
     &            VMAT, LDVF,
     &            WORK, LWORK,
     &            INFO2)
      IF (INFO2 .NE. 0) THEN
         INFO = 3
         RETURN
      END IF
C-----------------------------------------------------------------------
C Compute U = Q * Uw
C-----------------------------------------------------------------------
      CALL DGEMM('N','N',
     &           N, N, N,
     &           1.0D0,
     &           Q, N,
     &           W, N,
     &           0.0D0,
     &           UMAT, LDUF)

      RETURN
      END
