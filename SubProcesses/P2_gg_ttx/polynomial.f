      MODULE ML5_2_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INCLUDE 'coef_specs.inc'
      INCLUDE 'loop_max_coefs.inc'

C     Map associating a rank to each coefficient position
      INTEGER COEFTORANK_MAP(0:LOOPMAXCOEFS-1)
      DATA COEFTORANK_MAP(0:0)/1*0/
      DATA COEFTORANK_MAP(1:4)/4*1/
      DATA COEFTORANK_MAP(5:14)/10*2/
      DATA COEFTORANK_MAP(15:34)/20*3/

C     Map defining the number of coefficients for a symmetric tensor
C      of a given rank
      INTEGER NCOEF_R(0:3)
      DATA NCOEF_R/1,5,15,35/

C     Map defining the coef position resulting from the multiplication
C      of two lower rank coefs.
      INTEGER COMB_COEF_POS(0:LOOPMAXCOEFS-1,0:4)
      DATA COMB_COEF_POS(  0,  0:  4) /  0,  1,  2,  3,  4/
      DATA COMB_COEF_POS(  1,  0:  4) /  1,  5,  6,  8, 11/
      DATA COMB_COEF_POS(  2,  0:  4) /  2,  6,  7,  9, 12/
      DATA COMB_COEF_POS(  3,  0:  4) /  3,  8,  9, 10, 13/
      DATA COMB_COEF_POS(  4,  0:  4) /  4, 11, 12, 13, 14/
      DATA COMB_COEF_POS(  5,  0:  4) /  5, 15, 16, 19, 25/
      DATA COMB_COEF_POS(  6,  0:  4) /  6, 16, 17, 20, 26/
      DATA COMB_COEF_POS(  7,  0:  4) /  7, 17, 18, 21, 27/
      DATA COMB_COEF_POS(  8,  0:  4) /  8, 19, 20, 22, 28/
      DATA COMB_COEF_POS(  9,  0:  4) /  9, 20, 21, 23, 29/
      DATA COMB_COEF_POS( 10,  0:  4) / 10, 22, 23, 24, 30/
      DATA COMB_COEF_POS( 11,  0:  4) / 11, 25, 26, 28, 31/
      DATA COMB_COEF_POS( 12,  0:  4) / 12, 26, 27, 29, 32/
      DATA COMB_COEF_POS( 13,  0:  4) / 13, 28, 29, 30, 33/
      DATA COMB_COEF_POS( 14,  0:  4) / 14, 31, 32, 33, 34/
      DATA COMB_COEF_POS( 15,  0:  4) / 15, 35, 36, 40, 50/
      DATA COMB_COEF_POS( 16,  0:  4) / 16, 36, 37, 41, 51/
      DATA COMB_COEF_POS( 17,  0:  4) / 17, 37, 38, 42, 52/
      DATA COMB_COEF_POS( 18,  0:  4) / 18, 38, 39, 43, 53/
      DATA COMB_COEF_POS( 19,  0:  4) / 19, 40, 41, 44, 54/
      DATA COMB_COEF_POS( 20,  0:  4) / 20, 41, 42, 45, 55/
      DATA COMB_COEF_POS( 21,  0:  4) / 21, 42, 43, 46, 56/
      DATA COMB_COEF_POS( 22,  0:  4) / 22, 44, 45, 47, 57/
      DATA COMB_COEF_POS( 23,  0:  4) / 23, 45, 46, 48, 58/
      DATA COMB_COEF_POS( 24,  0:  4) / 24, 47, 48, 49, 59/
      DATA COMB_COEF_POS( 25,  0:  4) / 25, 50, 51, 54, 60/
      DATA COMB_COEF_POS( 26,  0:  4) / 26, 51, 52, 55, 61/
      DATA COMB_COEF_POS( 27,  0:  4) / 27, 52, 53, 56, 62/
      DATA COMB_COEF_POS( 28,  0:  4) / 28, 54, 55, 57, 63/
      DATA COMB_COEF_POS( 29,  0:  4) / 29, 55, 56, 58, 64/
      DATA COMB_COEF_POS( 30,  0:  4) / 30, 57, 58, 59, 65/
      DATA COMB_COEF_POS( 31,  0:  4) / 31, 60, 61, 63, 66/
      DATA COMB_COEF_POS( 32,  0:  4) / 32, 61, 62, 64, 67/
      DATA COMB_COEF_POS( 33,  0:  4) / 33, 63, 64, 65, 68/
      DATA COMB_COEF_POS( 34,  0:  4) / 34, 66, 67, 68, 69/

      END MODULE ML5_2_POLYNOMIAL_CONSTANTS


C     The subroutine to create the loop coefficients form the last
C      loop wf.
C     In this case of loop-induced process, the reduction is performed
C      at the loop
C     amplitude level so that no multiplication is performed.

      SUBROUTINE ML5_2_CREATE_LOOP_COEFS(LOOP_WF,RANK,LCUT_SIZE
     $ ,LOOP_GROUP_NUMBER,SYMFACT,MULTIPLIER)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      REAL*8 ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(ZERO,ONE))
      COMPLEX*16 CMPLX_ZERO
      PARAMETER (CMPLX_ZERO=(ZERO,ZERO))
      INTEGER    NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=41)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
C     
C     ARGUMENTS 
C     
      COMPLEX*16 LOOP_WF(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER RANK, SYMFACT, LCUT_SIZE, LOOP_GROUP_NUMBER, MULTIPLIER
C     
C     LOCAL VARIABLES 
C     
      COMPLEX*16 CONST
      INTEGER I,J
C     
C     GLOBAL VARIABLES
C     
      COMPLEX*16 LOOPCOEFS(0:LOOPMAXCOEFS-1,NLOOPGROUPS)
      COMMON/ML5_2_LCOEFS/LOOPCOEFS
C     
C     BEGIN CODE
C     
      CONST=CMPLX((ONE*MULTIPLIER)/SYMFACT,ZERO,KIND=16)

      CALL ML5_2_MERGE_WL(LOOP_WF,RANK,LCUT_SIZE,CONST,LOOPCOEFS(0
     $ ,LOOP_GROUP_NUMBER))

      END


C     Now the routines to update the wavefunctions



C     The subroutine to create the loop coefficients form the last
C      loop wf.
C     In this case of loop-induced process, the reduction is performed
C      at the loop
C     amplitude level so that no multiplication is performed.

      SUBROUTINE MP_ML5_2_CREATE_LOOP_COEFS(LOOP_WF,RANK,LCUT_SIZE
     $ ,LOOP_GROUP_NUMBER,SYMFACT,MULTIPLIER)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      REAL*16 ZERO,ONE
      PARAMETER (ZERO=0.0E0_16,ONE=1.0E0_16)
      COMPLEX*32 IMAG1
      PARAMETER (IMAG1=(ZERO,ONE))
      COMPLEX*32 CMPLX_ZERO
      PARAMETER (CMPLX_ZERO=(ZERO,ZERO))
      INTEGER    NLOOPGROUPS
      PARAMETER (NLOOPGROUPS=41)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)
C     
C     ARGUMENTS 
C     
      COMPLEX*32 LOOP_WF(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER RANK, SYMFACT, LCUT_SIZE, LOOP_GROUP_NUMBER, MULTIPLIER
C     
C     LOCAL VARIABLES 
C     
      COMPLEX*32 CONST
      INTEGER I,J
C     
C     GLOBAL VARIABLES
C     
      COMPLEX*32 LOOPCOEFS(0:LOOPMAXCOEFS-1,NLOOPGROUPS)
      COMMON/ML5_2_MP_LCOEFS/LOOPCOEFS
C     
C     BEGIN CODE
C     
      CONST=CMPLX((ONE*MULTIPLIER)/SYMFACT,ZERO,KIND=16)

      CALL MP_ML5_2_MERGE_WL(LOOP_WF,RANK,LCUT_SIZE,CONST,LOOPCOEFS(0
     $ ,LOOP_GROUP_NUMBER))

      END


C     Now the routines to update the wavefunctions



      SUBROUTINE ML5_2_EVAL_POLY(C,R,Q,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      COMPLEX*16 C(0:LOOPMAXCOEFS-1)
      INTEGER R
      COMPLEX*16 Q(0:3)
      COMPLEX*16 OUT

      OUT=C(0)
      IF (R.GE.1) THEN
        OUT=OUT+C(1)*Q(0)+C(2)*Q(1)+C(3)*Q(2)+C(4)*Q(3)
      ENDIF
      IF (R.GE.2) THEN
        OUT=OUT+C(5)*Q(0)*Q(0)+C(6)*Q(0)*Q(1)+C(7)*Q(1)*Q(1)+C(8)*Q(0)
     $   *Q(2)+C(9)*Q(1)*Q(2)+C(10)*Q(2)*Q(2)+C(11)*Q(0)*Q(3)+C(12)
     $   *Q(1)*Q(3)+C(13)*Q(2)*Q(3)+C(14)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.3) THEN
        OUT=OUT+C(15)*Q(0)*Q(0)*Q(0)+C(16)*Q(0)*Q(0)*Q(1)+C(17)*Q(0)
     $   *Q(1)*Q(1)+C(18)*Q(1)*Q(1)*Q(1)+C(19)*Q(0)*Q(0)*Q(2)+C(20)
     $   *Q(0)*Q(1)*Q(2)+C(21)*Q(1)*Q(1)*Q(2)+C(22)*Q(0)*Q(2)*Q(2)
     $   +C(23)*Q(1)*Q(2)*Q(2)+C(24)*Q(2)*Q(2)*Q(2)+C(25)*Q(0)*Q(0)
     $   *Q(3)+C(26)*Q(0)*Q(1)*Q(3)+C(27)*Q(1)*Q(1)*Q(3)+C(28)*Q(0)
     $   *Q(2)*Q(3)+C(29)*Q(1)*Q(2)*Q(3)+C(30)*Q(2)*Q(2)*Q(3)+C(31)
     $   *Q(0)*Q(3)*Q(3)+C(32)*Q(1)*Q(3)*Q(3)+C(33)*Q(2)*Q(3)*Q(3)
     $   +C(34)*Q(3)*Q(3)*Q(3)
      ENDIF
      END

      SUBROUTINE MP_ML5_2_EVAL_POLY(C,R,Q,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      COMPLEX*32 C(0:LOOPMAXCOEFS-1)
      INTEGER R
      COMPLEX*32 Q(0:3)
      COMPLEX*32 OUT

      OUT=C(0)
      IF (R.GE.1) THEN
        OUT=OUT+C(1)*Q(0)+C(2)*Q(1)+C(3)*Q(2)+C(4)*Q(3)
      ENDIF
      IF (R.GE.2) THEN
        OUT=OUT+C(5)*Q(0)*Q(0)+C(6)*Q(0)*Q(1)+C(7)*Q(1)*Q(1)+C(8)*Q(0)
     $   *Q(2)+C(9)*Q(1)*Q(2)+C(10)*Q(2)*Q(2)+C(11)*Q(0)*Q(3)+C(12)
     $   *Q(1)*Q(3)+C(13)*Q(2)*Q(3)+C(14)*Q(3)*Q(3)
      ENDIF
      IF (R.GE.3) THEN
        OUT=OUT+C(15)*Q(0)*Q(0)*Q(0)+C(16)*Q(0)*Q(0)*Q(1)+C(17)*Q(0)
     $   *Q(1)*Q(1)+C(18)*Q(1)*Q(1)*Q(1)+C(19)*Q(0)*Q(0)*Q(2)+C(20)
     $   *Q(0)*Q(1)*Q(2)+C(21)*Q(1)*Q(1)*Q(2)+C(22)*Q(0)*Q(2)*Q(2)
     $   +C(23)*Q(1)*Q(2)*Q(2)+C(24)*Q(2)*Q(2)*Q(2)+C(25)*Q(0)*Q(0)
     $   *Q(3)+C(26)*Q(0)*Q(1)*Q(3)+C(27)*Q(1)*Q(1)*Q(3)+C(28)*Q(0)
     $   *Q(2)*Q(3)+C(29)*Q(1)*Q(2)*Q(3)+C(30)*Q(2)*Q(2)*Q(3)+C(31)
     $   *Q(0)*Q(3)*Q(3)+C(32)*Q(1)*Q(3)*Q(3)+C(33)*Q(2)*Q(3)*Q(3)
     $   +C(34)*Q(3)*Q(3)*Q(3)
      ENDIF
      END

      SUBROUTINE ML5_2_ADD_COEFS(A,RA,B,RB)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I
      COMPLEX*16 A(0:LOOPMAXCOEFS-1),B(0:LOOPMAXCOEFS-1)
      INTEGER RA,RB

      DO I=0,NCOEF_R(RB)-1
        A(I)=A(I)+B(I)
      ENDDO
      END

      SUBROUTINE MP_ML5_2_ADD_COEFS(A,RA,B,RB)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I
      COMPLEX*32 A(0:LOOPMAXCOEFS-1),B(0:LOOPMAXCOEFS-1)
      INTEGER RA,RB

      DO I=0,NCOEF_R(RB)-1
        A(I)=A(I)+B(I)
      ENDDO
      END

      SUBROUTINE ML5_2_MERGE_WL(WL,R,LCUT_SIZE,CONST,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J
      COMPLEX*16 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER R,LCUT_SIZE
      COMPLEX*16 CONST
      COMPLEX*16 OUT(0:LOOPMAXCOEFS-1)

      DO I=1,LCUT_SIZE
        DO J=0,NCOEF_R(R)-1
          OUT(J)=OUT(J)+WL(I,J,I)*CONST
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_MERGE_WL(WL,R,LCUT_SIZE,CONST,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J
      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER R,LCUT_SIZE
      COMPLEX*32 CONST
      COMPLEX*32 OUT(0:LOOPMAXCOEFS-1)

      DO I=1,LCUT_SIZE
        DO J=0,NCOEF_R(R)-1
          OUT(J)=OUT(J)+WL(I,J,I)*CONST
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_2_UPDATE_WL_0_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_0_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_2_UPDATE_WL_1_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_1_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,0,I)*B(J,1,K)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,0,I)*B(J,2,K)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,0,I)*B(J,3,K)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,0,I)*B(J,4,K)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,1,I)*B(J,1,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,1,I)*B(J,2,K)+A(K,2,I)*B(J,1,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,2,I)*B(J,2,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,1,I)*B(J,3,K)+A(K,3,I)*B(J,1,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,2,I)*B(J,3,K)+A(K,3,I)*B(J,2,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,3,I)*B(J,3,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,1,I)*B(J,4,K)+A(K,4,I)*B(J,1,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,2,I)*B(J,4,K)+A(K,4,I)*B(J,2,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,3,I)*B(J,4,K)+A(K,4,I)*B(J,3,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,4,I)*B(J,4,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_2_UPDATE_WL_1_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_1_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,4
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_2_UPDATE_WL_0_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,0
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_0_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,0
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_2_UPDATE_WL_2_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_2_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,14
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE ML5_2_UPDATE_WL_2_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*16 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=(0.0D0,0.0D0)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.(0.0D0,0.0D0)) CYCLE
            DO L=0,14
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_2_1(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
      INTEGER I,J,K,L,M
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE
      INTEGER NEW_POSITION
      COMPLEX*32 UPDATER_COEF

C     Welcome to the computational heart of MadLoop...
      OUT(:,:,:)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
      DO J=1,OUT_SIZE
        DO M=0,4
          DO K=1,IN_SIZE
            UPDATER_COEF = B(J,M,K)
            IF (UPDATER_COEF.EQ.CMPLX(0.0E0_16,0.0E0_16,KIND=16)) CYCLE
            DO L=0,14
              NEW_POSITION = COMB_COEF_POS(L,M)
              DO I=1,LCUT_SIZE
                OUT(J,NEW_POSITION,I)=OUT(J,NEW_POSITION,I) + A(K,L,I)
     $           *UPDATER_COEF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE ML5_2_UPDATE_WL_3_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*16 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*16 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=(0.0D0,0.0D0)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,15,I)*B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,18,I)*B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,24,I)*B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,34,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END

      SUBROUTINE MP_ML5_2_UPDATE_WL_3_0(A,LCUT_SIZE,B,IN_SIZE,OUT_SIZE
     $ ,OUT)
      USE ML5_2_POLYNOMIAL_CONSTANTS
      INTEGER I,J,K
      COMPLEX*32 A(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 B(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      COMPLEX*32 OUT(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE)
      INTEGER LCUT_SIZE,IN_SIZE,OUT_SIZE

      DO I=1,LCUT_SIZE
        DO J=1,OUT_SIZE
          DO K=0,34
            OUT(J,K,I)=CMPLX(0.0E0_16,0.0E0_16,KIND=16)
          ENDDO
          DO K=1,IN_SIZE
            OUT(J,0,I)=OUT(J,0,I)+A(K,0,I)*B(J,0,K)
            OUT(J,1,I)=OUT(J,1,I)+A(K,1,I)*B(J,0,K)
            OUT(J,2,I)=OUT(J,2,I)+A(K,2,I)*B(J,0,K)
            OUT(J,3,I)=OUT(J,3,I)+A(K,3,I)*B(J,0,K)
            OUT(J,4,I)=OUT(J,4,I)+A(K,4,I)*B(J,0,K)
            OUT(J,5,I)=OUT(J,5,I)+A(K,5,I)*B(J,0,K)
            OUT(J,6,I)=OUT(J,6,I)+A(K,6,I)*B(J,0,K)
            OUT(J,7,I)=OUT(J,7,I)+A(K,7,I)*B(J,0,K)
            OUT(J,8,I)=OUT(J,8,I)+A(K,8,I)*B(J,0,K)
            OUT(J,9,I)=OUT(J,9,I)+A(K,9,I)*B(J,0,K)
            OUT(J,10,I)=OUT(J,10,I)+A(K,10,I)*B(J,0,K)
            OUT(J,11,I)=OUT(J,11,I)+A(K,11,I)*B(J,0,K)
            OUT(J,12,I)=OUT(J,12,I)+A(K,12,I)*B(J,0,K)
            OUT(J,13,I)=OUT(J,13,I)+A(K,13,I)*B(J,0,K)
            OUT(J,14,I)=OUT(J,14,I)+A(K,14,I)*B(J,0,K)
            OUT(J,15,I)=OUT(J,15,I)+A(K,15,I)*B(J,0,K)
            OUT(J,16,I)=OUT(J,16,I)+A(K,16,I)*B(J,0,K)
            OUT(J,17,I)=OUT(J,17,I)+A(K,17,I)*B(J,0,K)
            OUT(J,18,I)=OUT(J,18,I)+A(K,18,I)*B(J,0,K)
            OUT(J,19,I)=OUT(J,19,I)+A(K,19,I)*B(J,0,K)
            OUT(J,20,I)=OUT(J,20,I)+A(K,20,I)*B(J,0,K)
            OUT(J,21,I)=OUT(J,21,I)+A(K,21,I)*B(J,0,K)
            OUT(J,22,I)=OUT(J,22,I)+A(K,22,I)*B(J,0,K)
            OUT(J,23,I)=OUT(J,23,I)+A(K,23,I)*B(J,0,K)
            OUT(J,24,I)=OUT(J,24,I)+A(K,24,I)*B(J,0,K)
            OUT(J,25,I)=OUT(J,25,I)+A(K,25,I)*B(J,0,K)
            OUT(J,26,I)=OUT(J,26,I)+A(K,26,I)*B(J,0,K)
            OUT(J,27,I)=OUT(J,27,I)+A(K,27,I)*B(J,0,K)
            OUT(J,28,I)=OUT(J,28,I)+A(K,28,I)*B(J,0,K)
            OUT(J,29,I)=OUT(J,29,I)+A(K,29,I)*B(J,0,K)
            OUT(J,30,I)=OUT(J,30,I)+A(K,30,I)*B(J,0,K)
            OUT(J,31,I)=OUT(J,31,I)+A(K,31,I)*B(J,0,K)
            OUT(J,32,I)=OUT(J,32,I)+A(K,32,I)*B(J,0,K)
            OUT(J,33,I)=OUT(J,33,I)+A(K,33,I)*B(J,0,K)
            OUT(J,34,I)=OUT(J,34,I)+A(K,34,I)*B(J,0,K)
          ENDDO
        ENDDO
      ENDDO
      END
