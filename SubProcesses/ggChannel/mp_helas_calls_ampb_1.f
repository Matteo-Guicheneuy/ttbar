      SUBROUTINE ML5_2_MP_HELAS_CALLS_AMPB_1(P,NHEL,H,IC)
C     
      USE ML5_2_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=3)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=41, NLOOPGROUPS=41, NCTAMPS=79)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=120)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=10,NLOOPWAVEFUNCS=85)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_2_FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/ML5_2_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_2_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/ML5_2_MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/ML5_2_MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_2_MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NLOOPAMPS)
      COMMON/ML5_2_MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_CT_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

      CALL MP_VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))
      CALL MP_VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL MP_OXXXXX(P(0,3),MDL_MT,NHEL(3),+1*IC(3),W(1,3))
      CALL MP_IXXXXX(P(0,4),MDL_MT,NHEL(4),-1*IC(4),W(1,4))
      CALL MP_VVV1P0_1(W(1,1),W(1,2),GC_4,ZERO,ZERO,W(1,5))
C     Amplitude(s) for born diagram with ID 1
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),GC_5,AMP(1))
      CALL MP_FFV1_1(W(1,3),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,6))
C     Amplitude(s) for born diagram with ID 2
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),GC_5,AMP(2))
      CALL MP_FFV1_2(W(1,4),W(1,1),GC_5,MDL_MT,MDL_WT,W(1,7))
C     Amplitude(s) for born diagram with ID 3
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),GC_5,AMP(3))
      CALL MP_FFV1P0_3(W(1,4),W(1,3),GC_5,ZERO,ZERO,W(1,8))
C     Counter-term amplitude(s) for loop diagram number 4
      CALL MP_R2_GG_1_R2_GG_2_0(W(1,5),W(1,8),R2_GGG_1,R2_GGG_2,AMPL(1
     $ ,1))
C     Counter-term amplitude(s) for loop diagram number 5
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),R2_GQQ,AMPL(1,2))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQB_1EPS,AMPL(2,3))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQB_1EPS,AMPL(2,4))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQB_1EPS,AMPL(2,5))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQB_1EPS,AMPL(2,6))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQB_1EPS,AMPL(2,7))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQB_1EPS,AMPL(2,8))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQG_1EPS,AMPL(2,9))
      CALL MP_FFV1_0(W(1,4),W(1,3),W(1,5),UV_GQQT,AMPL(1,10))
      CALL MP_FFV1_2(W(1,4),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,9))
C     Counter-term amplitude(s) for loop diagram number 7
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,9),W(1,6),R2_QQQ,R2_QQT,AMPL(1,11))
      CALL MP_R2_QQ_2_0(W(1,9),W(1,6),UV_TMASS_1EPS,AMPL(2,12))
      CALL MP_R2_QQ_2_0(W(1,9),W(1,6),UV_TMASS,AMPL(1,13))
C     Counter-term amplitude(s) for loop diagram number 8
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),R2_GQQ,AMPL(1,14))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQB_1EPS,AMPL(2,15))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQB_1EPS,AMPL(2,16))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQB_1EPS,AMPL(2,17))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQB_1EPS,AMPL(2,18))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQB_1EPS,AMPL(2,19))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQB_1EPS,AMPL(2,20))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQG_1EPS,AMPL(2,21))
      CALL MP_FFV1_0(W(1,4),W(1,6),W(1,2),UV_GQQT,AMPL(1,22))
      CALL MP_FFV1_1(W(1,3),W(1,2),GC_5,MDL_MT,MDL_WT,W(1,10))
C     Counter-term amplitude(s) for loop diagram number 10
      CALL MP_R2_QQ_1_R2_QQ_2_0(W(1,7),W(1,10),R2_QQQ,R2_QQT,AMPL(1,23)
     $ )
      CALL MP_R2_QQ_2_0(W(1,7),W(1,10),UV_TMASS_1EPS,AMPL(2,24))
      CALL MP_R2_QQ_2_0(W(1,7),W(1,10),UV_TMASS,AMPL(1,25))
C     Counter-term amplitude(s) for loop diagram number 11
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),R2_GQQ,AMPL(1,26))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQB_1EPS,AMPL(2,27))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQB_1EPS,AMPL(2,28))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQB_1EPS,AMPL(2,29))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQB_1EPS,AMPL(2,30))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQB_1EPS,AMPL(2,31))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQB_1EPS,AMPL(2,32))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQG_1EPS,AMPL(2,33))
      CALL MP_FFV1_0(W(1,7),W(1,3),W(1,2),UV_GQQT,AMPL(1,34))
C     Counter-term amplitude(s) for loop diagram number 13
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),R2_GQQ,AMPL(1,35))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQB_1EPS,AMPL(2,36))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQB_1EPS,AMPL(2,37))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQB_1EPS,AMPL(2,38))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQB_1EPS,AMPL(2,39))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQB_1EPS,AMPL(2,40))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQB_1EPS,AMPL(2,41))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQG_1EPS,AMPL(2,42))
      CALL MP_FFV1_0(W(1,4),W(1,10),W(1,1),UV_GQQT,AMPL(1,43))
C     Counter-term amplitude(s) for loop diagram number 14
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),R2_GQQ,AMPL(1,44))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQB_1EPS,AMPL(2,45))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQB_1EPS,AMPL(2,46))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQB_1EPS,AMPL(2,47))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQB_1EPS,AMPL(2,48))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQB_1EPS,AMPL(2,49))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQB_1EPS,AMPL(2,50))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQG_1EPS,AMPL(2,51))
      CALL MP_FFV1_0(W(1,9),W(1,3),W(1,1),UV_GQQT,AMPL(1,52))
C     Counter-term amplitude(s) for loop diagram number 17
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GG,AMPL(1,53))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,54))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,55))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,56))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,57))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,58))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GB_1EPS,AMPL(2,59))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GG_1EPS,AMPL(2,60))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),UV_3GT,AMPL(1,61))
C     Counter-term amplitude(s) for loop diagram number 31
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,62))
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,63))
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,64))
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,65))
      CALL MP_R2_GG_1_0(W(1,5),W(1,8),R2_GGQ,AMPL(1,66))
C     Counter-term amplitude(s) for loop diagram number 32
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,67))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,68))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,69))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,70))
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,71))
C     Counter-term amplitude(s) for loop diagram number 34
      CALL MP_R2_GG_1_R2_GG_3_0(W(1,5),W(1,8),R2_GGQ,R2_GGT,AMPL(1,72))
C     Counter-term amplitude(s) for loop diagram number 35
      CALL MP_VVV1_0(W(1,1),W(1,2),W(1,8),R2_3GQ,AMPL(1,73))
C     At this point, all CT amps needed for (QCD=6), i.e. of split
C      order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 2000

      GOTO 1001
 2000 CONTINUE
      MP_CT_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END

