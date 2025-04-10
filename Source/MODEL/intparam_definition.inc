ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        G = 2 * DSQRT(AS*PI)  ! for the first init

        MDL_LHV = 1.000000D+00

        MDL_CONJG__CKM3X3 = 1.000000D+00

        MDL_CONJG__CKM22 = 1.000000D+00

        MDL_I4X33 = 0.000000D+00

        MDL_I1X33 = 0.000000D+00

        MDL_CKM3X3 = 1.000000D+00

        MDL_CKM33 = 1.000000D+00

        MDL_CKM22 = 1.000000D+00

        MDL_NCOL = 3.000000D+00

        MDL_CA = 3.000000D+00

        MDL_TF = 5.000000D-01

        MDL_CF = (4.000000D+00/3.000000D+00)

        MDL_COMPLEXI = DCMPLX(0.000000D+00,1.000000D+00)

        MDL_MZ__EXP__2 = MDL_MZ**2

        MDL_MZ__EXP__4 = MDL_MZ**4

        MDL_SQRT__2 = SQRT(DCMPLX(2.000000D+00))

        MDL_MH__EXP__2 = MDL_MH**2

        MDL_NCOL__EXP__2 = MDL_NCOL**2

        MDL_MT__EXP__2 = MDL_MT**2

        MDL_AEW = 1.000000D+00/AEWM1

        MDL_MW = SQRT(DCMPLX(MDL_MZ__EXP__2/2.000000D+00
     $   +SQRT(DCMPLX(MDL_MZ__EXP__4/4.000000D+00-(MDL_AEW*PI
     $   *MDL_MZ__EXP__2)/(MDL_GF*MDL_SQRT__2)))))

        MDL_SQRT__AEW = SQRT(DCMPLX(MDL_AEW))

        MDL_EE = 2.000000D+00*MDL_SQRT__AEW*SQRT(DCMPLX(PI))

        MDL_MW__EXP__2 = MDL_MW**2

        MDL_SW2 = 1.000000D+00-MDL_MW__EXP__2/MDL_MZ__EXP__2

        MDL_CW = SQRT(DCMPLX(1.000000D+00-MDL_SW2))

        MDL_SQRT__SW2 = SQRT(DCMPLX(MDL_SW2))

        MDL_SW = MDL_SQRT__SW2

        MDL_G1 = MDL_EE/MDL_CW

        MDL_GW = MDL_EE/MDL_SW

        MDL_V = (2.000000D+00*MDL_MW*MDL_SW)/MDL_EE

        MDL_V__EXP__2 = MDL_V**2

        MDL_LAM = MDL_MH__EXP__2/(2.000000D+00*MDL_V__EXP__2)

        MDL_YT = (MDL_YMT*MDL_SQRT__2)/MDL_V

        MDL_YTAU = (MDL_YMTAU*MDL_SQRT__2)/MDL_V

        MDL_MUH = SQRT(DCMPLX(MDL_LAM*MDL_V__EXP__2))

        MDL_AXIALZUP = (3.000000D+00/2.000000D+00)*(-(MDL_EE*MDL_SW)
     $   /(6.000000D+00*MDL_CW))-(1.000000D+00/2.000000D+00)*((MDL_CW
     $   *MDL_EE)/(2.000000D+00*MDL_SW))

        MDL_AXIALZDOWN = (-1.000000D+00/2.000000D+00)*(-(MDL_CW*MDL_EE)
     $   /(2.000000D+00*MDL_SW))+(-3.000000D+00/2.000000D+00)*(
     $   -(MDL_EE*MDL_SW)/(6.000000D+00*MDL_CW))

        MDL_VECTORZUP = (1.000000D+00/2.000000D+00)*((MDL_CW*MDL_EE)
     $   /(2.000000D+00*MDL_SW))+(5.000000D+00/2.000000D+00)*(-(MDL_EE
     $   *MDL_SW)/(6.000000D+00*MDL_CW))

        MDL_VECTORZDOWN = (1.000000D+00/2.000000D+00)*(-(MDL_CW*MDL_EE)
     $   /(2.000000D+00*MDL_SW))+(-1.000000D+00/2.000000D+00)*(
     $   -(MDL_EE*MDL_SW)/(6.000000D+00*MDL_CW))

        MDL_VECTORAUP = (2.000000D+00*MDL_EE)/3.000000D+00

        MDL_VECTORADOWN = -(MDL_EE)/3.000000D+00

        MDL_VECTORWMDXU = (1.000000D+00/2.000000D+00)*((MDL_EE)
     $   /(MDL_SW*MDL_SQRT__2))

        MDL_AXIALWMDXU = (-1.000000D+00/2.000000D+00)*((MDL_EE)
     $   /(MDL_SW*MDL_SQRT__2))

        MDL_VECTORWPUXD = (1.000000D+00/2.000000D+00)*((MDL_EE)
     $   /(MDL_SW*MDL_SQRT__2))

        MDL_AXIALWPUXD = -(1.000000D+00/2.000000D+00)*((MDL_EE)
     $   /(MDL_SW*MDL_SQRT__2))

        MDL_I2X33 = MDL_YT*MDL_CONJG__CKM3X3

        MDL_I3X33 = MDL_CKM3X3*MDL_YT

        MDL_VECTOR_TBGP = MDL_I1X33-MDL_I2X33

        MDL_AXIAL_TBGP = -MDL_I2X33-MDL_I1X33

        MDL_VECTOR_TBGM = MDL_I3X33-MDL_I4X33

        MDL_AXIAL_TBGM = -MDL_I4X33-MDL_I3X33

        MDL_GW__EXP__2 = MDL_GW**2

        MDL_CW__EXP__2 = MDL_CW**2

        MDL_EE__EXP__2 = MDL_EE**2

        MDL_SW__EXP__2 = MDL_SW**2

        MDL_YT__EXP__2 = MDL_YT**2

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      AS = G**2/4/PI

      MDL_SQRT__AS = SQRT(DCMPLX(AS))

      MDL_G__EXP__4 = G**4

      MDL_G__EXP__2 = G**2

      MDL_G__EXP__3 = G**3

      MDL_MU_R__EXP__2 = MU_R**2

C     
C     Parameters that should be updated for the loops.
C     
      IF (UPDATELOOP) THEN
        MDL_G_UVG_1EPS_ = -((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01
     $   *PI**2))*1.100000D+01*MDL_CA

        MDL_G_UVB_1EPS_ = ((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01
     $   *PI**2))*4.000000D+00*MDL_TF

        MDL_TMASS_UV_1EPS_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $   ,DCMPLX(MDL_COMPLEXI*((MDL_G__EXP__2)/(1.600000D+01*PI**2))
     $   *3.000000D+00*MDL_CF*MDL_MT))

        MDL_G_UVT_FIN_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $   ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01*PI**2))
     $   *4.000000D+00*MDL_TF*REGLOG(DCMPLX(MDL_MT__EXP__2
     $   /MDL_MU_R__EXP__2))))

        MDL_TMASS_UV_FIN_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $   ,DCMPLX(MDL_COMPLEXI*((MDL_G__EXP__2)/(1.600000D+01*PI**2))
     $   *MDL_CF*(4.000000D+00-3.000000D+00
     $   *REGLOG(DCMPLX(MDL_MT__EXP__2/MDL_MU_R__EXP__2)))*MDL_MT))

      ENDIF
C     
C     Definition of the EW coupling used in the write out of aqed
C     
      GAL(1) = 3.5449077018110318D0 / DSQRT(ABS(AEWM1))
      GAL(2) = 1D0

