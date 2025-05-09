C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,1)
C     
      SUBROUTINE MP_FFV1_1(F2, V3, COUP, M1, W1,F1)
      IMPLICIT NONE
      COMPLEX*32 CI
      PARAMETER (CI=(0Q0,1Q0))
      COMPLEX*32 COUP
      COMPLEX*32 F1(8)
      COMPLEX*32 F2(*)
      REAL*16 M1
      COMPLEX*32 P1(0:3)
      COMPLEX*32 V3(*)
      REAL*16 W1
      COMPLEX*32 DENOM
      F1(1) = +F2(1)+V3(1)
      F1(2) = +F2(2)+V3(2)
      F1(3) = +F2(3)+V3(3)
      F1(4) = +F2(4)+V3(4)
      P1(0) = -F1(1)
      P1(1) = -F1(2)
      P1(2) = -F1(3)
      P1(3) = -F1(4)
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      F1(5)= DENOM*CI*(F2(5)*(P1(0)*(-V3(5)+V3(8))+(P1(1)*(V3(6)-CI
     $ *(V3(7)))+(P1(2)*(+CI*(V3(6))+V3(7))+P1(3)*(-V3(5)+V3(8)))))
     $ +(F2(6)*(P1(0)*(V3(6)+CI*(V3(7)))+(P1(1)*(-1Q0)*(V3(5)+V3(8))
     $ +(P1(2)*(-1Q0)*(+CI*(V3(5)+V3(8)))+P1(3)*(V3(6)+CI*(V3(7))))))
     $ +M1*(F2(7)*(V3(5)+V3(8))+F2(8)*(V3(6)+CI*(V3(7))))))
      F1(6)= DENOM*(-CI)*(F2(5)*(P1(0)*(-V3(6)+CI*(V3(7)))+(P1(1)
     $ *(V3(5)-V3(8))+(P1(2)*(-CI*(V3(5))+CI*(V3(8)))+P1(3)*(V3(6)-CI
     $ *(V3(7))))))+(F2(6)*(P1(0)*(V3(5)+V3(8))+(P1(1)*(-1Q0)*(V3(6)
     $ +CI*(V3(7)))+(P1(2)*(+CI*(V3(6))-V3(7))-P1(3)*(V3(5)+V3(8)))))
     $ +M1*(F2(7)*(-V3(6)+CI*(V3(7)))+F2(8)*(-V3(5)+V3(8)))))
      F1(7)= DENOM*(-CI)*(F2(7)*(P1(0)*(V3(5)+V3(8))+(P1(1)*(-V3(6)+CI
     $ *(V3(7)))+(P1(2)*(-1Q0)*(+CI*(V3(6))+V3(7))-P1(3)*(V3(5)+V3(8)))
     $ ))+(F2(8)*(P1(0)*(V3(6)+CI*(V3(7)))+(P1(1)*(-V3(5)+V3(8))+(P1(2)
     $ *(-CI*(V3(5))+CI*(V3(8)))-P1(3)*(V3(6)+CI*(V3(7))))))+M1*(F2(5)
     $ *(-V3(5)+V3(8))+F2(6)*(V3(6)+CI*(V3(7))))))
      F1(8)= DENOM*CI*(F2(7)*(P1(0)*(-V3(6)+CI*(V3(7)))+(P1(1)*(V3(5)
     $ +V3(8))+(P1(2)*(-1Q0)*(+CI*(V3(5)+V3(8)))+P1(3)*(-V3(6)+CI
     $ *(V3(7))))))+(F2(8)*(P1(0)*(-V3(5)+V3(8))+(P1(1)*(V3(6)+CI
     $ *(V3(7)))+(P1(2)*(-CI*(V3(6))+V3(7))+P1(3)*(-V3(5)+V3(8)))))+M1
     $ *(F2(5)*(-V3(6)+CI*(V3(7)))+F2(6)*(V3(5)+V3(8)))))
      END


