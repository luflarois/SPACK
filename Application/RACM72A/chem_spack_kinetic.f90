 MODULE mod_chem_spack_kinetic
  
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: kinetic ! subroutine
 CONTAINS
   SUBROUTINE kinetic(jppj,Jphoto,rk,temp,xlw,Press,cosz,att,ijkbeg,ijkend,maxblock_size,nr)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the kinetic rates for the gas-phase.
!     This routine is automatically generated by SPACK.
!     Mechanism: ../Mechanism/RACM   
!     Species: ../Mechanism/ciRA72A
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     TEMP: temperature ([K]).
!     XLW: water massic fraction.
!     PRESS: pressure ([Pa]).
!     ATT: attenuation variable.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     RK: kinetic rates.
!
!------------------------------------------------------------------------
!
!     -- REMARKS
!
!------------------------------------------------------------------------
!
!     -- MODIFICATIONS
!
!------------------------------------------------------------------------
!
!     -- AUTHOR(S)
!
!     SPACK.
!
!------------------------------------------------------------------------
 
      IMPLICIT NONE
 
 
 
      INTEGER,INTENT(IN) :: jppj,ijkbeg,ijkend,maxblock_size,nr
      DOUBLE PRECISION,INTENT(IN) ::  xlw(maxblock_size),att(maxblock_size),cosz(maxblock_size)
      DOUBLE PRECISION,INTENT(IN) ::  temp(maxblock_size),Press(maxblock_size)
      DOUBLE PRECISION,INTENT(OUT) ::  rk(maxblock_size,nr)
      DOUBLE PRECISION,DIMENSION(maxblock_size) ::  Effko,Rapk,facteur
      DOUBLE PRECISION,DIMENSION(maxblock_size) :: YlH2O,SumM,azi
      DOUBLE PRECISION,INTENT(IN) ::  Jphoto(maxblock_size,jppj)
      INTEGER :: ijk
 
!     Compute third body.
!     Conversion = Avogadro*1d-6/Perfect gas constant.
!     PRESS in Pascal, SUMM in molecules/cm3, TEMP in Kelvin
 
      DO ijk=ijkbeg,ijkend
        SumM(ijk) = Press(ijk) * 7.243D16 / temp(ijk)
      END DO
 
!     Number of water molecules computed from the massic fraction
!     (absolute humidity)
 
      DO ijk=ijkbeg,ijkend
         YlH2O(ijk) = 29.d0*SumM(ijk)*xlw(ijk)/(18.d0+11.d0*xlw(ijk))
      END DO
 
!     For the zenithal angle at tropics
 
 
 
       DO ijk=ijkbeg,ijkend
       rk(ijk,  1) = Jphoto(ijk,  1)
       rk(ijk,  2) = Jphoto(ijk,  2)
       rk(ijk,  3) = Jphoto(ijk,  3)
       rk(ijk,  4) = Jphoto(ijk,  4)
       rk(ijk,  5) = Jphoto(ijk,  5)
       rk(ijk,  6) = Jphoto(ijk,  6)
       rk(ijk,  7) = Jphoto(ijk,  7)
       rk(ijk,  8) = Jphoto(ijk,  8)
       rk(ijk,  9) = Jphoto(ijk,  9)
       rk(ijk, 10) = Jphoto(ijk, 10)
       rk(ijk, 11) = Jphoto(ijk, 11)
       rk(ijk, 12) = Jphoto(ijk, 12)
       rk(ijk, 13) = Jphoto(ijk, 13)
       rk(ijk, 14) = Jphoto(ijk, 14)
       rk(ijk, 15) = Jphoto(ijk, 15)
       rk(ijk, 16) = Jphoto(ijk, 16)
       rk(ijk, 17) = Jphoto(ijk, 17)
       rk(ijk, 18) = Jphoto(ijk, 18)
       rk(ijk, 19) = Jphoto(ijk, 19)
       rk(ijk, 20) = Jphoto(ijk, 20)
       rk(ijk, 21) = Jphoto(ijk, 21)
       rk(ijk, 22) = Jphoto(ijk, 22)
       rk(ijk, 23) = Jphoto(ijk, 23)
      rk(ijk, 24) = SumM(ijk) * 6.0d-34 * (temp(ijk)/3.d2) ** (-2.3d0)
      rk(ijk, 24) = rk(ijk, 24) * SumM(ijk) * 0.2d0
      rk(ijk, 25) =  DEXP(-0.2555157957424871D+02   &
         - (  0.2060000000000000D+04 )/temp(ijk))
      rk(ijk, 26) =  DEXP(-0.2474064935803238D+02   &
         - ( -0.1100000000000000D+03 )/temp(ijk))
      rk(ijk, 26) = rk(ijk, 26) * SumM(ijk) * 0.8d0
      rk(ijk, 27) =  DEXP(-0.2416528521312882D+02   &
         - ( -0.7000000000000000D+02 )/temp(ijk))
      rk(ijk, 27) = rk(ijk, 27) * SumM(ijk) * 0.2d0
      rk(ijk, 28) =  0.2200000000000000D-09
      rk(ijk, 28) = rk(ijk, 28) * YlH2O(ijk)
      rk(ijk, 29) =  DEXP(-0.2716101748668281D+02   &
         - (  0.9400000000000000D+03 )/temp(ijk))
      rk(ijk, 30) =  DEXP(-0.3214088112211231D+02   &
         - (  0.5000000000000000D+03 )/temp(ijk))
      rk(ijk, 31) =  DEXP(-0.2375982010502066D+02   &
         - ( -0.2500000000000000D+03 )/temp(ijk))
      rk(ijk, 32) =  DEXP(-0.2656631037893612D+02   &
         - (  0.1600000000000000D+03 )/temp(ijk))
      rk(ijk, 33) = 2.2d-13 * dexp(620.0d0 / temp(ijk))   &
                    + 1.9d-33* SumM(ijk) * dexp(980.0d0 / temp(ijk))
      rk(ijk, 34) = 3.08d-34 * dexp(2820.0d0 / temp(ijk)) +    &
                    2.66d-54 * SumM(ijk) * dexp(3180.0d0 / temp(ijk))
      rk(ijk, 34) = rk(ijk, 34) * YlH2O(ijk)
      Effko(ijk) =  0.9000000000000000D-31* (temp(ijk) / 3.d2)   &
                   **(- ( 0.1500000000000000D+01))
      Rapk(ijk) =  0.3000000000000000D-10* (temp(ijk) / 3.d2)   &
                    **(- ( 0.0000000000000000D+00))
      rk(ijk, 35) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 36) =  DEXP(-0.2575921893902696D+02   &
         - ( -0.1200000000000000D+03 )/temp(ijk))
      Effko(ijk) =  0.9000000000000000D-31* (temp(ijk) / 3.d2)   &
                   **(- ( 0.2000000000000000D+01))
      Rapk(ijk) =  0.2200000000000000D-10* (temp(ijk) / 3.d2)   &
                    **(- ( 0.0000000000000000D+00))
      rk(ijk, 37) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      Effko(ijk) =  0.7000000000000000D-30* (temp(ijk) / 3.d2)   &
                   **(- ( 0.2600000000000000D+01))
      Rapk(ijk) =  0.1500000000000000D-10* (temp(ijk) / 3.d2)   &
                    **(- ( 0.5000000000000000D+00))
      rk(ijk, 38) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      Effko(ijk) =  0.2600000000000000D-29* (temp(ijk) / 3.d2)   &
                   **(- ( 0.3200000000000000D+01))
      Rapk(ijk) =  0.2400000000000000D-10* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1300000000000000D+01))
      rk(ijk, 39) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 40) =  0.2200000000000000D-10
      rk(ijk, 41) =  DEXP(-0.2632268829627837D+02   &
         - ( -0.2500000000000000D+03 )/temp(ijk))
      Effko(ijk) =  0.1800000000000000D-30* (temp(ijk) / 3.d2)   &
                   **(- ( 0.3200000000000000D+01))
      Rapk(ijk) =  0.4700000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1400000000000000D+01))
      rk(ijk, 42) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 43) =  DEXP( 0.6142746008608852D+02   &
         - (  0.1090000000000000D+05 )/temp(ijk))
      Effko(ijk) =  0.1800000000000000D-30* (temp(ijk) / 3.d2)   &
                   **(- ( 0.3200000000000000D+01))
      Rapk(ijk) =  0.4700000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1400000000000000D+01))
      facteur(ijk) = (Effko(ijk) * SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) /    &
                 Rapk(ijk))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 43) = facteur(ijk) * rk(ijk, 43)
      rk(ijk, 44) =  0.3500000000000000D-11
      rk(ijk, 45) =  DEXP(-0.2474064935803238D+02   &
         - (  0.3900000000000000D+03 )/temp(ijk))
      Effko(ijk) = 7.2d-15 * dexp(785.0d0 / temp(ijk))
      Rapk(ijk) = 4.1d-16 * dexp(1440.0d0 / temp(ijk))
      facteur(ijk) =1.9d-33 * dexp(725.0d0 / temp(ijk)) * SumM(ijk)
      rk(ijk, 46) = Effko(ijk) + facteur(ijk)/(1.0d0 + facteur(ijk) / Rapk(ijk))
      rk(ijk, 47) =  DEXP(-0.2736865685146106D+02   &
         - ( -0.3800000000000000D+03 )/temp(ijk))
      rk(ijk, 48) =  DEXP(-0.2693787393536860D+02   &
         - (  0.1400000000000000D+04 )/temp(ijk))
      rk(ijk, 49) =  DEXP(-0.2975128465212864D+02   &
         - (  0.2450000000000000D+04 )/temp(ijk))
      rk(ijk, 50) =  DEXP(-0.8860689615829534D+02   &
         - ( -0.5300000000000000D+03 )/temp(ijk))
      rk(ijk, 50) = rk(ijk, 50) * SumM(ijk) * 0.2d0
      rk(ijk, 51) =  DEXP(-0.2492297091482634D+02   &
         - ( -0.1700000000000000D+03 )/temp(ijk))
      rk(ijk, 52) =  DEXP(-0.3073211390514037D+02   &
         - (  0.1260000000000000D+04 )/temp(ijk))
      Effko(ijk) =  0.2200000000000000D-29* (temp(ijk) / 3.d2)   &
                   **(- ( 0.3900000000000000D+01))
      Rapk(ijk) =  0.1500000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.7000000000000000D+00))
      rk(ijk, 53) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 54) =  DEXP( 0.6117554523749536D+02   &
         - (  0.1100000000000000D+05 )/temp(ijk))
      Effko(ijk) =  0.2200000000000000D-29* (temp(ijk) / 3.d2)   &
                   **(- ( 0.3900000000000000D+01))
      Rapk(ijk) =  0.1500000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.7000000000000000D+00))
      facteur(ijk) = (Effko(ijk) * SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) /    &
                 Rapk(ijk))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 54) = facteur(ijk) * rk(ijk, 54)
      rk(ijk, 55) =  DEXP(-0.2779354004542632D+02   &
         - (  0.2450000000000000D+04 )/temp(ijk))
      rk(ijk, 56) =  DEXP(-0.2592627302369012D+02   &
         - (  0.2000000000000000D+04 )/temp(ijk))
      rk(ijk, 56) = rk(ijk, 56) * SumM(ijk) * 5.8d-7
      Effko(ijk) =  0.3000000000000000D-30* (temp(ijk) / 3.d2)   &
                   **(- ( 0.3300000000000000D+01))
      Rapk(ijk) =  0.1500000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.0000000000000000D+00))
      rk(ijk, 57) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk, 58) = 1.5d-13 * (1.0d0 + 2.439d-20 * SumM(ijk))
      rk(ijk, 59) =  0.6000000000000000D-10
      rk(ijk, 60) =  DEXP(-0.2486470200670236D+02   &
         - ( -0.1300000000000000D+02 )/temp(ijk))
      rk(ijk, 61) =  DEXP(-0.3943966082504782D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijk)))   &
        -  0.1361000000000000D+04/temp(ijk))
      rk(ijk, 62) =  DEXP(-0.3873183693007194D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijk)))   &
        -  0.4920000000000000D+03/temp(ijk))
      rk(ijk, 63) =  DEXP(-0.2597089008917893D+02   &
         - (  0.2600000000000000D+03 )/temp(ijk))
      rk(ijk, 64) =  DEXP(-0.2554908269405012D+02   &
         - (  0.1550000000000000D+03 )/temp(ijk))
      rk(ijk, 65) =  DEXP(-0.2483373978109839D+02   &
         - (  0.1250000000000000D+03 )/temp(ijk))
      rk(ijk, 66) =  DEXP(-0.2695807664268612D+02   &
         - ( -0.4380000000000000D+03 )/temp(ijk))
      rk(ijk, 67) =  DEXP(-0.2588705231053684D+02   &
         - ( -0.5000000000000000D+03 )/temp(ijk))
      rk(ijk, 68) =  DEXP(-0.2504325708070084D+02   &
         - ( -0.5000000000000000D+03 )/temp(ijk))
      rk(ijk, 69) =  DEXP(-0.2493639393515848D+02   &
         - ( -0.4480000000000000D+03 )/temp(ijk))
      rk(ijk, 70) =  DEXP(-0.2439627194190406D+02   &
         - ( -0.4100000000000000D+03 )/temp(ijk))
      rk(ijk, 71) =  DEXP(-0.2513781566332585D+02   &
         - ( -0.4440000000000000D+03 )/temp(ijk))
      rk(ijk, 72) =  0.1710000000000000D-09
      rk(ijk, 73) =  DEXP(-0.2703769427065081D+02   &
         - ( -0.3550000000000000D+03 )/temp(ijk))
      rk(ijk, 74) =  DEXP(-0.2564314676777420D+02   &
         - ( -0.3550000000000000D+03 )/temp(ijk))
      rk(ijk, 75) =  0.6000000000000000D-10
      rk(ijk, 76) =  0.9999999999999999D-11
      rk(ijk, 77) =  DEXP(-0.2591722318817020D+02   &
         - ( -0.3310000000000000D+03 )/temp(ijk))
      rk(ijk, 78) =  DEXP(-0.3970958044115977D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijk)))   &
        +  0.9200000000000000D+02/temp(ijk))
      rk(ijk, 79) =  0.3000000000000000D-11
      rk(ijk, 80) =  0.1140000000000000D-10
      rk(ijk, 81) =  0.1720000000000000D-10
      rk(ijk, 82) =  DEXP(-0.2470785953520939D+02   &
         - ( -0.1750000000000000D+03 )/temp(ijk))
      rk(ijk, 83) =  DEXP(-0.2429881660575334D+02   &
         - ( -0.1750000000000000D+03 )/temp(ijk))
      rk(ijk, 84) =  0.2700000000000000D-09
      rk(ijk, 85) =  DEXP(-0.2655601869289957D+02   &
         - ( -0.1900000000000000D+03 )/temp(ijk))
      rk(ijk, 86) =  DEXP(-0.2640724568430643D+02   &
         - ( -0.1900000000000000D+03 )/temp(ijk))
      rk(ijk, 87) =  DEXP(-0.2655601869289957D+02   &
         - ( -0.1900000000000000D+03 )/temp(ijk))
      rk(ijk, 88) =  0.4000000000000000D-13
      rk(ijk, 89) =  DEXP(-0.2875495121258095D+02   &
         - ( -0.5000000000000000D+03 )/temp(ijk))
      rk(ijk, 90) =  DEXP(-0.2596142928067470D+02   &
         - (  0.2600000000000000D+03 )/temp(ijk))
      rk(ijk, 91) =  DEXP(-0.2870983077730048D+02   &
         - (  0.1900000000000000D+04 )/temp(ijk))
      rk(ijk, 92) =  DEXP(-0.2729454887930734D+02   &
         - (  0.1900000000000000D+04 )/temp(ijk))
      rk(ijk, 93) =  DEXP(-0.2656631037893612D+02   &
         - (  0.1900000000000000D+04 )/temp(ijk))
      rk(ijk, 94) =  DEXP(-0.2729454887930734D+02   &
         - (  0.1900000000000000D+04 )/temp(ijk))
      rk(ijk, 95) =  DEXP(-0.3242614188587508D+02   &
         - (  0.1500000000000000D+03 )/temp(ijk))
      rk(ijk, 96) =  DEXP(-0.2887929417915106D+02   &
         - (  0.1000000000000000D+04 )/temp(ijk))
      rk(ijk, 97) =  0.2200000000000000D-10
      rk(ijk, 98) =  DEXP(-0.3986138645402777D+02   &
        + ( 0.2000000000000000D+01 * LOG(temp(ijk)))   &
        -  0.2282000000000000D+04/temp(ijk))
      rk(ijk, 99) =  DEXP(-0.2935139058906993D+02   &
         - (  0.4500000000000000D+03 )/temp(ijk))
      rk(ijk,100) =  DEXP(-0.2777720362610663D+02   &
         - ( -0.4500000000000000D+03 )/temp(ijk))
      rk(ijk,101) =  0.1000000000000000D-12
      rk(ijk,102) =  DEXP(-0.2624472675480866D+02   &
         - (  0.4460000000000000D+03 )/temp(ijk))
      rk(ijk,103) =  DEXP(-0.2745706780880511D+02   &
         - ( -0.4900000000000000D+03 )/temp(ijk))
      rk(ijk,104) =  0.1220000000000000D-10
      rk(ijk,105) =  DEXP(-0.3144773394155237D+02   &
         - (  0.5000000000000000D+03 )/temp(ijk))
      rk(ijk,106) =  DEXP(-0.3232611600944463D+02   &
         - (  0.2580000000000000D+04 )/temp(ijk))
      rk(ijk,107) =  DEXP(-0.3307320885289629D+02   &
         - (  0.1800000000000000D+04 )/temp(ijk))
      rk(ijk,108) =  DEXP(-0.3305717185398647D+02   &
         - (  0.8450000000000000D+03 )/temp(ijk))
      rk(ijk,109) =  DEXP(-0.3194352168795382D+02   &
         - (  0.2283000000000000D+04 )/temp(ijk))
      rk(ijk,110) =  DEXP(-0.3247698978846957D+02   &
         - (  0.1913000000000000D+04 )/temp(ijk))
      rk(ijk,111) =  DEXP(-0.3452882606405752D+02   &
         - (  0.7320000000000000D+03 )/temp(ijk))
      rk(ijk,112) =  0.2000000000000000D-15
      rk(ijk,113) =  DEXP(-0.3423129169516272D+02   &
         - (  0.2112000000000000D+04 )/temp(ijk))
      rk(ijk,114) =  0.2000000000000000D-17
      rk(ijk,115) =  DEXP(-0.3363861504496641D+02   &
         - (  0.1700000000000000D+04 )/temp(ijk))
      rk(ijk,116) =  0.2000000000000000D-10
      rk(ijk,117) =  0.9999999999999999D-11
      rk(ijk,118) =  0.3600000000000000D-10
      rk(ijk,119) =  DEXP(-0.3863712897853033D+02   &
         - ( -0.1044000000000000D+04 )/temp(ijk))
      rk(ijk,119) = rk(ijk,119) * SumM(ijk) * 0.2d0
      rk(ijk,120) =  0.5000000000000000D-10
      rk(ijk,121) =  0.3600000000000000D-10
      rk(ijk,122) =  DEXP(-0.3863712897853033D+02   &
         - ( -0.1044000000000000D+04 )/temp(ijk))
      rk(ijk,122) = rk(ijk,122) * SumM(ijk) * 0.2d0
      rk(ijk,123) =  0.9999999999999999D-11
      rk(ijk,124) =  0.3600000000000000D-10
      rk(ijk,125) =  DEXP(-0.3863712897853033D+02   &
         - ( -0.1044000000000000D+04 )/temp(ijk))
      rk(ijk,125) = rk(ijk,125) * SumM(ijk) * 0.2d0
      rk(ijk,126) =  0.5000000000000000D-10
      Effko(ijk) =  0.9700000000000000D-28* (temp(ijk) / 3.d2)   &
                   **(- ( 0.5600000000000000D+01))
      Rapk(ijk) =  0.9300000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1500000000000000D+01))
      rk(ijk,127) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk,128) =  DEXP( 0.6462080260895155D+02   &
         - (  0.1395400000000000D+05 )/temp(ijk))
      Effko(ijk) =  0.9700000000000000D-28* (temp(ijk) / 3.d2)   &
                   **(- ( 0.5600000000000000D+01))
      Rapk(ijk) =  0.9300000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1500000000000000D+01))
      facteur(ijk) = (Effko(ijk) * SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) /    &
                 Rapk(ijk))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk,128) = facteur(ijk) * rk(ijk,128)
      Effko(ijk) =  0.9700000000000000D-28* (temp(ijk) / 3.d2)   &
                   **(- ( 0.5600000000000000D+01))
      Rapk(ijk) =  0.9300000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1500000000000000D+01))
      rk(ijk,129) = (Effko(ijk) * SumM(ijk) /    &
                    ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *   &
                    0.6000D+00** (1.0d0 / (1.0d0 +    &
                     (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk,130) =  DEXP( 0.6462080260895155D+02   &
         - (  0.1395400000000000D+05 )/temp(ijk))
      Effko(ijk) =  0.9700000000000000D-28* (temp(ijk) / 3.d2)   &
                   **(- ( 0.5600000000000000D+01))
      Rapk(ijk) =  0.9300000000000000D-11* (temp(ijk) / 3.d2)   &
                    **(- ( 0.1500000000000000D+01))
      facteur(ijk) = (Effko(ijk) * SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) /    &
                 Rapk(ijk))) * 0.6d0 ** (1.0d0 / (1.0d0 +    &
                  (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))
      rk(ijk,130) = facteur(ijk) * rk(ijk,130)
      rk(ijk,131) =  DEXP(-0.2619593659063922D+02   &
         - ( -0.1800000000000000D+03 )/temp(ijk))
      rk(ijk,132) =  0.8700000000000000D-11
      rk(ijk,133) =  0.4000000000000000D-11
      rk(ijk,134) =  0.4000000000000000D-11
      rk(ijk,135) =  0.4000000000000000D-11
      rk(ijk,136) =  0.9000000000000000D-11
      rk(ijk,137) =  0.4000000000000000D-11
      rk(ijk,138) =  0.4000000000000000D-11
      rk(ijk,139) =  0.4000000000000000D-11
      rk(ijk,140) =  0.4000000000000000D-11
      rk(ijk,141) =  0.4000000000000000D-11
      rk(ijk,142) =  0.4000000000000000D-11
      rk(ijk,143) =  0.4000000000000000D-11
      rk(ijk,144) =  0.4000000000000000D-11
      rk(ijk,145) =  0.2000000000000000D-10
      rk(ijk,146) =  0.2000000000000000D-10
      rk(ijk,147) =  0.4000000000000000D-11
      rk(ijk,148) =  0.4000000000000000D-11
      rk(ijk,149) =  0.4000000000000000D-11
      rk(ijk,150) =  DEXP(-0.2859860514219026D+02   &
         - ( -0.8000000000000000D+03 )/temp(ijk))
      rk(ijk,151) =  DEXP(-0.2791870318838033D+02   &
         - ( -0.7000000000000000D+03 )/temp(ijk))
      rk(ijk,152) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,153) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,154) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,155) =  DEXP(-0.2929175232275020D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,156) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,157) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,158) =  DEXP(-0.2968674613099107D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,159) =  0.1500000000000000D-10
      rk(ijk,160) =  0.1500000000000000D-10
      rk(ijk,161) =  DEXP(-0.2861185036894027D+02   &
         - ( -0.9800000000000000D+03 )/temp(ijk))
      rk(ijk,162) =  DEXP(-0.2861185036894027D+02   &
         - ( -0.9800000000000000D+03 )/temp(ijk))
      rk(ijk,163) =  DEXP(-0.2861185036894027D+02   &
         - ( -0.9800000000000000D+03 )/temp(ijk))
      rk(ijk,164) =  DEXP(-0.2749125917355339D+02   &
         - ( -0.5500000000000000D+03 )/temp(ijk))
      rk(ijk,165) =  DEXP(-0.3549069430442799D+02   &
         - ( -0.2640000000000000D+04 )/temp(ijk))
      rk(ijk,166) =  DEXP(-0.2749125917355339D+02   &
         - ( -0.5500000000000000D+03 )/temp(ijk))
      rk(ijk,167) =  DEXP(-0.3549069430442799D+02   &
         - ( -0.2640000000000000D+04 )/temp(ijk))
      rk(ijk,168) =  DEXP(-0.2979384426654743D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,169) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,170) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,171) =  DEXP(-0.3002791688839384D+02   &
         - ( -0.4160000000000000D+03 )/temp(ijk))
      rk(ijk,172) =  DEXP(-0.2976809177044502D+02   &
         - ( -0.1580000000000000D+03 )/temp(ijk))
      rk(ijk,173) =  DEXP(-0.2998911891885285D+02   &
         - ( -0.4310000000000000D+03 )/temp(ijk))
      rk(ijk,174) =  DEXP(-0.2993360620892259D+02   &
         - ( -0.4670000000000000D+03 )/temp(ijk))
      rk(ijk,175) =  DEXP(-0.3076831695380433D+02   &
         - ( -0.6330000000000000D+03 )/temp(ijk))
      rk(ijk,176) =  DEXP(-0.2939711283840802D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,177) =  DEXP(-0.2955516977320235D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,178) =  DEXP(-0.3001916409728424D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,179) =  DEXP(-0.2962612150917463D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,180) =  DEXP(-0.3096643075705270D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,181) =  DEXP(-0.3096643075705270D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,182) =  DEXP(-0.3096643075705270D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,183) =  DEXP(-0.3096643075705270D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,184) =  DEXP(-0.3096643075705270D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,185) =  DEXP(-0.2416216508579258D+02   &
         - (  0.4400000000000000D+03 )/temp(ijk))
      rk(ijk,186) =  DEXP(-0.3585554469338197D+02   &
         - ( -0.2510000000000000D+04 )/temp(ijk))
      rk(ijk,187) =  DEXP(-0.2416216508579258D+02   &
         - (  0.4400000000000000D+03 )/temp(ijk))
      rk(ijk,188) =  DEXP(-0.3585554469338197D+02   &
         - ( -0.2510000000000000D+04 )/temp(ijk))
      rk(ijk,189) =  DEXP(-0.2800063657114302D+02   &
         - ( -0.5080000000000000D+03 )/temp(ijk))
      rk(ijk,190) =  DEXP(-0.2946360257967686D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,191) =  DEXP(-0.2996612940062816D+02   &
         - ( -0.7080000000000000D+03 )/temp(ijk))
      rk(ijk,192) =  DEXP(-0.2760146231368700D+02   &
         - ( -0.2110000000000000D+03 )/temp(ijk))
      rk(ijk,193) =  DEXP(-0.2800208479731938D+02   &
         - ( -0.4600000000000000D+03 )/temp(ijk))
      rk(ijk,194) =  DEXP(-0.2821262692175559D+02   &
         - ( -0.5220000000000000D+03 )/temp(ijk))
      rk(ijk,195) =  DEXP(-0.2902938805828271D+02   &
         - ( -0.6830000000000000D+03 )/temp(ijk))
      rk(ijk,196) =  DEXP(-0.2768442189265566D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,197) =  DEXP(-0.2784050834079527D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,198) =  DEXP(-0.2830632837836016D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,199) =  DEXP(-0.2790545796163031D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,200) =  DEXP(-0.2793212620871247D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,201) =  DEXP(-0.2793212620871247D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,202) =  DEXP(-0.2793212620871247D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,203) =  DEXP(-0.2793212620871247D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,204) =  DEXP(-0.2793212620871247D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,205) =  DEXP(-0.2660140169874739D+02   &
         - ( -0.5300000000000000D+03 )/temp(ijk))
      rk(ijk,206) =  DEXP(-0.2660140169874739D+02   &
         - ( -0.5300000000000000D+03 )/temp(ijk))
      rk(ijk,207) =  DEXP(-0.2791737074314655D+02   &
         - ( -0.5650000000000000D+03 )/temp(ijk))
      rk(ijk,208) =  DEXP(-0.2775318874990275D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,209) =  DEXP(-0.2825277830040182D+02   &
         - ( -0.7650000000000000D+03 )/temp(ijk))
      rk(ijk,210) =  DEXP(-0.3029028115286133D+02   &
         - ( -0.1000000000000000D+04 )/temp(ijk))
      rk(ijk,211) =  DEXP(-0.3078927231898031D+02   &
         - ( -0.1000000000000000D+04 )/temp(ijk))
      rk(ijk,212) =  DEXP(-0.3115100203358067D+02   &
         - ( -0.1000000000000000D+04 )/temp(ijk))
      rk(ijk,213) =  0.1200000000000000D-11
      rk(ijk,214) =  0.1200000000000000D-11
      rk(ijk,215) =  0.1200000000000000D-11
      rk(ijk,216) =  0.1200000000000000D-11
      rk(ijk,217) =  0.1200000000000000D-11
      rk(ijk,218) =  0.1200000000000000D-11
      rk(ijk,219) =  0.1200000000000000D-11
      rk(ijk,220) =  0.1200000000000000D-11
      rk(ijk,221) =  0.1200000000000000D-11
      rk(ijk,222) =  0.1200000000000000D-11
      rk(ijk,223) =  0.1200000000000000D-11
      rk(ijk,224) =  0.1200000000000000D-11
      rk(ijk,225) =  0.1200000000000000D-11
      rk(ijk,226) =  0.1200000000000000D-11
      rk(ijk,227) =  0.4000000000000000D-11
      rk(ijk,228) =  0.4000000000000000D-11
      rk(ijk,229) =  0.1200000000000000D-11
      rk(ijk,230) =  0.1200000000000000D-11
      rk(ijk,231) =  0.1200000000000000D-11
      rk(ijk,232) =  DEXP(-0.2942678860655414D+02   &
         - ( -0.1300000000000000D+04 )/temp(ijk))
      rk(ijk,233) =  DEXP(-0.3274868498278332D+02   &
         - ( -0.1510000000000000D+04 )/temp(ijk))
      rk(ijk,234) =  DEXP(-0.3101241587029453D+02   &
         - ( -0.1560000000000000D+04 )/temp(ijk))
      rk(ijk,235) =  DEXP(-0.3717963534647257D+02   &
         - ( -0.2950000000000000D+04 )/temp(ijk))
      rk(ijk,236) =  0.4000000000000000D-11
      rk(ijk,237) =  0.1200000000000000D-11
       END DO
 
   END SUBROUTINE kinetic
 
  END MODULE mod_chem_spack_kinetic
 
