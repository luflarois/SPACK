!------------------------------------------------------------------------
!     Subroutines for automatic generation of fortran files
!     needed for integration of gas-phase chemistry.
!     Authors: Bruno Sportisse and Pierre Plion, CEREA/ENPC
!     Date: February 2003.
!------------------------------------------------------------------------

SUBROUTINE wk1(nr,k)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetic rates case.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)         :: k



WRITE(nfick,10)nr,k

10   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16)

RETURN
END SUBROUTINE wk1
!------------------------------------------------------------------------

SUBROUTINE wk2(nr,k,tact)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetic rates case ARR2:
!     k(T) = b1 * exp(-b2/T)

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(OUT)            :: k
DOUBLE PRECISION, INTENT(IN OUT)         :: tact


!Modification, added (ABS(k) < petit)
IF((ABS(tact) < petit).OR.(ABS(k) < petit))THEN
!end of modification
  CALL wk1(nr,k)
ELSE
  k = LOG(k)
  WRITE(nfick,10)nr,k,tact
END IF

10   FORMAT(6X,'rk(ijk,',i3,') =  DEXP(',d23.16,3X,'&',/8X,  &
    ' - ( ',d23.16,' )/temp(ijk))')

RETURN
END SUBROUTINE wk2
!------------------------------------------------------------------------

SUBROUTINE wk3(nr,k,expt,tact)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetic rates case.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(OUT)            :: k
DOUBLE PRECISION, INTENT(IN OUT)         :: expt
DOUBLE PRECISION, INTENT(IN OUT)         :: tact



IF(ABS(expt) < petit) THEN
  CALL wk2(nr,k,tact)
ELSE
  k = LOG(k)
  IF(tact > petit)THEN
    WRITE(nfick,10)nr,k,expt,tact
  ELSE IF(tact < -petit) THEN
    WRITE(nfick,11)nr,k,expt,(-tact)
  ELSE
    WRITE(nfick,12)nr,k,expt
  END IF
END IF


10   FORMAT(6X,'rk(ijk,',i3,') =  DEXP(',d23.16,3X,'&',/,  &
    8X,'+ (',d23.16,' * LOG(temp(ijk)))',3X,'&',/, 8X,'- ',d23.16,'/temp(ijk))')


11   FORMAT(6X,'rk(ijk,',i3,') =  DEXP(',d23.16,3X,'&',/,  &
    8X,'+ (',d23.16,' * LOG(temp(ijk)))',3X,'&',/, 8X, '+ ',d23.16,'/temp(ijk))')
12   FORMAT(6X,'rk(ijk,',i3,') =  DEXP(',d23.16,3X,'&',/,  &
    8X,'+ (',d23.16,' * LOG(temp(ijk))) ) ')

RETURN
END SUBROUTINE wk3

!Modification, routine added
SUBROUTINE wkc9(nr,k1,expt,tact1,k2,tact2,k3,tact3,k4,tact4)
!------------------------------------------------------------------------
!     -- DESCRIPTION
! Arrhenius combinations
! of the general form ARR3+ARR2+ARR2+ARR2
! 9 parameters, 3 for the ARR3 and 2 for each following ARR2 (3*2=6):
! k1*(T^expt)*exp(-tact1/T)+k2*exp(-tact2/T)+k3*exp(-tact3/T)+k4*exp(-tact4/T)
! 6 combinations implemented (formats are more general):
!-----------------------------
! arr2+0+0+arr1 -> k1,0,tact1,0,0,0,0,k4,0
! arr2+arr2+0+0 -> k1,0,tact1,k2,tact2,0,0,0,0
! arr2+arr2+0+arr1 -> k1,0,tact1,k2,tact2,0,0,k4,0
! arr2+arr2+arr2+0 -> k1,0,tact1,k2,tact2,k3,tact3,0,0
! arr2+arr2+arr2+arr2 -> k1,0,tact1,k2,tact2,k3,tact3,k4,tact4
! arr3+arr2+arr2+arr1 -> k1,expt,tact1,k2,tact2,k3,tact3,k4,tact4
! ARR3 = ARR2 with T^0; ARR2 = ARR1 with e^(-0/T)
! olhio! aindã sem chequeo de erros!!!!
!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Madeleine Sánchez, 2009.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)            :: k1
DOUBLE PRECISION, INTENT(IN OUT)            :: k2
DOUBLE PRECISION, INTENT(IN OUT)            :: k3
DOUBLE PRECISION, INTENT(IN OUT)            :: k4
DOUBLE PRECISION, INTENT(IN OUT)         :: expt
DOUBLE PRECISION, INTENT(IN OUT)         :: tact1
DOUBLE PRECISION, INTENT(IN OUT)         :: tact2
DOUBLE PRECISION, INTENT(IN OUT)         :: tact3
DOUBLE PRECISION, INTENT(IN OUT)         :: tact4
IF(ABS(expt) < petit) THEN !arr2 + ..
  IF((ABS(tact2) < petit).AND.(ABS(k2) < petit))THEN !arr2 +..+
    IF((ABS(tact3) < petit).AND.(ABS(k3) < petit))THEN !arr2 +..+..+
      IF ((ABS(tact4) < petit).AND.(ABS(k4) > petit))THEN !arr2+arr1 *
!        WRITE(nfick,*) nr, k1, tact1, k4
        WRITE(nfick,10)nr,k1,tact1,k4
      END IF 
    END IF 
  ELSE !arr2+arr2 +..
    IF((ABS(tact3) < petit).AND.(ABS(k3) < petit))THEN !arr2+arr2+..+..
      IF ((ABS(tact4) < petit).AND.(ABS(k4) < petit))THEN !arr2+arr2  *
        WRITE(nfick,11)nr,k1,tact1,k2,tact2
      ELSE IF ((ABS(tact4) < petit).AND.(ABS(k4) > petit))THEN !arr2+arr2+arr1 *
        WRITE(nfick,12)nr,k1,tact1,k2,tact2,k4
      END IF 
    ELSE !arr2+arr2+arr2+..
      IF ((ABS(tact4) < petit).AND.(ABS(k4) < petit))THEN !arr2+arr2+arr2 *
        WRITE(nfick,13)nr,k1,tact1,k2,tact2,k3,tact3
      ELSE IF ((ABS(tact4) > petit))THEN !arr2+arr2+arr2+arr2 *
        WRITE(nfick,14)nr,k1,tact1,k2,tact2,k3,tact3,k4,tact4
      END IF 
    END IF
  END IF
ELSE !arr3+arr2+arr2+arr1 *
  IF (tact1 > petit)THEN
    WRITE(nfick,15)nr,k1,expt,tact1,k2,tact2,k3,tact3,k4
  ELSE IF (tact1 < -petit) THEN
    WRITE(nfick,16)nr,k1,expt,(-tact1),k2,tact2,k3,tact3,k4
  ELSE
    WRITE(nfick,17)nr,k1,expt,k2,tact2,k3,tact3,k4
  END IF
END IF

! arr2+arr1
10   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, & 
    8X, '* DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X, d23.16)

! arr2+arr2
11   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, & 
    8X, ' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16, 3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk))')

! arr2+arr2+arr1
12   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, & 
    8X, ' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X, d23.16)

! arr2+arr2+arr2
13   FORMAT(6X,'rk(ijk,',i3,') =  ',d23.16,3X,'&',/, & 
    8X, ' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk))')

! arr2+arr2+arr2+arr2 
14   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, & 
    8X, ' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16, 3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk))')

! arr3+arr2+arr2+arr1
! tact1 > 0
15   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, &
    8X, ' * DEXP(',d23.16,' * LOG(temp(ijk))',3X,'&',/,  &
    8X, '- ',d23.16,'/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X, d23.16)
! k1*(T^expt)*exp(-tact1/T)+k2*exp(-tact2/T)+k3*exp(-tact3/T)+k4*exp(-tact4/T)
! tact1 < 0
16   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, &
    8X, ' * DEXP(',d23.16,' *   LOG(temp(ijk))',3X,'&',/, &
    8X, '+ ',d23.16,'/temp(ijk)) + &',/,  &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X, d23.16)
! tact1 = 0
17   FORMAT(6X,'rk(ijk,',i3,') = ',d23.16,3X,'&',/, &
    8X, ' * DEXP(',d23.16,' * LOG(temp(ijk))) + &',/,  &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X,d23.16,3X,'&',/, & 
    8X,' * DEXP(-( ',d23.16,' )/temp(ijk)) + & ',/8X,   &
    8X, d23.16)

RETURN
END SUBROUTINE wkc9
! End of modification
!------------------------------------------------------------------------


!------------------------------------------------------------------------

SUBROUTINE wtroe(nr,b1,b2,b3,b4,b5)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetic rates case for Troe reactions:

!     k(T)=k0*M/(1+r)  * b5**(1/(1+[log10 r]**2))
!     with r= k0*M/kinf, k0=b1*(T/300)**(-b2), kinf=b3*(T/300)**(-b4)
!     b5=0.6 if not specified (case of KINETIC TROE4 ...).
!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion and Bruno Sportisse, 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)         :: b1
DOUBLE PRECISION, INTENT(IN OUT)         :: b2
DOUBLE PRECISION, INTENT(IN OUT)         :: b3
DOUBLE PRECISION, INTENT(IN OUT)         :: b4
DOUBLE PRECISION, INTENT(IN OUT)         :: b5



WRITE(nfick,10)b1,b2
WRITE(nfick,11)b3,b4
WRITE(nfick,12)nr,b5

!     Modification BS/KS 06/2002
!     Replace LOG10 by LOG and Effko*Rapk by Effko/Rapk

10   FORMAT(6X,'Effko(ijk) = ',d23.16,'* (temp(ijk) / 3.d2)',3X,'&',/  &
    8X,'           **(- (',d23.16,'))')
11   FORMAT(6X,'Rapk(ijk) = ',d23.16,'* (temp(ijk) / 3.d2)',3X,'&',/  &
    8X,'            **(- (',d23.16,'))')
12   FORMAT(6X,'rk(ijk,',i3,') = (Effko(ijk) * SumM(ijk) / ',3X,'&',/  &
    8X,'            ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *',3X,'&',/  &
    8X,'            ',d10.4,'** (1.0d0 / (1.0d0 + ',3X,'&',/  &
    8X,'             (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))')

RETURN
END SUBROUTINE wtroe
!------------------------------------------------------------------------

SUBROUTINE wtroe7(nr,b1,b2,b3,b4,b5,b6,b7)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetics for reactions computed from equilibria.
!     k(T)=k0*M/(1+r)  * b7**(1/(1+[log10 r]**2))
!     with r= k0*M/kinf
!     k0=b1*exp(-b2/T)*(T/300)**(-b3)
!     kinf=b4*exp(-b5/T)*(T/300)**(-b6)

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion and Bruno Sportisse 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)         :: b1
DOUBLE PRECISION, INTENT(IN OUT)         :: b2
DOUBLE PRECISION, INTENT(IN OUT)         :: b3
DOUBLE PRECISION, INTENT(IN OUT)         :: b4
DOUBLE PRECISION, INTENT(IN OUT)         :: b5
DOUBLE PRECISION, INTENT(IN OUT)         :: b6
DOUBLE PRECISION, INTENT(IN OUT)         :: b7



WRITE(nfick,10)b1,b3,b2
WRITE(nfick,11)b4,b6,b5
WRITE(nfick,12)nr,b7

10   FORMAT(6X,'Effko(ijk) = ',d23.16,'* (temp(ijk) / 3.d2)',3X,'&',/  &
    8X,'           **(- (',d23.16,'))*',3X,'&',/  &
    8X,'           dexp(-',d23.16,'/temp(ijk))')
11   FORMAT(6X,'Rapk(ijk) = ',d23.16,'* (temp(ijk) / 3.d2)',3X,'&',/  &
    8X,'            **(- (',d23.16,'))*',3X,'&',/  &
    8X,'           dexp(-',d23.16,'/temp(ijk))')
12   FORMAT(6X,'rk(ijk,',i3,') = (Effko(ijk)*SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) / ',3X,'&',/  &
    8X,'         Rapk(ijk))) * ',d10.4,'** (1.0d0 / (1.0d0 + ',3X,'&',/  &
    8X,'          (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))')
!     13    format(6x,'rk(',i3,') = facteur * (',D23.16,') * '/
!     &       5x,'&      dexp((',D23.16,')/temp)')

RETURN
END SUBROUTINE wtroe7
!------------------------------------------------------------------------

SUBROUTINE wtroe10(nr,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetic rates case for Troe reactions/MOCA.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion and Bruno Sportisse, 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)         :: b1
DOUBLE PRECISION, INTENT(IN OUT)         :: b2
DOUBLE PRECISION, INTENT(IN OUT)         :: b3
DOUBLE PRECISION, INTENT(IN OUT)         :: b4
DOUBLE PRECISION, INTENT(IN OUT)         :: b5
DOUBLE PRECISION, INTENT(IN OUT)         :: b6
DOUBLE PRECISION, INTENT(IN OUT)         :: b7
DOUBLE PRECISION, INTENT(IN OUT)         :: b8
DOUBLE PRECISION, INTENT(IN OUT)         :: b9
DOUBLE PRECISION, INTENT(IN OUT)         :: b10



!     ! ! ! MOCA computes QLN = log(k)
!     Modification BS/KS: 22/03/2002
!     Error for b5,b6 -> b6,b5

IF(b6 > petit)THEN
  WRITE(nfick,10)b4,b6,b5
ELSE IF(b6 < -petit)THEN
  WRITE(nfick,20)b4,(-b6),(b5)
ELSE
  WRITE(nfick,30)b4,b5
END IF

!     End modification.

!     Modification BS/KS: 20/03/2002
!     Error b1 --> log(b1)

IF(b3 > petit)THEN
  WRITE(nfick,11)LOG(b1),b2,b3
ELSE IF(b3 < -petit)THEN
  WRITE(nfick,21)LOG(b1),b2,(-b3)
ELSE
  WRITE(nfick,31)LOG(b1),b2
END IF

!     End modification.
WRITE(nfick,12)
IF(ABS(b8) < petit)THEN
  WRITE(nfick,14)nr,b7
ELSE
  WRITE(nfick,15)(1.d0-b7),b8,b7,b9,b10
  IF(ABS(b10) < petit) WRITE(nfick,16)
  WRITE(nfick,17)nr
END IF

10   FORMAT(6X,'Effko(ijk) = ',d23.16,'* SumM(ijk) * ',3X,'&',/  &
    8X,' DEXP(-(',d23.16,')/temp(ijk) + ',3X,'&',/  &
    8X, '(', d23.16,') * LOG(temp(ijk)/3.d2))')

20   FORMAT(6X,'Effko(ijk) = ',d23.16,'*DEXP(',d23.16,'/temp(ijk) + ',3X,'&',/  &
    8X,'(', d23.16,') * LOG(temp(ijk)/3.d2)) * SumM(ijk)')

30   FORMAT(6X,'Effko(ijk) = ',d23.16,'*DEXP(',/  &
    5X,'s   ',d23.16,' * LOG(temp(ijk)/3.d2)) * SumM(ijk)')

11   FORMAT(6X,'Rapk(ijk) = Effko(ijk) / DEXP(',d23.16,' + ',3X,'&',/  &
    8X,'     (',d23.16,')* LOG(temp(ijk)/3.d2) - ',3X,'&',/  &
    8X,'     (',d23.16,')/temp(ijk))')

21   FORMAT(6X,'Rapk(ijk) = Effko(ijk) / DEXP(',d23.16,' + ',3X,'&',/  &
    8X,'     (',d23.16,')* LOG(temp(ijk)/3.d2) + ',3X,'&',/  &
    8X,'     (',d23.16,')/temp(ijk))')

31   FORMAT(6X,'Rapk(ijk) = Effko(ijk) / DEXP(',d23.16,' + ',3X,'&',/,  &
    8X,'     (',d23.16,')* LOG(temp(ijk)/3.d2) ) ')

12   FORMAT(6X,'facteur(ijk) = 1.d0/(1.d0+LOG10(Rapk(ijk))**2)')
!     12    format(6x,'facteur = 0.d0',/
!     &   6x,'if(Rapk.gt.0.d0) facteur = 1.d0/(1.d0+LOG10(Rapk)**2)')
14   FORMAT(6X,'rk(ijk,',i3,') = DEXP( LOG(Effko(ijk)/(1.d0+Rapk(ijk))) + ',3X,'&',/  &
    8X,'             LOG(',d23.16,') * facteur(ijk))')

15   FORMAT(6X,'Fcent(ijk) = ',d23.16,' * DEXP(-temp(ijk)/',3X,'&',/  &
    8X,'(', d23.16,'))',/ 8X,'      + (',d23.16,') * DEXP(-temp(ijk)/',3X,'&', /  &
    8X,'(',d23.16,'))',3X,'&',/ 8X,'      +  DEXP (-(',d23.16,')/temp(ijk))')
16   FORMAT(6X,'Fcent(ijk) = Fcent(ijk) -1.d0')
17   FORMAT(6X,'rk(ijk,',i3,') = DEXP ( LOG(Effko(ijk)/(1.d0+Rapk(ijk))) ',3X,'&',/  &
    8X,'           + facteur(ijk)*LOG(Fcent(ijk)))')

RETURN
END SUBROUTINE wtroe10
!------------------------------------------------------------------------

SUBROUTINE wspec (nr,ispebp)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetic rates case for specific reactions.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion and Bruno Sportisse, 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'


INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN OUT)                  :: ispebp

! RACM
IF(ispebp == -1) THEN
  WRITE(nfick,11)nr
ELSE IF(ispebp == -2) THEN
  WRITE(nfick,12)nr
ELSE IF(ispebp == -3) THEN
  WRITE(nfick,13)nr
ELSE IF(ispebp == -4) THEN
  WRITE(nfick,14)nr
ELSE IF(ispebp == -5) THEN
  WRITE(nfick,15)nr
ELSE IF(ispebp == -6) THEN
  WRITE(nfick,16)nr
ELSE IF(ispebp == -7) THEN
  WRITE(nfick,17)nr
!CB2002
ELSE IF(ispebp == -8) THEN
  WRITE(nfick,18)nr
ELSE IF(ispebp == -9) THEN
  WRITE(nfick,19)nr
ELSE IF(ispebp == -10) THEN
  WRITE(nfick,110)nr
ELSE IF(ispebp == -11) THEN
  WRITE(nfick,111)nr
ELSE IF(ispebp == -12) THEN
  WRITE(nfick,112)nr
ELSE
  WRITE(*,*) 'ERROR: unknown specific reaction ',ispebp
  STOP
END IF

11   FORMAT(6X,'rk(ijk,',i3,') = SumM(ijk) * 6.0d-34 * (temp(ijk)/3.d2) ** (-2.3d0)')
!     BS 05/02/2003 values given by RACM
!     12   format(6x,'rk(',i3,') = 2.3d-13 * dexp(600.0d0 / temp)',/
!     &     5x,'&            + 1.73d-33* SumM * dexp(1000.0d0 / temp)')
!     13   format(6x,'rk(',i3,') = 3.22d-34 * dexp(2800.0d0 / temp) + ',/
!     &     5x,'&            2.38d-54 * SumM * dexp(3200.0d0 / temp)')
!     MODIF BS 06/06/2003 on the basis of CMAQ
12   FORMAT(6X,'rk(ijk,',i3,') = 2.2d-13 * dexp(620.0d0 / temp(ijk))',3X,'&',/  &
    8X,'            + 1.9d-33* SumM(ijk) * dexp(980.0d0 / temp(ijk))')
13   FORMAT(6X,'rk(ijk,',i3,') = 3.08d-34 * dexp(2820.0d0 / temp(ijk)) + ', 3X,'&',/  &
    8X,'            2.66d-54 * SumM(ijk) * dexp(3180.0d0 / temp(ijk))')
!     END MODIF
14   FORMAT(6X,'Effko(ijk) = 7.2d-15 * dexp(785.0d0 / temp(ijk))',/  &
    6X,'Rapk(ijk) = 4.1d-16 * dexp(1440.0d0 / temp(ijk))',/  &
    6X,'facteur(ijk) =1.9d-33 * dexp(725.0d0 / temp(ijk)) * SumM(ijk)',/  &
    6X,'rk(ijk,',i3,') = Effko(ijk) + facteur(ijk)/(1.0d0 + facteur(ijk) / Rapk(ijk))')
15   FORMAT(6X,'rk(ijk,',i3,') = 1.5d-13 * (1.0d0 + 2.439d-20 * SumM(ijk))')
16   FORMAT(6X,'Rapk(ijk) = 3.4d-30 * (300./temp(ijk))**(3.2D0)*SumM(ijk)',/  &
    6X,'Effko(ijk) = Rapk(ijk)/(4.77D-11*(300.D0/temp(ijk))**1.4D0)',/  &
    6X,'rk(ijk,',i3,')=(Rapk(ijk)/(1.D0+Effko(ijk)))*0.3D0**',3X,'&',/  &
    8X,'(1.0d0/(1.0d0+ ((LOG10(Effko(ijk))-0.12D0)/1.2D0)**2))')
17   FORMAT(6X,'rk(ijk,',i3,') = 2.0d-39 * YlH2O(ijk) * YlH2O(ijk)')
!
! CB2002
!     O + O2 -> O3 in NASA 2006
18   FORMAT(6X,'rk(ijk,',i3,') = SumM(ijk) * 6.0d-34 * (temp(ijk)/3.d2) ** (-2.4d0)')

!     N2O5 + H2O + H2O -> 2 HNO3 + H2O in CB2002
19   FORMAT(6X,'rk(ijk,',i3,') = 1.80d-39 * YlH2O(ijk) * YlH2O(ijk)')
!     HO2 + HO2 + H2O -> H2O2 + O2 + H2O in NASA 2006
110  FORMAT(6X,'rk(ijk,',i3,') = 1.7d-33 * SumM(ijk) * ', 3X,'&',/ &
    8X,'            dexp(1000.0d0 / temp(ijk)) * (1 + 1.4D-21 * ', 3X, '&',/ &
    8X,'            YlH2O(ijk) * dexp(2200.0d0 / temp(ijk)))')
!srf- 2 parenteses missing
!    8X,'            YlH2O(ijk) * dexp(2200.0d0 / temp(ijk)')

!     OH + HNO3 -> NO3 + H2O in NASA 2006
111 FORMAT(6X,'Effko(ijk) = 2.4d-14 * dexp(460.0d0 / temp(ijk))',/  &
    6X,'Rapk(ijk) = 2.7d-17 * dexp(2199.0d0 / temp(ijk))',/  &
    6X,'facteur(ijk) = 6.5d-34 * dexp(1335.0d0 / temp(ijk)) * SumM(ijk)',/  &
    6X,'rk(ijk,',i3,') = Effko(ijk) + facteur(ijk)/(1.0d0 + facteur(ijk) / Rapk(ijk))')
!     CO + OH -> H + CO2 then H + O2 -> HO2  in NASA 2006
112  FORMAT(6X,'Effko(ijk) = 1.5d-13  * (temp(ijk)/300.D0)**(0.6D0)',/  &

    6X,'Rapk(ijk) = 2.1d+09 * (temp(ijk)/300.D0)**(6.1D0)',/  &
!srf- paremteses missing
!    6X,'Rapk(ijk) = 2.1d+09 * (temp(ijk)/300.)**6.1)',/  &

    6X,'rk(ijk,',i3,')=(Effko(ijk)/(1.D0+Effko(ijk)/(Rapk(ijk)/SumM(ijk)))) * ',3X,'&',/  &
    8X,'0.6D0**(1.0d0+(DLOG10(Effko(ijk)/(Rapk(ijk)/SumM(ijk)))**2.0D0)**(-1.0D0))')

RETURN
END SUBROUTINE wspec
!------------------------------------------------------------------------

SUBROUTINE wtb(nr,ittb)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write third body kinetics.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion and Bruno Sportisse, 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'


INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN OUT)                  :: ittb


IF(ittb == 1) THEN
  WRITE(nfick,11)nr,nr
ELSE IF (ittb == 2) THEN
  WRITE(nfick,12)nr,nr
ELSE IF (ittb == 3) THEN
  WRITE(nfick,13)nr,nr
ELSE IF (ittb == 4) THEN
  WRITE(nfick,14)nr,nr
ELSE IF (ittb == 5) THEN
  WRITE(nfick,15)nr,nr
END IF

11   FORMAT(6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * SumM(ijk)')
!     Seinfeld pp 22: N2 0.78; O2 0.21
12   FORMAT(6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * SumM(ijk) * 0.2d0')
13   FORMAT(6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * SumM(ijk) * 0.8d0')
14   FORMAT(6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * YlH2O(ijk)')
15   FORMAT(6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * SumM(ijk) * 5.8d-7')

RETURN
END SUBROUTINE wtb
!------------------------------------------------------------------------

SUBROUTINE wcv(nr,b1,b2,b3,b4,b5,b6,b7)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetics for variable stoichiometry for MOCA.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)         :: b1
DOUBLE PRECISION, INTENT(IN OUT)         :: b2
DOUBLE PRECISION, INTENT(IN OUT)         :: b3
DOUBLE PRECISION, INTENT(IN OUT)         :: b4
DOUBLE PRECISION, INTENT(IN OUT)         :: b5
DOUBLE PRECISION, INTENT(IN OUT)         :: b6
DOUBLE PRECISION, INTENT(IN OUT)         :: b7


DOUBLE PRECISION :: d54,d65,d76

IF(ABS(b2) < petit)THEN
  CALL wk2(nr,b1,b3)
ELSE
  CALL wk3(nr,b1,b2,b3)
END IF
!     write(nficK,10)nr,b4,b5,b6,b7
d54 = (b5-b4)*.05D0
d65 = (b6-b5)*.05D0
d76 = (b7-b6)*.05D0
WRITE(nfick,12)
IF(ABS(b4) > petit)THEN
  WRITE(nfick,13)nr,nr,b4
ELSE
  WRITE(nfick,18)nr
END IF
WRITE(nfick,14)260.,280.
IF(ABS(b4) > petit .OR. ABS(d54) > petit)THEN
  WRITE(nfick,15)nr,nr,b4,260.,d54
ELSE
  WRITE(nfick,18)nr
END IF
WRITE(nfick,14)280.,300.
IF(ABS(b5) > petit .OR. ABS(d65) > petit)THEN
  WRITE(nfick,15)nr,nr,b5,280.,d65
ELSE
  WRITE(nfick,18)nr
END IF
WRITE(nfick,14)300.,320.
IF(ABS(b6) > petit .OR. ABS(d76) > petit)THEN
  WRITE(nfick,15)nr,nr,b6,300.,d76
ELSE
  WRITE(nfick,18)nr
END IF
WRITE(nfick,16)
IF(ABS(b7) > petit)THEN
  WRITE(nfick,13)nr,nr,b7
ELSE
  WRITE(nfick,18)nr
END IF
WRITE(nfick,17)

10   FORMAT('!',6X,i3,4(2X,d23.16))
12   FORMAT(6X,'if (temp(ijk).le.260.d0) then ')
13   FORMAT(7X,' rk(ijk,',i3,') = rk(ijk,',i3,') *  ',d23.16)
14   FORMAT(6X,'elseif(temp(ijk).gt.',d10.3,'.and.temp(ijk).le.', d10.3,') then')
15   FORMAT(7X,' rk(ijk,',i3,') = rk(ijk,',i3,') * (',  &
    d23.16,3X,'&',/8X,' + (temp(ijk)-',d10.3,') * (',d23.16,'))')
16   FORMAT(6X,'else')
17   FORMAT(6X,'endif',/'!')
18   FORMAT(7X,'rk(ijk,',i3,') = 0.d0')

RETURN
END SUBROUTINE wcv
!------------------------------------------------------------------------

SUBROUTINE wrcfe(nr,b1,b2,b3,b4,b5,b6)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetics for reactions computed from equilibria.
!     k(T)=k0*M/(1+r)  * 0.6**(1/(1+[log10 r]**2)) * b5* exp(-b6/T)
!     with r= k0*M/kinf, k0=b1*(T/300)**(-b2), kinf=b3*(T/300)**(-b4)

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion and Bruno Sportisse, 2002.

!------------------------------------------------------------------------
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
DOUBLE PRECISION, INTENT(IN OUT)         :: b1
DOUBLE PRECISION, INTENT(IN OUT)         :: b2
DOUBLE PRECISION, INTENT(IN OUT)         :: b3
DOUBLE PRECISION, INTENT(IN OUT)         :: b4
DOUBLE PRECISION, INTENT(IN OUT)         :: b5
DOUBLE PRECISION, INTENT(IN OUT)         :: b6




CALL wk2(nr,b5,b6)

WRITE(nfick,10)b1,b2
WRITE(nfick,11)b3,b4
WRITE(nfick,12)
WRITE(nfick,13)nr,nr

10   FORMAT(6X,'Effko(ijk) = ',d23.16,'* (temp(ijk) / 3.d2)',3X,'&',/  &
    8X,'           **(- (',d23.16,'))')
11   FORMAT(6X,'Rapk(ijk) = ',d23.16,'* (temp(ijk) / 3.d2)',3X,'&',/  &
    8X,'            **(- (',d23.16,'))')
12   FORMAT(6X,'facteur(ijk) = (Effko(ijk) * SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) / ',  &
    3X,'&',/ 8X,'         Rapk(ijk))) * 0.6d0 ** (1.0d0 / (1.0d0 + ',3X,'&',/  &
    8X,'          (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))')
13   FORMAT(6X,'rk(ijk,',i3,') = facteur(ijk) * rk(ijk,',i3,')')
!     13    format(6x,'rk(',i3,') = facteur * (',D23.16,') * '/
!     &       5x,'&      dexp((',D23.16,')/temp)')


RETURN
END SUBROUTINE wrcfe
!------------------------------------------------------------------------

SUBROUTINE wphot(nr,ntabphot,b,tabphot,it)

!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write kinetics for tabulated photolysis at angles TABPHOT(I)
!     for 1<=I<=NTABPHOT.
!     The corresponding values for photolysis rates are BP(I,NR)

!     at 90: no photolysis.
!     at 0 and 90: first derivative = 0.

!     Interpolation with standard cubic spline (second derivative =0
!     at boundaries).

!     Correction factor according to IT (IT =1 : O3 >2 OH).

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, 2003.

!------------------------------------------------------------------------
INCLUDE 'parametre'
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN)                      :: ntabphot
DOUBLE PRECISION, INTENT(IN)             :: b(ntabphotmax)
DOUBLE PRECISION, INTENT(IN)             :: tabphot(ntabphotmax)
INTEGER, INTENT(IN OUT)                  :: it
INTEGER :: i,j
DOUBLE PRECISION :: a(ntabphotmax), c(4,nintphotmax)



!      print*,'nr= ',nr
!      print*,'ntabphot= ',ntabphot
!      print*,'b= ',b
!      print*,'tabphot= ',tabphot
!      print*,'it= ',it

DO i=1,ntabphot
  a(i)=tabphot(i)
END DO


CALL spl3(ntabphot,a,b,c)


DO i = 1,ntabphot-1
  WRITE(6,100)b(i),(c(j,i),j=1,4),b(i+1)
END DO

WRITE(nfick,12)

DO i = 1,ntabphot-1
  WRITE(nfick,11)a(i),a(i+1)
  WRITE(nfick,14)nr,c(4,i)
  WRITE(nfick,15)nr,c(3,i),a(i),nr
  WRITE(nfick,15)nr,c(2,i),a(i),nr
  WRITE(nfick,15)nr,c(1,i),a(i),nr
  
  
END DO
WRITE(nfick,16)nr,b(ntabphot)
WRITE(nfick,13)
WRITE(nfick,10)nr,nr


IF (it == 1) THEN
  WRITE(nfick,103)nr,nr
END IF

!srf
!10   FORMAT(6X,'if(att.lt.0.99999) rk(ijk,',i3,') = rk(ijk,',i3,') * att',/'!')
10   FORMAT(6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * att(ijk)',/'!')
11   FORMAT(6X,'elseif(azi(ijk).ge.',d9.2, ' .and. azi(ijk).lt.',d9.2,') then')
14   FORMAT(7X,'rk(ijk,',i3,')=',d23.16)
15   FORMAT(7X,'rk(ijk,',i3,')=',d23.16,'+(azi(ijk)-',d9.2,') * rk(ijk,',i3,')')
12   FORMAT('!',/6X,'if(azi(ijk).lt.0.D0)then',/,7X,'STOP')
16   FORMAT(6X,'elseif(azi(ijk).ge.90.D0)then',/,7X,'rk(ijk,',i3,')=',d23.16)
13   FORMAT(6X,'endif')
100  FORMAT(6(2X,d23.16))

103  FORMAT(6X,'VO2  = 3.2D-11 * exp(70.D0 /temp(ijk)) * (0.21D0*SumM(ijk))',/  &
    6X,'VN2  = 1.8D-11 * exp(110.D0/temp(ijk)) * (0.79D0*SumM(ijk))',/  &
    6X,'VH2O = 2.2D-10 * YlH2O(ijk)'/  &
    6X,'rk(ijk,',i3,') = rk(ijk,',i3,') * VH2O / (VH2O + VN2 + VO2) ')

RETURN
END SUBROUTINE wphot


!KML--------------------------------------------------------------------

SUBROUTINE wphot_tuvonline(nr)

INCLUDE 'parametre'
INCLUDE 'nficfort'


WRITE(nfick,10)nr,nr


10      FORMAT(7X,'rk(ijk,',i3,') = Jphoto(ijk,',i3,')')
RETURN
END SUBROUTINE wphot_tuvonline
!KML--------------------------------------------------------------------


!------------------------------------------------------------------------

SUBROUTINE ww(nr,NE,ie1,ie2,ie3)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write reaction rates.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN OUT)                  :: NE
INTEGER, INTENT(IN OUT)                  :: ie1
INTEGER, INTENT(IN OUT)                  :: ie2
INTEGER, INTENT(IN OUT)                  :: ie3


IF    (NE == 1)THEN
  WRITE(nficw,10)nr,nr,ie1
!     write(nficFF,100)nr,ie1
ELSE IF(NE == 2)THEN
  WRITE(nficw,11)nr,nr,ie1,ie2
!     write(nficFF,110)nr,ie1,ie2
ELSE IF(NE == 3)THEN
  WRITE(nficw,12)nr,nr,ie1,ie2,ie3
!     write(nficFF,120)nr,ie1,ie2,ie3
END IF

10   FORMAT(6X,'w(ijk,',i3,') =  rk(ijk,',i3,') * Y(ijk,',i3,')')
11   FORMAT(6X,'w(ijk,',i3,') =  rk(ijk,',i3,') * Y(ijk,',i3,')', ' * Y(ijk,',i3,')' )
12   FORMAT(6X,'w(ijk,',i3,') =  rk(ijk,',i3,') * Y(ijk,',i3,')', ' * Y(ijk,',i3,')',  &
    ' * Y(ijk,',i3,')' )

RETURN
END SUBROUTINE ww
!------------------------------------------------------------------------

SUBROUTINE dw(nr,NE,ie1,ie2,ie3)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Write derivative of reaction rates.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'

INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN OUT)                  :: NE
INTEGER, INTENT(IN OUT)                  :: ie1
INTEGER, INTENT(IN OUT)                  :: ie2
INTEGER, INTENT(IN OUT)                  :: ie3


IF    (NE == 1)THEN
  WRITE(nficdw,10)nr,ie1,nr
!     write(nficJJ,100)nr
ELSE IF(NE == 2)THEN
  WRITE(nficdw,11)nr,ie1,nr,ie2
!     write(nficJJ,110)nr,ie2
  WRITE(nficdw,11)nr,ie2,nr,ie1
!     write(nficJJ,111)nr,ie1
ELSE IF(NE == 3)THEN
  WRITE(nficdw,12)nr,ie1,nr,ie2,ie3
!     write(nficJJ,120)nr,ie2,ie3
  WRITE(nficdw,12)nr,ie2,nr,ie1,ie3
!     write(nficJJ,121)nr,ie1,ie3
  WRITE(nficdw,12)nr,ie3,nr,ie1,ie2
!     write(nficJJ,122)nr,ie1,ie2
END IF

10   FORMAT(6X,'dw(ijk,',i3,',',i3,') =  rk(ijk,',i3,')')
11   FORMAT(6X,'dw(ijk,',i3,',',i3,') =  rk(ijk,',i3,') * Y(ijk,',i3,')' )
12   FORMAT(6X,'dw(ijk,',i3,',',i3,') =  rk(ijk,',i3,') * Y(ijk,',i3,')', ' * Y(ijk,',i3,')' )

RETURN
END SUBROUTINE dw
!------------------------------------------------------------------------

SUBROUTINE wfj(s,nr,jer)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Update chemical production term and Jacobian Matrix.

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Pierre Plion, 2002.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'
INCLUDE 'parametre'
INCLUDE 'auxnom'

DOUBLE PRECISION, INTENT(IN OUT)         :: s(nespmax,nrmax)
INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN OUT)                  :: jer(3,nrmax)


INTEGER :: ie,je

DO ie = 1,nespmax
  IF(ABS(s(ie,nr)) > petit)THEN
!     Products
    IF(s(ie,nr) > 0.d0)THEN
!     Products with  stoichiometry = 1
      IF(ABS(s(ie,nr)-1.d0) < petit)THEN
        WRITE(nficf,30) ie,ie,nr
        DO je = 1,3
          IF(jer(je,nr) > 0) THEN
            WRITE(nficj,31) ie,jer(je,nr),ie,jer(je,nr),nr,jer(je,nr)
          END IF
        END DO
      ELSE
        WRITE(nficf,10) ie,ie,s(ie,nr),nr
        DO je = 1,3
          IF(jer(je,nr) > 0) WRITE(nficj,11) ie,jer(je,nr),ie,jer(je,nr),  &
              s(ie,nr),nr,jer(je,nr)
        END DO
      END IF
    ELSE
!     Reactants
      IF(ABS(s(ie,nr)+1.d0) < petit)THEN
!     Reactants with stochiometry = 1
        WRITE(nficf,40) ie,ie,nr
        DO je = 1,3
          IF(jer(je,nr) > 0) WRITE(nficj,41)  &
              ie,jer(je,nr),ie,jer(je,nr),nr,jer(je,nr)
        END DO
      ELSE
        WRITE(nficf,20) ie,ie,(-s(ie,nr)),nr
        DO je = 1,3
          IF(jer(je,nr) > 0) WRITE(nficj,21) ie,jer(je,nr),ie,jer(je,nr),  &
              (-s(ie,nr)),nr,jer(je,nr)
        END DO
      END IF
    END IF
  END IF
END DO

10   FORMAT(6X,'chem(ijk,',I3,') = chem(ijk,',I3,') + ',d23.16,' * w(ijk,',i3,')')
11   FORMAT(6X,'JacC(ijk,',I3,',',i3,') = JacC(ijk,',I3,',',i3,')+',  &
    d23.16,'*dw(ijk,',i3,',',i3,')')

20   FORMAT(6X,'chem(ijk,',I3,') = chem(ijk,',I3,') - ',d23.16,' * w(ijk,',i3,')')
21   FORMAT(6X,'JacC(ijk,',I3,',',i3,') = JacC(ijk,',I3,',',i3,')-',  &
    d23.16,'*dw(ijk,',i3,',',i3,')')

30   FORMAT(6X,'chem(ijk,',I3,') = chem(ijk,',I3,') + w(ijk,',i3,')')
40   FORMAT(6X,'chem(ijk,',I3,') = chem(ijk,',I3,') - w(ijk,',i3,')')
31   FORMAT(6X,'JacC(ijk,',I3,',',i3,') = JacC(ijk,',I3,',',i3,  &
    ') + dw(ijk,',i3,',',i3,')')
41   FORMAT(6X,'JacC(ijk,',I3,',',i3,') = JacC(ijk,',I3,',',i3,  &
    ') - dw(ijk,',i3,',',i3,')')

!~ 10   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') + ',d23.16,' * w(ijk,',i3,')')
!~ 11   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,')+',  &
    !~ d23.16,'*dw(ijk,',i3,',',i3,')')

!~ 20   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') - ',d23.16,' * w(ijk,',i3,')')
!~ 21   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,')-',  &
    !~ d23.16,'*dw(ijk,',i3,',',i3,')')

!~ 30   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') + w(ijk,',i3,')')
!~ 40   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') - w(ijk,',i3,')')
!~ 31   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,  &
    !~ ') + dw(ijk,',i3,',',i3,')')
!~ 41   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,  &
    !~ ') - dw(ijk,',i3,',',i3,')')
    
RETURN
END SUBROUTINE wfj
!------------------------------------------------------------------------

SUBROUTINE wpl(s,nr,jer)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Update production and loss terms (P-Lc formulation).

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, 2003.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'
INCLUDE 'parametre'
INCLUDE 'auxnom'

DOUBLE PRECISION, INTENT(IN OUT)         :: s(nespmax,nrmax)
INTEGER, INTENT(IN OUT)                  :: nr
INTEGER, INTENT(IN OUT)                  :: jer(3,nrmax)


INTEGER :: ie

DO ie = 1,nespmax
  IF(ABS(s(ie,nr)) > petit)THEN
!     Products
    IF(s(ie,nr) > 0.d0)THEN
!     Products with  stoichiometry = 1
      IF(ABS(s(ie,nr)-1.d0) < petit)THEN
        WRITE(nficprod,30) ie,ie,nr
!        WRITE(nficprod,30)nom_aux(ie),nom_aux(ie),nr
      ELSE
        WRITE(nficprod,10) ie,ie,s(ie,nr),nr
!        WRITE(nficprod,10)nom_aux(ie),nom_aux(ie),s(ie,nr),nr
      END IF
    ELSE
!     Reactants
      IF(ABS(s(ie,nr)+1.d0) < petit)THEN
!     Reactants with stochiometry = 1
        WRITE(nficloss,40) ie,ie,nr,ie
!        WRITE(nficloss,40)nom_aux(ie),nom_aux(ie),nr,nom_aux(ie)
      ELSE
        WRITE(nficloss,20) ie,ie,(-s(ie,nr)),nr,ie
!        WRITE(nficloss,20)nom_aux(ie),nom_aux(ie),(-s(ie,nr)),nr,nom_aux(ie)
      END IF
!     Correction if bimolecular reactant: none
!     Because d(c*c)/dc=c and not 2*c as computed !
!     nm=0
!     do je=1,3
!     if (jer(je,nr).eq.ie) nm=nm+1
!     enddo
!     if (nm.gt.1) write(nficloss,50)ie,ie,nm*1.
    END IF
  END IF
END DO

10   FORMAT(6X,'prod(ijk,',I3,') = prod(ijk,',I3,') + ',d23.16,' * w(ijk,',i3,')')
30   FORMAT(6X,'prod(ijk,',I3,') = prod(ijk,',I3,') + w(ijk,',i3,')')

20   FORMAT(6X,'loss(ijk,',I3,') = loss(ijk,',I3,') + ',d23.16,' * dw(ijk,',i3,',',  &
    I3,')')
40   FORMAT(6X,'loss(ijk,',I3,') = loss(ijk,',I3,') + dw(ijk,',i3,',',I3,')')

50   FORMAT(6X,'loss(ijk,',i3,') = loss(ijk,',i3,')/',d23.16)

!~ 10   FORMAT(6X,'prod(ijk,',A,') = prod(ijk,',A,') + ',d23.16,' * w(ijk,',i3,')')
!~ 30   FORMAT(6X,'prod(ijk,',A,') = prod(ijk,',A,') + w(ijk,',i3,')')

!~ 20   FORMAT(6X,'loss(ijk,',A,') = loss(ijk,',A,') + ',d23.16,' * dw(ijk,',i3,',',  &
    !~ A,')')
!~ 40   FORMAT(6X,'loss(ijk,',A,') = loss(ijk,',A,') + dw(ijk,',i3,',',A,')')

!~ 50   FORMAT(6X,'loss(ijk,',i3,') = loss(ijk,',i3,')/',d23.16)

RETURN
END SUBROUTINE wpl
!------------------------------------------------------------------------

SUBROUTINE wnonzero(s,nrtot,jer)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Create file for nonzero entries of Jacobian matrix:
!     non_zero.dat

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, 2003.

!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'nficfort'
INCLUDE 'parametre'

DOUBLE PRECISION, INTENT(IN OUT)         :: s(nespmax,nrmax)
INTEGER, INTENT(IN)                      :: nrtot
INTEGER, INTENT(IN)                     :: jer(3,nrmax)
DOUBLE PRECISION :: jnz(nespmax,nespmax)
INTEGER :: nr
INTEGER :: ie,je,nztot

DO ie=1,nespmax
  DO je=1,nespmax
    jnz(ie,je)=0
  END DO
END DO

DO nr=1,nrtot
  DO ie = 1,nespmax
    IF(ABS(s(ie,nr)) > petit)THEN
!     Products
      IF(s(ie,nr) > 0.d0)THEN
        DO je = 1,3
          IF(jer(je,nr) > 0) THEN
            jnz(ie,jer(je,nr))=1
          END IF
        END DO
      ELSE
!     Reactants
        DO je = 1,3
          IF(jer(je,nr) > 0) THEN
            jnz(ie,jer(je,nr))=1
          END IF
        END DO
      END IF
    END IF
  END DO
END DO

nztot=0
DO ie=1,nespmax
  DO je=1,nespmax
    IF (jnz(ie,je) == 1) nztot=nztot+1
  END DO
END DO

WRITE(nficnz,*)nztot
DO ie=1,nespmax
  DO je=1,nespmax
    IF (jnz(ie,je) == 1) WRITE(nficnz,*)ie,je
  END DO
END DO

RETURN
END SUBROUTINE wnonzero
!------------------------------------------------------------------------

SUBROUTINE spl3(n1,a,b,c)
!------------------------------------------------------------------------
!     Determination des coefficients des arcs de cubiques permettant
!     d'interpoler b(a) en respectant la condition de derivee premiere
!     nulle au bord

!     n1 nombre de points, n0 nombre de segments

!     Pour a(i) < aa < a(i+1)
!     On ecrira
!     bb = c(1,i) + c(2,i)*aa + c(3,i)*aa**2 + c(4,i)*aa**3
!     ou aa = angle -a(i)

!     c(1,i) = b(i)

!     bp = c(2,i) + 2*c(3,i)*aa + 3*c(4,i)*aa**2
!     bs = 2*c(3,i) + 6*c(4,i)*aa

!     l'annulation de la derivee premiere au premier point fournit
!     c(2,1) = 0

!     1<i<n0
!     la continuite de la fonction s'ecrit
!     bb(i,a(i+1) = bb(i+1,a(i+1)) = b(i+1)
!     c(1,i)+c(2,i)*d(i)+c(3,i)*d(i)**2+c(4,i)*d(i)**3 = b(i+1)
!     ou d(i) = a(i+1)-a(i)

!     la continuite des derivees premieres s'ecrit

!     1<i<(n0-1)
!     bp(i,a(i+1)) = bp(i+1,a(i+1))
!     c(2,i) +2*c(3,i)*d(i)+3*c(4,i)*d(i)**2 = c2(i+1)

!     la continuite des derivees secondes s'ecrit

!     1<i<(n0-1)
!     bs(i,a(i+1)) = bs(i+1,a(i+1))
!     c(3,i)+3.*c(4,i)*d(i) = c(3,i+1)

!     L'annulation de la derivee premiere au dernier point s'ecrit

!     bp(n0,a(n1)) = 0
!     c(2,n0)+2.*c(3,n0)*d(n0)+3.*c(4,n0)*d(n0)**2 = 0

!     i/ La continuite de la fonction sert de relation de recurrence
!     pour calculer c(4,i) en fonction de c(3,i)
!     ii/ La continuite de la derivee premiere sert de relation de recurrence
!     pour calculer c(2,i+1) en fonction de c(2,i),c(3,i),c(4,i)
!     iii/ La continuite de la derivee seconde sert de relation de recurrence
!     pour calculer c(3,i+1) en fonction de c(3,i),c(4,i)

!     On utilise donc comme inconnue auxilliaire c(3,1) qui sera determine
!     par la condition sur la derivee premiere au dernier point
!------------------------------------------------------------------------
IMPLICIT NONE
INCLUDE 'parametre'

INTEGER, INTENT(IN)                      :: n1
DOUBLE PRECISION, INTENT(IN)             :: a(ntabphotmax)
DOUBLE PRECISION, INTENT(IN)             :: b(ntabphotmax)
DOUBLE PRECISION, INTENT(OUT)            :: c(4,nintphotmax)
INTEGER :: i, n0

DOUBLE PRECISION :: d(ntabphotmax)
DOUBLE PRECISION :: c20(nintphotmax),c21(nintphotmax)
DOUBLE PRECISION :: c30(nintphotmax)
DOUBLE PRECISION :: c31(nintphotmax),c40(nintphotmax)
DOUBLE PRECISION :: c41(nintphotmax),cc31

n0=n1-1
DO i = 1,n0
  c(1,i) = b(i)
  d(i) = a(i+1)-a(i)
END DO
c20(1) = 0.d0
c21(1) = 0.d0
c30(1) = 0.d0
c31(1) = 1.d0
c40(1) = (b(2)-c(1,1)-c20(1)*d(1)-c30(1)* d(1)**2)/d(1)**3
c41(1) = (             -c21(1)*d(1)-c31(1)*d(1)**2)/d(1)**3
DO i = 2,n0
!     c(2,i) +2*c(3,i)*d(i)+3*c(4,i)*d(i)**2 = c2(i+1)
  c20(i) = c20(i-1)+2.d0*c30(i-1)*d(i-1)+3.d0*c40(i-1)*d(i-1)**2
  c21(i) = c21(i-1)+2.d0*c31(i-1)*d(i-1)+3.d0*c41(i-1)*d(i-1)**2
!     c(3,i)+3.*c(4,i)*d(i) = c(3,i+1)
  c30(i) = c30(i-1) +3.d0*c40(i-1)*d(i-1)
  c31(i) = c31(i-1) +3.d0*c41(i-1)*d(i-1)
!     c(1,i)+c(2,i)*d(i)+c(3,i)*d(i)**2+c(4,i)*d(i)**3 = b(i+1)
!     c4(i) = (b(i+1)-c(1,i)-c(2,i)*d(i)-c(3,i)*d(i)**2)/d(i)**3
  c40(i) = (b(i+1)-c(1,i)-c20(i)*d(i)-c30(i)* d(i)**2)/d(i)**3
  c41(i) = (             -c21(i)*d(i)-c31(i)*d(i)**2)/d(i)**3
END DO
!     c(2,n0)+2.*c(3,n0)*d(n0)+3.*c(4,n0)*d(n0)**2 = 0
!     (c20(n0)+c21(n0)*cc31) + 2.*(c30(n0)+c31(n0)*cc31)....
cc31 = c20(n0)+2.d0*c30(n0)*d(n0)+3.d0*c40(n0)*d(n0)**2
cc31 = -cc31 /(c21(n0)+2.d0*c31(n0)*d(n0)+3.d0*c41(n0)*d(n0)**2)

DO i = 1,n0
  c(2,i) = c20(i)+cc31*c21(i)
  c(3,i) = c30(i)+cc31*c31(i)
  c(4,i) = c40(i)+cc31*c41(i)
END DO

RETURN
END SUBROUTINE spl3













