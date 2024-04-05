SUBROUTINE  creac(mot,imot,nmot)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-28  Time: 12:54:58
 
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization for gas-phase reactions.

!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     NR: number of reactions.
!     S: stoichiometric matrix

!     -- OUTPUT VARIABLES

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN)          :: mot(nbmot)
INTEGER, INTENT(IN)                      :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot


CHARACTER (LEN=12) :: nam(5)
DIMENSION  inam(5)

nr=nr+1

IF (nr > nrmax) THEN
  WRITE(*,*)'ERROR: bad dimension for nr>nrmax.'
  STOP
END IF

!     Check reactants

molec(nr)=1
nam(1)(1:imot(1)) = mot(1)(1:imot(1))
inam(1)=imot(1)
icurseur=1
100  icurseur = icurseur+1
IF (mot(icurseur)(1:imot(icurseur)) == '+') THEN
  icurseur = icurseur+1
  molec(nr)=molec(nr)+1
  nam(molec(nr))(1:imot(icurseur))= mot(icurseur)(1:imot(icurseur))
  inam(molec(nr))=imot(icurseur)
  GO TO 100
END IF

IF (molec(nr) > 3) THEN
  WRITE(*,*)'ERROR: more than 3 reactants.'
  STOP
END IF

DO je = 1 ,molec(nr)
  indic=0
  DO ie = 1 ,nesp(1)
    IF (inam(je) == inom(ie))  THEN
      IF (nam(je)(1:inam(je)) == nom(ie)(1:inom(ie))) THEN
        jer(je,nr)=ie
        s(ie,nr) = s(ie,nr) -1.d0
        indic=1
      END IF
    END IF
  END DO
  
  IF (indic == 0) THEN
    WRITE(6,*)'ERROR: the following reactant is unknown ', nam(je)(1:inam(je))
    STOP
  END IF
END DO

CALL ww(nr,molec(nr),jer(1,nr),jer(2,nr),jer(3,nr))
CALL dw(nr,molec(nr),jer(1,nr),jer(2,nr),jer(3,nr))

!     Check products

IF ((mot(icurseur)(1:imot(icurseur)) /= '>').AND.  &
      (mot(icurseur)(1:imot(icurseur)) /= '->')) THEN
  WRITE(*,*)'ERROR: > or -> expected '
  STOP
END IF
icurseur = icurseur+1

200  stoieir=1.d0
IF (icurseur < nmot) THEN
  IF ((mot(icurseur+1)(1:imot(icurseur+1)) /= '+').AND.  &
        (mot(icurseur+1)(1:imot(icurseur+1)) /= '-')) THEN
    CALL reel(stoieir,mot(icurseur),imot(icurseur))
    icurseur = icurseur+1
  END IF
END IF
nam(1)(1:imot(icurseur))= mot(icurseur)(1:imot(icurseur))

indic=0
DO ie = 1 ,nesp(1)
  IF ((imot(icurseur) == inom(ie)).AND.  &
        (nam(1)(1:imot(icurseur)) == nom(ie)(1:inom(ie)))) THEN
    s(ie,nr)=s(ie,nr)+stoieir
    indic=1
  END IF
END DO
IF (indic == 0) THEN
  WRITE(6,*)'WARNING: product unknown ',nam(1)(1:imot(icurseur))
END IF

icurseur = icurseur+1
IF (icurseur <= nmot) THEN
  IF ((mot(icurseur)(1:imot(icurseur)) == '+').OR.  &
        (mot(icurseur)(1:imot(icurseur)) == '-')) THEN
    icurseur = icurseur+1
    GO TO 200
  ELSE
    WRITE(*,*)'ERROR: + or - expected.'
    STOP
  END IF
END IF

!     Check gas/liquid

DO i=1,nesp(1)
  IF ((s(i,nr) /= 0.).AND.(indaq(i) /= indaqr(nr))) THEN
    WRITE(*,*)'ERROR: the phases are not coherent'
    WRITE(*,*)'reaction ',nr,' species ',nom(i)
    STOP
  END IF
END DO
RETURN
END SUBROUTINE  creac

!------------------------------------------------------------------------

SUBROUTINE  kinreac(mot,imot,nmot,ntuvonline)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization for kinetics of gas-phase reactions.
!     The different routines associated to the kinetic laws are
!     called.
!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!     BP(.,NR): coefficient for kinetic rates.

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)

INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN OUT)      :: mot(nbmot)
INTEGER, INTENT(IN OUT)                  :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot
INTEGER, INTENT(IN OUT)                  :: ntuvonline




DIMENSION b(ntabphotmax)

!     Arrhenius' law
!     Gas-phase only
i=2

100  CONTINUE
IF (mot(i)(1:4) == 'ARR1') THEN
  nb(nr)=1
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL wk1(nr,bp(1,nr))
ELSE IF (mot(i)(1:4) == 'ARR2') THEN
  nb(nr)=2
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL wk2(nr,bp(1,nr),bp(2,nr))
ELSE IF (mot(i)(1:4) == 'ARR3') THEN
  nb(nr)=3
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL wk3(nr,bp(1,nr),bp(2,nr),bp(3,nr))
ELSE IF (mot(i)(1:5) == 'ARRC2') THEN
  nb(nr)=2
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  iprecalc(nr)=2
ELSE IF (mot(i)(1:5) == 'ARRC3') THEN
  nb(nr)=3
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  iprecalc(nr)=3
 
! Modification to allow
! Arrhenius combinations of the general form
! ARR3+ARR2+ARR2+ARR2
ELSE IF (mot(i)(1:5) == 'ARRC9') THEN
  nb(nr)=3   !tentativamente, revisar
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL reel(bp(5,nr),mot(i+5),imot(i+5))
  CALL reel(bp(6,nr),mot(i+6),imot(i+6))
  CALL reel(bp(7,nr),mot(i+7),imot(i+7))
  CALL reel(bp(8,nr),mot(i+8),imot(i+8))
  CALL reel(bp(9,nr),mot(i+9),imot(i+9))
  CALL wkc9 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr),  &
      bp(6,nr),bp(7,nr),bp(8,nr),bp(9,nr))
!end of modification
 
!     PHOT: Photolysis.
  
ELSE IF (mot(i)(1:4) == 'PHOT') THEN
  
  nb(nr)=10
  nrphot = nrphot + 1
  IF (ntabphot == 0) THEN
    WRITE(*,*) 'ERROR: tabulated angles not defined for photolysis'
    WRITE(*,*)'Needs to be defined by SET TABULATION ...'
    STOP
  END IF
!     Read according increasing or decreasing sequence
  IF (ireversetab == 0) THEN
    DO j=1,ntabphot
      CALL reel(bp(j,nr),mot(i+j),imot(i+j))
      b(j)=bp(j,nr)
    END DO
  ELSE
    DO j=1,ntabphot
      jj=ntabphot+1-j
      CALL reel(bp(jj,nr),mot(i+j),imot(i+j))
      b(jj)=bp(jj,nr)
    END DO
  END IF
  
  
  IF(ntuvonline == 0) THEN
    CALL wphot(nr,ntabphot,b,tabphot,0)
  ELSE
    CALL wphot_tuvonline(nr)
  END IF
!     TROE: TROE/Fall-off for the general case.
ELSE IF (mot(i)(1:5) == 'TROE4') THEN
  nb(nr)=4
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL wtroe (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),0.6D0)
  
ELSE IF (mot(i)(1:5) == 'TROE5') THEN
  nb(nr)=4
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL reel(bp(5,nr),mot(i+5),imot(i+5))
  CALL wtroe (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr))
  
ELSE IF (mot(i)(1:5) == 'TROE7') THEN
  nb(nr)=4
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL reel(bp(5,nr),mot(i+5),imot(i+5))
  CALL reel(bp(6,nr),mot(i+6),imot(i+6))
  CALL reel(bp(7,nr),mot(i+7),imot(i+7))
  CALL wtroe7 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr),  &
      bp(6,nr),bp(7,nr))
  
!     TROE12: TROE/Fall-off for MOCA.
ELSE IF (mot(i)(1:6) == 'TROE10') THEN
  nb(nr)=4
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL reel(bp(5,nr),mot(i+5),imot(i+5))
  CALL reel(bp(6,nr),mot(i+6),imot(i+6))
  CALL reel(bp(7,nr),mot(i+7),imot(i+7))
  CALL reel(bp(8,nr),mot(i+8),imot(i+8))
  CALL reel(bp(9,nr),mot(i+9),imot(i+9))
  CALL reel(bp(10,nr),mot(i+10),imot(i+10))
  CALL wtroe10 (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr),  &
      bp(6,nr),bp(7,nr),bp(8,nr),bp(9,nr),bp(10,nr))
  
!     CVAR/MOCA: temperature-dependent stoichiometry.
ELSE IF (mot(i)(1:5) == 'CVAR') THEN
  nb(nr)=5
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL reel(bp(5,nr),mot(i+5),imot(i+5))
  CALL reel(bp(6,nr),mot(i+6),imot(i+6))
  CALL reel(bp(7,nr),mot(i+7),imot(i+7))
  CALL wcv (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr),  &
      bp(6,nr),bp(7,nr))
  
!     RCFE: Reactions Calculated From Equilibria.
ELSE IF (mot(i)(1:4) == 'RCFE') THEN
  nb(nr)=8
  CALL reel(bp(1,nr),mot(i+1),imot(i+1))
  CALL reel(bp(2,nr),mot(i+2),imot(i+2))
  CALL reel(bp(3,nr),mot(i+3),imot(i+3))
  CALL reel(bp(4,nr),mot(i+4),imot(i+4))
  CALL reel(bp(5,nr),mot(i+5),imot(i+5))
  CALL reel(bp(6,nr),mot(i+6),imot(i+6))
  CALL wrcfe (nr,bp(1,nr),bp(2,nr),bp(3,nr),bp(4,nr),bp(5,nr), bp(6,nr))
  
!     SPEC: specific reactions.
ELSE IF (mot(i)(1:4) == 'SPEC') THEN
  nb(nr)=5
  CALL entier(ispebp(nr),mot(i+1),imot(i+1))
  CALL wspec (nr,ispebp(nr))
  
!     EXTRA: specific reaction with corrected factors
!     O3 -> 2. OH with corrected photolysis
ELSE IF (mot(i)(1:5) == 'EXTRA') THEN
  nb(nr)=10
  IF (ntabphot == 0) THEN
    WRITE(*,*)'ERROR: tabulation not given for photolysis'
    WRITE(*,*)'Needs to be defined by SET TABULATION ...'
    STOP
  END IF
!     Read according increasing or decreasing sequence
  IF (ireversetab == 0) THEN
    DO j=1,ntabphot
      CALL reel(bp(j,nr),mot(i+j),imot(i+j))
      b(j)=bp(j,nr)
    END DO
  ELSE
    DO j=1,ntabphot
      jj=ntabphot+1-j
      CALL reel(bp(jj,nr),mot(i+j),imot(i+j))
      b(jj)=bp(jj,nr)
    END DO
  END IF
  CALL wphot(nr,ntabphot,b,tabphot,1)
  
!     Third body.
ELSE IF (mot(i)(1:2) == 'TB') THEN
  IF (mot(i+1)(1:1) == 'M') THEN
    ittb(nr)=1
  ELSE IF (mot(i+1)(1:2) == 'O2') THEN
    ittb(nr)=2
  ELSE IF (mot(i+1)(1:2) == 'N2') THEN
    ittb(nr)=3
  ELSE IF (mot(i+1)(1:3) == 'H2O') THEN
    ittb(nr)=4
  ELSE IF ((mot(i+1)(1:2) == 'H2').AND.(imot(i+1) == 2)) THEN
    ittb(nr)=5
  ELSE
    WRITE(*,*)'ERROR: syntax for Third Body'
    WRITE(*,*)'M, O2, N2, H20 or H2 expected.'
    STOP
  END IF
  i=i+2
ELSE
  WRITE(*,*) 'ERROR: unknown syntax for KINETIC definition ',mot(2)
  STOP
END IF

!     Modification BS/KS 21/05/2002
!     Case of a third body reaction: need for kinetics.

IF ((ittb(nr) /= 0).AND.(nb(nr) == 0)) GO TO 100
IF (ittb(nr) /= 0) CALL wtb(nr,ittb(nr))

!     Update the chemical production term and the Jacobian matrix.
CALL wfj(s,nr,jer)

!     Update the production and loss terms (P-Lc formulation)
CALL wpl(s,nr,jer)

!     Conversion mol/l -> molec/cm3.
IF (indaqr(nr) == 1) THEN
  IF(molec(nr) == 2) bp(1,nr)=bp(1,nr)*1.0D3/av
  IF(molec(nr) == 3) bp(1,nr)=bp(1,nr)*1.0D6/av**2
END IF

RETURN
END SUBROUTINE  kinreac
!------------------------------------------------------------------------

SUBROUTINE  cdis(mot,imot,nmot)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization for ionic dissociations.

!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     NR: number of reactions.
!     S: stoichiometric matrix

!     -- OUTPUT VARIABLES

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN)          :: mot(nbmot)
INTEGER, INTENT(IN)                      :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot



CHARACTER (LEN=12) :: nam(5)
DIMENSION  inam(5)

nequil=nequil+1
IF (nequil > nionx) THEN
  WRITE(*,*)'ERROR: bad dimension for nionx=',nionx
  STOP
END IF

!     Reactant.
nam(1)(1:imot(1))=mot(1)(1:imot(1))
inam(1)=imot(1)

indic=0
DO ie = 1 ,nesp(1)
  IF ((inam(1) == inom(ie)).AND.  &
        (nam(1)(1:inam(1)) == nom(ie)(1:inom(ie)))) THEN
    indic=ie
    iesp(nequil,1) = ie
  END IF
END DO
IF (indic == 0) THEN
  WRITE(6,*)'ERROR: unknown species ',nam(1)
  STOP
END IF
IF (indaq(indic) == 0) THEN
  WRITE(6,*)'ERROR: this is a gaseous species ',nam(1)
  STOP
END IF

icurseur=1
icurseur = icurseur+1
IF (mot(icurseur)(1:imot(icurseur)) /= '=') THEN
  WRITE(*,*)'ERROR: syntax = expected.'
  STOP
END IF

icurseur = icurseur+1

!     Products.
nneq(nequil)=2
inam(nneq(nequil))=imot(icurseur)
nam(nneq(nequil))(1:imot(icurseur))= mot(icurseur)(1:imot(icurseur))

200  stoeq=1.d0

!     BS 2003: to be checked !
IF (icurseur < nmot) THEN
  IF ((mot(icurseur+1)(1:imot(icurseur+1)) /= '+').AND.  &
        (mot(icurseur+1)(1:imot(icurseur+1)) /= '-')) THEN
    CALL reel(stoeq,mot(icurseur),imot(icurseur))
    icurseur = icurseur+1
  END IF
END IF

nam(nneq(nequil))(1:imot(icurseur))= mot(icurseur)(1:imot(icurseur))
inam(nneq(nequil))=imot(icurseur)


WRITE(*,*)'Ionic species: ', nam(nneq(nequil))(1:inam(nneq(nequil)))

indic=0
DO ie = 1 ,nesp(1)
  IF ((inam(nneq(nequil)) == inom(ie)).AND.  &
      (nam(nneq(nequil))(1:imot(icurseur)) == nom(ie)(1:inom(ie))))  &
      THEN
  seqion(ie,nequil) = seqion(ie,nequil) + stoeq
  iesp(nequil,nneq(nequil)) = ie
  indic=1
  IF (indaq(ie) == 0) THEN
    WRITE(*,*)'ERROR: gaseous species ',nom(ie), ' in dissociation ',nequil
    STOP
  END IF
END IF
END DO
IF (indic == 0) THEN
  WRITE(6,*)'ERROR: unknown species in dissociation ', nam(nneq(nequil))
  STOP
END IF

icurseur = icurseur+1

IF (icurseur <= nmot) THEN
  IF ((mot(icurseur)(1:imot(icurseur)) == '+').OR.  &
        (mot(icurseur)(1:imot(icurseur)) == '-'))  THEN
    icurseur = icurseur + 1
    nneq(nequil) = nneq(nequil) + 1
    GO TO 200
  ELSE
    WRITE(*,*)'ERROR: syntax + expected.'
    STOP
  END IF
END IF

!     Check lumping for more than 4 ions.
IF (nneq(nequil) > 3) THEN
  WRITE(*,*)'ERROR: more than 2 products in dissociation. See routine "lump"'
  STOP
END IF

RETURN
END SUBROUTINE  cdis
!------------------------------------------------------------------------

SUBROUTINE  kindis(mot,imot,nmot)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization of kinetics for ionic dissociations.

!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN OUT)      :: mot(nbmot)
INTEGER, INTENT(IN OUT)                  :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot





IF (mot(2)(1:5) == 'ARRC2') THEN
  CALL reel(xk1(nequil),mot(3),imot(3))
  CALL reel(xk2(nequil),mot(4),imot(4))
ELSE
  WRITE(*,*)nom(2)(1:4),' ','ERROR: syntax.'
  STOP
END IF
RETURN
END SUBROUTINE  kindis
!------------------------------------------------------------------------

SUBROUTINE  chenry(mot,imot,nmot,ieq)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization for Henry's equilibrium.

!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN)          :: mot(nbmot)
INTEGER, INTENT(IN)                      :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot
INTEGER, INTENT(IN)                      :: ieq

CHARACTER (LEN=12) :: nam(5)

DIMENSION  inam(5)

!     Henry gas --> liq
!     Creation of Henry liq --> gas.

nr=nr+1

IF (nr > nrmax) THEN
  WRITE(*,*)'ERROR: Bad dimension for nr>nrmax'
  STOP
END IF

!     Reactant.

molec(nr)=1

nam(1)(1:imot(1)) = mot(1)(1:imot(1))
inam(1)=imot(1)
indic=0
DO ie = 1 ,nesp(1)
  IF (inam(1) == inom(ie))  THEN
    IF (nam(1)(1:inam(1)) == nom(ie)(1:inom(ie))) THEN
      indr=ie
      jer(1,nr) = ie
      s(ie,nr) = -1.d0
      indic=1
    END IF
  END IF
END DO

IF (indic == 0) THEN
  WRITE(6,*)'ERROR: unknown species ',nam(1)
  STOP
END IF

!     Product.

nam(1)(1:imot(3)) = mot(3)(1:imot(3))
inam(1)=imot(3)
indic=0
DO ie = 1 ,nesp(1)
  IF (inam(1) == inom(ie))  THEN
    IF (nam(1)(1:inam(1)) == nom(ie)(1:inom(ie))) THEN
      indp=ie
      s(ie,nr) = 1.d0
      indic=1
    END IF
  END IF
END DO

IF (indic == 0) THEN
  WRITE(6,*)'ERROR: unknown species ',nam(1)
  STOP
END IF
indaqr(nr)=indaq(jer(1,nr))

!     Creation reaction L->G

nr=nr+1
IF (nr > nrmax) THEN
  WRITE(*,*)'ERROR: Bad dimension for nr>nrmax.'
  STOP
END IF

molec(nr)=1

s(indp,nr)=-1.
s(indr,nr)=1.d0
jep(nr)=indr
jer(1,nr)=indp

!     Lumping.

iheq(nr-1)=ieq
iheq(nr)=ieq
j1=jer(1,nr)
j2=jer(1,nr-1)
ihreac(j1)=nr
ihreac(j2)=nr-1

IF (indaq(jer(1,nr)) /= 1) THEN
  WRITE(*,*)'ERROR: phases are not coherent for'
  WRITE(*,*)'Henry ',nr
  WRITE(*,*)'Species ',nom(i),' has to be aqueous.'
  STOP
END IF
IF (indaq(jep(nr)) /= 0) THEN
  WRITE(*,*)'ERROR phases are not coherent for'
  WRITE(*,*)'Henry ',nr
  WRITE(*,*)'Species ',nom(i),' has to be gaseous.'
  STOP
END IF
indaqr(nr)=1

RETURN
END SUBROUTINE  chenry
!------------------------------------------------------------------------

SUBROUTINE  kinhenry(mot,imot,nmot)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization for kinetics of Henry's equilibrium.

!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN OUT)      :: mot(nbmot)
INTEGER, INTENT(IN OUT)                  :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot





!     Reaction G ->L

nr=nr-1

nb(nr)=6
i=jer(1,nr)

IF (rmol(i) == 0.) THEN
  ilettre=1
  
  DO WHILE (ilettre <= inom(i))
    ncomp=0
    IF (nom(i)(ilettre:ilettre) == 'H') THEN
      ncomp=1
    ELSE IF (nom(i)(ilettre:ilettre) == 'C') THEN
      ncomp=12
    ELSE IF (nom(i)(ilettre:ilettre) == 'O') THEN
      ncomp=16
    ELSE IF (nom(i)(ilettre:ilettre) == 'N') THEN
      ncomp=14
    ELSE IF (nom(i)(ilettre:ilettre) == 'S') THEN
      ncomp=32
    ELSE
      WRITE(*,*)'ERROR: species ',nom(i),' unknown molar mass.'
      STOP
    END IF
    nn=1
    ilet=ilettre+1
    IF (inom(i) >= ilet) THEN
      IF (nom(i)(ilet:ilet) == '2') nn=2
      IF (nom(i)(ilet:ilet) == '3') nn=3
      IF (nom(i)(ilet:ilet) == '4') nn=4
      IF (nom(i)(ilet:ilet) == '5') nn=5
      IF (nom(i)(ilet:ilet) == '6') nn=6
    END IF
    rmol(i)=rmol(i)+ncomp*nn
    IF (nn > 1) ilettre=ilettre+1
    ilettre=ilettre+1
  END DO
  WRITE(*,*)'Computed molar mass for ',nom(i),'=',rmol(i)
ELSE
  WRITE(*,*)'Read molar mass for ',nom(i),'=',rmol(i)
END IF

IF (rmol(i) == 0.) THEN
  WRITE(*,*)'ERROR: Species ',nom(i),' unknown molar mass.'
  STOP
END IF

!     Reaction L->G

nr=nr+1
nb(nr)=7
IF (mot(2)(1:5) == 'ARRC2') THEN
  CALL reel(bp(1,nr),mot(3),imot(3))
  CALL reel(bp(2,nr),mot(4),imot(4))
  CALL reel(bp(3,nr),mot(5),imot(5))
  iprecalc(nr)=2
ELSE IF (mot(2)(1:4) == 'ARR1') THEN
  CALL reel(bp(1,nr),mot(3),imot(3))
  bp(2,nr)=0.d0
ELSE IF (mot(2)(1:4) == 'ARR2') THEN
  CALL reel(bp(1,nr),mot(3),imot(3))
  CALL reel(bp(2,nr),mot(4),imot(4))
ELSE
  WRITE(*,*)'ERROR: syntax kinetics Henry ',nr-1
  STOP
END IF

RETURN
END SUBROUTINE  kinhenry
!------------------------------------------------------------------------

SUBROUTINE initphase(ntuvonline,filename)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization for multiphase description.
!     IPHASE = 1 : multiphase description.
!     IPHASE = 2 : gas-phase.
!     IPHASE = 3 : aqueous-phase.

!     IRMONODI: index reaction in one phase --> multiphase.

!     Notice that Henry's reactions are only to be taken into account for
!     the multiphase case.
!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'

INTEGER,INTENT(IN) :: ntuvonline
CHARACTER(LEN=*),INTENT(IN) :: fileName

DO ir=1,nr
  
  ip=indaqr(ir)+2
  nrp(1)=nrp(1)+1
  IF ((nb(ir) /= 6).AND.(nb(ir) /= 7)) THEN
    nrp(ip)=nrp(ip)+1
    irmonodi(nrp(ip),ip)=ir
  END IF
  
  IF (molec(ir) == 1) THEN
    nrmol1(1)=nrmol1(1)+1
    imolec1(nrmol1(1),1)=nrp(1)
    IF ((nb(ir) /= 6).AND.(nb(ir) /= 7)) THEN
      nrmol1(ip)=nrmol1(ip)+1
      imolec1(nrmol1(ip),ip)=ir
    END IF
    
  ELSE IF (molec(ir) == 2) THEN
    nrmol2(1)=nrmol2(1)+1
    imolec2(nrmol2(1),1)=nrp(1)
    nrmol2(ip)=nrmol2(ip)+1
    imolec2(nrmol2(ip),ip)=ir
    
  ELSE IF (molec(ir) == 3) THEN
    nrmol3(1)=nrmol3(1)+1
    imolec3(nrmol3(1),1)=nrp(1)
    nrmol3(ip)=nrmol3(ip)+1
    imolec3(nrmol3(ip),ip)=ir
  END IF
  
  IF (nb(ir) == 1) THEN
    narr1(1)=narr1(1)+1
    iarr1(narr1(1),1)=nrp(1)
    narr1(ip)=narr1(ip)+1
    iarr1(narr1(ip),ip)=ir
  ELSE IF (nb(ir) == 2) THEN
    narr2(1)=narr2(1)+1
    iarr2(narr2(1),1)=nrp(1)
    narr2(ip)=narr2(ip)+1
    iarr2(narr2(ip),ip)=ir
  ELSE IF (nb(ir) == 3) THEN
    narr3(1)=narr3(1)+1
    iarr3(narr3(1),1)=nrp(1)
    narr3(ip)=narr3(ip)+1
    iarr3(narr3(ip),ip)=ir
  ELSE IF (nb(ir) == 4) THEN
    narr4(1)=narr4(1)+1
    iarr4(narr4(1),1)=nrp(1)
    narr4(ip)=narr4(ip)+1
    iarr4(narr4(ip),ip)=ir
  ELSE IF (nb(ir) == 5) THEN
    narr5(1)=narr5(1)+1
    iarr5(narr5(1),1)=nrp(1)
    narr5(ip)=narr5(ip)+1
    iarr5(narr5(ip),ip)=ir
  ELSE IF (nb(ir) == 6) THEN
    narr6(1)=narr6(1)+1
    iarr6(narr6(1),1)=nrp(1)
  ELSE IF (nb(ir) == 7) THEN
    narr7(1)=narr7(1)+1
    iarr7(narr7(1),1)=nrp(1)
  ELSE IF (nb(ir) == 8) THEN
    narr8(1)=narr8(1)+1
    iarr8(narr8(1),1)=nrp(1)
    narr8(ip)=narr8(ip)+1
    iarr8(narr8(ip),ip)=ir
  END IF
  
  IF (indaqr(ir) == 1) THEN
    IF (molec(ir) == 2) THEN
      naq2(1)=naq2(1)+1
      inaq2(naq2(1),1)=ir
      naq2(3)=naq2(3)+1
      inaq2(naq2(3),1)=ir
    END IF
  END IF
  
  IF (ittb(ir) /= 0) THEN
    nthird(1)=nthird(1)+1
    indthird(nthird(1),1)=nrp(1)
    nthird(ip)=nthird(ip)+1
    indthird(nthird(ip),ip)=ir
  END IF
  
END DO

!     write(*,*)' Total number of reactions         =',nrp(1)
WRITE(*,*)' Number of gas-phase reactions     =',nrp(2)
WRITE(*,*)' Number of aqueous-phase reactions =',nrp(3)
WRITE(*,*) ' Nbr of Henry reversible reactions =',nrp(1)-nrp(2)-nrp(3)
WRITE(*,*)' Third body reactions              =',nthird(1)
WRITE(*,*)'#############################################'

WRITE(UNIT=78,FMT='(A,I3.3,A)') '  INTEGER,PARAMETER :: nr   =',nrp(2),'!Number of gas-phase reactions'
WRITE(UNIT=78,FMT='(A,I3.3,A)') '  INTEGER,PARAMETER :: nrt  =',nrp(1),'!Total Number of reactions'
WRITE(UNIT=78,FMT='(A,I3.3,A)') '  INTEGER,PARAMETER :: nrh2o=',nrp(3),'!Number of aqueous-phase reactions'

 WRITE(UNIT=78,FMT='(A)') '  '

 !IF Fastjx then call the map
 PRINT *,'Cinet->Ntuvonline=',ntuvonline
 
 !IF(ntuvonline==1) CALL JxSpack(fileName,nrphot,78)
 CALL JxSpack(fileName,ntuvonline,78)


 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') 'END MODULE chem1_list' 
 
 CLOSE(UNIT=78)


RETURN
END SUBROUTINE initphase

