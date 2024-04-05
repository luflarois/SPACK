SUBROUTINE initconst
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Initialization of constants.
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
!     Michel Pirre, 28/01/2008 : change to create chem1_lisaq

!------------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'

av = 6.022D23

RETURN
END SUBROUTINE initconst

!------------------------------------------------------------------------

SUBROUTINE lectdata(y0,neq,indicaq)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Read data.
!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!     SLUMP: lumped stoichiometric matrix.
!     XLPHY: physical lumping.
!     INDPUR: ii=indpur(i,j) true label of J-th  species in lumping I.
!     Y0: initial conditions.
!     S: stoichiometrix matrix.
!     NALG: physical algebraic onstraints.

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
INCLUDE 'nficfort'
INCLUDE 'auxnom'

DOUBLE PRECISION, INTENT(IN OUT)         :: y0(nespmax)
INTEGER, INTENT(OUT)                     :: neq
INTEGER, INTENT(OUT)                     :: indicaq
INTEGER :: ntuvonline
CHARACTER (LEN=20) :: fd
CHARACTER(LEN=100) :: fileName



ipiste=15
nmaster=ipiste
ipiste=ipiste+1
!     Initialization of constants and main arrays.
CALL initconst
CALL initcinet

nalg=0

!     Read master file.
OPEN(nmaster,FILE='inSPACK',STATUS='old')
WRITE(*,*)'-------------Input files-----------'

!     Read photolysis rates option
READ(nmaster,*)
READ(nmaster,*) ntuvonline

READ(nmaster,*)
iunit = 0
READ(nmaster,*)

!     Chemical mechanism
READ(nmaster,*)fd
WRITE(6,*)'File for chemical mechanism: ',fd
filemeca=fd
print*,'filemeca=',filemeca
do i=len_trim(filemeca),1,-1
  if(filemeca(i:i)=='/') exit
enddo
chemical_mechanism=trim(filemeca(i+1:len_trim(filemeca)))
print*, 'chemical_mechanism=',chemical_mechanism

CALL openfic(ipiste,ifdth,fd,0)
READ(nmaster,*)
READ(nmaster,*)
READ(nmaster,*)fd
!      read(nmaster,*)
READ(nmaster,*)
READ(nmaster,*)
READ(nmaster,*)fileName !To FastJX
print *, 'FastJx adapter Mechanism: ',fileName

CLOSE(nmaster)
!     Chemical species
WRITE(6,*)'File for chemical species: ',fd
filespecies=fd
CALL openfic(ipiste,ifdin,fd,0)

!     Modif BS multiphase: units molec. cm-3 (ippb=0) ou ppb (1)

READ(ifdin,*)
READ(ifdin,*)
!     read(Ifdin,*)ippb
!     read(Ifdin,*)
!     End modif BS

ippb = 0
READ(ifdin,*)nesp(2),nesp(3)

nesp(1)=nesp(2)+nesp(3)
indicaq=0
IF (nesp(3) > 0) indicaq=1
indicaqcom=indicaq
WRITE(6,*)'Number of multiphase species (nesp): ',nesp(1)
WRITE(6,*)'Number of gas-phase species: ',nesp(2)
WRITE(6,*)'Number of aqueous-phase species: ',nesp(3)
WRITE(6,*)'Max number of species (nespmax): ',nespmax
WRITE(6,*)'Max number of reactions (nrmax): ',nrmax
IF (nesp(1) > nespmax) THEN
  WRITE(*,*)'ERROR: dimension, nesp>nespmax'
  CALL halte
END IF
CALL lectci(ifdin)

!     ..v.7..x....v....x....v....x....v....x....v....x....v....x....v....x.I
!     PP 12 02 2002
!     Check species name
DO ie = 1,nesp(1)
  WRITE(6,777)ie,nom(ie)
END DO
777  FORMAT(2X,i3,2X,a10)

! change MP 29/01/08 to create chem1_listaq for aqueous species
IF (nesp(3) > 0)THEN
CALL lectciaq(ifdin)
!DO ie = 1,nesp(1)
!  WRITE(6,777)ie,nom(ie)
!END DO
END IF
!end change MP

!     ..v.7..x....v....x....v....x....v....x....v....x....v....x....v....x.I
!     Read chemical mechanism

WRITE(*,*)
WRITE(*,*)'-----------Chemical mechanism---------'
WRITE(*,*)
CALL lectcinet(ifdth,indicaq,ntuvonline,fileName)
!     Modif BS for gas-phase only
!     call convcinet(iunit,indicaq)


!     Modif BS for gas-phase only
!     call inition
!     call initlphy
!     call initphot
!     End modif gas-phase only.

!     Dimension

neq=ndiff(1)*nbrem

100  FORMAT(a10)

RETURN
END SUBROUTINE lectdata

SUBROUTINE openfic(ipiste,ifd,fd,inew)



INTEGER, INTENT(IN OUT)                  :: ipiste
INTEGER, INTENT(OUT)                     :: ifd
CHARACTER (LEN=20), INTENT(IN)           :: fd
INTEGER, INTENT(IN)                  :: inew

ifd =ipiste
ipiste=ipiste+1
IF (inew == 0) OPEN(ifd,FILE = fd,STATUS = 'old')
IF (inew == 1) OPEN(ifd,FILE = fd,STATUS = 'new')
RETURN
END SUBROUTINE openfic

SUBROUTINE convcinet(iunit,indicaq)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Conversion of units: define conversion factors.
!     Change units for concentrations (kinetic rates): debug(1) ((2),(3))
!     Initial conditions in molec.cm-3

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



INTEGER, INTENT(IN OUT)                  :: iunit
INTEGER, INTENT(IN OUT)                  :: indicaq

!     Conversion factors ppm <-> molec/cm3

debug(1) = 2.46D10
debug(2) = 2.46D13
debug(3) = 6.05D26

IF (iunitgas == -1) THEN
  WRITE(*,*)'ERROR:  gas units undefined.'
  STOP
END IF

IF ((indicaq == 1).AND.(iunitaq == -1)) THEN
  WRITE(*,*)'ERROR: liquid units undefined'
  STOP
END IF

IF (iunitgas /= iunit) THEN
  IF ((iunitgas == 0).AND.(iunit == 1)) THEN
    debug2=1.d0/debug(1)
    debug3=1.d0/debug(1)**2
  ELSE IF ((iunitgas == 1).AND.(iunit == 0)) THEN
    debug2=debug(1)
    debug3=debug(1)**2
  ELSE
    WRITE(*,*)'ERROR: gas units conversion'
    STOP
  END IF
END IF

IF ((iunit == 1).AND.(iunitaq /= -1)) THEN
  WRITE(*,*)'ERROR: liquid units conversion'
  STOP
END IF

RETURN
END SUBROUTINE convcinet
!------------------------------------------------------------------------

SUBROUTINE convci(ippb,iunit,nn,y0)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Conversion of units for initial conditions.
!     Change units for concentrations (kinetic rates): debug(1) ((2),(3))

!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     IPPB: 1 (0) if initial conditions in ppb (mplec/cm3)
!     IUNIT: 1(2) if computation in molec/cm3 (ppm)
!     NN: number of multiphase species

!     -- INPUT/OUTPUT VARIABLES

!     Y0: initial conditions

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

INTEGER, INTENT(IN OUT)                  :: ippb
INTEGER, INTENT(IN OUT)                  :: iunit
INTEGER, INTENT(IN)                      :: nn
DOUBLE PRECISION, INTENT(OUT)            :: y0(nespmax)



!     Conversion factor for concentrations: debug(1)
!     Conversion factor for kinetic rates: debug(2) and debug(3)
!     Initial condtions in molec/cm3

!     Conversion ppm <-> molec/cm3

debug(1) = 2.46D10
debug(2) = 2.46D13
debug(3) = 6.05D26

IF (iunit /= ippb) THEN
  WRITE(*,*)'Initial concentrations without unit conversion'
  IF (iunit == 0) THEN
    IF (ippb == 1) THEN
      WRITE(*,*) 'concentrations initiales en ppb converties mol/cm3'
      
      DO i=1,nn
        y0(i) = y0(i) * debug(1)
      END DO
    ELSE
      WRITE(*,*)'ERROR: initial concentrations w/o conversion'
      STOP
    END IF
  ELSE IF (iunit == 1) THEN
    IF (ippb == 0) THEN
      WRITE(*,*)'Initial concentrations in mol/cm3 -> ppb'
      DO i=1,nn
        y0(i) = y0(i) / debug(1)
      END DO
    ELSE
      WRITE(*,*) 'ERROR: no conversion for initial concentrations'
      STOP
    END IF
  ELSE
    WRITE(*,*)'ERROR: no conversion for initial concentrations'
    STOP
  END IF
END IF

RETURN
END SUBROUTINE convci
!--------------------------------------

SUBROUTINE halte
IMPLICIT DOUBLE PRECISION (a-h,o-z)
STOP
RETURN
END SUBROUTINE halte
!--------------------------------------
!     subroutine precalcul
!cccccccccccccccccccccCCCcccccccccccccc
!     routine de precalcul des parametres
!     a partir des donnees d'entree
!cccccccccccccccccccccccccccccccccccccc
!     implicit REAL (a-h,o-z)
!     include 'parametre'
!     include 'ficcom'
!     common/comprec/bpsave(4,nrmax)

!cccc
!     preprocess des ctes cinetiques
!cccc
!     do i=1,nr
!     if (iprecalc(i).eq.2) then
!     bp(1,i)=bpsave(1,i)*dexp(-bp(2,i)/bp(3,i))
!     elseif (iprecalc(i).eq.3) then
!     bp(1,i)=bpsave(1,i)*dexp(-bp(3,i)/bp(4,i))
!     endif
!     enddo

!     return
!     end
