SUBROUTINE lectcinet(ifdth,indicaq,ntuvonline,fileName)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Read chemical mechanism.
!------------------------------------------------------------------------

!     -- INPUT VARIABLES

!     IFDTH: input file
!     FPAR: physical parameters file

!     -- INPUT/OUTPUT VARIABLES

!     -- OUTPUT VARIABLES

!     NESP: number of species.
!     NR: number of reactions.
!     XLPHY: physical lumping.
!     INDPUR: ii=indpur(i,j) true label of j-th species in i-th lumping.
!     NALG: number of algebraic constraints.

!------------------------------------------------------------------------

!     -- REMARKS

!------------------------------------------------------------------------

!     -- MODIFICATIONS

!------------------------------------------------------------------------

!     -- AUTHOR(S)

!     Bruno Sportisse, CEREA, 2003.

!------------------------------------------------------------------------
IMPLICIT REAL (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'

INTEGER, INTENT(IN OUT)                  :: ifdth
INTEGER, INTENT(IN OUT)                  :: indicaq
INTEGER, INTENT(IN OUT)                  :: ntuvonline
CHARACTER(LEN=*) :: FileName
INTEGER, PARAMETER :: nbmot=100

COMMON/nblanc/nblanc

CHARACTER (LEN=500) :: chdon
CHARACTER (LEN=500) :: mot(nbmot)
DIMENSION imot(nbmot)


imp=6
iattente=0
iarret=0

!     Loop for reading the file.
!     Write the first lines of the output routines
CALL write_header(ntuvonline)

100  READ(ifdth,'(a)')chdon
WRITE(imp,'(a)')chdon

!     Build the sequence and decompose in words.


nblanc=nbmot
CALL part(chdon,mot,imot,nmot)



!     Case of long sequences (two lines).

102  CONTINUE
IF(mot(nmot)(1:2) == '//')THEN
  nchold=nmot
  READ(ifdth,'(a)')chdon
  WRITE(imp,'(a)')chdon
  nblanc=nbmot-nchold+1
  
  
  CALL part(chdon,mot(nchold),imot(nchold),nchplus)
  nmot=nchold-1+nchplus
END IF

IF (mot(nmot)(1:2) == '//') GO TO 102

!     Call the different routines according to the keyword.

IF(mot(1)(1:3) == 'KIN')  THEN
  IF (iattente == 0) THEN
    WRITE(*,*)'ERROR: kinetics before reaction ',nr
    STOP
  ELSE IF (iattente == 1) THEN
    CALL kinreac(mot,imot,nmot,ntuvonline)
  ELSE IF (iattente == 2) THEN
    CALL kindis(mot,imot,nmot)
  ELSE IF (iattente == 3) THEN
    CALL kinhenry(mot,imot,nmot)
  END IF
  iattente=0
  
ELSE IF (mot(1)(1:3) == 'SET') THEN
  
  IF (mot(2)(1:4) == 'UNIT')THEN
    CALL lectunit(mot,imot,nmot)
    
  ELSE IF (mot(2)(1:10) == 'TABULATION') THEN
    CALL lect_tabulation(mot,imot,nmot)
    
  ELSE
    WRITE(*,*)'ERROR: UNKNOWN SET FUNCTIONS'
    STOP
  END IF
!     Symbols for commented lines: %, !,
ELSE IF (mot(1)(1:1) == '%')THEN
ELSE IF (mot(1)(1:1) == '!')THEN
  
!     END.
ELSE IF (mot(1)(1:3) == 'END')THEN
  iarret=1
ELSE IF (iattente /= 0) THEN
  WRITE(*,*)'ERROR: I wait for kinetics'
  STOP
ELSE
  CALL reaction(mot,imot,nmot,iattente)
END IF

IF (iarret == 0) GO TO 100


CALL write_end

!     Write file non_zero.dat

CALL wnonzero(s,nr,jer)

WRITE(*,*)'########################################'
WRITE(*,*)'########################################'
WRITE(6,*)'Summary for the kinetic scheme'
WRITE(*,*)'Total number of reactions =',nr
WRITE(*,*)'Number of photolytic reactions =',nrphot
WRITE(*,*)'Number of dissociation equilibria =',nequil

IF (indicaq == 0) THEN
  WRITE(*,*)'Gas-phase chemistry'
END IF

IF (indicaq == 1) THEN
  WRITE(*,*)'Multiphase (gas-phase and aqueous-phase) chemistry'
END IF

CALL initphase(ntuvonline,fileName)

RETURN
END SUBROUTINE lectcinet

!---------------------------------------------------

SUBROUTINE lectunit(mot,imot,nmot)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Read units.
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
IMPLICIT REAL (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN OUT)      :: mot(nbmot)
INTEGER, INTENT(IN OUT)                  :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot





iunitgas = 0
IF (mot(2)(1:imot(2)) == 'GAS') THEN
  IF (mot(3)(1:imot(3)) == 'MOLCM3') THEN
    iunitgas=0
  ELSE IF (mot(3)(1:imot(3)) == 'PPB') THEN
    iunitgas=1
  ELSE
    WRITE(*,*)'ERROR: unknown units for gas-phase kinetics'
    STOP
  END IF
ELSE IF (mot(2)(1:imot(2)) == 'AQ') THEN
  IF (mot(3)(1:imot(3)) == 'MOLL') THEN
    iunitaq=2
  ELSE
    WRITE(*,*)'ERROR: unknown units for aqueous-phase kinetics'
    STOP
  END IF
END IF
RETURN
END SUBROUTINE lectunit
!---------------------------------------------------

SUBROUTINE lect_tabulation(mot,imot,nmot)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Read tabulation for photolysis.
!     The format is: TABULATION N DEGREES D1 D2 ... DN
!     N is the number of tabulated angles of values Di (in degrees).
!     The sequence may be increasing or decreasing.
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
IMPLICIT REAL (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN OUT)      :: mot(nbmot)
INTEGER, INTENT(IN OUT)                  :: imot(nbmot)
INTEGER, INTENT(IN OUT)                  :: nmot



DIMENSION t(ntabphotmax)

ireversetab=0

CALL entier(ntabphot,mot(3),imot(3))
IF (ntabphot == 0) THEN
  WRITE(*,*)'ERROR: number of tabulated angles to be checked in'
  WRITE(*,*)'subroutine lect-tabulation'
  STOP
END IF
IF (ntabphot > ntabphotmax) THEN
  WRITE(*,*)'ERROR: ntabphot>ntabphotmax'
  STOP
END IF

DO i=1,ntabphot
  CALL reel(t(i),mot(i+4),imot(i+4))
END DO
!     Check increasing order
IF (t(1) > t(2)) THEN
  ireversetab=1
  DO j=1,ntabphot
    tabphot(j)=t(ntabphot+1-j)
  END DO
ELSE
  DO j=1,ntabphot
    tabphot(j)=t(j)
  END DO
END IF
DO j=1,ntabphot-1
  IF (tabphot(j) >= tabphot(j+1)) THEN
    WRITE(*,*) 'ERROR: the tabulation has to be strictly monotonic'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE lect_tabulation
!---------------------------------------------------

SUBROUTINE reaction(mot,imot,nmot,iattente)

!     -- DESCRIPTION

!     Read reactions.
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
!     implicit REAL (a-h,o-z)
INCLUDE 'parametre'
INTEGER, PARAMETER :: nbmot=100

CHARACTER (LEN=500), INTENT(IN OUT)      :: mot(nbmot)
INTEGER, INTENT(IN OUT)                  :: imot(nbmot)
INTEGER, INTENT(IN)                      :: nmot
INTEGER, INTENT(OUT)                     :: iattente




DO i=1,nmot
  IF (mot(i)(1:imot(i)) == '>') THEN
    IF (iattente /= 0) THEN
      WRITE(*,*)'ERROR: too many symbols >'
      STOP
    ELSE
      iattente=1
    END IF
  ELSE IF (mot(i)(1:imot(i)) == '->') THEN
    IF (iattente /= 0) THEN
      WRITE(*,*)'ERROR: too many symbols ->'
      STOP
    ELSE
      iattente=1
    END IF
  ELSE IF (mot(i)(1:imot(i)) == '=') THEN
    IF (iattente /= 0) THEN
      WRITE(*,*)'ERROR: too many symbols ='
      STOP
    ELSE
      iattente=2
    END IF
  ELSE IF (mot(i)(1:imot(i)) == '<H>') THEN
    IF (iattente /= 0) THEN
      WRITE(*,*)'ERROR: too many symbols <H>'
      STOP
    ELSE
      iattente=3
    END IF
  ELSE IF (mot(i)(1:imot(i)) == '=H=') THEN
    IF (iattente /= 0) THEN
      WRITE(*,*)'ERROR: too many symbols =H='
      STOP
    ELSE
      iattente=4
    END IF
  END IF
END DO

IF (iattente == 0) THEN
  WRITE(*,*)'ERROR: lack of symbol for reaction'
  STOP
ELSE IF (iattente == 1) THEN
  CALL creac(mot,imot,nmot)
ELSE IF (iattente == 2) THEN
  CALL cdis(mot,imot,nmot)
ELSE IF (iattente == 3) THEN
  CALL chenry(mot,imot,nmot,0)
ELSE IF (iattente == 4) THEN
  CALL chenry(mot,imot,nmot,1)
  iattente=3
END IF

RETURN
END SUBROUTINE reaction














