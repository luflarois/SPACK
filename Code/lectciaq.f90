SUBROUTINE lectciaq(ifdin)
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Read chemical species.
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

!     Michel Pirre, 28/01/2008  change to create chem1aq_list
!
!------------------------------------------------------------------------

!     CHDON:  string variable corresponding to one line of the file.
!     MOT: array composed of the words in CHDON.
!     IMOT: array of sizes for words of MOT.
!     NMOT: number of words.
IMPLICIT DOUBLE PRECISION (a-h,o-z)
INCLUDE 'parametre'
INCLUDE 'ficcom'
INCLUDE 'auxnom' !LFR


INTEGER, INTENT(IN OUT)                   :: ifdin
INTEGER, PARAMETER :: nbmot=100

COMMON/nblanc/nblanc

CHARACTER (LEN=500) :: chdon
CHARACTER (LEN=500) :: mot(nbmot)
CHARACTER(LEN=3) :: cie, cieg
CHARACTER(LEN=2) :: caq
DIMENSION imot(nbmot)
INTEGER :: ii
DIMENSION igas(nespmax), inomaq(nespmax)
CHARACTER(LEN=500) :: nomaq(nespmax)
CHARACTER(LEN=20) :: cacco(nespmax)
imp=6
!     Loop for reading the input file.

READ(ifdin,*)
nomaq=''
OPEN(UNIT=178,FILE='chem1aq_list.f90',STATUS='replace')
WRITE(UNIT=178,FMT='(A)') 'MODULE chem1aq_list'
WRITE(UNIT=178,FMT='(A)') '  IMPLICIT NONE'
 WRITE(UNIT=178,FMT='(A)') '  '
 WRITE(UNIT=178,FMT='(A)') '  '
 
!WRITE(UNIT=178,FMT='(A,A,A,A)') '  CHARACTER(LEN=24),PARAMETER :: chemical_mechanismaq='&
!          ,char(39),chemical_mechanism(1:len_trim(chemical_mechanism)),char(39)
 
WRITE(UNIT=178,FMT='(A)') '  INTEGER,PARAMETER :: maxnspeciesaq= 200'
WRITE(UNIT=178,FMT='(A,I3.3)') '  INTEGER,PARAMETER :: nspeciesaq=',nesp(3)
 WRITE(UNIT=178,FMT='(A)') '  '
 WRITE(UNIT=178,FMT='(A)') '  '
 WRITE(UNIT=178,FMT='(A)') '  !Name of speciesaq '
WRITE(UNIT=178,FMT='(A,I3.3,A)') '  CHARACTER(LEN=8),PARAMETER,'//&
                'DIMENSION(nspeciesaq) :: spcaq_name=(/ &'
caq='aq'
DO ie=1,nesp(3)
  READ(ifdin,'(a)')chdon
  print*,ie,chdon
  nblanc=nbmot
  CALL part(chdon,mot,imot,nmot)
  indaq(ie)=0
  nomaq(ie)(1:imot(1))=mot(1)(1:imot(1))
!  nomaq_aux(ie)=nomaq(ie)
  cacco(ie)=mot(2)
  inomaq(ie)=imot(1)
  ii=4-len(trim(nomaq(ie)))
  IF(ie==1) THEN
    WRITE(UNIT=178,FMT='(A,A,A)') '   '//char(39),trim(nomaq(ie))//caq//repeat(' ',ii),char(39)//' & !' 
  ELSE
    WRITE(UNIT=178,FMT='(A,A,A)') '     ,'//char(39),trim(nomaq(ie))//caq//repeat(' ',ii),char(39)//' & !' 
  END IF
  DO k=1,ie-1
    IF (nomaq(k)(1:inomaq(k)) == nomaq(ie)(1:inomaq(ie))) THEN
      WRITE(*,*)'ERROR: aqueous species already initialized: ',nomaq(ie)
      STOP
    END IF
  END DO
!  END IF
  itesaq=0
  DO k=1,nesp(2)
!     IF (nom(k)(1:inom(k)) == nomaq(ie)(1:inomaq(ie))) THEN
      IF (nom_aux(k) == nomaq(ie)(1:inomaq(ie))) THEN
      itesaq=itesaq+1
      igas(ie)=k
     END IF
  END DO
  IF (itesaq==0)THEN
     WRITE(*,*)'ERROR: aqueous species has no corresponding gas species: ',nomaq(ie)
     STOP
  END IF
  IF (itesaq > 1)THEN
     WRITE(*,*)'ERROR: gas species already initialized: ',nomaq(ie)
     STOP
  END IF  
END DO
 WRITE(UNIT=178,FMT='(A)') '   /)' 
 WRITE(UNIT=178,FMT='(A)') '  '
 WRITE(UNIT=178,FMT='(A)') '  '

 WRITE(UNIT=178,FMT='(A)') '  !Number of each specie   '
DO ie=1,nesp(3)
  WRITE(cie,FMT='(I3.3)') ie
  ii=4-len(trim(nomaq(ie)))
  WRITE(UNIT=178,FMT='(A)') '  INTEGER,PARAMETER :: '//trim(nomaq(ie))//caq//&
              repeat(' ',ii)//'='//cie
END DO

 WRITE(UNIT=178,FMT='(A)') '  '
 WRITE(UNIT=178,FMT='(A)') '  '

 WRITE(UNIT=178,FMT='(A)') '!     number of the corresponding gaseous species'
 WRITE(UNIT=178,FMT='(A)') '  INTEGER,PARAMETER,DIMENSION(nspeciesaq) :: ind_gas=(/&'
 DO ie=1,nesp(3)
    WRITE(cie,FMT='(I3.3)') ie
    WRITE(cieg,FMT='(I3.3)')igas(ie)
    IF(ie<(nesp(3))) THEN
      WRITE(UNIT=178,FMT='(A)') &
      '    '//cieg//' ,   & ! '//trim(nomaq(ie))//caq//' - '//cie
    ELSE
      WRITE(UNIT=178,FMT='(A)') &
      '    '//cieg//'     & ! '//trim(nomaq(ie))//caq//' - '//cie
    END IF  
 END DO
 WRITE(UNIT=178,FMT='(A)') '    /)'
 WRITE(UNIT=178,FMT='(A)') '  '
 
 WRITE(UNIT=178,FMT='(A)') '! accomodation coefficient'
 WRITE(UNIT=178,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspeciesaq) :: acco=(/&'
 DO ie=1,nesp(3)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<(nesp(3))) THEN
      WRITE(UNIT=178,FMT='(A)') &
      '    '//cacco(ie)//' ,   & ! '//trim(nomaq(ie))//caq//' - '//cie
    ELSE
      WRITE(UNIT=178,FMT='(A)') &
      '    '//cacco(ie)//'     & ! '//trim(nomaq(ie))//caq//' - '//cie
    END IF  
 END DO
 WRITE(UNIT=178,FMT='(A)') '    /)'
    
 WRITE(UNIT=178,FMT='(A)') '  '
 WRITE(UNIT=178,FMT='(A)') 'END MODULE chem1aq_list' 
 CLOSE(UNIT=178)
RETURN
END SUBROUTINE lectciaq
