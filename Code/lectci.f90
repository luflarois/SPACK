SUBROUTINE lectci(ifdin)
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
CHARACTER(LEN=3) :: cie
CHARACTER(LEN=3) :: esp
CHARACTER(LEN=1) :: caer
DIMENSION imot(nbmot)
INTEGER :: ii,ia
LOGICAL :: dontHaveCO2=.true.

imp=6
!     Loop for reading the input file.

READ(ifdin,*)
nom=''
OPEN(UNIT=78,FILE='chem1_list.f90',STATUS='replace')
WRITE(UNIT=78,FMT='(A)') 'MODULE chem1_list'
WRITE(UNIT=78,FMT='(A)') '  IMPLICIT NONE'
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
 
WRITE(UNIT=78,FMT='(A,A,A,A)') '  CHARACTER(LEN=24),PARAMETER :: chemical_mechanism='&
          ,char(39),chemical_mechanism(1:len_trim(chemical_mechanism)),char(39)
 
WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: maxnspecies= 200'
WRITE(UNIT=78,FMT='(A,I3.3)') '  INTEGER,PARAMETER :: nspecies=',nesp(2)
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  !Name of species '
WRITE(UNIT=78,FMT='(A,I3.3,A)') '  CHARACTER(LEN=8),PARAMETER,'//&
                'DIMENSION(nspecies) :: spc_name=(/ &'
DO ie=1,nesp(2)
  READ(ifdin,'(a)')chdon
  print*,ie,chdon
  nblanc=nbmot
  CALL part(chdon,mot,imot,nmot)
  indaq(ie)=0
  nom(ie)(1:imot(1))=mot(1)(1:imot(1))
  nom_aux(ie)=nom(ie)
  !LFR
  !LSPR
  cwei(ie)=mot(2)
  csou(ie)=mot(3)
  cdry(ie)=mot(4)
  cwet(ie)=mot(5)
  !cpas(ie)=mot(6)
  cfdd(ie)=mot(6)
  !cfut(ie)=mot(7)
  coff(ie)=mot(7)
  !coff(ie)=mot(8)
  ctra(ie)=mot(8)
  hstar(ie)=mot(9)
  f0(ie)=mot(10)
  difrat(ie)=mot(11)
  uv_eq(ie)=mot(12)
  dvj(ie)=mot(13)
  dhr(ie)=mot(14)
  ak0(ie)=mot(15)
  dak(ie)=mot(16)
  aer_eq(ie)=mot(17) !matrix


  !LFR
  inom(ie)=imot(1)
  ii=4-len(trim(nom(ie)))
  IF(ie==1) THEN
    WRITE(UNIT=78,FMT='(A,A,A)') '   '//char(39),trim(nom(ie))//repeat(' ',ii),char(39)//' & !' 
  ELSE
    WRITE(UNIT=78,FMT='(A,A,A)') '     ,'//char(39),trim(nom(ie))//repeat(' ',ii),char(39)//' & !' 
  END IF
  DO k=1,ie-1
    IF (nom(k)(1:inom(k)) == nom(ie)(1:inom(ie))) THEN
      WRITE(*,*)'ERROR: species already initialized: ',nom(ie)
      STOP
    END IF
  END DO
END DO
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A)') '  !Number of each specie   '
DO ie=1,nesp(2)
  WRITE(cie,FMT='(I3.3)') ie
  ii=4-len(trim(nom(ie)))
  WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: '//trim(nom(ie))//&
              repeat(' ',ii)//'='//cie
  IF(trim(nom(ie))=='CO2') dontHaveCO2=.false.
END DO
  ie=ie+1
  WRITE(cie,FMT='(I3.3)') ie
  IF (dontHaveCO2)  WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: '//'CO2'//&
              repeat(' ',1)//'='//cie

 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A)') '  !for memory allocattion: '
 WRITE(UNIT=78,FMT='(A)') '  !This parameters are use for documentation only. '
 WRITE(UNIT=78,FMT='(A)') '  !Use them in a program in substitution of numerical terms.'
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: src = 1 ! source term '
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: ddp = 2 ! dry deposition '
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: wdp = 3 ! wet deposition '
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: fdda = 4! four-dim assimilation '
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: offline = 5! ! off-line emissions: '
 WRITE(UNIT=78,FMT='(A)') '                                  !'//&
                '=1, emission will be read from file'
 WRITE(UNIT=78,FMT='(A)') '                                  '//&
                 '!=0, emission will be calculated during the model simulation (on-line emission)'
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: transport = 6! ! off-line emissions: '
 WRITE(UNIT=78,FMT='(A)') '                                  !'//&
                '=1, yes'
 WRITE(UNIT=78,FMT='(A)') '                                  '//&
                 '!=0, no transport'
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: on = 1'
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER :: off = 0'


 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
 
 WRITE(UNIT=78,FMT='(A)') '  ! spaction(specie,[1=source,2=drydep,3=wetdep,4=fdda,5=offline emission,6=transport]) ]) '
 WRITE(UNIT=78,FMT='(A)') '  INTEGER,PARAMETER,DIMENSION(6,nspecies) :: spc_alloc='//&
              'RESHAPE((/ &'
	       
  DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//csou(ie)//' , '//cdry(ie)//' , '//cwet(ie)//' , ' &
      //cfdd(ie)//' , '//coff(ie)//' , '//ctra(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//csou(ie)//' , '//cdry(ie)//' , '//cwet(ie)//' , ' &
      //cfdd(ie)//' , '//coff(ie)//' , '//ctra(ie)//'    & ! '//trim(nom(ie))//' - '//cie 
   END IF  
  END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /),(/6,nspecies/))'
    
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A,I3.3,A)') '  INTEGER,PARAMETER,DIMENSION(nspecies) :: spc_uveq=(/ &'
 DO ie=1,nesp(2)
  IF(ie/=nesp(2)) THEN
    WRITE(UNIT=78,FMT='(A)') '   '//uv_eq(ie)//', & !'//trim(nom(ie))
  ELSE
    WRITE(UNIT=78,FMT='(A)') '   '//uv_eq(ie)//'  & !'//trim(nom(ie))
  END IF
 END DO
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
 
 !Crossreference between Mechanism and Matrix
 WRITE(UNIT=78,FMT='(A,I3.3,A)') '  INTEGER,PARAMETER,DIMENSION(5) :: spc_matEqv=(/ &'
 mat_eq='0'
 
 DO ie=1,nesp(2)
   DO ia=1,5
      write(caer,FMT='(I1)') ia
      PRINT *,'LFR->ia,caer,ie,aer_eqv: ',ia,'|'//caer//'|',ie,'|'//aer_eq(ie)//'|'
      IF (trim(aer_eq(ie))==trim(caer)) THEN
         write(esp,FMT='(I3.3)') ie
         PRINT *,'LFR->ie,esp:',ie,esp
         mat_eq(ia)=esp
         exit
      END IF
   END DO
 END DO
 
 DO ia=1,4
   WRITE(UNIT=78,FMT='(A)') '   '//mat_eq(ia)//', &'
 END DO 
 WRITE(UNIT=78,FMT='(A)') '   '//mat_eq(5)//' &'
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A)') '!     HENRYS LAW COEFFICIENTS'
 WRITE(UNIT=78,FMT='(A)') '!     Henrys law coefficient'
 WRITE(UNIT=78,FMT='(A)') '!     [KH298]=mole/(l atm)'
 WRITE(UNIT=78,FMT='(A)') '!     Referencias em R. Sander (1999)'
 WRITE(UNIT=78,FMT='(A)') '!     Compilation of Henry Law Constants '
 WRITE(UNIT=78,FMT='(A)') '!     for Inorganic and Organic Species '
 WRITE(UNIT=78,FMT='(A)') '!     of Potential Importance in '
 WRITE(UNIT=78,FMT='(A)') '!     Environmental Chemistry (Version 3) '
 WRITE(UNIT=78,FMT='(A)') '!     http://www.henrys-law.org '
 WRITE(UNIT=78,FMT='(A)') '!     * indica artigos nao encontrados nesse endereço eletronico'
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: hstar=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//hstar(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//hstar(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /)'
    
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '! [1] Noziere B. et al. The uptake of methyl vinyl ketone, methacrolein,'
 WRITE(UNIT=78,FMT='(A)') '! and 2-methyl-3-butene-2-olonto sulfuric acid solutions,Journal of Physical'
 WRITE(UNIT=78,FMT='(A)') '! Chemistry A, Vol.110, No.7, 2387-2395, 2006.'
 WRITE(UNIT=78,FMT='(A)') '! [2] Abraham M. H. et al. Partition of compounds from gas to water and'
 WRITE(UNIT=78,FMT='(A)') '! from gas to physiological saline at 310K: Linear free energy relationships,'
 WRITE(UNIT=78,FMT='(A)') '! elsevier, 2006.'
 

 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: f0=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//f0(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//f0(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /)'
    
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: difrat=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//difrat(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//difrat(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /)'


 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A)') '!     DIFFUSION COEFFICIENTS'
 WRITE(UNIT=78,FMT='(A)') '!     [DV]=cm2/s (assumed: 1/SQRT(molar mass) when not known)'
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: dvj=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//dvj(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//dvj(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /)'



 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '!     -DH/R (for temperature correction)'
 WRITE(UNIT=78,FMT='(A)') '!     [-DH/R]=K'
 WRITE(UNIT=78,FMT='(A)') '!     Referencias em R. Sander (1999)'
 WRITE(UNIT=78,FMT='(A)') '!     Compilation of Henry Law Constants'
 WRITE(UNIT=78,FMT='(A)') '!     for Inorganic and Organic Species '
 WRITE(UNIT=78,FMT='(A)') '!     of Potential Importance in '
 WRITE(UNIT=78,FMT='(A)') '!     Environmental Chemistry (Version 3)'
 WRITE(UNIT=78,FMT='(A)') '!     http://www.henrys-law.org '
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: dhr=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//dhr(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//dhr(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /)'

 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
!----------
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: weight=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//cwei(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//cwei(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
!----------
!----------srf
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: init_ajust=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//'1.0'//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//'1.0'//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
!----------
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: emiss_ajust=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//'1.0'//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//'1.0'//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '
!----------



 WRITE(UNIT=78,FMT='(A)') '!    ACID DISSOCIATION CONSTANT AT 298K '
 WRITE(UNIT=78,FMT='(A)') '!     [mole/liter of liquid water]'
 WRITE(UNIT=78,FMT='(A)') '!     Referencias: Barth et al. JGR 112, D13310 2007'
 WRITE(UNIT=78,FMT='(A)') '!     Martell and Smith, 1976, Critical stability'
 WRITE(UNIT=78,FMT='(A)') '!     vol1-4 Plenum Press New York'
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: ak0=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//ak0(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//ak0(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '   /)' 
 WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '

 WRITE(UNIT=78,FMT='(A)') '!     Temperature correction factor for'
 WRITE(UNIT=78,FMT='(A)') '!     acid dissociation constants'
 WRITE(UNIT=78,FMT='(A)') '!     [K]'
 WRITE(UNIT=78,FMT='(A)') '!     Referencias: Barth et al. JGR 112, D13310 2007'
 WRITE(UNIT=78,FMT='(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: dak=(/&'
 DO ie=1,nesp(2)
    WRITE(cie,FMT='(I3.3)') ie
    IF(ie<nesp(2)) THEN
      WRITE(UNIT=78,FMT='(A)') &
      '    '//dak(ie)//' ,   & ! '//trim(nom(ie))//' - '//cie
    ELSE
      WRITE(UNIT=78,FMT='(A)') &
      '    '//dak(ie)//'     & ! '//trim(nom(ie))//' - '//cie
    END IF  
 END DO
 
 WRITE(UNIT=78,FMT='(A)') '    /)'

WRITE(UNIT=78,FMT='(A)') '  '
 WRITE(UNIT=78,FMT='(A)') '  '



! change MP 29/01/2008: the following lines can be removed
!READ(ifdin,*)
!DO ie=nesp(2)+1,nesp(2)+nesp(3)
! READ(ifdin,'(a)')chdon
!     write(imp,'(a)')chdon
!  nblanc=nbmot
!  CALL part(chdon,mot,imot,nmot)
!  indaq(ie)=1
!  nom(ie)(1:imot(1))=mot(1)(1:imot(1))
!  inom(ie)=imot(1)
!  DO k=1,ie-1
!    IF (nom(k)(1:inom(k)) == nom(ie)(1:inom(ie))) THEN
!      WRITE(*,*)'ERROR: species already initialized 1: ',nom(ie)
!      STOP
!    END IF
!  END DO
!END DO
! end change MP
RETURN
END SUBROUTINE lectci
