SUBROUTINE JxSpack(fileName,ntuvonline,nOutFile)
   IMPLICIT NONE
   CHARACTER(LEN=*),INTENT(IN) :: fileName !Input file name
   INTEGER  :: nReactions !Total of reactions 
   INTEGER,INTENT(IN) :: nOutFile !File to write the output reactions
   INTEGER,INTENT(IN) :: ntuvonline
   INTEGER,PARAMETER :: nFile=22 !#file to read

   INTEGER,PARAMETER :: maxJcomb=5 !Maximum of J reactions combinations
   INTEGER,PARAMETER :: sizeOfSpecieName=7 !Size of specie name
   INTEGER,PARAMETER :: maxnReactions=200!Maximum of reactions
   CHARACTER(LEN=80) :: header
   CHARACTER(LEN=20) :: title
   CHARACTER(LEN=3) :: sp,ss
   CHARACTER(LEN=80) :: fmt_txt
   INTEGER :: nOfReaction !Number od reaction
   INTEGER,DIMENSION(maxnReactions) :: specProd !Total of species production by reaction
   INTEGER :: i,n

   TYPE tf
     double precision :: factor !Factor of production
     CHARACTER(LEN=sizeOfSpecieName) :: specie !Name of specie    
   END TYPE tf
   TYPE(tf),DIMENSION(maxnReactions,maxJcomb) :: spcData

 
  IF(ntuvonline > 0) then 
      OPEN(UNIT=nFile,FILE=fileName)
      !Reading the 5 first lines - header
      DO i=1,5
 	READ(nFile,FMT='(A80)') header
      END DO
      !Reading the reactions and species' production factors
      spcData(:,:)%specie='NONE'
      DO i=1,maxnReactions
 	READ(nFile,FMT=*,end=110) nOfReaction,specProd(i),(spcData(i,n),n=1,specProd(i))
 	!WRITE (*,*) i,nOfReaction,specProd(i),(spcData(i,n),n=1,specProd(i))
 	!pause
      END DO
   END IF
   110 continue 
   
   IF(ntuvonline == 0) then 
      nReactions=i-1  
      CLOSE(UNIT=22)
      title='FAST-JX'
   ELSEIF(ntuvonline == 0) then 
      nReactions=1
      spcData(:,:)%specie='NONE'
      specProd(:) = 0
      spcData(:,:)%factor = 0.
      title='LUT'
   ELSEIF(ntuvonline == 2) then 
      nReactions=1
      spcData(:,:)%specie='NONE'
      specProd(:) = 0
      spcData(:,:)%factor = 0.
      title='FAST-TUV'
   ELSE
      print*,'Invalid photojmethod - ntuvonline=',ntuvonline
      stop 100
   
   ENDIF
   
   !not oficial bug fix
   SELECT CASE(ntuvonline)
   CASE(0)
        title = 'LUT'
   CASE(1)
        title = 'FAST-JX'
   CASE(2)
        title = 'FAST-TUV'
   END SELECT

   !Writing the adapt fastJX to Spack routine

   !First part: An array with number of species from Jx which makes a reaction   
   WRITE(nOutFile,FMT='(A)') ' !------------------------------------------------------------------------------  '
   WRITE(nOutFile,FMT='(A)') '  '
   WRITE(nOutFile,FMT='(A)') ' ! Photolysis Rate Calculation: method used (LUT=look_up_table, FAST-JX= on-line)'
   WRITE(nOutFile,FMT='(A)') '  '

   WRITE(nOutFile,FMT='(A,A,A,A)') '  CHARACTER(LEN=10),PARAMETER :: PhotojMethod= ',char(39),trim(title),char(39)
   WRITE(nOutFile,FMT='(A,i5,A,i5)') '  INTEGER,PARAMETER :: maxJcomb=',maxJcomb, ', nr_photo=',nReactions
   WRITE(nOutFile,FMT='(A)') '  INTEGER,PARAMETER,DIMENSION(nr_photo) :: nfactors=(/ &'
   DO i=1,nReactions-1      
   WRITE(nOutFile,FMT='(A,I2,A,I3)') '                           ',specProd(i),', &!',i
   END DO
   WRITE(nOutFile,FMT='(A,I2,A,I3)') '                           ',specProd(nReactions),'/)!',i
   WRITE(nOutFile,FMT='(A)') ' '
   WRITE(sp,FMT='(I2.2)') maxJcomb  

  !Second Part: An array with multiplication factor for each species from Jx
  WRITE(nOutFile,FMT='(A,I2.2,A)') '  DOUBLE PRECISION,PARAMETER,DIMENSION(maxJcomb,nr_photo) :: factor=RESHAPE((/ &'
   DO i=1,nReactions-1
      WRITE(nOutFile,FMT='(A,$)')   '                                               '    
      DO n=1,maxJcomb   
         WRITE(nOutFile,FMT='(D10.3,",",$)')  spcData(i,n)%factor
      END DO
      WRITE(nOutFile,FMT='(A)') ' &'
   END DO
   WRITE(nOutFile,FMT='(A,$)')   '                                               '    
   DO n=1,maxJcomb-1   
      WRITE(nOutFile,FMT='(D10.3,",",$)')  spcData(i,n)%factor
   END DO
   WRITE(nOutFile,FMT='(D10.3,A)') spcData(i,maxJcomb)%factor,'/),(/maxJcomb,nr_photo/))'

   WRITE(ss,FMT='(I2.2)') sizeOfSpecieName

   
   WRITE(nOutFile,FMT='(A,I2.2,A)') '  CHARACTER(LEN='//ss//'),PARAMETER,DIMENSION(maxJcomb,nr_photo) :: JReactionComp=RESHAPE((/ &'    
    DO i=1,nReactions-1
      WRITE(nOutFile,FMT='(A,$)')   '                                               '    
      DO n=1,maxJcomb   
         WRITE(nOutFile,FMT='(A,",",$)')  '"'//spcData(i,n)%specie//'"'
      END DO
      WRITE(nOutFile,FMT='(A)') ' &'
   END DO
   WRITE(nOutFile,FMT='(A,$)')   '                                               '    
   DO n=1,maxJcomb-1   
      WRITE(nOutFile,FMT='(A,",",$)')  '"'//spcData(i,n)%specie//'"'
   END DO
   WRITE(nOutFile,FMT='(A,A)')  '"'//spcData(i,maxJcomb)%specie,'"/),(/maxJcomb,nr_photo/))'

   WRITE(nOutFile,FMT='(A)') ''
    
         
END SUBROUTINE JxSpack

!PROGRAM teste
!
!   CALL JxSpack('RACM_to_FASTJX.dat',23, 50)
!
!END PROGRAM teste
