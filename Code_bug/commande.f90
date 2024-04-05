PROGRAM commande
 
! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-28  Time: 12:55:30
 
!------------------------------------------------------------------------

!     -- DESCRIPTION

!     Master code for SPACK:
!     Simplified Preprocessor for Atmospheric Chemical Kinetics.

!     The input files are:
!     - a file describing the mechanism in a symbolic way,
!     - a file of chemical species.

!     The output files are F77 routines describing the chemical
!     production terms.
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
  COMMON/comprec/bpsave(4,nrmax)


  DIMENSION y0(nespmax)


  !     Initialization of data
  !     NEQ has the right dimension.

  CALL lectdata(y0,neq,indicaq)

  DO i=1,nrmax
    bpsave(1,i) =bp(1,i)
  END DO

  !     call precalcul

!LFR>   OPEN(45,FILE='Chemical.f90',STATUS='new')
!LFR>   nwrite=45
!LFR>   WRITE(nwrite,FMT='(A)') 'MODULE Chemical'
!LFR>   WRITE(nwrite,FMT='(A)') '  IMPLICIT NONE'
!LFR>   WRITE(nwrite,FMT='(A)') '  '
!LFR>   CLOSE(UNIT=45)
!LFR>   CALL SYSTEM('cat SpecName.f90 >> Chemical.f90')
!LFR> 
!LFR>   OPEN(45,FILE='Chemical.f90',STATUS='old',POSITION='append')
!LFR>   WRITE(nwrite,FMT='(A)') '  '
!LFR>   WRITE(nwrite,FMT='(A)') '  CONTAINS'
!LFR>   WRITE(nwrite,FMT='(A)') '  '
!LFR>   
!LFR>   WRITE(nwrite,FMT='(A)') ' SUBROUTINE Make_interp_spec(ns,name_sp,iintersp)'
!LFR>   WRITE(nwrite,FMT='(A)') '   IMPLICIT NONE'
!LFR>   WRITE(nwrite,FMT='(A)') '   INTEGER,INTENT(IN) :: ns'
!LFR>   WRITE(nwrite,FMT='(A)') '   CHARACTER(LEN=*),DIMENSION(ns),INTENT(IN) ::  name_sp'
!LFR>   WRITE(nwrite,FMT='(A)') '   INTEGER,DIMENSION(ns),INTENT(OUT) :: iintersp'
!LFR>   WRITE(nwrite,FMT='(A)') '   INTEGER :: i,j'
!LFR>   WRITE(nwrite,FMT='(A)') '   '
!LFR>   WRITE(nwrite,FMT='(A)') '   iintersp=0'
!LFR>   WRITE(nwrite,FMT='(A)') '   '
!LFR>   WRITE(nwrite,FMT='(A)') '   DO i=1,ns'
!LFR>   WRITE(nwrite,FMT='(A)') '     DO j=1,nspec'
!LFR>   WRITE(nwrite,FMT='(A)') '       IF(TRIM(name_sp(i))==TRIM(spname(j))) THEN'
!LFR>   WRITE(nwrite,FMT='(A)') '	   iintersp(i)=j'
!LFR>   WRITE(nwrite,FMT='(A)') '	   EXIT'
!LFR>   WRITE(nwrite,FMT='(A)') '	 END IF'
!LFR>   WRITE(nwrite,FMT='(A)') '      END DO'
!LFR>   WRITE(nwrite,FMT='(A)') '    END DO'
!LFR>   WRITE(nwrite,FMT='(A)') '         '
!LFR>   WRITE(nwrite,FMT='(A)') '   END SUBROUTINE Make_interp_spec'
!LFR>   WRITE(nwrite,FMT='(A)') '         '
!LFR>   WRITE(nwrite,FMT='(A)') '         '
!LFR> 
!LFR>   CLOSE(UNIT=45)

!LFR>   CALL SYSTEM('cat kinetic.f90 >> Chemical.f90')
!LFR>   CALL SYSTEM('cat fexchem.f90 >> Chemical.f90')
!LFR>   CALL SYSTEM('cat jacdchemdc.f90 >> Chemical.f90')
!LFR> !  CALL SYSTEM('cat fexloss.f90 >> Chemical.f90')
!LFR> !  CALL SYSTEM('cat fexprod.f90 >> Chemical.f90')
!LFR>   CALL SYSTEM('cat rates.f90 >> Chemical.f90')
!LFR>   CALL SYSTEM('cat dratedc.f90 >> Chemical.f90')
!LFR>   CALL SYSTEM('echo "END MODULE Chemical" >> Chemical.f90')
!LFR>   CALL SYSTEM('rm -f kinetic.f90 fexchem.f90 jacdchemdc.f90 SpecName.f90')
!LFR>   CALL SYSTEM('rm -f fexloss.f90 fexprod.f90 rates.f90 dratedc.f90')
!LFR>   CALL SYSTEM('rm -f fexloss.f90 fexprod.f90')

END PROGRAM commande







