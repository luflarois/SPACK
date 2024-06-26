MODULE chem1_list
  IMPLICIT NONE
  
  
  CHARACTER(LEN=24),PARAMETER :: chemical_mechanism='ANALY'
  INTEGER,PARAMETER :: maxnspecies= 200
  INTEGER,PARAMETER :: nspecies=002
  
  
  !Name of species 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies) :: spc_name=(/ &
   'O3  ' & !
     ,'CO  ' & !
   /)
  
  
  !Number of each specie   
  INTEGER,PARAMETER :: O3  =001
  INTEGER,PARAMETER :: CO  =002
  
  
  !for memory allocattion: 
  !This parameters are use for documentation only. 
  !Use them in a program in substitution of numerical terms.
  INTEGER,PARAMETER :: src = 1 ! source term 
  INTEGER,PARAMETER :: ddp = 2 ! dry deposition 
  INTEGER,PARAMETER :: wdp = 3 ! wet deposition 
  INTEGER,PARAMETER :: fdda = 4! four-dim assimilation 
  INTEGER,PARAMETER :: offline = 5! ! off-line emissions: 
                                  !=1, emission will be read from file
                                  !=0, emission will be calculated during the model simulation (on-line emission)
  INTEGER,PARAMETER :: transport = 6! ! off-line emissions: 
                                  !=1, yes
                                  !=0, no transport
  INTEGER,PARAMETER :: on = 1
  INTEGER,PARAMETER :: off = 0
  
  
  ! spaction(specie,[1=source,2=drydep,3=wetdep,4=fdda,5=offline emission,6=transport]) ]) 
  INTEGER,PARAMETER,DIMENSION(6,nspecies) :: spc_alloc=RESHAPE((/ &
    0 , 1 , 1 , 1 , 0 , 1 ,   & ! O3 - 001
    1 , 1 , 1 , 1 , 1 , 1    & ! CO - 002
    /),(/6,nspecies/))
  
  
  INTEGER,PARAMETER,DIMENSION(nspecies) :: spc_uveq=(/ &
   09, & !O3
   00  & !CO
   /)
  
  
!     HENRYS LAW COEFFICIENTS
!     Henrys law coefficient
!     [KH298]=mole/(l atm)
!     Referencias em http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!     * indica artigos nao encontrados nesse endere�o eletronico
  REAL,PARAMETER,DIMENSION(nspecies) :: hstar=(/&
    1.20E-2              ,   & ! O3 - 001
    9.90E-4                  & ! CO - 002
    /)
  
! [1] Noziere B. et al. The uptake of methyl vinyl ketone, methacrolein,
! and 2-methyl-3-butene-2-olonto sulfuric acid solutions,Journal of Physical
! Chemistry A, Vol.110, No.7, 2387-2395, 2006.
! [2] Abraham M. H. et al. Partition of compounds from gas to water and
! from gas to physiological saline at 310K: Linear free energy relationships,
! elsevier, 2006.
  REAL,PARAMETER,DIMENSION(nspecies) :: f0=(/&
    0.1 ,   & ! O3 - 001
    0.0     & ! CO - 002
    /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: difrat=(/&
    1.6 ,   & ! O3 - 001
    0.0     & ! CO - 002
    /)
  
  
!     DIFFUSION COEFFICIENTS
!     [DV]=cm2/s (assumed: 1/SQRT(molar mass) when not known)
  REAL,PARAMETER,DIMENSION(nspecies) :: dvj=(/&
    0.1470000	    ,   & ! O3 - 001
    0.1890000	        & ! CO - 002
    /)
  
  
!     -DH/R (for temperature correction)
!     [-DH/R]=K
  REAL,PARAMETER,DIMENSION(nspecies) :: dhr=(/&
    2500.         ,   & ! O3 - 001
    1300.             & ! CO - 002
    /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: weight=(/&
    48.  ,   & ! O3 - 001
    28.      & ! CO - 002
   /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: init_ajust=(/&
    1.0 ,   & ! O3 - 001
    1.0     & ! CO - 002
   /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: emiss_ajust=(/&
    1.0 ,   & ! O3 - 001
    1.0     & ! CO - 002
   /)
  
  
!    ACID DISSOCIATION CONSTANT AT 298K 
!     [mole/liter of liquid water]
!     Referencias: Barth et al. JGR 112, D13310 2007
!     Martell and Smith, 1976, Critical stability
!     vol1-4 Plenum Press New York
  REAL,PARAMETER,DIMENSION(nspecies) :: ak0=(/&
    0.00E+00	    ,   & ! O3 - 001
    0.00E+00	        & ! CO - 002
   /)
  
  
!     Temperature correction factor for
!     acid dissociation constants
!     [K]
!     Referencias: Barth et al. JGR 112, D13310 2007
  REAL,PARAMETER,DIMENSION(nspecies) :: dak=(/&
    0.         ,   & ! O3 - 001
    0.             & ! CO - 002
    /)
  
  
  INTEGER,PARAMETER :: nr   =001!Number of gas-phase reactions
  INTEGER,PARAMETER :: nrt  =001!Total Number of reactions
  INTEGER,PARAMETER :: nrh2o=000!Number of aqueous-phase reactions
  
 !------------------------------------------------------------------------------  
  
 ! Photolysis Rate Calculation: method used (LUT=look_up_table, FAST-JX= on-line)
  
  CHARACTER(LEN=10),PARAMETER :: PhotojMethod= 'LUT'
  INTEGER,PARAMETER :: maxJcomb=    5, nr_photo=    1
  INTEGER,PARAMETER,DIMENSION(nr_photo) :: nfactors=(/ &
                            0/)!  1
 
  DOUBLE PRECISION,PARAMETER,DIMENSION(maxJcomb,nr_photo) :: factor=RESHAPE((/ &
                                                0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00/),(/maxJcomb,nr_photo/))
  CHARACTER(LEN=07 ),PARAMETER,DIMENSION(maxJcomb,nr_photo) :: JReactionComp=RESHAPE((/ &
                                               "NONE   ","NONE   ","NONE   ","NONE   ","NONE   "/),(/maxJcomb,nr_photo/))

  
END MODULE chem1_list
