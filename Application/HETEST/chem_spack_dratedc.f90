  
 MODULE mod_chem_spack_dratedc
  
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dratedc ! subroutine
 CONTAINS
  
   SUBROUTINE dratedc(rk,y,dw,ngas,ijkbeg,ijkend,maxblock_size,nr)
 
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This routine computes the derivative of reaction  rates.
!     This routine is automatically generated by SPACK.
!     Mechanism: ../Mechanism/HETEST 
!     Species: ../Mechanism/ciHETR 
!
!------------------------------------------------------------------------
!
!     -- INPUT VARIABLES
!
!     RK: kinetic rates.
!     Y: chemical concentrations.
!
!     -- INPUT/OUTPUT VARIABLES
!
!     -- OUTPUT VARIABLES
!
!     DW: derivative of reaction rates wrt Y.
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
 
 
 
     INTEGER	       , INTENT(IN)  :: ngas                      
     INTEGER	       , INTENT(IN)  :: ijkbeg			  
     INTEGER	       , INTENT(IN)  :: ijkend			  
     INTEGER	       , INTENT(IN)  :: maxblock_size		  
     INTEGER	       , INTENT(IN)  :: nr			  
     DOUBLE PRECISION , INTENT(IN)  :: rk(maxblock_size,nr)	  
     DOUBLE PRECISION , INTENT(IN)  :: y(maxblock_size,NGAS) 	  
     DOUBLE PRECISION , INTENT(OUT) :: dw(maxblock_size,nr,NGAS) 
     INTEGER :: ijk						  
 
     DO ijk=ijkbeg,ijkend
  
  
      dw(ijk,  1,  5) =  rk(ijk,  1)
      dw(ijk,  2,  1) =  rk(ijk,  2)
      dw(ijk,  3,  1) =  rk(ijk,  3)
      dw(ijk,  4,  8) =  rk(ijk,  4)
      dw(ijk,  5,  9) =  rk(ijk,  5)
      dw(ijk,  6, 10) =  rk(ijk,  6)
      dw(ijk,  7,  6) =  rk(ijk,  7)
      dw(ijk,  8,  6) =  rk(ijk,  8)
      dw(ijk,  9,  3) =  rk(ijk,  9)
      dw(ijk, 10, 17) =  rk(ijk, 10)
      dw(ijk, 11, 17) =  rk(ijk, 11)
      dw(ijk, 12, 18) =  rk(ijk, 12)
      dw(ijk, 13, 23) =  rk(ijk, 13)
      dw(ijk, 14, 23) =  rk(ijk, 14)
      dw(ijk, 15, 25) =  rk(ijk, 15)
      dw(ijk, 16, 26) =  rk(ijk, 16)
      dw(ijk, 17, 27) =  rk(ijk, 17)
      dw(ijk, 18, 30) =  rk(ijk, 18)
      dw(ijk, 19, 31) =  rk(ijk, 19)
      dw(ijk, 20, 12) =  rk(ijk, 20)
      dw(ijk, 21, 12) =  rk(ijk, 21) * Y(ijk,  1)
      dw(ijk, 21,  1) =  rk(ijk, 21) * Y(ijk, 12)
      dw(ijk, 22, 13) =  rk(ijk, 22)
      dw(ijk, 23, 13) =  rk(ijk, 23)
      dw(ijk, 24, 13) =  rk(ijk, 24)
      dw(ijk, 25, 13) =  rk(ijk, 25) * Y(ijk,  2)
      dw(ijk, 25,  2) =  rk(ijk, 25) * Y(ijk, 13)
      dw(ijk, 26,  1) =  rk(ijk, 26) * Y(ijk, 14)
      dw(ijk, 26, 14) =  rk(ijk, 26) * Y(ijk,  1)
      dw(ijk, 27,  1) =  rk(ijk, 27) * Y(ijk, 15)
      dw(ijk, 27, 15) =  rk(ijk, 27) * Y(ijk,  1)
      dw(ijk, 28, 14) =  rk(ijk, 28) * Y(ijk, 15)
      dw(ijk, 28, 15) =  rk(ijk, 28) * Y(ijk, 14)
      dw(ijk, 29,  3) =  rk(ijk, 29) * Y(ijk, 14)
      dw(ijk, 29, 14) =  rk(ijk, 29) * Y(ijk,  3)
      dw(ijk, 30, 15) =  rk(ijk, 30) * Y(ijk, 15)
      dw(ijk, 30, 15) =  rk(ijk, 30) * Y(ijk, 15)
      dw(ijk, 31, 15) =  rk(ijk, 31) * Y(ijk, 15)
      dw(ijk, 31, 15) =  rk(ijk, 31) * Y(ijk, 15)
      dw(ijk, 32, 12) =  rk(ijk, 32) * Y(ijk,  4)
      dw(ijk, 32,  4) =  rk(ijk, 32) * Y(ijk, 12)
      dw(ijk, 33, 12) =  rk(ijk, 33) * Y(ijk,  5)
      dw(ijk, 33,  5) =  rk(ijk, 33) * Y(ijk, 12)
      dw(ijk, 34, 12) =  rk(ijk, 34) * Y(ijk,  5)
      dw(ijk, 34,  5) =  rk(ijk, 34) * Y(ijk, 12)
      dw(ijk, 35, 14) =  rk(ijk, 35) * Y(ijk,  4)
      dw(ijk, 35,  4) =  rk(ijk, 35) * Y(ijk, 14)
      dw(ijk, 36, 14) =  rk(ijk, 36) * Y(ijk,  5)
      dw(ijk, 36,  5) =  rk(ijk, 36) * Y(ijk, 14)
      dw(ijk, 37, 14) =  rk(ijk, 37) * Y(ijk,  6)
      dw(ijk, 37,  6) =  rk(ijk, 37) * Y(ijk, 14)
      dw(ijk, 38, 15) =  rk(ijk, 38) * Y(ijk,  4)
      dw(ijk, 38,  4) =  rk(ijk, 38) * Y(ijk, 15)
      dw(ijk, 39, 15) =  rk(ijk, 39) * Y(ijk,  5)
      dw(ijk, 39,  5) =  rk(ijk, 39) * Y(ijk, 15)
      dw(ijk, 40, 10) =  rk(ijk, 40)
      dw(ijk, 41, 15) =  rk(ijk, 41) * Y(ijk,  6)
      dw(ijk, 41,  6) =  rk(ijk, 41) * Y(ijk, 15)
      dw(ijk, 42, 14) =  rk(ijk, 42) * Y(ijk,  8)
      dw(ijk, 42,  8) =  rk(ijk, 42) * Y(ijk, 14)
      dw(ijk, 43, 14) =  rk(ijk, 43) * Y(ijk,  9)
      dw(ijk, 43,  9) =  rk(ijk, 43) * Y(ijk, 14)
      dw(ijk, 44, 14) =  rk(ijk, 44) * Y(ijk, 10)
      dw(ijk, 44, 10) =  rk(ijk, 44) * Y(ijk, 14)
      dw(ijk, 45,  1) =  rk(ijk, 45) * Y(ijk,  4)
      dw(ijk, 45,  4) =  rk(ijk, 45) * Y(ijk,  1)
      dw(ijk, 46,  1) =  rk(ijk, 46) * Y(ijk,  5)
      dw(ijk, 46,  5) =  rk(ijk, 46) * Y(ijk,  1)
      dw(ijk, 47,  4) =  rk(ijk, 47) * Y(ijk,  4)
      dw(ijk, 47,  4) =  rk(ijk, 47) * Y(ijk,  4)
      dw(ijk, 48,  6) =  rk(ijk, 48) * Y(ijk,  4)
      dw(ijk, 48,  4) =  rk(ijk, 48) * Y(ijk,  6)
      dw(ijk, 49,  6) =  rk(ijk, 49) * Y(ijk,  5)
      dw(ijk, 49,  5) =  rk(ijk, 49) * Y(ijk,  6)
      dw(ijk, 50,  6) =  rk(ijk, 50) * Y(ijk,  5)
      dw(ijk, 50,  5) =  rk(ijk, 50) * Y(ijk,  6)
      dw(ijk, 51,  7) =  rk(ijk, 51)
      dw(ijk, 52,  6) =  rk(ijk, 52) * Y(ijk,  6)
      dw(ijk, 52,  6) =  rk(ijk, 52) * Y(ijk,  6)
      dw(ijk, 53, 14) =  rk(ijk, 53) * Y(ijk,  2)
      dw(ijk, 53,  2) =  rk(ijk, 53) * Y(ijk, 14)
      dw(ijk, 54, 11) =  rk(ijk, 54) * Y(ijk, 14)
      dw(ijk, 54, 14) =  rk(ijk, 54) * Y(ijk, 11)
      dw(ijk, 55, 16) =  rk(ijk, 55) * Y(ijk, 14)
      dw(ijk, 55, 14) =  rk(ijk, 55) * Y(ijk, 16)
      dw(ijk, 56, 17) =  rk(ijk, 56) * Y(ijk, 14)
      dw(ijk, 56, 14) =  rk(ijk, 56) * Y(ijk, 17)
      dw(ijk, 57, 18) =  rk(ijk, 57) * Y(ijk, 14)
      dw(ijk, 57, 14) =  rk(ijk, 57) * Y(ijk, 18)
      dw(ijk, 58, 17) =  rk(ijk, 58) * Y(ijk,  6)
      dw(ijk, 58,  6) =  rk(ijk, 58) * Y(ijk, 17)
      dw(ijk, 59, 19) =  rk(ijk, 59) * Y(ijk,  4)
      dw(ijk, 59,  4) =  rk(ijk, 59) * Y(ijk, 19)
      dw(ijk, 60, 19) =  rk(ijk, 60) * Y(ijk, 15)
      dw(ijk, 60, 15) =  rk(ijk, 60) * Y(ijk, 19)
      dw(ijk, 61, 19) =  rk(ijk, 61) * Y(ijk, 19)
      dw(ijk, 61, 19) =  rk(ijk, 61) * Y(ijk, 19)
      dw(ijk, 62, 19) =  rk(ijk, 62) * Y(ijk,  6)
      dw(ijk, 62,  6) =  rk(ijk, 62) * Y(ijk, 19)
      dw(ijk, 63, 17) =  rk(ijk, 63) * Y(ijk, 20)
      dw(ijk, 63, 20) =  rk(ijk, 63) * Y(ijk, 17)
      dw(ijk, 64, 19) =  rk(ijk, 64) * Y(ijk, 22)
      dw(ijk, 64, 22) =  rk(ijk, 64) * Y(ijk, 19)
      dw(ijk, 65, 20) =  rk(ijk, 65) * Y(ijk,  1)
      dw(ijk, 65,  1) =  rk(ijk, 65) * Y(ijk, 20)
      dw(ijk, 66, 22) =  rk(ijk, 66) * Y(ijk,  4)
      dw(ijk, 66,  4) =  rk(ijk, 66) * Y(ijk, 22)
      dw(ijk, 67, 22) =  rk(ijk, 67) * Y(ijk, 12)
      dw(ijk, 67, 12) =  rk(ijk, 67) * Y(ijk, 22)
      dw(ijk, 68, 20) =  rk(ijk, 68) * Y(ijk, 16)
      dw(ijk, 68, 16) =  rk(ijk, 68) * Y(ijk, 20)
      dw(ijk, 69, 20) =  rk(ijk, 69)
      dw(ijk, 70, 20) =  rk(ijk, 70) * Y(ijk, 15)
      dw(ijk, 70, 15) =  rk(ijk, 70) * Y(ijk, 20)
      dw(ijk, 71, 20) =  rk(ijk, 71) * Y(ijk, 15)
      dw(ijk, 71, 15) =  rk(ijk, 71) * Y(ijk, 20)
      dw(ijk, 72, 22) =  rk(ijk, 72) * Y(ijk, 14)
      dw(ijk, 72, 14) =  rk(ijk, 72) * Y(ijk, 22)
      dw(ijk, 73, 22) =  rk(ijk, 73) * Y(ijk, 14)
      dw(ijk, 73, 14) =  rk(ijk, 73) * Y(ijk, 22)
      dw(ijk, 74, 14) =  rk(ijk, 74) * Y(ijk, 21)
      dw(ijk, 74, 21) =  rk(ijk, 74) * Y(ijk, 14)
      dw(ijk, 75, 22) =  rk(ijk, 75) * Y(ijk,  5)
      dw(ijk, 75,  5) =  rk(ijk, 75) * Y(ijk, 22)
      dw(ijk, 76, 22) =  rk(ijk, 76) * Y(ijk, 15)
      dw(ijk, 76, 15) =  rk(ijk, 76) * Y(ijk, 22)
      dw(ijk, 77, 14) =  rk(ijk, 77) * Y(ijk, 24)
      dw(ijk, 77, 24) =  rk(ijk, 77) * Y(ijk, 14)
      dw(ijk, 78, 12) =  rk(ijk, 78) * Y(ijk, 24)
      dw(ijk, 78, 24) =  rk(ijk, 78) * Y(ijk, 12)
      dw(ijk, 79, 20) =  rk(ijk, 79) * Y(ijk, 24)
      dw(ijk, 79, 24) =  rk(ijk, 79) * Y(ijk, 20)
      dw(ijk, 80, 22) =  rk(ijk, 80) * Y(ijk, 22)
      dw(ijk, 80, 22) =  rk(ijk, 80) * Y(ijk, 22)
      dw(ijk, 81, 25) =  rk(ijk, 81)
      dw(ijk, 82, 22) =  rk(ijk, 82) * Y(ijk, 22)
      dw(ijk, 82, 22) =  rk(ijk, 82) * Y(ijk, 22)
      dw(ijk, 83, 22) =  rk(ijk, 83) * Y(ijk, 22)
      dw(ijk, 83, 22) =  rk(ijk, 83) * Y(ijk, 22)
      dw(ijk, 84, 12) =  rk(ijk, 84) * Y(ijk, 23)
      dw(ijk, 84, 23) =  rk(ijk, 84) * Y(ijk, 12)
      dw(ijk, 85, 17) =  rk(ijk, 85) * Y(ijk, 28)
      dw(ijk, 85, 28) =  rk(ijk, 85) * Y(ijk, 17)
      dw(ijk, 86, 28) =  rk(ijk, 86) * Y(ijk,  1)
      dw(ijk, 86,  1) =  rk(ijk, 86) * Y(ijk, 28)
      dw(ijk, 87, 27) =  rk(ijk, 87) * Y(ijk, 12)
      dw(ijk, 87, 12) =  rk(ijk, 87) * Y(ijk, 27)
      dw(ijk, 88, 27) =  rk(ijk, 88) * Y(ijk,  4)
      dw(ijk, 88,  4) =  rk(ijk, 88) * Y(ijk, 27)
      dw(ijk, 89, 27) =  rk(ijk, 89) * Y(ijk, 27)
      dw(ijk, 89, 27) =  rk(ijk, 89) * Y(ijk, 27)
      dw(ijk, 90, 28) =  rk(ijk, 90) * Y(ijk, 15)
      dw(ijk, 90, 15) =  rk(ijk, 90) * Y(ijk, 28)
      dw(ijk, 91, 13) =  rk(ijk, 91) * Y(ijk, 29)
      dw(ijk, 91, 29) =  rk(ijk, 91) * Y(ijk, 13)
      dw(ijk, 92, 27) =  rk(ijk, 92) * Y(ijk, 14)
      dw(ijk, 92, 14) =  rk(ijk, 92) * Y(ijk, 27)
      dw(ijk, 93, 14) =  rk(ijk, 93) * Y(ijk, 29)
      dw(ijk, 93, 29) =  rk(ijk, 93) * Y(ijk, 14)
      dw(ijk, 94, 27) =  rk(ijk, 94) * Y(ijk,  5)
      dw(ijk, 94,  5) =  rk(ijk, 94) * Y(ijk, 27)
      dw(ijk, 95, 27) =  rk(ijk, 95) * Y(ijk, 15)
      dw(ijk, 95, 15) =  rk(ijk, 95) * Y(ijk, 27)
      dw(ijk, 96, 31) =  rk(ijk, 96) * Y(ijk, 12)
      dw(ijk, 96, 12) =  rk(ijk, 96) * Y(ijk, 31)
       END DO
 
   END SUBROUTINE dratedc
 
  END MODULE mod_chem_spack_dratedc
 
