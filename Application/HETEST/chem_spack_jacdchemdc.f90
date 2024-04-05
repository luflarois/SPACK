                                                                                                                                                      
 MODULE mod_chem_spack_jacdchemdc                                                                                                                     
                                                                                                                                                      
   USE mod_chem_spack_dratedc, ONLY: dratedc  ! subroutine                                                                                            
   IMPLICIT NONE                                                                                                                                      
   PRIVATE                                                                                                                                            
   PUBLIC :: jacdchemdc ! subroutine                                                                                                                  
 CONTAINS                                                                                                                                             
                                                                                                                                                      
   SUBROUTINE jacdchemdc(y,rk,JacC,ngas,ijkbeg,ijkend,maxblock_size,nr)                                                                               
                                                                                                                                                      
!------------------------------------------------------------------------                                                                             
!                                                                                                                                                     
!     -- DESCRIPTION                                                                                                                                  
!                                                                                                                                                     
!     This routine computes the Jacobian matrix for the gas-phase.                                                                                    
!     This routine is automatically generated by SPACK.                                                                                               
!     Mechanism: ../Mechanism/HETEST                                                                                                                  
!     Species: ../Mechanism/ciHETR                                                                                                                    
!                                                                                                                                                     
!------------------------------------------------------------------------                                                                             
!                                                                                                                                                     
!     -- INPUT VARIABLES                                                                                                                              
!                                                                                                                                                     
!     Y: chemical concentrations.                                                                                                                     
!     RK: kinetic rates.                                                                                                                              
!                                                                                                                                                     
!     -- INPUT/OUTPUT VARIABLES                                                                                                                       
!                                                                                                                                                     
!     -- OUTPUT VARIABLES                                                                                                                             
!                                                                                                                                                     
!     JACC: Jacobian matrix.                                                                                                                          
!                                                                                                                                                     
!------------------------------------------------------------------------                                                                             
!                                                                                                                                                     
!     -- REMARKS                                                                                                                                      
!                                                                                                                                                     
!     The matrix JACC could be stored in a low-dimensional vector.                                                                                    
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
                                                                                                                                                      
                                                                                                                                                      
       INTEGER 	 , INTENT(IN)  :: ngas                                                                                                                
       INTEGER 	 , INTENT(IN)  :: ijkbeg			                                                                                                           
       INTEGER 	 , INTENT(IN)  :: ijkend			                                                                                                           
       INTEGER 	 , INTENT(IN)  :: maxblock_size 		                                                                                                    
       INTEGER 	 , INTENT(IN)  :: nr 				                                                                                                             
       DOUBLE PRECISION , INTENT(IN)  :: rk(maxblock_size,nr)		                                                                                       
       DOUBLE PRECISION , INTENT(IN)  :: y(maxblock_size,NGAS) 	                                                                                      
       DOUBLE PRECISION , INTENT(OUT) :: JacC(maxblock_size,NGAS,NGAS)                                                                                
       								                                                                                                                                       
       DOUBLE PRECISION :: dw(maxblock_size,nr,NGAS)			                                                                                               
       INTEGER :: ijk							                                                                                                                          
                                                                                                                                                      
                                                                                                                                                      
      CALL dratedc(rk,y,dw,ngas,ijkbeg,ijkend,maxblock_size,nr)                                                                                       
                                                                                                                                                      
      DO ijk=ijkbeg,ijkend                                                                                                                            
      JacC(ijk,  4,  5) =  + dw(ijk,  1,  5) &
                           + dw(ijk, 33,  5) &
                           + dw(ijk, 49,  5)
      JacC(ijk,  5,  5) =  - dw(ijk,  1,  5) &
                           - dw(ijk, 33,  5) &
                           - dw(ijk, 34,  5) &
                           - dw(ijk, 36,  5) &
                           - dw(ijk, 39,  5) &
                           - dw(ijk, 46,  5) &
                           - dw(ijk, 50,  5) &
                           - dw(ijk, 75,  5) &
                           - dw(ijk, 94,  5)
      JacC(ijk, 12,  5) =  + dw(ijk,  1,  5) &
                           - dw(ijk, 33,  5) &
                           - dw(ijk, 34,  5)
      JacC(ijk,  1,  1) =  - dw(ijk,  2,  1) &
                           - dw(ijk,  3,  1) &
                           - dw(ijk, 21,  1) &
                           - dw(ijk, 26,  1) &
                           - dw(ijk, 27,  1) &
                           - dw(ijk, 45,  1) &
                           - dw(ijk, 46,  1) &
                           - dw(ijk, 65,  1) &
                           - dw(ijk, 86,  1)
      JacC(ijk, 13,  1) =  + dw(ijk,  2,  1)
      JacC(ijk, 12,  1) =  + dw(ijk,  3,  1) &
                           - dw(ijk, 21,  1)
      JacC(ijk,  4,  8) =  + dw(ijk,  4,  8)
      JacC(ijk,  8,  8) =  - dw(ijk,  4,  8) &
                           - dw(ijk, 42,  8)
      JacC(ijk, 14,  8) =  + dw(ijk,  4,  8) &
                           - dw(ijk, 42,  8)
      JacC(ijk,  5,  9) =  + dw(ijk,  5,  9)
      JacC(ijk,  9,  9) =  - dw(ijk,  5,  9) &
                           - dw(ijk, 43,  9)
      JacC(ijk, 14,  9) =  + dw(ijk,  5,  9) &
                           - dw(ijk, 43,  9)
      JacC(ijk,  5, 10) =  + dw(ijk,  6, 10) &
                           + dw(ijk, 40, 10) &
                           + dw(ijk, 44, 10)
      JacC(ijk, 10, 10) =  - dw(ijk,  6, 10) &
                           - dw(ijk, 40, 10) &
                           - dw(ijk, 44, 10)
      JacC(ijk, 15, 10) =  + dw(ijk,  6, 10) &
                           + dw(ijk, 40, 10)
      JacC(ijk,  4,  6) =  + dw(ijk,  7,  6) &
                           - dw(ijk, 48,  6) &
                           + dw(ijk, 49,  6)
      JacC(ijk,  6,  6) =  - dw(ijk,  7,  6) &
                           - dw(ijk,  8,  6) &
                           - dw(ijk, 37,  6) &
                           - dw(ijk, 41,  6) &
                           - dw(ijk, 48,  6) &
                           - dw(ijk, 49,  6) &
                           - dw(ijk, 50,  6) &
                          - 0.2000000000000000D+01*dw(ijk, 52,  6) &
                          - 0.2000000000000000D+01*dw(ijk, 52,  6) &
                           - dw(ijk, 58,  6) &
                           - dw(ijk, 62,  6)
      JacC(ijk,  5,  6) =  + dw(ijk,  8,  6) &
                           + dw(ijk, 37,  6) &
                          + 0.7000000000000000D+00*dw(ijk, 41,  6) &
                          + 0.2000000000000000D+01*dw(ijk, 48,  6) &
                           - dw(ijk, 50,  6) &
                          + 0.2000000000000000D+01*dw(ijk, 52,  6) &
                          + 0.2000000000000000D+01*dw(ijk, 52,  6) &
                           + dw(ijk, 62,  6)
      JacC(ijk, 12,  6) =  + dw(ijk,  8,  6)
      JacC(ijk,  3,  3) =  - dw(ijk,  9,  3) &
                           - dw(ijk, 29,  3)
      JacC(ijk, 14,  3) = + 0.2000000000000000D+01*dw(ijk,  9,  3) &
                           - dw(ijk, 29,  3)
      JacC(ijk,  2, 17) =  + dw(ijk, 10, 17)
      JacC(ijk, 11, 17) =  + dw(ijk, 10, 17) &
                           + dw(ijk, 11, 17) &
                           + dw(ijk, 56, 17) &
                           + dw(ijk, 58, 17) &
                           + dw(ijk, 63, 17) &
                           + dw(ijk, 85, 17)
      JacC(ijk, 17, 17) =  - dw(ijk, 10, 17) &
                           - dw(ijk, 11, 17) &
                           - dw(ijk, 56, 17) &
                           - dw(ijk, 58, 17) &
                           - dw(ijk, 63, 17) &
                           - dw(ijk, 85, 17)
      JacC(ijk, 15, 17) = + 0.2000000000000000D+01*dw(ijk, 11, 17) &
                           + dw(ijk, 56, 17) &
                           + dw(ijk, 58, 17) &
                           + dw(ijk, 63, 17) &
                           + dw(ijk, 85, 17)
      JacC(ijk, 14, 18) =  + dw(ijk, 12, 18) &
                          - 0.6500000000000000D+00*dw(ijk, 57, 18)
      JacC(ijk, 15, 18) =  + dw(ijk, 12, 18)
      JacC(ijk, 17, 18) =  + dw(ijk, 12, 18) &
                          + 0.3500000000000000D+00*dw(ijk, 57, 18)
      JacC(ijk, 18, 18) =  - dw(ijk, 12, 18) &
                           - dw(ijk, 57, 18)
      JacC(ijk,  6, 23) =  + dw(ijk, 13, 23) &
                           + dw(ijk, 84, 23)
      JacC(ijk, 20, 23) =  + dw(ijk, 13, 23)
      JacC(ijk, 23, 23) =  - dw(ijk, 13, 23) &
                           - dw(ijk, 14, 23) &
                           - dw(ijk, 84, 23)
      JacC(ijk,  5, 23) =  + dw(ijk, 14, 23)
      JacC(ijk, 22, 23) =  + dw(ijk, 14, 23) &
                           + dw(ijk, 84, 23)
      JacC(ijk, 20, 25) = + 0.2000000000000000D+01*dw(ijk, 15, 25)
      JacC(ijk, 25, 25) =  - dw(ijk, 15, 25) &
                           - dw(ijk, 81, 25)
      JacC(ijk, 20, 26) = + 0.2000000000000000D+01*dw(ijk, 16, 26)
      JacC(ijk, 26, 26) =  - dw(ijk, 16, 26)
      JacC(ijk, 12, 27) =  + dw(ijk, 17, 27) &
                           - dw(ijk, 87, 27)
      JacC(ijk, 27, 27) =  - dw(ijk, 17, 27) &
                           - dw(ijk, 87, 27) &
                           - dw(ijk, 88, 27) &
                          - 0.2000000000000000D+01*dw(ijk, 89, 27) &
                          - 0.2000000000000000D+01*dw(ijk, 89, 27) &
                           - dw(ijk, 92, 27) &
                           - dw(ijk, 94, 27) &
                           - dw(ijk, 95, 27)
      JacC(ijk, 28, 27) =  + dw(ijk, 17, 27) &
                           + dw(ijk, 87, 27) &
                           + dw(ijk, 88, 27) &
                          + 0.2000000000000000D+01*dw(ijk, 89, 27) &
                          + 0.2000000000000000D+01*dw(ijk, 89, 27) &
                           + dw(ijk, 92, 27)
      JacC(ijk,  6, 30) =  + dw(ijk, 18, 30)
      JacC(ijk, 28, 30) =  + dw(ijk, 18, 30)
      JacC(ijk, 30, 30) =  - dw(ijk, 18, 30)
      JacC(ijk, 14, 31) =  + dw(ijk, 19, 31) &
                           + dw(ijk, 96, 31)
      JacC(ijk, 28, 31) =  + dw(ijk, 19, 31)
      JacC(ijk, 31, 31) =  - dw(ijk, 19, 31) &
                           - dw(ijk, 96, 31)
      JacC(ijk,  1, 12) =  + dw(ijk, 20, 12) &
                           - dw(ijk, 21, 12)
      JacC(ijk, 12, 12) =  - dw(ijk, 20, 12) &
                           - dw(ijk, 21, 12) &
                           - dw(ijk, 32, 12) &
                           - dw(ijk, 33, 12) &
                           - dw(ijk, 34, 12) &
                           - dw(ijk, 67, 12) &
                           - dw(ijk, 78, 12) &
                           - dw(ijk, 84, 12) &
                           - dw(ijk, 87, 12) &
                           - dw(ijk, 96, 12)
      JacC(ijk, 12, 13) =  + dw(ijk, 22, 13) &
                           + dw(ijk, 23, 13)
      JacC(ijk, 13, 13) =  - dw(ijk, 22, 13) &
                           - dw(ijk, 23, 13) &
                           - dw(ijk, 24, 13) &
                           - dw(ijk, 25, 13) &
                           - dw(ijk, 91, 13)
      JacC(ijk, 14, 13) = + 0.2000000000000000D+01*dw(ijk, 24, 13) &
                           + dw(ijk, 25, 13) &
                           + dw(ijk, 91, 13)
      JacC(ijk,  2, 13) =  - dw(ijk, 25, 13)
      JacC(ijk,  2,  2) =  - dw(ijk, 25,  2) &
                           - dw(ijk, 53,  2)
      JacC(ijk, 13,  2) =  - dw(ijk, 25,  2)
      JacC(ijk, 14,  2) =  + dw(ijk, 25,  2) &
                           - dw(ijk, 53,  2)
      JacC(ijk, 15, 13) =  + dw(ijk, 25, 13)
      JacC(ijk, 15,  2) =  + dw(ijk, 25,  2) &
                           + dw(ijk, 53,  2)
      JacC(ijk,  1, 14) =  - dw(ijk, 26, 14)
      JacC(ijk, 14,  1) =  - dw(ijk, 26,  1) &
                           + dw(ijk, 27,  1)
      JacC(ijk, 14, 14) =  - dw(ijk, 26, 14) &
                           - dw(ijk, 28, 14) &
                           - dw(ijk, 29, 14) &
                           - dw(ijk, 35, 14) &
                           - dw(ijk, 36, 14) &
                           - dw(ijk, 37, 14) &
                           - dw(ijk, 42, 14) &
                           - dw(ijk, 43, 14) &
                           - dw(ijk, 44, 14) &
                           - dw(ijk, 53, 14) &
                           - dw(ijk, 54, 14) &
                           - dw(ijk, 55, 14) &
                           - dw(ijk, 56, 14) &
                          - 0.6500000000000000D+00*dw(ijk, 57, 14) &
                           - dw(ijk, 72, 14) &
                           - dw(ijk, 73, 14) &
                           - dw(ijk, 74, 14) &
                           - dw(ijk, 77, 14) &
                           - dw(ijk, 92, 14) &
                           - dw(ijk, 93, 14)
      JacC(ijk, 15,  1) =  + dw(ijk, 26,  1) &
                           - dw(ijk, 27,  1)
      JacC(ijk, 15, 14) =  + dw(ijk, 26, 14) &
                           - dw(ijk, 28, 14) &
                           + dw(ijk, 29, 14) &
                           + dw(ijk, 37, 14) &
                           + dw(ijk, 53, 14) &
                           + dw(ijk, 54, 14) &
                           + dw(ijk, 56, 14) &
                           + dw(ijk, 72, 14) &
                           + dw(ijk, 92, 14)
      JacC(ijk,  1, 15) =  - dw(ijk, 27, 15)
      JacC(ijk, 14, 15) =  + dw(ijk, 27, 15) &
                           - dw(ijk, 28, 15) &
                           + dw(ijk, 38, 15) &
                          + 0.7000000000000000D+00*dw(ijk, 41, 15) &
                           + dw(ijk, 71, 15)
      JacC(ijk, 15, 15) =  - dw(ijk, 27, 15) &
                           - dw(ijk, 28, 15) &
                          - 0.2000000000000000D+01*dw(ijk, 30, 15) &
                          - 0.2000000000000000D+01*dw(ijk, 30, 15) &
                          - 0.2000000000000000D+01*dw(ijk, 31, 15) &
                          - 0.2000000000000000D+01*dw(ijk, 31, 15) &
                           - dw(ijk, 38, 15) &
                           - dw(ijk, 39, 15) &
                           - dw(ijk, 41, 15) &
                           - dw(ijk, 60, 15) &
                           - dw(ijk, 70, 15) &
                           - dw(ijk, 71, 15) &
                           - dw(ijk, 76, 15) &
                           - dw(ijk, 90, 15) &
                           - dw(ijk, 95, 15)
      JacC(ijk,  3, 14) =  - dw(ijk, 29, 14)
      JacC(ijk, 15,  3) =  + dw(ijk, 29,  3)
      JacC(ijk,  3, 15) =  + dw(ijk, 30, 15) &
                           + dw(ijk, 30, 15) &
                           + dw(ijk, 31, 15) &
                           + dw(ijk, 31, 15)
      JacC(ijk,  4, 12) =  - dw(ijk, 32, 12) &
                           + dw(ijk, 33, 12)
      JacC(ijk,  4,  4) =  - dw(ijk, 32,  4) &
                           - dw(ijk, 35,  4) &
                           - dw(ijk, 38,  4) &
                           - dw(ijk, 45,  4) &
                          - 0.2000000000000000D+01*dw(ijk, 47,  4) &
                          - 0.2000000000000000D+01*dw(ijk, 47,  4) &
                           - dw(ijk, 48,  4) &
                           - dw(ijk, 59,  4) &
                           - dw(ijk, 66,  4) &
                           - dw(ijk, 88,  4)
      JacC(ijk,  5, 12) =  + dw(ijk, 32, 12) &
                           - dw(ijk, 33, 12) &
                           - dw(ijk, 34, 12)
      JacC(ijk,  5,  4) =  + dw(ijk, 32,  4) &
                           + dw(ijk, 38,  4) &
                           + dw(ijk, 45,  4) &
                          + 0.2000000000000000D+01*dw(ijk, 47,  4) &
                          + 0.2000000000000000D+01*dw(ijk, 47,  4) &
                          + 0.2000000000000000D+01*dw(ijk, 48,  4) &
                           + dw(ijk, 59,  4) &
                           + dw(ijk, 66,  4) &
                           + dw(ijk, 88,  4)
      JacC(ijk, 12,  4) =  - dw(ijk, 32,  4)
      JacC(ijk,  6, 12) =  + dw(ijk, 34, 12) &
                           + dw(ijk, 84, 12)
      JacC(ijk,  6,  5) =  + dw(ijk, 34,  5) &
                           + dw(ijk, 46,  5) &
                           - dw(ijk, 49,  5) &
                           - dw(ijk, 50,  5)
      JacC(ijk,  4, 14) =  - dw(ijk, 35, 14)
      JacC(ijk,  8, 14) =  + dw(ijk, 35, 14) &
                           - dw(ijk, 42, 14)
      JacC(ijk,  8,  4) =  + dw(ijk, 35,  4)
      JacC(ijk, 14,  4) =  - dw(ijk, 35,  4) &
                           + dw(ijk, 38,  4)
      JacC(ijk,  5, 14) =  - dw(ijk, 36, 14) &
                           + dw(ijk, 37, 14) &
                           + dw(ijk, 42, 14) &
                           + dw(ijk, 44, 14)
      JacC(ijk,  9, 14) =  + dw(ijk, 36, 14) &
                           - dw(ijk, 43, 14)
      JacC(ijk,  9,  5) =  + dw(ijk, 36,  5)
      JacC(ijk, 14,  5) =  - dw(ijk, 36,  5)
      JacC(ijk,  6, 14) =  - dw(ijk, 37, 14) &
                           + dw(ijk, 43, 14)
      JacC(ijk, 14,  6) =  - dw(ijk, 37,  6) &
                          + 0.7000000000000000D+00*dw(ijk, 41,  6)
      JacC(ijk, 15,  6) =  + dw(ijk, 37,  6) &
                           - dw(ijk, 41,  6) &
                           + dw(ijk, 58,  6) &
                           + dw(ijk, 62,  6)
      JacC(ijk,  4, 15) =  - dw(ijk, 38, 15)
      JacC(ijk,  5, 15) =  + dw(ijk, 38, 15) &
                           - dw(ijk, 39, 15) &
                          + 0.7000000000000000D+00*dw(ijk, 41, 15)
      JacC(ijk, 15,  4) =  - dw(ijk, 38,  4) &
                           + dw(ijk, 59,  4)
      JacC(ijk, 10, 15) =  + dw(ijk, 39, 15)
      JacC(ijk, 10,  5) =  + dw(ijk, 39,  5)
      JacC(ijk, 15,  5) =  - dw(ijk, 39,  5)
      JacC(ijk,  6, 15) =  - dw(ijk, 41, 15)
      JacC(ijk,  9, 15) = + 0.3000000000000000D+00*dw(ijk, 41, 15)
      JacC(ijk,  9,  6) = + 0.3000000000000000D+00*dw(ijk, 41,  6) &
                           + dw(ijk, 58,  6)
      JacC(ijk,  5,  8) =  + dw(ijk, 42,  8)
      JacC(ijk,  6,  9) =  + dw(ijk, 43,  9)
      JacC(ijk, 10, 14) =  - dw(ijk, 44, 14)
      JacC(ijk, 14, 10) =  - dw(ijk, 44, 10)
      JacC(ijk,  1,  4) =  - dw(ijk, 45,  4)
      JacC(ijk,  4,  1) =  - dw(ijk, 45,  1)
      JacC(ijk,  5,  1) =  + dw(ijk, 45,  1) &
                           - dw(ijk, 46,  1)
      JacC(ijk,  1,  5) =  - dw(ijk, 46,  5)
      JacC(ijk,  6,  1) =  + dw(ijk, 46,  1)
      JacC(ijk,  6,  4) =  - dw(ijk, 48,  4)
      JacC(ijk,  7,  6) =  + dw(ijk, 50,  6)
      JacC(ijk,  7,  5) =  + dw(ijk, 50,  5)
      JacC(ijk,  5,  7) =  + dw(ijk, 51,  7)
      JacC(ijk,  6,  7) =  + dw(ijk, 51,  7)
      JacC(ijk,  7,  7) =  - dw(ijk, 51,  7)
      JacC(ijk,  2, 14) =  - dw(ijk, 53, 14)
      JacC(ijk, 11, 11) =  - dw(ijk, 54, 11)
      JacC(ijk, 11, 14) =  - dw(ijk, 54, 14) &
                           + dw(ijk, 56, 14)
      JacC(ijk, 14, 11) =  - dw(ijk, 54, 11)
      JacC(ijk, 15, 11) =  + dw(ijk, 54, 11)
      JacC(ijk, 14, 16) =  - dw(ijk, 55, 16)
      JacC(ijk, 16, 16) =  - dw(ijk, 55, 16) &
                           - dw(ijk, 68, 16)
      JacC(ijk, 16, 14) =  - dw(ijk, 55, 14)
      JacC(ijk, 19, 16) =  + dw(ijk, 55, 16) &
                           + dw(ijk, 68, 16)
      JacC(ijk, 19, 14) =  + dw(ijk, 55, 14) &
                          + 0.6500000000000000D+00*dw(ijk, 57, 14)
      JacC(ijk, 14, 17) =  - dw(ijk, 56, 17)
      JacC(ijk, 17, 14) =  - dw(ijk, 56, 14) &
                          + 0.3500000000000000D+00*dw(ijk, 57, 14)
      JacC(ijk, 18, 14) =  - dw(ijk, 57, 14)
      JacC(ijk, 19, 18) = + 0.6500000000000000D+00*dw(ijk, 57, 18)
      JacC(ijk,  6, 17) =  - dw(ijk, 58, 17)
      JacC(ijk,  9, 17) =  + dw(ijk, 58, 17)
      JacC(ijk, 11,  6) =  + dw(ijk, 58,  6)
      JacC(ijk, 17,  6) =  - dw(ijk, 58,  6) &
                           + dw(ijk, 62,  6)
      JacC(ijk,  4, 19) =  - dw(ijk, 59, 19)
      JacC(ijk,  5, 19) =  + dw(ijk, 59, 19) &
                           + dw(ijk, 62, 19) &
                           + dw(ijk, 64, 19)
      JacC(ijk, 15, 19) =  + dw(ijk, 59, 19) &
                           - dw(ijk, 60, 19) &
                          + 0.6600000000000000D+00*dw(ijk, 61, 19) &
                          + 0.6600000000000000D+00*dw(ijk, 61, 19) &
                           + dw(ijk, 62, 19)
      JacC(ijk, 17, 19) =  + dw(ijk, 59, 19) &
                          + 0.1330000000000000D+01*dw(ijk, 61, 19) &
                          + 0.1330000000000000D+01*dw(ijk, 61, 19) &
                           + dw(ijk, 62, 19) &
                           + dw(ijk, 64, 19)
      JacC(ijk, 17,  4) =  + dw(ijk, 59,  4)
      JacC(ijk, 19, 19) =  - dw(ijk, 59, 19) &
                           - dw(ijk, 60, 19) &
                          - 0.2000000000000000D+01*dw(ijk, 61, 19) &
                          - 0.2000000000000000D+01*dw(ijk, 61, 19) &
                           - dw(ijk, 62, 19) &
                           - dw(ijk, 64, 19)
      JacC(ijk, 19,  4) =  - dw(ijk, 59,  4)
      JacC(ijk, 18, 19) =  + dw(ijk, 60, 19)
      JacC(ijk, 18, 15) =  + dw(ijk, 60, 15)
      JacC(ijk, 19, 15) =  - dw(ijk, 60, 15)
      JacC(ijk,  6, 19) =  - dw(ijk, 62, 19)
      JacC(ijk, 19,  6) =  - dw(ijk, 62,  6)
      JacC(ijk, 11, 20) =  + dw(ijk, 63, 20)
      JacC(ijk, 15, 20) =  + dw(ijk, 63, 20) &
                           - dw(ijk, 70, 20) &
                           - dw(ijk, 71, 20)
      JacC(ijk, 17, 20) =  - dw(ijk, 63, 20)
      JacC(ijk, 20, 17) =  - dw(ijk, 63, 17)
      JacC(ijk, 20, 20) =  - dw(ijk, 63, 20) &
                           - dw(ijk, 65, 20) &
                           - dw(ijk, 68, 20) &
                           - dw(ijk, 69, 20) &
                           - dw(ijk, 70, 20) &
                           - dw(ijk, 71, 20) &
                           - dw(ijk, 79, 20)
      JacC(ijk, 21, 17) =  + dw(ijk, 63, 17)
      JacC(ijk, 21, 20) =  + dw(ijk, 63, 20) &
                           + dw(ijk, 68, 20) &
                           + dw(ijk, 69, 20) &
                           + dw(ijk, 70, 20)
      JacC(ijk,  5, 22) =  + dw(ijk, 64, 22) &
                           + dw(ijk, 66, 22) &
                           - dw(ijk, 75, 22)
      JacC(ijk, 17, 22) =  + dw(ijk, 64, 22)
      JacC(ijk, 19, 22) =  - dw(ijk, 64, 22)
      JacC(ijk, 20, 19) =  + dw(ijk, 64, 19)
      JacC(ijk, 20, 22) =  + dw(ijk, 64, 22) &
                           + dw(ijk, 66, 22) &
                           + dw(ijk, 67, 22) &
                           + dw(ijk, 72, 22)
      JacC(ijk, 22, 19) =  - dw(ijk, 64, 19)
      JacC(ijk, 22, 22) =  - dw(ijk, 64, 22) &
                           - dw(ijk, 66, 22) &
                           - dw(ijk, 67, 22) &
                           - dw(ijk, 72, 22) &
                           - dw(ijk, 73, 22) &
                           - dw(ijk, 75, 22) &
                           - dw(ijk, 76, 22) &
                          - 0.2000000000000000D+01*dw(ijk, 80, 22) &
                          - 0.2000000000000000D+01*dw(ijk, 80, 22) &
                          - 0.2000000000000000D+01*dw(ijk, 82, 22) &
                          - 0.2000000000000000D+01*dw(ijk, 82, 22) &
                          - 0.2000000000000000D+01*dw(ijk, 83, 22) &
                          - 0.2000000000000000D+01*dw(ijk, 83, 22)
      JacC(ijk,  1, 20) =  - dw(ijk, 65, 20)
      JacC(ijk, 20,  1) =  - dw(ijk, 65,  1)
      JacC(ijk, 22, 20) =  + dw(ijk, 65, 20) &
                           + dw(ijk, 71, 20)
      JacC(ijk, 22,  1) =  + dw(ijk, 65,  1)
      JacC(ijk,  4, 22) =  - dw(ijk, 66, 22)
      JacC(ijk, 20,  4) =  + dw(ijk, 66,  4)
      JacC(ijk, 22,  4) =  - dw(ijk, 66,  4)
      JacC(ijk, 12, 22) =  - dw(ijk, 67, 22)
      JacC(ijk, 20, 12) =  + dw(ijk, 67, 12)
      JacC(ijk, 22, 12) =  - dw(ijk, 67, 12) &
                           + dw(ijk, 78, 12) &
                           + dw(ijk, 84, 12)
      JacC(ijk, 16, 20) =  - dw(ijk, 68, 20)
      JacC(ijk, 19, 20) =  + dw(ijk, 68, 20)
      JacC(ijk, 20, 16) =  - dw(ijk, 68, 16)
      JacC(ijk, 21, 16) =  + dw(ijk, 68, 16)
      JacC(ijk, 20, 15) =  - dw(ijk, 70, 15) &
                           - dw(ijk, 71, 15)
      JacC(ijk, 21, 15) =  + dw(ijk, 70, 15)
      JacC(ijk, 14, 20) =  + dw(ijk, 71, 20) &
                           + dw(ijk, 79, 20)
      JacC(ijk, 22, 15) =  + dw(ijk, 71, 15) &
                           - dw(ijk, 76, 15)
      JacC(ijk, 14, 22) =  - dw(ijk, 72, 22) &
                           - dw(ijk, 73, 22)
      JacC(ijk, 15, 22) =  + dw(ijk, 72, 22) &
                           - dw(ijk, 76, 22)
      JacC(ijk, 20, 14) =  + dw(ijk, 72, 14) &
                           + dw(ijk, 74, 14)
      JacC(ijk, 22, 14) =  - dw(ijk, 72, 14) &
                           - dw(ijk, 73, 14) &
                           + dw(ijk, 77, 14)
      JacC(ijk, 21, 22) =  + dw(ijk, 73, 22)
      JacC(ijk, 21, 14) =  + dw(ijk, 73, 14) &
                           - dw(ijk, 74, 14)
      JacC(ijk, 14, 21) =  - dw(ijk, 74, 21)
      JacC(ijk, 20, 21) =  + dw(ijk, 74, 21)
      JacC(ijk, 21, 21) =  - dw(ijk, 74, 21)
      JacC(ijk, 22,  5) =  - dw(ijk, 75,  5)
      JacC(ijk, 23, 22) =  + dw(ijk, 75, 22)
      JacC(ijk, 23,  5) =  + dw(ijk, 75,  5)
      JacC(ijk, 24, 22) =  + dw(ijk, 76, 22)
      JacC(ijk, 24, 15) =  + dw(ijk, 76, 15)
      JacC(ijk, 14, 24) =  - dw(ijk, 77, 24) &
                           + dw(ijk, 78, 24) &
                           + dw(ijk, 79, 24)
      JacC(ijk, 22, 24) =  + dw(ijk, 77, 24) &
                           + dw(ijk, 78, 24)
      JacC(ijk, 24, 14) =  - dw(ijk, 77, 14)
      JacC(ijk, 24, 24) =  - dw(ijk, 77, 24) &
                           - dw(ijk, 78, 24) &
                           - dw(ijk, 79, 24)
      JacC(ijk, 12, 24) =  - dw(ijk, 78, 24)
      JacC(ijk, 14, 12) =  + dw(ijk, 78, 12) &
                           + dw(ijk, 96, 12)
      JacC(ijk, 24, 12) =  - dw(ijk, 78, 12)
      JacC(ijk, 20, 24) =  - dw(ijk, 79, 24)
      JacC(ijk, 24, 20) =  - dw(ijk, 79, 20)
      JacC(ijk, 26, 20) =  + dw(ijk, 79, 20)
      JacC(ijk, 26, 24) =  + dw(ijk, 79, 24)
      JacC(ijk, 25, 22) =  + dw(ijk, 80, 22) &
                           + dw(ijk, 80, 22)
      JacC(ijk, 22, 25) = + 0.2000000000000000D+01*dw(ijk, 81, 25)
      JacC(ijk, 26, 22) =  + dw(ijk, 82, 22) &
                           + dw(ijk, 82, 22) &
                           + dw(ijk, 83, 22) &
                           + dw(ijk, 83, 22)
      JacC(ijk, 12, 23) =  - dw(ijk, 84, 23)
      JacC(ijk, 23, 12) =  - dw(ijk, 84, 12)
      JacC(ijk, 11, 28) =  + dw(ijk, 85, 28)
      JacC(ijk, 15, 28) =  + dw(ijk, 85, 28) &
                           - dw(ijk, 90, 28)
      JacC(ijk, 17, 28) =  - dw(ijk, 85, 28)
      JacC(ijk, 28, 17) =  - dw(ijk, 85, 17)
      JacC(ijk, 28, 28) =  - dw(ijk, 85, 28) &
                           - dw(ijk, 86, 28) &
                           - dw(ijk, 90, 28)
      JacC(ijk, 29, 17) =  + dw(ijk, 85, 17)
      JacC(ijk, 29, 28) =  + dw(ijk, 85, 28) &
                           + dw(ijk, 90, 28)
      JacC(ijk,  1, 28) =  - dw(ijk, 86, 28)
      JacC(ijk, 27, 28) =  + dw(ijk, 86, 28)
      JacC(ijk, 27,  1) =  + dw(ijk, 86,  1)
      JacC(ijk, 28,  1) =  - dw(ijk, 86,  1)
      JacC(ijk, 27, 12) =  - dw(ijk, 87, 12) &
                           + dw(ijk, 96, 12)
      JacC(ijk, 28, 12) =  + dw(ijk, 87, 12)
      JacC(ijk,  4, 27) =  - dw(ijk, 88, 27)
      JacC(ijk,  5, 27) =  + dw(ijk, 88, 27) &
                           - dw(ijk, 94, 27)
      JacC(ijk, 27,  4) =  - dw(ijk, 88,  4)
      JacC(ijk, 28,  4) =  + dw(ijk, 88,  4)
      JacC(ijk, 28, 15) =  - dw(ijk, 90, 15)
      JacC(ijk, 29, 15) =  + dw(ijk, 90, 15)
      JacC(ijk, 13, 29) =  - dw(ijk, 91, 29)
      JacC(ijk, 14, 29) =  + dw(ijk, 91, 29) &
                           - dw(ijk, 93, 29)
      JacC(ijk, 28, 13) =  + dw(ijk, 91, 13)
      JacC(ijk, 28, 29) =  + dw(ijk, 91, 29) &
                           + dw(ijk, 93, 29)
      JacC(ijk, 29, 13) =  - dw(ijk, 91, 13)
      JacC(ijk, 29, 29) =  - dw(ijk, 91, 29) &
                           - dw(ijk, 93, 29)
      JacC(ijk, 14, 27) =  - dw(ijk, 92, 27)
      JacC(ijk, 15, 27) =  + dw(ijk, 92, 27) &
                           - dw(ijk, 95, 27)
      JacC(ijk, 27, 14) =  - dw(ijk, 92, 14)
      JacC(ijk, 28, 14) =  + dw(ijk, 92, 14) &
                           + dw(ijk, 93, 14)
      JacC(ijk, 29, 14) =  - dw(ijk, 93, 14)
      JacC(ijk, 27,  5) =  - dw(ijk, 94,  5)
      JacC(ijk, 30, 27) =  + dw(ijk, 94, 27)
      JacC(ijk, 30,  5) =  + dw(ijk, 94,  5)
      JacC(ijk, 27, 15) =  - dw(ijk, 95, 15)
      JacC(ijk, 31, 27) =  + dw(ijk, 95, 27)
      JacC(ijk, 31, 15) =  + dw(ijk, 95, 15)
      JacC(ijk, 12, 31) =  - dw(ijk, 96, 31)
      JacC(ijk, 27, 31) =  + dw(ijk, 96, 31)
      JacC(ijk, 31, 12) =  - dw(ijk, 96, 12)
       END DO                                                                                                                                         
                                                                                                                                                      
   END SUBROUTINE jacdchemdc                                                                                                                          
                                                                                                                                                      
  END MODULE mod_chem_spack_jacdchemdc                                                                                                                
                                                                                                                                                      
