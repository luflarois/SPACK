                                                                                
 MODULE mod_chem_spack_fexloss                                                  
                                                                                
   IMPLICIT NONE                                                                
   PRIVATE                                                                      
   PUBLIC :: fexloss ! subroutine                                               
 CONTAINS                                                                       
                                                                                
   SUBROUTINE fexloss(dw,loss,ngas,ijkbegin,ijkend,maxblock_size,nr)            
                                                                                
!------------------------------------------------------------------------       
!                                                                               
!     -- DESCRIPTION                                                            
!                                                                               
!     This routine computes the chemical loss  term L in a P-Lc formulation.    
!     This routine is automatically generated by SPACK.                         
!     Mechanism: ../Mechanism/HETEST                                            
!     Species: ../Mechanism/ciHETR                                              
!                                                                               
!------------------------------------------------------------------------       
!                                                                               
!     -- INPUT VARIABLES                                                        
!                                                                               
!     DW: derivative of reaction rates wrt Y.                                   
!                                                                               
!     -- INPUT/OUTPUT VARIABLES                                                 
!                                                                               
!     -- OUTPUT VARIABLES                                                       
!                                                                               
!     LOSS: array of chemical loss terms.                                       
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
                                                                                
     INTEGER	       , INTENT(IN)  :: NGAS	                                      
     INTEGER	       , INTENT(IN)  :: ijkbegin		                                 
     INTEGER	       , INTENT(IN)  :: ijkend			                                  
     INTEGER	       , INTENT(IN)  :: maxblock_size		                            
     INTEGER	       , INTENT(IN)  :: nr		 	                                     
     DOUBLE PRECISION , INTENT(IN)  :: dw(maxblock_size,nr,NGAS)                
     DOUBLE PRECISION , INTENT(OUT) :: loss(maxblock_size,NGAS)                 
     INTEGER :: ijk						                                                       
                                                                                
                                                                                
!     Chemical loss terms.                                                      
                                                                                
      DO ijk=ijkbegin,ijkend                                                    
      loss(ijk,  1) = + dw(ijk,  2,  1) &
          + dw(ijk,  3,  1) &
          + dw(ijk, 21,  1) &
          + dw(ijk, 26,  1) &
          + dw(ijk, 27,  1) &
          + dw(ijk, 45,  1) &
          + dw(ijk, 46,  1) &
          + dw(ijk, 65,  1) &
          + dw(ijk, 86,  1)
      loss(ijk,  2) = + dw(ijk, 25,  2) &
          + dw(ijk, 53,  2)
      loss(ijk,  3) = + dw(ijk,  9,  3) &
          + dw(ijk, 29,  3)
      loss(ijk,  4) = + dw(ijk, 32,  4) &
          + dw(ijk, 35,  4) &
          + dw(ijk, 38,  4) &
          + dw(ijk, 45,  4) &
          +  0.2000000000000000D+01 * dw(ijk, 47,  4) &
          + dw(ijk, 48,  4) &
          + dw(ijk, 59,  4) &
          + dw(ijk, 66,  4) &
          + dw(ijk, 88,  4)
      loss(ijk,  5) = + dw(ijk,  1,  5) &
          + dw(ijk, 33,  5) &
          + dw(ijk, 34,  5) &
          + dw(ijk, 36,  5) &
          + dw(ijk, 39,  5) &
          + dw(ijk, 46,  5) &
          + dw(ijk, 50,  5) &
          + dw(ijk, 75,  5) &
          + dw(ijk, 94,  5)
      loss(ijk,  6) = + dw(ijk,  7,  6) &
          + dw(ijk,  8,  6) &
          + dw(ijk, 37,  6) &
          + dw(ijk, 41,  6) &
          + dw(ijk, 48,  6) &
          + dw(ijk, 49,  6) &
          + dw(ijk, 50,  6) &
          +  0.2000000000000000D+01 * dw(ijk, 52,  6) &
          + dw(ijk, 58,  6) &
          + dw(ijk, 62,  6)
      loss(ijk,  7) = + dw(ijk, 51,  7)
      loss(ijk,  8) = + dw(ijk,  4,  8) &
          + dw(ijk, 42,  8)
      loss(ijk,  9) = + dw(ijk,  5,  9) &
          + dw(ijk, 43,  9)
      loss(ijk, 10) = + dw(ijk,  6, 10) &
          + dw(ijk, 40, 10) &
          + dw(ijk, 44, 10)
      loss(ijk, 11) = + dw(ijk, 54, 11)
      loss(ijk, 12) = + dw(ijk, 20, 12) &
          + dw(ijk, 21, 12) &
          + dw(ijk, 32, 12) &
          + dw(ijk, 33, 12) &
          + dw(ijk, 34, 12) &
          + dw(ijk, 67, 12) &
          + dw(ijk, 78, 12) &
          + dw(ijk, 84, 12) &
          + dw(ijk, 87, 12) &
          + dw(ijk, 96, 12)
      loss(ijk, 13) = + dw(ijk, 22, 13) &
          + dw(ijk, 23, 13) &
          + dw(ijk, 24, 13) &
          + dw(ijk, 25, 13) &
          + dw(ijk, 91, 13)
      loss(ijk, 14) = + dw(ijk, 26, 14) &
          + dw(ijk, 28, 14) &
          + dw(ijk, 29, 14) &
          + dw(ijk, 35, 14) &
          + dw(ijk, 36, 14) &
          + dw(ijk, 37, 14) &
          + dw(ijk, 42, 14) &
          + dw(ijk, 43, 14) &
          + dw(ijk, 44, 14) &
          + dw(ijk, 53, 14) &
          + dw(ijk, 54, 14) &
          + dw(ijk, 55, 14) &
          + dw(ijk, 56, 14) &
          +  0.6500000000000000D+00 * dw(ijk, 57, 14) &
          + dw(ijk, 72, 14) &
          + dw(ijk, 73, 14) &
          + dw(ijk, 74, 14) &
          + dw(ijk, 77, 14) &
          + dw(ijk, 92, 14) &
          + dw(ijk, 93, 14)
      loss(ijk, 15) = + dw(ijk, 27, 15) &
          + dw(ijk, 28, 15) &
          +  0.2000000000000000D+01 * dw(ijk, 30, 15) &
          +  0.2000000000000000D+01 * dw(ijk, 31, 15) &
          + dw(ijk, 38, 15) &
          + dw(ijk, 39, 15) &
          + dw(ijk, 41, 15) &
          + dw(ijk, 60, 15) &
          + dw(ijk, 70, 15) &
          + dw(ijk, 71, 15) &
          + dw(ijk, 76, 15) &
          + dw(ijk, 90, 15) &
          + dw(ijk, 95, 15)
      loss(ijk, 16) = + dw(ijk, 55, 16) &
          + dw(ijk, 68, 16)
      loss(ijk, 17) = + dw(ijk, 10, 17) &
          + dw(ijk, 11, 17) &
          + dw(ijk, 56, 17) &
          + dw(ijk, 58, 17) &
          + dw(ijk, 63, 17) &
          + dw(ijk, 85, 17)
      loss(ijk, 18) = + dw(ijk, 12, 18) &
          + dw(ijk, 57, 18)
      loss(ijk, 19) = + dw(ijk, 59, 19) &
          + dw(ijk, 60, 19) &
          +  0.2000000000000000D+01 * dw(ijk, 61, 19) &
          + dw(ijk, 62, 19) &
          + dw(ijk, 64, 19)
      loss(ijk, 20) = + dw(ijk, 63, 20) &
          + dw(ijk, 65, 20) &
          + dw(ijk, 68, 20) &
          + dw(ijk, 69, 20) &
          + dw(ijk, 70, 20) &
          + dw(ijk, 71, 20) &
          + dw(ijk, 79, 20)
      loss(ijk, 21) = + dw(ijk, 74, 21)
      loss(ijk, 22) = + dw(ijk, 64, 22) &
          + dw(ijk, 66, 22) &
          + dw(ijk, 67, 22) &
          + dw(ijk, 72, 22) &
          + dw(ijk, 73, 22) &
          + dw(ijk, 75, 22) &
          + dw(ijk, 76, 22) &
          +  0.2000000000000000D+01 * dw(ijk, 80, 22) &
          +  0.2000000000000000D+01 * dw(ijk, 82, 22) &
          +  0.2000000000000000D+01 * dw(ijk, 83, 22)
      loss(ijk, 23) = + dw(ijk, 13, 23) &
          + dw(ijk, 14, 23) &
          + dw(ijk, 84, 23)
      loss(ijk, 24) = + dw(ijk, 77, 24) &
          + dw(ijk, 78, 24) &
          + dw(ijk, 79, 24)
      loss(ijk, 25) = + dw(ijk, 15, 25) &
          + dw(ijk, 81, 25)
      loss(ijk, 26) = + dw(ijk, 16, 26)
      loss(ijk, 27) = + dw(ijk, 17, 27) &
          + dw(ijk, 87, 27) &
          + dw(ijk, 88, 27) &
          +  0.2000000000000000D+01 * dw(ijk, 89, 27) &
          + dw(ijk, 92, 27) &
          + dw(ijk, 94, 27) &
          + dw(ijk, 95, 27)
      loss(ijk, 28) = + dw(ijk, 85, 28) &
          + dw(ijk, 86, 28) &
          + dw(ijk, 90, 28)
      loss(ijk, 29) = + dw(ijk, 91, 29) &
          + dw(ijk, 93, 29)
      loss(ijk, 30) = + dw(ijk, 18, 30)
      loss(ijk, 31) = + dw(ijk, 19, 31) &
          + dw(ijk, 96, 31)
      END DO   
                                                                                
   END SUBROUTINE fexloss                                                       
                                                                                
  END MODULE mod_chem_spack_fexloss                                             
                                                                                
