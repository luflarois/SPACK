                                                                                
 MODULE mod_chem_spack_fexprod                                                  
                                                                                
   IMPLICIT NONE                                                                
   PRIVATE                                                                      
   PUBLIC :: fexprod ! subroutine                                               
 CONTAINS                                                                       
                                                                                
   SUBROUTINE fexprod(w,prod,ngas,ijkbegin,ijkend,maxblock_size,nr)             
                                                                                
!------------------------------------------------------------------------       
!                                                                               
!     -- DESCRIPTION                                                            
!                                                                               
!     This routine computes the production  term P in a P-Lc formulation.       
!     This routine is automatically generated by SPACK.                         
!     Mechanism: ../Mechanism/RACM                                              
!     Species: ../Mechanism/ciRA72A                                             
!                                                                               
!------------------------------------------------------------------------       
!                                                                               
!     -- INPUT VARIABLES                                                        
!                                                                               
!     W: reaction rates.                                                        
!                                                                               
!     -- INPUT/OUTPUT VARIABLES                                                 
!                                                                               
!     -- OUTPUT VARIABLES                                                       
!                                                                               
!     PROD: array of chemical production terms.                                 
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
                                                                                
      INTEGER  	, INTENT(IN)  :: ngas                                           
      INTEGER  	, INTENT(IN)  :: ijkbegin		                                     
      INTEGER  	, INTENT(IN)  :: ijkend 		                                      
      INTEGER  	, INTENT(IN)  :: maxblock_size		                                
      INTEGER  	, INTENT(IN)  :: nr		      	                                    
      DOUBLE PRECISION , INTENT(IN)  :: w(maxblock_size,nr)	                    
      DOUBLE PRECISION , INTENT(OUT) :: prod(maxblock_size,NGAS)                
      INTEGER :: ijk						                                                      
                                                                                
                                                                                
!     Chemical production terms.                                                
                                                                                
      DO ijk=ijkbegin,ijkend                                                    
      prod(ijk,  1) = + w(ijk, 24) &
          + w(ijk,165) &
          + w(ijk,167)
      prod(ijk,  2) = + w(ijk, 33) &
          + w(ijk, 34) &
          +  0.6000000000000000D-02 * w(ijk,107) &
          +  0.1100000000000000D-01 * w(ijk,108) &
          +  0.1000000000000000D-02 * w(ijk,109) &
          +  0.1000000000000000D-02 * w(ijk,110) &
          +  0.2000000000000000D-01 * w(ijk,111) &
          +  0.2000000000000000D-01 * w(ijk,112)
      prod(ijk,  3) = + w(ijk,  1) &
          + w(ijk,  4) &
          + w(ijk,  7) &
          + w(ijk, 36) &
          + w(ijk, 52)
      prod(ijk,  4) = + w(ijk,  5) &
          +  0.6500000000000000D+00 * w(ijk,  6) &
          + w(ijk,  8) &
          + w(ijk, 21) &
          + w(ijk, 35) &
          + w(ijk, 40) &
          + w(ijk, 41) &
          + w(ijk, 43) &
          +  0.7000000000000000D+00 * w(ijk, 44) &
          + w(ijk, 45) &
          + w(ijk, 47) &
          + w(ijk, 48) &
          +  0.2000000000000000D+01 * w(ijk, 50) &
          +  0.2000000000000000D+01 * w(ijk, 51) &
          + w(ijk, 54) &
          +  0.2000000000000000D+01 * w(ijk, 55) &
          + w(ijk, 90) &
          +  0.5000000000000000D+00 * w(ijk, 96) &
          +  0.4000000000000000D+00 * w(ijk,105) &
          +  0.7000000000000000D+00 * w(ijk,115) &
          + w(ijk,128) &
          + w(ijk,130) &
          + w(ijk,131) &
          + w(ijk,132) &
          +  0.9409999999999999D+00 * w(ijk,133) &
          +  0.8760000000000000D+00 * w(ijk,134) &
          +  0.7390000000000000D+00 * w(ijk,135) &
          + w(ijk,136) &
          + w(ijk,137) &
          + w(ijk,138) &
          +  0.8470000000000000D+00 * w(ijk,139) &
          +  0.8000000000000000D+00 * w(ijk,140) &
          +  0.6500000000000000D+00 * w(ijk,141) &
          +  0.9500000000000000D+00 * w(ijk,142) &
          +  0.9500000000000000D+00 * w(ijk,143) &
          + w(ijk,144) &
          + w(ijk,145) &
          + w(ijk,146) &
          + w(ijk,147) &
          + w(ijk,148) &
          +  0.2000000000000000D+01 * w(ijk,149) &
          +  0.5000000000000000D+00 * w(ijk,191) &
          +  0.5160000000000000D+00 * w(ijk,209) &
          +  0.5000000000000000D+00 * w(ijk,211) &
          + w(ijk,212) &
          + w(ijk,213) &
          + w(ijk,214) &
          + w(ijk,215) &
          + w(ijk,216) &
          + w(ijk,217) &
          + w(ijk,218) &
          + w(ijk,219) &
          + w(ijk,220) &
          + w(ijk,221) &
          + w(ijk,222) &
          + w(ijk,223) &
          + w(ijk,224) &
          + w(ijk,225) &
          + w(ijk,226) &
          + w(ijk,227) &
          + w(ijk,228) &
          + w(ijk,229) &
          + w(ijk,230) &
          +  0.2000000000000000D+01 * w(ijk,231) &
          + w(ijk,236) &
          + w(ijk,237)
      prod(ijk,  5) = +  0.3500000000000000D+00 * w(ijk,  6) &
          + w(ijk, 37) &
          + w(ijk, 46) &
          + w(ijk, 49) &
          + w(ijk, 54) &
          + w(ijk, 88) &
          +  0.6000000000000000D+00 * w(ijk, 89)
      prod(ijk,  6) = + w(ijk, 53)
      prod(ijk,  7) = + w(ijk, 38) &
          + w(ijk,118) &
          + w(ijk,121) &
          + w(ijk,124)
      prod(ijk,  8) = + w(ijk, 39) &
          +  0.3000000000000000D+00 * w(ijk, 44) &
          + w(ijk, 91) &
          + w(ijk, 92) &
          + w(ijk, 93) &
          + w(ijk, 94) &
          +  0.2000000000000000D+00 * w(ijk, 95) &
          +  0.5000000000000000D+00 * w(ijk, 96) &
          + w(ijk, 97)
      prod(ijk,  9) = + w(ijk, 42)
      prod(ijk, 11) = + w(ijk, 57)
      prod(ijk, 12) = + w(ijk, 10) &
          + w(ijk, 11) &
          + w(ijk, 12) &
          +  0.1870000000000000D+01 * w(ijk, 17) &
          +  0.1550000000000000D+01 * w(ijk, 18) &
          + w(ijk, 19) &
          + w(ijk, 22) &
          +  0.1000000000000000D-01 * w(ijk, 59) &
          +  0.3600000000000000D-01 * w(ijk, 63) &
          + w(ijk, 76) &
          +  0.2000000000000000D+01 * w(ijk, 80) &
          + w(ijk, 81) &
          +  0.4100000000000000D+00 * w(ijk, 82) &
          + w(ijk, 91) &
          +  0.2000000000000000D+01 * w(ijk, 93) &
          + w(ijk, 94) &
          +  0.8000000000000000D+00 * w(ijk, 95) &
          +  0.4300000000000000D+00 * w(ijk,106) &
          +  0.3700000000000000D+00 * w(ijk,107) &
          +  0.3000000000000000D+00 * w(ijk,108) &
          +  0.3600000000000000D+00 * w(ijk,109) &
          +  0.3600000000000000D+00 * w(ijk,110) &
          +  0.1400000000000000D+00 * w(ijk,111) &
          +  0.1400000000000000D+00 * w(ijk,112) &
          +  0.5400000000000000D+00 * w(ijk,113) &
          +  0.6600000000000000D+00 * w(ijk,114) &
          +  0.1300000000000000D+00 * w(ijk,115)
      prod(ijk, 13) = + w(ijk,  1) &
          + w(ijk,  3) &
          + w(ijk,  8) &
          + w(ijk, 26) &
          + w(ijk, 27) &
          +  0.9000000000000000D-01 * w(ijk,109) &
          +  0.9000000000000000D-01 * w(ijk,110)
      prod(ijk, 14) = + w(ijk,  2)
      prod(ijk, 15) = + w(ijk,  4) &
          + w(ijk,  5) &
          +  0.3500000000000000D+00 * w(ijk,  6) &
          +  0.2000000000000000D+01 * w(ijk,  9) &
          + w(ijk, 13) &
          + w(ijk, 14) &
          + w(ijk, 15) &
          +  0.2000000000000000D+01 * w(ijk, 28) &
          + w(ijk, 30) &
          + w(ijk, 41) &
          +  0.7000000000000000D+00 * w(ijk, 44) &
          +  0.2000000000000000D-01 * w(ijk, 59) &
          +  0.1200000000000000D+00 * w(ijk,106) &
          +  0.4000000000000000D+00 * w(ijk,107) &
          +  0.6300000000000000D+00 * w(ijk,108) &
          +  0.2800000000000000D+00 * w(ijk,109) &
          +  0.2800000000000000D+00 * w(ijk,110) &
          +  0.8500000000000000D+00 * w(ijk,111) &
          +  0.8500000000000000D+00 * w(ijk,112) &
          +  0.7000000000000001D-01 * w(ijk,113) &
          +  0.2100000000000000D+00 * w(ijk,114) &
          +  0.3600000000000000D-01 * w(ijk,115) &
          + w(ijk,120) &
          + w(ijk,123) &
          + w(ijk,126)
      prod(ijk, 16) = +  0.6500000000000000D+00 * w(ijk,  6) &
          +  0.2000000000000000D+01 * w(ijk, 11) &
          + w(ijk, 12) &
          + w(ijk, 13) &
          + w(ijk, 14) &
          +  0.8000000000000000D+00 * w(ijk, 18) &
          + w(ijk, 19) &
          + w(ijk, 20) &
          + w(ijk, 21) &
          + w(ijk, 22) &
          + w(ijk, 23) &
          + w(ijk, 29) &
          + w(ijk, 32) &
          + w(ijk, 40) &
          + w(ijk, 43) &
          + w(ijk, 56) &
          + w(ijk, 57) &
          + w(ijk, 58) &
          +  0.2800000000000000D+00 * w(ijk, 59) &
          +  0.3810000000000000D+00 * w(ijk, 63) &
          +  0.2500000000000000D+00 * w(ijk, 64) &
          +  0.4900000000000000D-01 * w(ijk, 65) &
          +  0.1000000000000000D+00 * w(ijk, 73) &
          +  0.1000000000000000D+00 * w(ijk, 74) &
          +  0.5000000000000000D-01 * w(ijk, 75) &
          + w(ijk, 76) &
          + w(ijk, 79) &
          + w(ijk, 80) &
          +  0.4900000000000000D+00 * w(ijk, 82) &
          +  0.5000000000000000D+00 * w(ijk, 83) &
          + w(ijk, 84) &
          +  0.3500000000000000D+00 * w(ijk, 87) &
          +  0.4000000000000000D+00 * w(ijk, 89) &
          + w(ijk, 91) &
          + w(ijk, 93) &
          +  0.5000000000000000D+00 * w(ijk, 96) &
          +  0.2600000000000000D+00 * w(ijk,106) &
          +  0.2500000000000000D+00 * w(ijk,107) &
          +  0.2200000000000000D+00 * w(ijk,108) &
          +  0.3000000000000000D+00 * w(ijk,109) &
          +  0.3000000000000000D+00 * w(ijk,110) &
          +  0.1000000000000000D+00 * w(ijk,111) &
          +  0.1000000000000000D+00 * w(ijk,112) &
          +  0.2900000000000000D+00 * w(ijk,113) &
          +  0.2900000000000000D+00 * w(ijk,114) &
          +  0.8000000000000000D-01 * w(ijk,115) &
          +  0.2000000000000000D-01 * w(ijk,119) &
          +  0.2000000000000000D-01 * w(ijk,122) &
          +  0.2000000000000000D-01 * w(ijk,125) &
          + w(ijk,131) &
          + w(ijk,132) &
          +  0.7420000000000000D+00 * w(ijk,133) &
          +  0.5990000000000000D+00 * w(ijk,134) &
          +  0.6060000000000000D+00 * w(ijk,135) &
          + w(ijk,136) &
          + w(ijk,137) &
          + w(ijk,138) &
          +  0.8470000000000000D+00 * w(ijk,139) &
          +  0.8000000000000000D+00 * w(ijk,140) &
          +  0.6500000000000000D+00 * w(ijk,141) &
          +  0.9500000000000000D+00 * w(ijk,142) &
          +  0.9500000000000000D+00 * w(ijk,143) &
          + w(ijk,144) &
          +  0.7700000000000000D+00 * w(ijk,147) &
          + w(ijk,148) &
          +  0.6600000000000000D+00 * w(ijk,171) &
          + w(ijk,172) &
          +  0.9920000000000000D+00 * w(ijk,173) &
          +  0.9460000000000000D+00 * w(ijk,174) &
          +  0.9930000000000000D+00 * w(ijk,175) &
          + w(ijk,176) &
          + w(ijk,177) &
          + w(ijk,178) &
          + w(ijk,179) &
          +  0.2000000000000000D+01 * w(ijk,180) &
          +  0.2000000000000000D+01 * w(ijk,181) &
          + w(ijk,182) &
          + w(ijk,183) &
          +  0.2000000000000000D+01 * w(ijk,184) &
          + w(ijk,185) &
          + w(ijk,187) &
          +  0.8800000000000000D+00 * w(ijk,189) &
          + w(ijk,190) &
          +  0.5000000000000000D+00 * w(ijk,191) &
          +  0.5000000000000000D+00 * w(ijk,192) &
          +  0.4880000000000000D+00 * w(ijk,193) &
          +  0.4380000000000000D+00 * w(ijk,194) &
          +  0.4890000000000000D+00 * w(ijk,195) &
          +  0.5000000000000000D+00 * w(ijk,196) &
          +  0.5010000000000000D+00 * w(ijk,197) &
          +  0.5100000000000000D+00 * w(ijk,198) &
          +  0.5060000000000000D+00 * w(ijk,199) &
          + w(ijk,200) &
          + w(ijk,201) &
          + w(ijk,202) &
          + w(ijk,203) &
          + w(ijk,204) &
          +  0.3800000000000000D+00 * w(ijk,207) &
          +  0.5000000000000000D+00 * w(ijk,208) &
          + w(ijk,210) &
          +  0.5000000000000000D+00 * w(ijk,211) &
          + w(ijk,213) &
          + w(ijk,214) &
          +  0.7920000000000000D+00 * w(ijk,215) &
          +  0.6990000000000000D+00 * w(ijk,216) &
          +  0.8450000000000000D+00 * w(ijk,217) &
          + w(ijk,218) &
          + w(ijk,219) &
          + w(ijk,220) &
          + w(ijk,221) &
          + w(ijk,222) &
          + w(ijk,223) &
          + w(ijk,224) &
          + w(ijk,225) &
          + w(ijk,226) &
          +  0.7700000000000000D+00 * w(ijk,229) &
          + w(ijk,230) &
          + w(ijk,233)
      prod(ijk, 17) = +  0.6000000000000000D-01 * w(ijk,107) &
          +  0.7000000000000001D-01 * w(ijk,108)
      prod(ijk, 18) = +  0.3000000000000000D-01 * w(ijk,107) &
          +  0.6000000000000000D-01 * w(ijk,108)
      prod(ijk, 23) = +  0.8600000000000000D+00 * w(ijk, 59) &
          +  0.3500000000000000D+00 * w(ijk,109) &
          +  0.3500000000000000D+00 * w(ijk,110) &
          +  0.4600000000000000D+00 * w(ijk,112) &
          +  0.3540000000000000D+00 * w(ijk,139) &
          +  0.3700000000000000D+00 * w(ijk,179) &
          +  0.2290000000000000D+00 * w(ijk,199) &
          +  0.4000000000000000D+00 * w(ijk,221)
      prod(ijk, 24) = +  0.2500000000000000D+00 * w(ijk,141) &
          +  0.8000000000000000D-01 * w(ijk,179) &
          +  0.4000000000000000D+00 * w(ijk,181) &
          +  0.4000000000000000D+00 * w(ijk,201) &
          +  0.4000000000000000D+00 * w(ijk,223)
      prod(ijk, 31) = +  0.1000000000000000D+00 * w(ijk,116) &
          + w(ijk,117) &
          + w(ijk,118) &
          +  0.2000000000000000D-01 * w(ijk,119) &
          + w(ijk,120) &
          + w(ijk,121) &
          +  0.2000000000000000D-01 * w(ijk,122) &
          + w(ijk,123) &
          + w(ijk,124) &
          +  0.2000000000000000D-01 * w(ijk,125) &
          + w(ijk,126)
      prod(ijk, 32) = + w(ijk, 13) &
          +  0.1300000000000000D+00 * w(ijk, 17) &
          +  0.4500000000000000D+00 * w(ijk, 18) &
          + w(ijk, 22) &
          + w(ijk, 23) &
          +  0.5000000000000000D-01 * w(ijk, 59) &
          +  0.1000000000000000D-01 * w(ijk, 63) &
          +  0.8000000000000000D-01 * w(ijk, 82) &
          +  0.3500000000000000D+00 * w(ijk, 85) &
          +  0.3500000000000000D+00 * w(ijk, 87) &
          + w(ijk, 88) &
          +  0.4000000000000000D+00 * w(ijk, 89) &
          +  0.4000000000000000D+00 * w(ijk,105) &
          + w(ijk,106) &
          +  0.6400000000000000D+00 * w(ijk,107) &
          +  0.2000000000000000D-01 * w(ijk,108) &
          +  0.9000000000000000D+00 * w(ijk,109) &
          +  0.9000000000000000D+00 * w(ijk,110) &
          +  0.4000000000000000D-01 * w(ijk,112) &
          +  0.4000000000000000D+00 * w(ijk,113) &
          +  0.7000000000000000D+00 * w(ijk,115) &
          + w(ijk,131) &
          +  0.4700000000000000D-01 * w(ijk,133) &
          +  0.2100000000000000D-01 * w(ijk,134) &
          +  0.1600000000000000D+01 * w(ijk,136) &
          + w(ijk,137) &
          +  0.6060000000000000D+00 * w(ijk,139) &
          +  0.2500000000000000D+00 * w(ijk,141) &
          + w(ijk,146) &
          +  0.2870000000000000D+00 * w(ijk,149) &
          +  0.1330000000000000D+01 * w(ijk,171) &
          +  0.7500000000000000D+00 * w(ijk,172) &
          +  0.8100000000000001D+00 * w(ijk,173) &
          +  0.8290000000000000D+00 * w(ijk,174) &
          +  0.7530000000000000D+00 * w(ijk,175) &
          +  0.1550000000000000D+01 * w(ijk,176) &
          +  0.1250000000000000D+01 * w(ijk,177) &
          +  0.7550000000000000D+00 * w(ijk,178) &
          +  0.1090000000000000D+01 * w(ijk,179) &
          + w(ijk,180) &
          +  0.1400000000000000D+01 * w(ijk,181) &
          + w(ijk,182) &
          + w(ijk,183) &
          + w(ijk,184) &
          + w(ijk,185) &
          + w(ijk,186) &
          +  0.2000000000000000D+01 * w(ijk,187) &
          + w(ijk,188) &
          +  0.7500000000000000D+00 * w(ijk,189) &
          +  0.7500000000000000D+00 * w(ijk,190) &
          +  0.9600000000000000D+00 * w(ijk,191) &
          +  0.9100000000000000D-01 * w(ijk,193) &
          +  0.7600000000000000D-01 * w(ijk,194) &
          +  0.8000000000000000D+00 * w(ijk,196) &
          +  0.5010000000000000D+00 * w(ijk,197) &
          +  0.3400000000000000D+00 * w(ijk,199) &
          +  0.4000000000000000D+00 * w(ijk,201) &
          + w(ijk,206) &
          +  0.2070000000000000D+00 * w(ijk,209) &
          +  0.2020000000000000D+00 * w(ijk,211) &
          +  0.5040000000000000D+00 * w(ijk,212) &
          + w(ijk,213) &
          +  0.4800000000000000D-01 * w(ijk,215) &
          +  0.2100000000000000D-01 * w(ijk,216) &
          +  0.1600000000000000D+01 * w(ijk,218) &
          + w(ijk,219) &
          +  0.6860000000000001D+00 * w(ijk,221) &
          +  0.4000000000000000D+00 * w(ijk,223) &
          + w(ijk,228) &
          +  0.2800000000000000D+00 * w(ijk,231) &
          + w(ijk,233)
      prod(ijk, 33) = + w(ijk, 14) &
          +  0.2000000000000000D+00 * w(ijk, 21) &
          + w(ijk, 60) &
          +  0.3350000000000000D+00 * w(ijk, 63) &
          +  0.2500000000000000D-01 * w(ijk, 65) &
          +  0.8800000000000000D+00 * w(ijk, 84) &
          +  0.8000000000000000D-01 * w(ijk, 86) &
          +  0.2500000000000000D+00 * w(ijk, 96) &
          +  0.4400000000000000D+00 * w(ijk,107) &
          +  0.9900000000000000D+00 * w(ijk,108) &
          +  0.6500000000000000D+00 * w(ijk,111) &
          +  0.1600000000000000D+00 * w(ijk,114) &
          + w(ijk,132) &
          +  0.2330000000000000D+00 * w(ijk,133) &
          +  0.2110000000000000D+00 * w(ijk,134) &
          +  0.1500000000000000D+00 * w(ijk,135) &
          +  0.2000000000000000D+00 * w(ijk,136) &
          +  0.9399999999999999D+00 * w(ijk,137) &
          +  0.1710000000000000D+01 * w(ijk,138) &
          +  0.8000000000000000D+00 * w(ijk,140) &
          +  0.4600000000000000D+00 * w(ijk,147) &
          +  0.1240000000000000D+01 * w(ijk,149) &
          +  0.7500000000000000D+00 * w(ijk,172) &
          +  0.5800000000000000D+00 * w(ijk,173) &
          +  0.5230000000000000D+00 * w(ijk,174) &
          +  0.4110000000000000D+00 * w(ijk,175) &
          +  0.3500000000000000D+00 * w(ijk,176) &
          +  0.6690000000000000D+00 * w(ijk,177) &
          +  0.9320000000000001D+00 * w(ijk,178) &
          + w(ijk,180) &
          +  0.3000000000000000D+00 * w(ijk,189) &
          +  0.6400000000000000D+00 * w(ijk,191) &
          + w(ijk,192) &
          +  0.7240000000000000D+00 * w(ijk,193) &
          +  0.6770000000000000D+00 * w(ijk,194) &
          +  0.4970000000000000D+00 * w(ijk,195) &
          +  0.6000000000000000D+00 * w(ijk,196) &
          +  0.8590000000000000D+00 * w(ijk,197) &
          +  0.9409999999999999D+00 * w(ijk,198) &
          + w(ijk,200) &
          +  0.3500000000000000D+00 * w(ijk,207) &
          +  0.6500000000000000D+00 * w(ijk,209) &
          +  0.6400000000000000D+00 * w(ijk,211) &
          +  0.1210000000000000D+01 * w(ijk,212) &
          + w(ijk,214) &
          +  0.2430000000000000D+00 * w(ijk,215) &
          +  0.2390000000000000D+00 * w(ijk,216) &
          +  0.1870000000000000D+00 * w(ijk,217) &
          +  0.2000000000000000D+00 * w(ijk,218) &
          +  0.9399999999999999D+00 * w(ijk,219) &
          +  0.1710000000000000D+01 * w(ijk,220) &
          + w(ijk,222) &
          +  0.4600000000000000D+00 * w(ijk,229) &
          +  0.1240000000000000D+01 * w(ijk,231)
      prod(ijk, 34) = +  0.8000000000000000D+00 * w(ijk, 21) &
          +  0.2500000000000000D+00 * w(ijk, 64) &
          +  0.1200000000000000D+00 * w(ijk, 84) &
          +  0.4100000000000000D+00 * w(ijk, 86) &
          +  0.3000000000000000D-01 * w(ijk, 96) &
          +  0.3000000000000000D-01 * w(ijk,107) &
          +  0.1600000000000000D+00 * w(ijk,108) &
          +  0.5300000000000000D+00 * w(ijk,111) &
          +  0.6230000000000000D+00 * w(ijk,133) &
          +  0.7220000000000000D+00 * w(ijk,134) &
          +  0.6420000000000000D+00 * w(ijk,135) &
          +  0.6000000000000000D-01 * w(ijk,137) &
          +  0.2900000000000000D+00 * w(ijk,138) &
          +  0.8000000000000000D+00 * w(ijk,140) &
          +  0.4640000000000000D+00 * w(ijk,149) &
          +  0.1800000000000000D-01 * w(ijk,173) &
          +  0.2400000000000000D+00 * w(ijk,174) &
          +  0.4190000000000000D+00 * w(ijk,175) &
          +  0.8100000000000000D-01 * w(ijk,177) &
          +  0.3130000000000000D+00 * w(ijk,178) &
          + w(ijk,180) &
          +  0.1490000000000000D+00 * w(ijk,191) &
          +  0.1270000000000000D+00 * w(ijk,193) &
          +  0.3300000000000000D+00 * w(ijk,194) &
          +  0.5810000000000000D+00 * w(ijk,195) &
          +  0.1410000000000000D+00 * w(ijk,197) &
          +  0.5690000000000000D+00 * w(ijk,198) &
          + w(ijk,200) &
          +  0.1100000000000000D+00 * w(ijk,207) &
          +  0.1670000000000000D+00 * w(ijk,209) &
          +  0.1490000000000000D+00 * w(ijk,211) &
          +  0.2850000000000000D+00 * w(ijk,212) &
          +  0.6700000000000000D+00 * w(ijk,215) &
          +  0.8280000000000000D+00 * w(ijk,216) &
          +  0.8800000000000000D+00 * w(ijk,217) &
          +  0.6000000000000000D-01 * w(ijk,219) &
          +  0.2900000000000000D+00 * w(ijk,220) &
          + w(ijk,222) &
          +  0.4690000000000000D+00 * w(ijk,231)
      prod(ijk, 35) = +  0.3600000000000000D-01 * w(ijk, 63) &
          +  0.1500000000000000D+00 * w(ijk, 83) &
          +  0.2500000000000000D+00 * w(ijk, 96) &
          +  0.5000000000000000D+00 * w(ijk,114) &
          +  0.6300000000000000D-01 * w(ijk,133) &
          +  0.1200000000000000D+01 * w(ijk,142) &
          +  0.3500000000000000D+00 * w(ijk,143) &
          + w(ijk,144) &
          +  0.1190000000000000D+00 * w(ijk,173) &
          +  0.6500000000000000D+00 * w(ijk,182) &
          +  0.3700000000000000D+00 * w(ijk,183) &
          + w(ijk,184) &
          +  0.1000000000000000D+00 * w(ijk,193) &
          +  0.6500000000000000D+00 * w(ijk,202) &
          +  0.3700000000000000D+00 * w(ijk,203) &
          + w(ijk,204) &
          +  0.6300000000000000D-01 * w(ijk,215) &
          +  0.1300000000000000D+01 * w(ijk,224) &
          +  0.7400000000000000D+00 * w(ijk,225) &
          + w(ijk,226)
      prod(ijk, 36) = + w(ijk, 79) &
          +  0.8000000000000000D-01 * w(ijk, 82) &
          +  0.1500000000000000D+00 * w(ijk, 83) &
          +  0.2500000000000000D+00 * w(ijk, 96) &
          +  0.6000000000000000D+00 * w(ijk,113) &
          +  0.6200000000000000D+00 * w(ijk,114) &
          +  0.6500000000000000D+00 * w(ijk,142) &
          +  0.6000000000000000D+00 * w(ijk,143) &
          + w(ijk,144) &
          +  0.5400000000000000D+00 * w(ijk,147) &
          +  0.5000000000000000D-02 * w(ijk,173) &
          +  0.3500000000000000D+00 * w(ijk,182) &
          +  0.6300000000000000D+00 * w(ijk,183) &
          + w(ijk,184) &
          +  0.4000000000000000D+00 * w(ijk,189) &
          +  0.4000000000000000D-02 * w(ijk,193) &
          +  0.3500000000000000D+00 * w(ijk,202) &
          +  0.6300000000000000D+00 * w(ijk,203) &
          + w(ijk,204) &
          +  0.5400000000000000D+00 * w(ijk,207) &
          +  0.7000000000000000D+00 * w(ijk,224) &
          +  0.1260000000000000D+01 * w(ijk,225) &
          + w(ijk,226) &
          +  0.5400000000000000D+00 * w(ijk,229)
      prod(ijk, 37) = +  0.1300000000000000D+00 * w(ijk, 59) &
          +  0.5000000000000000D+00 * w(ijk,142) &
          +  0.9500000000000000D+00 * w(ijk,143) &
          + w(ijk,182) &
          + w(ijk,183) &
          + w(ijk,202) &
          + w(ijk,203) &
          +  0.5000000000000000D+00 * w(ijk,224) &
          + w(ijk,225)
      prod(ijk, 38) = +  0.9000000000000000D+00 * w(ijk,101) &
          +  0.9000000000000000D+00 * w(ijk,102) &
          +  0.3900000000000000D+00 * w(ijk,109) &
          +  0.3900000000000000D+00 * w(ijk,110) &
          +  0.7900000000000000D+00 * w(ijk,112) &
          +  0.4460000000000000D+00 * w(ijk,139) &
          +  0.4000000000000000D+00 * w(ijk,141) &
          +  0.5500000000000000D+00 * w(ijk,179) &
          +  0.6000000000000000D+00 * w(ijk,181) &
          +  0.7710000000000000D+00 * w(ijk,199) &
          +  0.6000000000000000D+00 * w(ijk,201) &
          +  0.6000000000000000D+00 * w(ijk,221) &
          +  0.6000000000000000D+00 * w(ijk,223)
      prod(ijk, 39) = +  0.3500000000000000D+00 * w(ijk, 83)
      prod(ijk, 40) = +  0.2400000000000000D-01 * w(ijk, 65) &
          +  0.4100000000000000D+00 * w(ijk, 82) &
          +  0.6000000000000000D+00 * w(ijk, 89) &
          +  0.3000000000000000D+00 * w(ijk,189)
      prod(ijk, 41) = +  0.6000000000000000D+00 * w(ijk,105) &
          + w(ijk,116) &
          +  0.5900000000000000D-01 * w(ijk,133) &
          +  0.1240000000000000D+00 * w(ijk,134) &
          +  0.2610000000000000D+00 * w(ijk,135) &
          +  0.1530000000000000D+00 * w(ijk,139) &
          +  0.2000000000000000D+00 * w(ijk,140) &
          +  0.3500000000000000D+00 * w(ijk,141) &
          +  0.5000000000000000D-01 * w(ijk,142) &
          +  0.5000000000000000D-01 * w(ijk,143) &
          + w(ijk,148) &
          + w(ijk,169) &
          + w(ijk,170) &
          + w(ijk,190) &
          +  0.5000000000000000D+00 * w(ijk,191) &
          + w(ijk,208) &
          +  0.4840000000000000D+00 * w(ijk,209) &
          +  0.2000000000000000D+01 * w(ijk,210) &
          +  0.1500000000000000D+01 * w(ijk,211) &
          + w(ijk,212) &
          + w(ijk,230)
      prod(ijk, 42) = +  0.4000000000000000D+00 * w(ijk, 89) &
          +  0.4000000000000000D+00 * w(ijk,105) &
          +  0.3000000000000000D+00 * w(ijk,115) &
          + w(ijk,127)
      prod(ijk, 43) = + w(ijk,129)
      prod(ijk, 44) = + w(ijk,150)
      prod(ijk, 45) = +  0.1300000000000000D+00 * w(ijk,113) &
          + w(ijk,151) &
          + w(ijk,152) &
          + w(ijk,153) &
          + w(ijk,154) &
          + w(ijk,155) &
          + w(ijk,156) &
          + w(ijk,157) &
          + w(ijk,158) &
          + w(ijk,159) &
          + w(ijk,160) &
          + w(ijk,161) &
          + w(ijk,162) &
          + w(ijk,163) &
          + w(ijk,166) &
          + w(ijk,168) &
          + w(ijk,232)
      prod(ijk, 46) = +  0.1100000000000000D+00 * w(ijk,114) &
          + w(ijk,164)
      prod(ijk, 47) = +  0.3600000000000000D-01 * w(ijk, 63) &
          +  0.3700000000000000D+00 * w(ijk,106) &
          +  0.1400000000000000D+00 * w(ijk,107) &
          +  0.1500000000000000D+00 * w(ijk,109) &
          +  0.1500000000000000D+00 * w(ijk,110) &
          +  0.1000000000000000D-01 * w(ijk,112) &
          +  0.2200000000000000D+00 * w(ijk,113) &
          +  0.1100000000000000D+00 * w(ijk,114) &
          +  0.1100000000000000D+00 * w(ijk,115)
      prod(ijk, 48) = +  0.1000000000000000D+00 * w(ijk,107) &
          +  0.1400000000000000D+00 * w(ijk,108) &
          +  0.7000000000000001D-01 * w(ijk,112) &
          +  0.1300000000000000D+00 * w(ijk,113) &
          +  0.2100000000000000D+00 * w(ijk,114) &
          + w(ijk,165) &
          + w(ijk,167) &
          + w(ijk,186) &
          + w(ijk,188) &
          +  0.5000000000000000D+00 * w(ijk,192) &
          +  0.4990000000000000D+00 * w(ijk,193) &
          +  0.4950000000000000D+00 * w(ijk,194) &
          +  0.4950000000000000D+00 * w(ijk,195) &
          +  0.5000000000000000D+00 * w(ijk,196) &
          +  0.4990000000000000D+00 * w(ijk,197) &
          +  0.4900000000000000D+00 * w(ijk,198) &
          +  0.4940000000000000D+00 * w(ijk,199) &
          +  0.5000000000000000D+00 * w(ijk,207) &
          +  0.5000000000000000D+00 * w(ijk,208) &
          +  0.4840000000000000D+00 * w(ijk,209)
      prod(ijk, 49) = + w(ijk, 12) &
          + w(ijk, 15) &
          + w(ijk, 61) &
          +  0.6500000000000000D+00 * w(ijk, 85) &
          +  0.1900000000000000D+00 * w(ijk,107) &
          +  0.2300000000000000D+00 * w(ijk,108) &
          +  0.3000000000000000D-01 * w(ijk,109) &
          +  0.3000000000000000D-01 * w(ijk,110) &
          +  0.1500000000000000D+00 * w(ijk,133) &
          +  0.3100000000000000D-01 * w(ijk,134) &
          + w(ijk,145) &
          +  0.5000000000000000D+00 * w(ijk,192) &
          +  0.5080000000000000D+00 * w(ijk,193) &
          +  0.5540000000000000D+00 * w(ijk,194) &
          +  0.5070000000000000D+00 * w(ijk,195) &
          +  0.5000000000000000D+00 * w(ijk,196) &
          +  0.5010000000000000D+00 * w(ijk,197) &
          +  0.5100000000000000D+00 * w(ijk,198) &
          +  0.5060000000000000D+00 * w(ijk,199) &
          + w(ijk,200) &
          + w(ijk,201) &
          + w(ijk,202) &
          + w(ijk,203) &
          + w(ijk,204) &
          +  0.2000000000000000D+01 * w(ijk,205) &
          + w(ijk,206) &
          +  0.5000000000000000D+00 * w(ijk,207) &
          +  0.5000000000000000D+00 * w(ijk,208) &
          +  0.5160000000000000D+00 * w(ijk,209) &
          +  0.1550000000000000D+00 * w(ijk,215) &
          +  0.4000000000000000D-01 * w(ijk,216) &
          + w(ijk,227) &
          + w(ijk,234)
      prod(ijk, 50) = + w(ijk, 16) &
          + w(ijk, 62) &
          +  0.1000000000000000D+00 * w(ijk,107) &
          +  0.1800000000000000D+00 * w(ijk,108) &
          +  0.2000000000000000D+00 * w(ijk,111) &
          +  0.1600000000000000D+00 * w(ijk,112) &
          +  0.4800000000000000D-01 * w(ijk,133) &
          +  0.2450000000000000D+00 * w(ijk,134) &
          +  0.1330000000000000D+00 * w(ijk,135) &
          +  0.1400000000000000D-01 * w(ijk,174) &
          +  0.1300000000000000D-01 * w(ijk,175) &
          +  0.6000000000000000D-02 * w(ijk,193) &
          +  0.1800000000000000D-01 * w(ijk,194) &
          +  0.1500000000000000D-01 * w(ijk,195) &
          +  0.5300000000000000D-01 * w(ijk,215) &
          +  0.2620000000000000D+00 * w(ijk,216) &
          +  0.1550000000000000D+00 * w(ijk,217)
      prod(ijk, 51) = +  0.5830000000000000D+00 * w(ijk, 63) &
          +  0.4400000000000000D+00 * w(ijk, 86) &
          + w(ijk, 90)
      prod(ijk, 52) = +  0.7500000000000000D+00 * w(ijk, 64)
      prod(ijk, 53) = +  0.9510000000000000D+00 * w(ijk, 65)
      prod(ijk, 54) = + w(ijk, 66)
      prod(ijk, 55) = + w(ijk, 67)
      prod(ijk, 56) = + w(ijk, 68)
      prod(ijk, 57) = + w(ijk, 69) &
          + w(ijk, 70)
      prod(ijk, 58) = + w(ijk, 71)
      prod(ijk, 59) = + w(ijk, 72)
      prod(ijk, 60) = +  0.1000000000000000D+00 * w(ijk, 75) &
          + w(ijk, 97)
      prod(ijk, 61) = +  0.9000000000000000D+00 * w(ijk, 73)
      prod(ijk, 62) = +  0.9000000000000000D+00 * w(ijk, 74)
      prod(ijk, 63) = +  0.8500000000000000D+00 * w(ijk, 75)
      prod(ijk, 64) = +  0.9800000000000000D+00 * w(ijk,119)
      prod(ijk, 65) = +  0.9800000000000000D+00 * w(ijk,122)
      prod(ijk, 66) = +  0.9800000000000000D+00 * w(ijk,125)
      prod(ijk, 67) = + w(ijk, 16) &
          + w(ijk, 19) &
          + w(ijk, 22) &
          + w(ijk, 23) &
          + w(ijk, 77) &
          + w(ijk, 81) &
          +  0.6500000000000000D+00 * w(ijk, 87) &
          + w(ijk, 92) &
          + w(ijk, 94) &
          +  0.1500000000000000D+00 * w(ijk,109) &
          +  0.1500000000000000D+00 * w(ijk,110) &
          +  0.1300000000000000D+00 * w(ijk,113) &
          +  0.2800000000000000D+00 * w(ijk,114) &
          +  0.7000000000000000D+00 * w(ijk,115) &
          + w(ijk,128) &
          + w(ijk,146) &
          +  0.2300000000000000D+00 * w(ijk,147) &
          + w(ijk,187) &
          +  0.1200000000000000D+00 * w(ijk,189) &
          + w(ijk,228) &
          +  0.2300000000000000D+00 * w(ijk,229)
      prod(ijk, 68) = + w(ijk, 20) &
          +  0.5100000000000000D+00 * w(ijk, 82) &
          +  0.5000000000000000D+00 * w(ijk, 83) &
          +  0.2000000000000000D+00 * w(ijk, 95) &
          +  0.5000000000000000D+00 * w(ijk, 96) &
          + w(ijk,130)
      prod(ijk, 69) = + w(ijk, 78) &
          +  0.3000000000000000D-01 * w(ijk,107) &
          +  0.1200000000000000D+00 * w(ijk,108) &
          +  0.2000000000000000D-01 * w(ijk,109) &
          +  0.2000000000000000D-01 * w(ijk,110) &
          +  0.4200000000000000D+00 * w(ijk,111) &
          +  0.4200000000000000D+00 * w(ijk,112)
      prod(ijk, 70) = +  0.8000000000000000D+00 * w(ijk, 95) &
          +  0.8000000000000000D+00 * w(ijk, 98) &
          +  0.4300000000000000D+00 * w(ijk, 99) &
          +  0.1100000000000000D+00 * w(ijk,100) &
          +  0.9000000000000000D+00 * w(ijk,101) &
          +  0.9000000000000000D+00 * w(ijk,102) &
          +  0.1000000000000000D+00 * w(ijk,103) &
          +  0.1300000000000000D+00 * w(ijk,104)
      prod(ijk, 71) = +  0.2000000000000000D+00 * w(ijk, 98) &
          +  0.5700000000000000D+00 * w(ijk, 99) &
          +  0.8900000000000000D+00 * w(ijk,100) &
          +  0.1000000000000000D+00 * w(ijk,101) &
          +  0.1000000000000000D+00 * w(ijk,102) &
          +  0.9000000000000000D+00 * w(ijk,103) &
          +  0.8700000000000000D+00 * w(ijk,104)
      prod(ijk, 72) = +  0.1500000000000000D+00 * w(ijk, 59) &
          +  0.1000000000000000D+00 * w(ijk, 73) &
          +  0.1000000000000000D+00 * w(ijk, 74) &
          +  0.5000000000000000D-01 * w(ijk, 75) &
          +  0.4900000000000000D+00 * w(ijk, 82) &
          +  0.5000000000000000D+00 * w(ijk, 83) &
          +  0.7000000000000001D-01 * w(ijk, 86) &
          +  0.3500000000000000D+00 * w(ijk, 87) &
          + w(ijk, 88) &
          + w(ijk, 89) &
          +  0.5000000000000000D+00 * w(ijk, 96) &
          + w(ijk,105) &
          +  0.1300000000000000D+00 * w(ijk,109) &
          +  0.1300000000000000D+00 * w(ijk,110) &
          +  0.4800000000000000D-01 * w(ijk,133) &
          +  0.3340000000000000D+00 * w(ijk,134) &
          +  0.4160000000000000D+00 * w(ijk,135) &
          +  0.1600000000000000D+00 * w(ijk,147) &
          +  0.8500000000000001D-01 * w(ijk,173) &
          +  0.2450000000000000D+00 * w(ijk,174) &
          +  0.3220000000000000D+00 * w(ijk,175) &
          +  0.8000000000000000D-01 * w(ijk,189) &
          +  0.7099999999999999D-01 * w(ijk,193) &
          +  0.2370000000000000D+00 * w(ijk,194) &
          +  0.3180000000000000D+00 * w(ijk,195) &
          +  0.8000000000000000D-01 * w(ijk,207) &
          +  0.5100000000000000D-01 * w(ijk,215) &
          +  0.3910000000000000D+00 * w(ijk,216) &
          +  0.5870000000000000D+00 * w(ijk,217) &
          +  0.1600000000000000D+00 * w(ijk,229)
      END DO   
                                                                                
   END SUBROUTINE fexprod                                                       
                                                                                
  END MODULE mod_chem_spack_fexprod                                             
                                                                                
