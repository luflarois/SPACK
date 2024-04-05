SUBROUTINE entier(i,ch,ich)
IMPLICIT DOUBLE PRECISION (a-h,o-z)


INTEGER, INTENT(IN OUT)                  :: i
CHARACTER (LEN=500), INTENT(IN OUT)      :: ch
INTEGER, INTENT(IN OUT)                  :: ich

!     Read an integer in a string


CHARACTER (LEN=1) :: chtemp
CHARACTER (LEN=20) :: chtemp2
WRITE(chtemp,'(i1)')ich

chtemp2='(i'//chtemp//')'
READ(ch,chtemp2)i
RETURN
END SUBROUTINE entier



SUBROUTINE part(chdon,mot,imot,nmot)
IMPLICIT REAL (a-h,o-z)
INTEGER, PARAMETER :: nbmot=100


CHARACTER (LEN=500), INTENT(IN)          :: chdon
CHARACTER (LEN=500), INTENT(OUT)         :: mot(nbmot)
INTEGER, INTENT(OUT)                     :: imot(nbmot)
INTEGER, INTENT(OUT)                     :: nmot
!     Decomposition of a sequence (CHDON) in NMOT words (MOT)

COMMON/nblanc/nblanc





DO  i0=1,nblanc
  mot(i0)='                                                    '
END DO
nmot=0
iblanc=1
ilettre=0
DO  nl=1,500
  IF(chdon(nl:nl) /= ' ')THEN
    IF(iblanc == 1)THEN
      iblanc=0
      nmot=nmot+1
      ilettre=1
      mot(nmot)(ilettre:ilettre)=chdon(nl:nl)
    ELSE
      ilettre=ilettre+1
      mot(nmot)(ilettre:ilettre)=chdon(nl:nl)
    END IF
  ELSE
    IF(iblanc == 0)imot(nmot)=ilettre
    iblanc=1
  END IF
END DO

RETURN
END SUBROUTINE part


!     *****************************

SUBROUTINE reel(r,ch,ich)
IMPLICIT DOUBLE PRECISION (a-h,o-z)

DOUBLE PRECISION, INTENT(IN OUT)         :: r
CHARACTER (LEN=500), INTENT(IN OUT)      :: ch
INTEGER, INTENT(IN OUT)                  :: ich

!     Read a real in a string


CHARACTER (LEN=2) :: chtemp
CHARACTER (LEN=20) :: chtemp2
WRITE(chtemp,'(i2)')ich

chtemp2='(e'//chtemp//'.0)'
READ(ch,chtemp2)r
RETURN
END SUBROUTINE reel





