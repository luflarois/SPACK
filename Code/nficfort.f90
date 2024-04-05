  DOUBLE PRECISION  :: petit
  PARAMETER (petit = 1.d-30)

  ! Labels for files kinetic.f, fexchem.f and jacdchemdc.f
  ! Label for non_zero.dat
  ! Labels for fexprod.f, fexloss.f,rates.f,dratedc.f

  INTEGER :: nfick,nficf,nficj,nficnz,nficloss,nficprod
  INTEGER :: nficw,nficdw
  COMMON/labfiles/nfick,nficf,nficj,nficnz,nficloss,nficprod, nficw,nficdw

  !   Files for species and mechanism
  CHARACTER (LEN=20) :: filemeca,filespecies
  COMMON/namefiles/filemeca, filespecies

  ! Current free label
  INTEGER :: ipiste
  COMMON/piste/ipiste

