DIMENSION nrp(nphase)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2006-04-28  Time: 12:59:13
 
DIMENSION nb(nrmax)
DIMENSION molec(nrmax)

DIMENSION iemonodi(nespmax,nphase)
DIMENSION iedimono(nespmax,nphase)
DIMENSION irmonodi(nrmax,nphase)

DIMENSION jer(3,nrmax)

DIMENSION bp(nbpmax,nrmax)
DIMENSION ispebp(nrmax)

DIMENSION tabphot(ntabphotmax)
DIMENSION iphot(nrmax),iphotinv(nrmax), ntab(nphotmax)
DIMENSION xtab(nphotmax,nangl)
DIMENSION ytab(nphotmax,nangl)
DIMENSION cpolg(nphotmax,nlo)
DIMENSION slump(nespmax,nrmax),seff(nespmax,nrmax)
CHARACTER (LEN=12) :: nom
!
DIMENSION nom(nespmax),inom(nespmax)

!-----------------------------------------------------
COMMON/cinetique1/bp,debug(5), s(nespmax,nrmax),slump,seff,  &
    nr,nrp,nb, jer,molec,ispebp,  &
    imolec1(nrmax,nphase),imolec2(nrmax,nphase), imolec3(nrmax,nphase),  &
    nrmol1(nphase),nrmol2(nphase),nrmol3(nphase), indaqr(nrmax),  &
    narr1(nphase),narr2(nphase),narr3(nphase),  &
    narr4(nphase),narr5(nphase),narr6(nphase), narr7(nphase),narr8(nphase),  &
    iarr1(nrmax,nphase),iarr2(nrmax,nphase),  &
    iarr3(nrmax,nphase),iarr4(nrmax,nphase),  &
    iarr5(nrmax,nphase),iarr6(nrmax,nphase),  &
    iarr7(nrmax,nphase),iarr8(nrmax,nphase), naq2(nphase),inaq2(nrmax,nphase),  &
    iemonodi,irmonodi,iedimono,iphasecom,  &
    indthird(nrmax,nphase),nthird(nphase),ittb(nrmax),  &
    iprecalc(nrmax),iunitaq,iunitgas

COMMON/photolyse/xtab,ytab,cpolg,tabphot,  &
    iphot,iphotinv,ntab,ntabphot,ireversetab,nrphot

COMMON/especes/nom,inom
COMMON/dimens/nesp(nphase),nalg,ndiff(nphase)


COMMON/reducphys/xlphy(nespmax,5),indpur(nespmax,5),  &
    idlump(nespmax),indlump(nespmax), idifford(nespmax)

COMMON/indicateurs/indicaqcom,ireductcom
COMMON/indic1/ixl
COMMON/aqueous/alpha(nespmax),theta,dg(nespmax),diam,xl0,xlmax,  &
    t0,t1cloud,rmol(nespmax), indmod,iswich,jep(nrmax),  &
    indaq(nespmax), ihreac(nespmax),iheq(nrmax),ieqhion(nrmax),  &
    jhplus,johmoin
COMMON/comph/ph,iph,nion,iion(nespmax),ival(nespmax),nitph
COMMON/const/av
COMMON/equil/xk1(nionx),xk2(nionx), seqion(nespmax,nionx),  &
    iesp(nequilx,10), nneq(nrmax),  &
    jion1(nequilx),jion2(nequilx), jaq(nequilx+nionx),  &
    jgaz(nequilx+nionx), nr1(nequilx),nr2(nequilx),nequil,  &
    nequil11,nequil12,nequil13, nequil21,nequil22,nequil23,nequil3,  &
    nequil41,nequil42,nequil51,nequil52, iequil11(neqmax),iequil12(neqmax),  &
    iequil13(neqmax),iequil21(neqmax), iequil22(neqmax),iequil23(neqmax),  &
    iequil3(neqmax),iequil41(neqmax,2), iequil51(neqmax,5),jhpoh(nionx),  &
    iequil52(neqmax,5),iequil42(neqmax,2)



