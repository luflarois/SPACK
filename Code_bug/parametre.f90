INTEGER :: nrmax,nphotmax,nespmax,nphase
INTEGER :: ntabphotmax,nintphotmax,nbpmax
INTEGER :: nequilx,nionx,nangl,nlo,ntab1,neqmax
PARAMETER(nrmax=350)
PARAMETER(nphotmax=2)
PARAMETER(ntabphotmax=12)
PARAMETER(nintphotmax=ntabphotmax-1)
PARAMETER(nbpmax=12+ntabphotmax)
PARAMETER(nespmax=135)
PARAMETER(nphase=3)
PARAMETER(nequilx=25)
PARAMETER(nionx=25)
PARAMETER(nangl=20,nlo=3*nangl,ntab1=11)
PARAMETER(neqmax = nionx+nequilx)
