PRO mrs_alm_convol, Alm, Filter
vs = size( Filter)
Lmax = MIN([ vs[1]-1, Alm.lmax])
for l=1, Lmax do Alm.alm[l,*,*] = Alm.alm[l,*,*] * Filter[l]
for l=Lmax+1, Alm.lmax do Alm.alm[l,*,*] = 0
END
