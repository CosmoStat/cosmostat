pro mkcmp, nuTheo, gTheo, nu, g, Data, sig=sig, Dim=Dim, nuindex=nuindex, tabg, tabnu,tabgm,tabgs

if not keyword_set(sig) then sig=3.
if not keyword_set(Dim) then Dim=128
if not keyword_set(nuindex) then nuindex=-1.

s=100.
NbrPix=512
MaxNu=3.
theogenus, nuTheo, gTheo, sig=sig, NuIndex=NuIndex, NbrPix=NbrPix, MaxNu=MaxNu
gTheo = -gTheo
; gTheo = gTheo / max(gTheo)
plot, nuTheo, gTheo

;cmd = ' -D ' + STRCOMPRESS(STRING(Dim), /REMOVE_ALL) + $
;      ' -g -c '  + STRCOMPRESS(STRING(sig), /REMOVE_ALL)
;  mkgenus, Data, nu, g, opt=cmd, /simu

Nsim=50
Np = 49
tabg = fltarr(Np, Nsim)
tabnu = fltarr(Np, Nsim)

for i = 0,Nsim-1 do begin 
  print, "SIMU ", i+1
  InitV = 10L*(i+1)
  cmd =' -D ' + STRCOMPRESS(STRING(Dim), /REMOVE_ALL) + $
    ' -c ' + STRCOMPRESS(STRING(sig), /REMOVE_ALL) + $
    ' -G ' + STRCOMPRESS(STRING(nuindex), /REMOVE_ALL)+ $
    ' -I ' + STRCOMPRESS(STRING(InitV), /REMOVE_ALL)
  ; mkgauss, Data, OPT=cmd
  
  cmd = ' -D ' + STRCOMPRESS(STRING(Dim), /REMOVE_ALL) + $
      ' -c '  + STRCOMPRESS(STRING(sig), /REMOVE_ALL)+ $
       ' -I ' + STRCOMPRESS(STRING(InitV), /REMOVE_ALL) + $
      ' -g -G ' + STRCOMPRESS(STRING(nuindex), /REMOVE_ALL) 
  mkgenus, Data, nu, g, opt=cmd, /simu
  g = -g
  g = g / max(g)
  tabg(*,i) = g
  tabnu(*,i) = nu
  oplot,  nu, g, psym=1
end
writefits, 'tabg.fits', tabg
writefits, 'tabnu.fits', tabnu
;  tabg = readfits('tabg.fits')
info, tabg

g1 = gTheo / max(gTheo)
g2 = g / max(g)
plot, nuTheo, g1
oplot,  nu, g2, psym=1

tabgm = fltarr(Np)
tabgs = fltarr(Np)
for i = 0,Np-1 do begin
  tabgm(i) = mean(tabg(i,*))
  tabgs(i) = sigma(tabg(i,*))
end
g = tabgm  
info, g
; save, filename='RESGAUSS_C3_SIM50.xdr', nuTheo, gTheo, nu, g, tabgs, tabgm, tabg

; plotg, nuTheo, gTheo, nu, g, tabgs
end

pro makeallcmp
mkcmp, nuTheo, gTheo, nu, g, Data, sig=1., tabg, tabnu,tabgm,tabgs
save, filename='RESGAUSS_C1_SIM50.xdr', nuTheo, gTheo, nu, g, tabgs, tabgm, tabg

mkcmp, nuTheo, gTheo, nu, g, Data, sig=2., tabg, tabnu,tabgm,tabgs
save, filename='RESGAUSS_C2_SIM50.xdr', nuTheo, gTheo, nu, g, tabgs, tabgm, tabg

mkcmp, nuTheo, gTheo, nu, g, Data, sig=4., tabg, tabnu,tabgm,tabgs
save, filename='RESGAUSS_C4_SIM50.xdr', nuTheo, gTheo, nu, g, tabgs, tabgm, tabg

mkcmp, nuTheo, gTheo, nu, g, Data, sig=8., tabg, tabnu,tabgm,tabgs
save, filename='RESGAUSS_C8_SIM50.xdr', nuTheo, gTheo, nu, g, tabgs, tabgm, tabg

end


pro plotg, nuTheo, gTheo, nu, g, tabgs
if keyword_set(PS) then setps, filename='fig_genus_theo_sig3.ps', /portrait
g1 = gTheo / max(gTheo)
g2 = g / max(g)
yrang=[-1.,1.5]
;plot, nuTheo, g1
;oplot,  nu, g2, psym=1
ploterror, nu, g2, tabgs/max(g), psym=3, yticks=5, xthick=3,xrange=xrang, yrange=yrang,  $
         ticklen=0.05,yminor=2,xminor=2, xcharsize=1.5, ycharsize = 1.5, $
         thick=1, xticklen=0.1, ythick=2,/ynozero, title=title
oplot, nuTheo, g1
if keyword_set(PS) then endps
end


pro readdat
g1 = readfits('g_virgo128_c0p1.fits')
g2 = readfits('g_virgo128_c0p5.fits')
g3 = readfits('g_virgo128_c1.fits')
g4 = readfits('g_virgo128_c1p2.fits')
g5 = readfits('g_virgo128_c1p4.fits')
g6 = readfits('g_virgo128_c1p5.fits')
g7 = readfits('g_virgo128_c1p7.fits')
g8 = readfits('g_virgo128_c2.fits')
g9 = readfits('g_virgo128_c3.fits')
g10 = readfits('g_virgo128_c5.fits')

n10 = readfits('g_n51p3_c5.fits')

end
