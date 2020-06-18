;+
; NAME:
;        mrs_sparse_master
;
; PURPOSE:
;    Compute using sparsity both the power spectrum and the denoised-inpainted map  from a masked noisy spherical image.
;    The noise is assumed to Gaussian, but not stationary. A Root Mean Square (RMS) map must be given as an input parameter.
;
; CALLING:
;       mrs_sparse_master, Map, Mask, RMS, InpMap, InpPS, opt=opt, plot=plot, Clin=Clin, niter=niter,  MatMask=MatMask, NoisePS=NoisePS, Nrea=Nrea
;
; INPUTS:
;     Map  -- IDL array of healpix map: Noisy and masked input image 
;     Mask -- IDL array of healpix map: Mask image, which should contain only 1 or 0. Zero pixels indicates missing values.
;     RMS -- IDL array of healpix map:  RMS image
; 
; OUTPUTS:
;     InpMap Map  -- IDL array of healpix map: Inpainted and denoised image.
;     InpPS-- 1D IDL fltarr:  Power Spectrum of the inpainted map
;
; INPUT KEYWORDS:
;               niter -- scalar: number of iterations. Default is 40.
;               opt  -- string:  additional option for the mrs_alm_inpainting program
;               Nrea -- scalar: number of realizations for the noise power spectrum estimation. Default is 50. Only used if Clin and NoisePS are not set.
;
; INPUT/OUTPUT KEYWORDS:
;               Clin  -- 1D IDL fltarr:  Power Spectrum of the solution.
;         MatMask -- 2D IDL fltarr: matrix related to the Mask.
;
; EXTERNAL CALLS:
;       mrs_alm_inpainting (C++ program)
;
; EXAMPLE:
;       Compute the  power spectrum and the denoised map of an healpix image 
;               mrs_sparse_master, Map, Mask, RMS,  InpMap, InpPS
;         
; HISTORY:
;       Written: Jean-Luc Starck, 2010
;       June, 2010 File creation
;--------------------------------------------------------------------------------------------------------
 
;=====================================================================


;================================

function mrs_get_master_ps, PDat, Mask, MatMask=MatMask,  niter=niter, TrueCl=TrueCL, NoisePS=NoisePS, Tmat=Tmat

if not keyword_set(niter) then niter=20
if not keyword_set(NoisePS) then NoisePS=0.

Fsky = float( N_ELEMENTS(Mask) / float ( total(mask)) )
if not keyword_set(MatMask) then  MatMask  = mrs_matmask(Mask, lmax=lmax)

; Firt initialisation: Fsky correction and noise subtraction
PRes = PDat*Fsky  

Tmat = transpose(MatMask)

for i=0,Niter-1 do begin
   resi = PDat - MatMask # PRes
   PRes += Tmat # resi
   
   ; Positivity constraint
   ind = where(Pres LT 0, c)
   if c GT 0 then Pres[ind] = 0
   
   if keyword_set(TrueCl) then begin
      plotcl, TrueCl, thick=2
      oplotcl, Pres
   end
end
Pres = Pres-NoisePS
; Positivity constraint
   ind = where(Pres LT 0, c)
   if c GT 0 then Pres[ind] = 0
   
return, Pres
end
 
;================================

pro mrs_sparse_master, Map, Mask, RMS, InpMap, InpPS, opt=opt, plot=plot, Clin=Clin, niter=niter,  MatMask=MatMask, NoisePS=NoisePS, Nrea=Nrea, m1=m1

if N_PARAMS() LT 5  then begin 
        print, 'CALLING SEQUENCE: mrs_sparse_master, Map, Mask, RMS, InpMap, InpPS, opt=opt, plot=plot, Clin=Clin, niter=niter,  MatMask=MatMask, NoisePS=NoisePS, Nrea=Nrea'
        Ret=-1
        goto, DONE
end

if N_ELEMENTS(Map) NE N_ELEMENTS(Mask) or N_ELEMENTS(Map) NE N_ELEMENTS(RMS)  then begin
      print, "ERROR:  input maps don't have the same size .... "
      help, Map
      help, Mask
      help, RMS
       goto, DONE
end

ind = where(RMS eq 0 and mask EQ 1, c)
if c GT 0 then RMS[ind] = 1

if not keyword_set(Clin)  then BEGIN
   PDat  = mrs_powspec(Map)
   nside = gnside(Map)
   if not keyword_set(NoisePS) then begin
		if not keyword_set(Nrea) then Nrea=50
		PN=0
		for i=0, Nrea-1 do begin N = randomn(seed, nside^2*12)*rms &  Pn = pn  + mrs_powspec(n) & end
		Pn = Pn / float(Nrea)
		NoisePS = Pn
    end
   PSSol = mrs_get_master_ps( PDat, Mask, MatMask=MatMask, niter=20, NoisePS=NoisePS)
   Clin = PSSol
end

MapFN = gettmpfilename()  
RMSFN = gettmpfilename() 
ClFN = gettmpfilename() 
MaskFN = gettmpfilename() 

mrs_write, MapFN, Map
mrs_write, RMSFN, RMS
mrs_write, MaskFN, Mask
if keyword_set(Clin) then writefits, ClFN, Clin
if keyword_set(niter) then opti = ' -i ' + STRC(niter) + ' '  else opti=' ' 

InpFN = gettmpfilename() 
PSFN =  gettmpfilename() 

if not keyword_set(m1) then cmd = 'mrs_alm_inpainting  -m2 -W ' + opti $
else cmd = 'mrs_alm_inpainting  -m1 ' + opti 

 if keyword_set(opt) then cmd = cmd + ' ' + opt + ' '
 
if keyword_set(clin) then cmd = cmd + ' -r ' +  RMSFN + ' -S ' + ClFN  + ' '  +   MapFN  + ' ' + MaskFN + ' ' + InpFN + ' ' + PSFN $
else  cmd = cmd + ' -r ' +  RMSFN +   ' '  +   MapFN  + ' ' + MaskFN + ' ' + InpFN + ' ' + PSFN

print, cmd
spawn, cmd

InpMap = mrs_read(InpFN)
InpPS = readfits(PSFN)

if keyword_set(plot) then begin
winbs, win=0
plotcl, InpPS
if keyword_set(Clin) then oplotcl, Clin, line=2, thick=2
 tvs, InpMap
 end

delete, MapFN
delete, RMSFN
delete, MaskFN
if keyword_set(Clin) then  delete, ClFN
delete, InpFN
delete, PSFN

DONE:

end

;================================

; THIS DOES NOT WORK
pro bad_sparse_master,  Data,  Mask,  RMSMap,  InpMap, InpPowSpec, inpiter=inpiter, niter=niter, MatMask=MatMask, TrueCl=TrueCl, PDat=PDat, plot=plot, lmax=lmax, optinp=optinp
; Data ;Healpix Masked data
; Mask = Healpix mask
COMMON C_PLANCK

nside=gnside(Data)
; print, nside
if not keyword_set(lmax) then lmax =3*nside
if not keyword_set(MatMask) then  MatMask  = mrs_matmask(Mask, lmax=lmax)
if not keyword_set(inpiter) then inpiter =40


PDat = mrs_powspec(Data, normval=PSNormVal)

master_niter=50
Npix = nside^2*12L
MinRMS = min(RMSMap)

InpMap = Data
PMap = PDat

Rea=0
if Rea EQ 1 then begin
Prea=0
Ni=20
for i=0,Ni-1 do begin
  Rea = randomn(seed, Npix) * RMSMap * Mask
  Prea = Prea  + mrs_powspec(Rea, normval=PSNormVal)
  end
  Prea = Prea / float(Ni)
 end
 
  PNoise = PDat
  PNoise[*] = MinRMS^2. / PSNormVal^2.
  PNoise = MatMask # PNoise
 
 
; print, MinRMS
Resi = 0
Map =0 

PSSol = mrs_get_master_ps( PDat, Mask, MatMask=MatMask, niter=master_niter, NoisePS=Prea,Tmat=Tmat)


for i=0,niter-1 do begin
    
    mrs_filter_inpainting, Data, Mask, RMSMap, clin=PSSol,  InpMap, InpPowSpec, niter=inpiter , opt=optinp
    ; mrs_filter_inpainting, Data, Mask, RMSMap,  InpMap, InpPowSpec, opt='-i40   '
    PSSol = InpPowSpec
    
    if Rea EQ 0 then Resi = (Data - InpMap) / RMSMap*MinRMS * Mask  $
    else Resi = (Data - InpMap)  * Mask
    
    Presi = mrs_powspec(Resi, normval=PSNormVal)
    
    PSSol = PSSol + Tmat # (Presi - PNoise)
    ind = where(PSSol LT 0, c)
    if c GT 0 then PSSol[ind] = 0
       
   if keyword_set(plot) then begin
      wset, 1
      plotcl, PNoise - Presi, tit='Resi'
      ; oplotcl, PSSol, line=2
   end
   
   FN = 'xx'+ STRC(i+1)+'.xdr'
   save, filename=FN, Resi, Presi, PNoise, Map, PMap, PSSol, InpMap, InpPowSpec 

   if keyword_set(plot) then begin
   mollview, window=2, /nest, InpMap, tit=STRC(i+1)
   wset, 0
   plotcl, PSSol, tit='Master + Inp'
   ; oplotcl, InpPowSpec, thick=3, line=4
  end
  
   ; info, InpPowSpec-PSSol
   if keyword_set(TrueCl) then begin
   		    if keyword_set(plot) then  oplotcl, TrueCl, thick=2, line=2
  		    DiffPS = InpPowSpec-TrueCl
  		    DiffPS = DiffPS[0:2* nside]
     	    Err = sigma(DiffPS)
  		    print, "Iter ", i+1, ",  Error = ", Err, ' ' , sigma(PSSol-TrueCl)
          end
 end

end
;================================

;================================
 
 pro test128
 cmb = rims('cmb128.fits')
 pcmb = mrs_powspec(cmb)
 rms = rims('rms128_snr_5_20.fits')*10.
 m = rims('wmap_mask128.fits')

r10 = rms
rmsmap=r10  
 Noise = randomn(seed, 128L^2*12) * r10  

 data = cmb + Noise
 PDat = mrs_powspec(data)

 datam = data*m

PDat = mrs_powspec(datam)


 Pxn = mrs_crosspowspec(datam, Noise)
 Pnx = mrs_crosspowspec(Noise , datam)

 
sparse_master,  datam,  M,  R10,  x1, x2, niter=20, MatMask=MatMask, TrueCl=pcmb, /plot
restore, /verb, "xx19.xdr"
data = datam
mask = m
Resi = (Data - InpMap)  * Mask
pr = mrs_powspec(resi)
pxr = mrs_crosspowspec(InpMap*m , resi)
px = mrs_powspec(InpMap*m)
pdat = mrs_powspec(datam)
pok = mrs_powspec(cmb*m)

nside=gnside(Data)
Npix = nside^2*12L
MinRMS = min(RMSMap)

 Rea = randomn(seed, Npix) * RMSMap * Mask
 PRea = mrs_powspec(Rea, normval=PSNormVal)
 Presi =  mrs_powspec(resi)

 PL = PDat
 PL[*] = MinRMS^2. / PSNormVal^2.
 PL = MatMask # PL
 
ResiNorm  = (Data - InpMap) / RMSMap*MinRMS * Mask
PResiNorm = mrs_powspec(ResiNorm)

Pz =   mrs_powspec(randomn(seed, Npix) * MinRMS * Mask)
 
Cm = cmb*m
PCm = mrs_powspec(cm)





Pmaster=mrs_get_master_ps( PCm, m, MatMask=MatMask,  niter=50)

Pmaster=mrs_get_master_ps( PDat, m, MatMask=MatMask,  niter=50)
plotcl, Pmaster

mrs_cole_inpainting, datam, m, r10, clin=pcmb,  x1, x2, opt=' -v -W  -m2 ', /plot
mrs_cole_inpainting, datam, m, r10, clin=Pmaster,  x1, x2, opt=' -v -W  -m2 ', /plot
sparse_master,  datam,  M,  R10,  x1, x2, niter=5, MatMask=MatMask, TrueCl=pcmb


 mrs_write, 'data128.fits', datam
 mrs_write, 'rms.fits', rms * 10
 writefits, 'ps_cmb128.fits', pcmb

pd = mrs_powspec(data*m)
 cmd = 'mrs_alm_inpainting -v -m1 -r rms.fits  -S ps_cmb128.fits  data128.fits wmap_mask128.fits tt ttp'
 spawn, cmd
 t = rims('tt.fits')
tvs, t
 
 
cmd = 'mrs_alm_inpainting -v  -m1  -r rms.fits   data128.fits wmap_mask128.fits tt ttp'
spawn, cmd
t = rims('tt.fits')
tvs, t

I2 = mrs_alm_inpainting(datam, mask=m, Niter=20,  OutPowSpec=OutPowSpec,  lmax=lmax,  Opt=Opt)


 end
 
 
 ;================================

pro test512
restore, /verb, "xxwmap.xdr"
m = rims('wmap-sap_gmca-june11-finalmask_galps.fits')
mask=m
rms = readfits('wmap-sap_gmca-june11-compsep_estim_src_std.fits')
rms = rms[*,0]
rmsmap=rms
nside=512L
cl=0
c= getcmb(Cl=cl, nside=nside)
cmb = c * 1.e-3
Cl = cl * 1.e-6
pcmb = mrs_powspec(cmb)
Noise = randomn(seed, nside^2*12)*rms  
data = cmb + Noise
datam = data*m
PDat = mrs_powspec(datam)

pn =readfits('pn.fits')
runrea = 0
if keyword_set(runrea) then begin
Nrea=50
PN=0
for i=0, Nrea-1 do begin N = randomn(seed, nside^2*12)*rms &  Pn = pn  + mrs_powspec(n) & end
Pn = Pn / float(Nrea)
writefits, 'pn.fits', pn
end


PSSol = mrs_get_master_ps( PDat, Mask, MatMask=MatMask, niter=20, NoisePS=Pn)

 mrs_sparse_master, datam, Mask, RMS, InpMap, InpPS,  /plot,  MatMask=MatMask, NoisePS=pn
 
 
;  sparse_master,  datam,  M,  rms,  x1, x2, niter=3, MatMask=MatMask, TrueCl=cl, /plot
save, filename='xxwmap1.xdr',  datam,  M,  rms,  cmb,  i4, psi4, pn, nside, data, MatMask

end


