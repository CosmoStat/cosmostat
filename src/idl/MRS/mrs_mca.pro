;+
; NAME:
;      mrs_mca
;	
; PURPOSE:
;	Apply the sparse component analysis method called Morphological Component Analysis, including 
;	a hard thresholding with linear decreasing threhold, on a spherical map in Healpix representation 
;	(NESTED format) using several basis decompositions and there transforms on the sphere selected 
;	in the following list:
;
;		1: Isotropic Undecimated Wavelet
;		2: Pyramidal Wavelet
;		3: Othogonal Wavelet Transform (on each face)
;		4: ALM
;		5: Dirac
;		6: Curvelet
;		7: DCT (on each face)
;		8: A Trou Wavelet (on each face)
;		Default transforms are Pyramidal Wavelet and ALM
;
;	With the selection of only one transform and the use of a mask, MCA will make an inpainting of the input data.
;		
; CALLING:
;		mrs_mca, data_in, data_out, Bounded=Bounded, residual=residual, CstSigma=CstSigma, SelectTrans=SelectTrans, Positivity=Positivity, niter=niter, 
;				mad=mad, mom=mom, expo=expo, LastThreshold=LastThreshold, FirstThreshold=FirstThreshold, SigmaNoise=SigmaNoise, NbrScale=NbrScale, lmax=lmax, 
;				Mask=Mask, tabNameTrans=tabNameTrans, FirstWTDetectScale=FirstWTDetectScale, soft=soft, DCTblocksize=DCTblocksize, fit=fit, nomean=nomean
;
; INPUTS:
;		data_in : Input 1D IDL array of a Healpix map, image to be analysed.
;
; OUTPUTS:
;		data_out : Output 2D IDL array[*, NbTrans] of Healpix maps, components estimated from data_in. NbTrans is the number 
;					of selected transforms (via SelectTrans keyword), by default there are 2 transforms.
;         
; INPUT KEYWORDS:
;		SelectTrans : 1D int array with the code number of the selected transforms. SelectTrans[i] must be >=1 and <=9. Default value: SelectTrans = [2,4]
;		niter : int, iteration number of the MCA algorithm. Default value is 10.
;		Mask : 1D IDL array of Healpix map, mask applied to data_in. Inpainting on the masked areas.
;		expo : scalar, if set use an exponential decreasing thresholding instead of linear decreasing thresholding.
;		mom : scalar, if set use a linear decreasing thresholding with MOM as a first threshold.
;		mad : scalar, if set use a linear decreasing thresholding with MAD as a first threshold.
;		fit : scalar, if set fit the threshold levels to ALM decomposition of data_in.
;		soft : scalar, if set use soft thresholding instead of hard thresholding.
;		SigmaNoise : float, standard deviation of the noise, assumed gaussian. Default value is 1.
;		NbrScale : int, number of scale decompositions for wavelet transfroms. Default value is 5.
;		lmax : int, maximum l number of spherical harmonics. Default value is 3*nside, max value is 3000.
;		Bounded : scalar, if set constraints the reconstructed components of data_out to be bounded by the min and max of data_in.
;		Positivity : scalar, if set constraints the reconstructed components of data_out to be positive.
;		CstSigma : scalar, if set and if a mask is applied, constraints the decompositions coefficients to have the same standard deviation inside and outside the masked area.
;		nomean : scalar, if set remove the mean of the reconstructed components of data_out. Work only with keywords mask and CstSigma.
;		DCTblocksize : int, size of the blocks for DCT transform (if selected). Default value is the nside parameter of data_in.
;		FirstWTDetectScale : int, for isotropic wavelet only, set all wavelet coefficients from Scale < FirstWTDetectScale to 0.
;
; INPUT/OUTPUT KEYWORDS:
;		LastThreshold : float, last threshold level. Default value is 0.
;		FirstThreshold : float, first threshold level. Default is automatically estimated.
;
; OUTPUT KEYWORDS:
;		residual : 1D IDL array of a Healpix map, final residual.
;		tabNameTrans : string array, list of the possible transforms. 
;						tabNameTrans = ['Unknown', 'Isotropic Undecimated Wavelet', 'Pyramidal Wavelet', 'Orthogonal Wavelet Transform', 'ALM', 'DIRAC', 'Curvelet', 'DCT', 'a trous WT']
;
; EXAMPLE:
;       Compute the MCA on a Healpix image data, considering 4 components: Curvelet, ALM, Pyramidal Wavelet and Dirac
;          mrs_mca, Data, Components, SelectTrans=[2,4,5,6]
;         
; HISTORY:
;	Written: Jerome Bobin, 2005
;	2005 File creation
;--------------------------------------------------------------------------------------------------------


;=====================================================

function maxalm, s, lmax=lmax
mrs_almtrans, s, S_ALM, /norm, lmax=lmax
 ; info, S_ALM.alm[1:*,*]

m = max( ABS(double(S_ALM.ALM[1:*,*])))*sqrt(double(2.))
return, m
end

;=====================================================
function get_lambda_alm, s, lmax=lmax,niter=niter
mrs_almtrans, s, S_ALM, /norm, lmax=lmax
 ; info, S_ALM.alm[1:*,*]
ind = reverse(sort(abs(S_ALM.alm[*,*])))
ntot = n_elements(ind)
el = findgen(niter)*ntot/niter
el1 = ind(el)
lambda_threshold  = abs(S_ALM.alm(el1))
return, lambda_threshold
end


function maxdct, s, blocksize=blocksize
mrs_dcttrans, s, TransDct, blockSize=BlockSize
TransDct[0,0,*]=0
m = max( ABS(TransDct))
return, m
end

;=====================================================

function maxpwt, s, NbrScale=NbrScale
m=0
mrs_pwttrans, s, Trans, NbrScale=NbrScale
for j =0,NbrScale-2 do begin
    Scale = mrs_wtget(Trans,j,NormVal=NormVal)
    ms = mrs_absmax(Scale) / NormVal
    if m lt ms then m = ms
    end
return, m
end

;=====================================================

function maxwt, s, NbrScale=NbrScale
m=0
mrs_wttrans, s, Trans, NbrScale=NbrScale
for j =0,NbrScale-2 do begin
    Scale = mrs_wtget(Trans,j,NormVal=NormVal)
    ms = mrs_absmax(Scale) / NormVal
    if m lt ms then m = ms
    end
return, m
end
;=====================================================

function maxowt, s, NbrScale=NbrScale
mrs_owttrans, s,Trans,  NbrScale=NbrScale
m = mrs_absmax(Trans.coef)
return, m
end

;=====================================================

function maxawt, s, NbrScale=NbrScale
mrs_attrans, s, Trans,  NbrScale=NbrScale, /modif
m = mrs_absmax(Trans.coef)
return, m
end

;=====================================================

function maxcur, s, NbrScale=NbrScale
m=0
mrs_curtrans, s, Trans, NbrScale=NbrScale, /overlap, opt=' -O ',/undec
;for j =0,NbrScale-2 do begin
;    nb_ech_rid = Trans.TABNBRSCALERID[j]
;    for k=0,nb_ech_rid-1 do begin
   FirstScale=1   
for s2d = 0,Trans.NbrScale-2 do begin
  for s1d = 0,Trans.TabNbrScaleRid[s2d]-1 do begin
    Band = mrs_curget(Trans, s2d, s1d)
    
    vs = size(Trans.TabNorm)
    NxTNorm = vs[1]
    NyTNorm = vs[2]
    
    if s2d LT FirstScale-1 then Band[*] = 0 $
    else BEGIN    
      s2d1 = s2d
      s1d1 = s1d
      if s2d1 GE NyTNorm then s2d1 = NyTNorm - 1
      if s1d1 GE NxTNorm then s1d1 = NxTNorm - 1
      if s1d1 EQ Trans.TabNbrScaleRid[s2d]-1 and Trans.TabNbrScaleRid[s2d] lt NxTNorm then s1d1 = s1d1-1
      NormVal = Trans.TabNorm[s1d1,s2d1] * sqrt(Trans.TabBlockSize[s2d])
      
    ;Scale = mrs_curget(Trans,j,k);ImaMad=NormVal)
    ms = mrs_absmax(Band) / NormVal
    if m lt ms then m = ms
     endelse
    endfor
endfor
return, m
end

;=====================================================

pro mca_proj_cur, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, lmax=lmax, SigmaNoise=SigmaNoise, Soft=Soft

mrs_curtrans, Ima, Trans, NbrScale=NbrScale,  /overlap, opt=' -O ', /undec
MaxCoef=0.
; Noise estimation
if not  keyword_set(MAD)  and not keyword_set(SigmaNoise) then begin
   Scale = mrs_wtget(Trans.WT, 0,NormVal=NormVal)
   SigmaNoise=get_noise(Scale) / NormVal
   print, 'Estimated SigmaNoise  = ', SigmaNoise
end
;if not keyword_set(NbrScale) then NbrScale=4
;if not keyword_set(NSigma) then NSigma=3.
if not keyword_set(FirstScale) then FirstScale=1
if not keyword_set(NSigma) then NSigma=40.

vs = size(Trans.TabNorm)
NxTNorm = vs[1]
NyTNorm = vs[2]

TabThreshold = fltarr(NxTNorm, NyTNorm)
TabNoise = fltarr(NxTNorm, NyTNorm)

NSigma = Lambda

for s2d = 0,Trans.NbrScale-2 do begin
  for s1d = 0,Trans.TabNbrScaleRid[s2d]-1 do begin
    Band = mrs_curget(Trans, s2d, s1d)
    if s2d LT FirstScale-1 then Band[*] = 0 $
    else BEGIN
      s2d1 = s2d
      s1d1 = s1d
      if s2d1 GE NyTNorm then s2d1 = NyTNorm - 1
      if s1d1 GE NxTNorm then s1d1 = NxTNorm - 1
      if s1d1 EQ Trans.TabNbrScaleRid[s2d]-1 and Trans.TabNbrScaleRid[s2d] lt NxTNorm then s1d1 = s1d1-1
      NormVal = Trans.TabNorm[s1d1,s2d1] * sqrt(Trans.TabBlockSize[s2d])
      
      ; Threshold level  
      if s2d EQ 0 and s1d EQ 0 then NS=NSigma+0.1 else NS=NSigma
      ;if not keyword_set(MAD) then 
      ThresholdLevel = NS*SigmaNoise*NormVal 
      ;else 
      ;ThresholdLevel = NS * median( abs(Band)) / 0.6745
      TabThreshold[s1d1,s2d1] = ThresholdLevel
      TabNoise[s1d1,s2d1] = ThresholdLevel / NS
      
      Band2 = mrs_absthreshold(Band, ThresholdLevel)
      mrs_curput, Trans, Band2, s2d, s1d
      Band = mrs_diff(Band, Band2)
      ms = mrs_absmax(Band)  / NormVal
      if MaxCoef lt ms then MaxCoef = ms
      END
    
  end
end

print, 'Reconstruction '
mrs_currec, Trans, Ima 





end

;=====================================================

function get_nbr_trans, TabC
vs = size(TabC)
if type_code(TabC[0]) EQ 8 then  NbrTrans = vs[1] else begin
  if vs[0] EQ 1 then  NbrTrans = 1 else NbrTrans = vs[2]
end
return, NbrTrans
end

;=====================================================

function mca_add, TabC
Data = TabC[0]
NbrTrans = get_nbr_trans(TabC)
for i=1,NbrTrans-1 do begin
   if type_code(Data) EQ 8 then Data.T_Sky =   Data.T_Sky + TabC[i].T_Sky $
   else Data = Data + TabC[i]
end
return, Data
end

;=====================================================

function mca_residual, Data, TabC, mask=mask
Residual = Data
NbrTrans = get_nbr_trans(TabC)
; print, 'NbrTrans=',NbrTrans
if  type_code(Data) EQ 8 then begin
  if NbrTrans EQ 1 then  Residual = mrs_diff(Residual, TabC[0]) $
  else for i=0,NbrTrans-1 do  Residual = mrs_diff(Residual, TabC[i])
end else begin
  if NbrTrans EQ 1 then Residual =  Residual - TabC[*] $
  else for i=0,NbrTrans-1 do  Residual =  Residual - TabC[*,i]
end

if keyword_set(mask) then begin
   if type_code(Data) EQ 8 then Residual.T_Sky = Residual.T_Sky * Mask.T_Sky $
   else Residual =  Residual * Mask
end

return, Residual
end

;=====================================================
 
pro mca_proj_dirac, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, SigmaNoise=SigmaNoise 
Scale = Ima 
NormVal = 1.
if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
else ThresholdLevel = 5. * mrs_mad(Scale) 
TScale = mrs_absthreshold(Scale, ThresholdLevel)
Scale = mrs_diff(Scale, TScale)
MaxCoef = mrs_absmax(Scale)  / NormVal
; print, 'PROJ DIRAC:  Count = ', Count, ' MaxCoef = ', MaxCoef 
end

;=====================================================

pro mca_proj_dct, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef,  SigmaNoise=SigmaNoise , soft=soft, BlockSize=BlockSize, l2=l2
; mrs_info, ima, mes='Proj DCT' 
mrs_dcttrans, Ima, Trans, blocksize=BlockSize
NormVal = 1.
if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
else ThresholdLevel = 5. * mrs_mad(Trans) 
ZeroDCT = Trans[0,0,*]
Trans[0,0,*] = 0
TScale = mrs_absthreshold(Trans, ThresholdLevel, l2=l2)
TScale[0,0,*] = ZeroDCT
mrs_dctrec, TScale , Ima,  blocksize=BlockSize
TScale[0,0,*] = 0
Trans =  Trans - TScale
MaxCoef = mrs_absmax(Trans )  / NormVal
; mrs_info, ima, mes='  ==> Rec DCT' 
 end

;=====================================================

pro mca_proj_wave, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, lmax=lmax, SigmaNoise=SigmaNoise, pyr=pyr, FirstWTDetectScale=FirstWTDetectScale, soft=soft
; if not keyword_set(Pyr) then mrs_wttrans, Ima, Trans, NbrScale=NbrScale, lmax=lmax, /Healpix_with_Glesp   $
if not keyword_set(Pyr) then mrs_wttrans, Ima, Trans, NbrScale=NbrScale, lmax=lmax     $
else mrs_pwttrans, Ima, Trans, NbrScale=NbrScale, lmax=lmax
MaxCoef=0.
if not keyword_set(FirstWTDetectScale) then FirstWTDetectScale=0

for j =0,NbrScale-1 do begin
    Scale = mrs_wtget(Trans, j, NormVal=NormVal)
    if j LT FirstWTDetectScale then begin
       mrs_set, Scale, 0. 
       mrs_wtput, Trans, Scale, j
    end else begin
       if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
       else ThresholdLevel = 5. * mrs_mad(Scale)
       ; print, " WT scale ", j+1, "   ThresholdLevel = ", ThresholdLevel 
        TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
        mrs_wtput, Trans, TScale, j
        Scale = mrs_diff(Scale, TScale)
        ms = mrs_absmax(Scale)  / NormVal
        if MaxCoef lt ms then MaxCoef = ms
       end
    end
if not keyword_set(Pyr) then mrs_wtrec,  Trans, Ima $
else mrs_pwtrec, Trans, Ima 
end


;=====================================================

pro mca_proj_owt, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, SigmaNoise=SigmaNoise , soft=soft, l2=l2

mrs_owttrans, Ima, Trans,  NbrScale=NbrScale
MaxCoef=0.
EndX = Trans.nx
EndY = Trans.ny
BorderSize = 0
NormVal = 1.
for j =0,NbrScale-2 do begin
    HalfX = EndX/2
    HalfY = EndY/2
    for d=0,2 do begin
      if d EQ 0 then Scale = Trans.coef[HalfX+BorderSize:EndX-1-BorderSize, BorderSize:HalfY-1-BorderSize, *] $
      else if d EQ 1 then Scale = Trans.coef[BorderSize:HalfX-1-BorderSize,HalfY+BorderSize:EndY-1-BorderSize, *] $
      else Scale = Trans.coef[HalfX+BorderSize:EndX-1-BorderSize, HalfY+BorderSize:EndY-1-BorderSize, *]
   
      if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
      else ThresholdLevel = 5. * mrs_mad(Scale)
      TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
      if d EQ 0 then Trans.coef[HalfX+BorderSize:EndX-1-BorderSize, BorderSize:HalfY-1-BorderSize, *] = TScale $
      else if d EQ 1 then  Trans.coef[BorderSize:HalfX-1-BorderSize,HalfY+BorderSize:EndY-1-BorderSize, *] = TScale $
      else Trans.coef[HalfX+BorderSize:EndX-1-BorderSize, HalfY+BorderSize:EndY-1-BorderSize, *]= TScale 
      ms = mrs_absmax(Scale)  / NormVal
      if MaxCoef lt ms then MaxCoef = ms
   end
    EndX = EndX / 2
    EndY = EndY / 2
end
j = NbrScale-1
Scale = Trans.coef[BorderSize:EndX-1-BorderSize, BorderSize:EndY-1-BorderSize, *]
if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
else ThresholdLevel = 5. * mrs_mad(Scale)
TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
Trans.coef[BorderSize:EndX-1-BorderSize, BorderSize:EndY-1-BorderSize, *] = TScale
ms = mrs_absmax(Scale)  / NormVal
if MaxCoef lt ms then MaxCoef = ms
    

mrs_owtrec, Trans, Ima
end
;=====================================================
; mrs_mca, map128, tabc, mask=mask128, niter=30,residual=residual,selecttrans=[8], /expo, NbrScale=5
pro mca_proj_at, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, SigmaNoise=SigmaNoise, KillLastScale=KillLastScale, soft=soft, l2=l2
; mrs_attrans, Ima, Trans,  NbrScale=NbrScale, /modif
mrs_attrans, Ima, Trans,  NbrScale=3, /modif
MaxCoef=0.
BorderSize = 0
NormVal = 1.
vs = size(Ima)
npix = vs[1]
nside = npix2nside(npix)
Nx = nside
Ny = nside

for j =0,Trans.NbrScale-2 do begin
      Scale = Trans.coef[BorderSize:Nx-1-BorderSize, BorderSize:Ny-1-BorderSize, j, *]
      NormVal= Trans.TabNorm[j]
      if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
      else ThresholdLevel = 5. * mrs_mad(Scale)
      TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
      ;ind = where(TScale LT 0, c)
      ;if c GT 0 then TScale[ind] = 0
      Trans.coef[BorderSize:Nx-1-BorderSize, BorderSize:Ny-1-BorderSize, j, *] = TScale 
     
     ms = mrs_absmax(Scale)  / NormVal
     if MaxCoef lt ms then MaxCoef = ms
 end
j = Trans.NbrScale-1
NormVal= Trans.TabNorm[j]
Scale = Trans.coef[BorderSize:Nx-1-BorderSize, BorderSize:Ny-1-BorderSize, j, *]
if keyword_set(KillLastScale) then Scale[*]=0 ; $
;else begin
;  if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
;  else ThresholdLevel = 5. * mrs_mad(Scale)
;  TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
;  Trans.coef[BorderSize:Nx-1-BorderSize, BorderSize:Ny-1-BorderSize, j, *] = TScale
;  ms = mrs_absmax(Scale)  / NormVal
;  if MaxCoef lt ms then MaxCoef = ms
;end
;info, trans.coef 
mrs_atrec, Trans, Ima
end


;=====================================================
;  mrs_mca, gd64, tc1, selecttrans=[4], /cmb, niter=200,  mask=gm64, /expo
pro mca_proj_alm, Ima, S_ALM=S_ALM, Lambda, Mad=Mad, MaxCoef=MaxCoef,  lmax=lmax, SigmaNoise=SigmaNoise, cmb=cmb, soft=soft, LambdaCst=LambdaCst

if not keyword_set(CMB) then begin
  mrs_almtrans, Ima, S_ALM, /norm, lmax=lmax
  ; info, S_ALM.alm[1:*,*]
  Scale = S_ALM.alm[1:*,*]
  TScale = Scale
  NormVal = double(1.) / sqrt(double(2.))
  if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
  else ThresholdLevel = 5. * median( abs(Scale)) / 0.6745
 
  ind = where (ABS(Scale) GT ThresholdLevel, c)
  print, " MaxAbs(Scale) = ", max(ABS(Scale)), " Threshold = ", ThresholdLevel,  ' > = ', c
 
  TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
  S_ALM.alm[1:*,*] = TScale
  Scale = Scale - TScale
  MaxCoef = max(ABS(Scale)) / NormVal
end else begin
  mrs_almtrans, Ima, S_ALM, /norm, lmax=lmax, /tab
  Zero = S_ALM.alm[0,0,*]
  Scale = S_ALM.alm
  TScale = Scale
  NormVal = 1. / sqrt(2.)
  ThresholdLevel = Lambda*SigmaNoise*NormVal
 
  TScale = mrs_absthreshold(Scale, 10*LambdaCst, /l2)
  ; TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
  TScale[0,0,*] = Zero
  S_ALM.alm = TScale
  Scale = Scale - TScale
  MaxCoef = max(ABS(Scale)) / NormVal
  
  P = mrs_alm2spec(S_ALM)
  Sig= 1.
  mrs_pswtfil, P, Pfil, Nscale=Nscale, SigmaNoise=Sig, NSigma=5., NbrIter=1, L1Regul=L1Regul 
  Pfil = Pfil + Sig
 
  if cmb EQ 2 then begin 
  for l=3,S_ALM.lmax do begin
     SigNoise = sqrt(Pfil[l]/2.)
     SigRe = sigma(S_ALM.alm[l, 0:S_ALM.TabNbrM[l]-1, 0])
     SigIM = sigma(S_ALM.alm[l, 0:S_ALM.TabNbrM[l]-1, 1])
     if SigRe GT  0 then  CoefRe = 1 + LambdaCst *( SigNoise / SigRe - 1.) else CoefRe =  0
     if SigIm GT  0 then  CoefIM = 1 + LambdaCst *( SigNoise / SigIm - 1.) else CoefIM =  0
     S_ALM.alm[l, 0:S_ALM.TabNbrM[l]-1, 0] = S_ALM.alm[l, 0:S_ALM.TabNbrM[l]-1, 0] * CoefRe
     S_ALM.alm[l, 0:S_ALM.TabNbrM[l]-1, 1] = S_ALM.alm[l, 0:S_ALM.TabNbrM[l]-1, 1] * CoefIM
  end
  end
  
  S_ALM.alm[0,0,*] = Zero
end



; S_ALM.alm[0,*] = 0
mrs_almrec, S_ALM, Ima

;mrs_almtrans, Ima, Z, /norm, lmax=lmax, /tab, /PSP
;loadct, 15
;load, alog(z.alm[*,*,0]+1)
; P = mrs_alm2spec(S_ALM)

end

;=====================================================

function maxcmblet, s, lmax=lmax, wtp=wtp
MaxCoef=0.
NormVal = 0.43
cmblet_trans, s, S_CMB, lmax=lmax, optwt=optwt, /reim, /NoWT
MaxCoef = max(ABS(S_CMB.ReIMA)) / NormVal
ms = max(ABS(S_CMB.ImIMA)) / NormVal
if MaxCoef lt ms then MaxCoef = ms
print, MaxCoef

if keyword_set(wtp) then begin
P = sqrt(S_CMB.spec1d)
ms = sqrt(max(ABS(P)))
if MaxCoef lt ms then MaxCoef = ms
end

return, MaxCoef
end

;=====================================================

pro mca_proj_cmblet, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef,  lmax=lmax, SigmaNoise=SigmaNoise, optwt=optwt , T=T, wtp=wtp
; S_CMB = {WtSpec: WtSpec, PowSpecIma: PowSpecIma, Phase: Phase, log: log, ALM : Alm, mu: mu, psi1: psi1, NoWt:NoWt} 
Soft=1
MaxCoef=0.
NormVal = 0.43
cmblet_trans, Ima, S_CMB, lmax=lmax, optwt=optwt, /reim, /NoWT

; REAL part
Scale =  S_CMB.ReIMA
TScale = Scale
if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
else ThresholdLevel = 5. * median( abs(Scale)) / 0.6745
TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
S_CMB.ReIMA = TScale
Scale = Scale - TScale
MaxCoef = max(ABS(Scale)) / NormVal

;Imaginary part
Scale = S_CMB.ImIMA
TScale = Scale
if not keyword_set(MAD) then ThresholdLevel = Lambda*SigmaNoise*NormVal $
else ThresholdLevel = 5. * median( abs(Scale)) / 0.6745
TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
S_CMB.ImIMA = TScale
Scale = Scale - TScale
Ms = max(ABS(Scale)) / NormVal
if MaxCoef lt ms then MaxCoef = ms

; Wavelet part
if keyword_set(wtp) and Lambda NE 0 then begin
P = S_CMB.spec1d
; NoiseSpectrum = P - P / (Lambda + 1.)
mrs_pswtfil, P, Pfil, Nscale=Nscale, SigmaNoise=Lambda, NSigma=Lamda, NbrIter=1, L1Regul=L1Regul
; Pfil = Pfil + 1.
Pfil[0]=0 
S_CMB.spec1d = Pfil
 
;Seuil = Lambda^2
;S_CMB.spec1d = mrs_absthreshold(P, Seuil, /soft)
end

cmblet_rec, S_CMB, Ima,/zeromean, biais=biais

end

;=====================================================


function old_maxcmblet, Ima, lmax=lmax
mrs_almtrans, Ima, ALM, lmax=lmax, /tab,  /norm
P = mrs_alm2spec(alm)
MaxCoef = sqrt(max(P[1:*]))
return, MaxCoef
end
;=====================================================

pro old_mca_proj_cmblet, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef,  lmax=lmax, SigmaNoise=SigmaNoise, optwt=optwt , T=T
; S_CMB = {WtSpec: WtSpec, PowSpecIma: PowSpecIma, Phase: Phase, log: log, ALM : Alm, mu: mu, psi1: psi1, NoWt:NoWt} 
MaxCoef=0.
NormVal = 1.

mrs_almtrans, Ima, ALM, lmax=lmax, /tab,  /norm 
P = mrs_alm2spec(alm)
NoiseSpectrum = P - P / (Lambda + 1.)
Sig= sqrt(NoiseSpectrum)
mrs_pswtfil, P, Pfil, Nscale=Nscale, SigmaNoise=Sig, NSigma=Lamda, NbrIter=1, L1Regul=L1Regul 

; P = P - NoiseSpectrum
ind = where(Pfil LT 0, c)
if c GT 0 then Pfil[ind] = 0
mrs_wiener, Ima, NoiseSpectrum, Ima, SignalPrior=Pfil, alm=ALM, Spec1D=Spec1D
; plotcl, p, line=2
; oplotcl, pfil, line=1
mrs_info, pfil, mes='PFIL'
mrs_tv, ima
end

;=====================================================
 
function  mrs_get_mom, NbrTrans, TabMaxTrans
if NbrTrans EQ 1 then FirstThreshold=TabMaxTrans[0] $
   else begin
      ind = sort(TabMaxTrans)
      ind = reverse(ind)
      F1 = TabMaxTrans[ ind[0]] 
      F2 = TabMaxTrans[ ind[1]]
      S1 = F2 + (F1-F2)*0.05
      S2 = F2 + F2 * 0.05
      FirstThreshold = MIN([S1,S2])
   end
   print, '       MOM: ', TabMaxTrans, ' == ', FirstThreshold
return, FirstThreshold
end
;=====================================================

; mrs_mca, map350, tc, selecttrans=[7,8], nbrscale=6, niter=30, /expo, DCTblocksize=32, positivity=[1,1], SigmaNoise=0.0123198, /KillLastScale
pro mrs_mca, Bounded=Bounded, residual=residual, data, tabC, CstSigma=CstSigma, SelectTrans=SelectTrans, Positivity=Positivity, niter=niter, mad=mad, disp=disp, mom=mom, expo=expo, LastThreshold=LastThreshold, FirstThreshold=FirstThreshold, SigmaNoise=SigmaNoise, NbrScale=NbrScale, lmax=lmax, Mask=Mask, cmb=cmb, tabNameTrans=tabNameTrans, FirstWTDetectScale=FirstWTDetectScale, wtp=wtp, soft=soft, log=log, l2=l2, IterResi=IterResi, film=film,  DCTblocksize=DCTblocksize, KillLastScale=KillLastScale,fit=fit, nomean=nomean, FirstIterConstr=FirstIterConstr, TabSupport=TabSupport, FirstGuess=FirstGuess

if not keyword_set(FirstIterConstr)  then FirstIterConstr=0
if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_mca, Data, TabC, SelectTrans=SelectTrans'
	print, '                     SelectTrans=array with values: '
	print,  '                       1: Wavelet '
	print,  '                       2: Pyramidal Wavelet '
	print,  '                       3: Othogonal Wavelet Transform (on each face)' 
	print,  '                       4: ALM '
  	print,  '                       5: Dirac'
 	print,  '                       6: Curvelet'
	print,  '                       7: DCT (on each face)'
	print,  '                       8: Isotropic WT (on each face)'
	print,  '                      Default transforms are [2,4]'
        goto, DONE
end

nside=0
if type_code(Data) EQ 8 then  GLESP=1 else begin
   GLESP=0
   npix = (size(Data))[1]
   nside = npix2nside(npix)
end

if keyword_set(TabSupport)  then print, 'SUPPORT'

if keyword_set(Mask) and keyword_set(CstSigma)  then begin
     if GLESP EQ 1 then begin
        indMaskOK = where (Mask.t_sky NE 0)
	SigMaskOK = sigma( data.t_sky[indMaskOK] )
	indMaskZero =  where (Mask.t_sky EQ 0)
      end else begin
        indMaskOK = where (Mask NE 0)
	SigMaskOK = sigma( data[indMaskOK] )
	indMaskZero =  where (Mask EQ 0)
      end
end

if keyword_set(Bounded) then begin
     if GLESP EQ 1 then begin
 	BoundMin = min(data.t_sky)
	BoundMax = max(data.t_sky)
     end else begin
 	BoundMin = min(data)
	BoundMax = max(data)
     end
end


NbrTotalTrans=8
tabNameTrans=strarr(NbrTotalTrans+1)
tabNameTrans = ['Unknown', 'Wavelet',  'Pyramidal Wavelet', 'Orthogonal Wavelet Transform', 'ALM', 'DIRAC', 'Curvelet', 'DCT', 'a trous WT', 'CMBLET']
T_Wave = 1
T_PWave = 2
T_OWT = 3
T_ALM = 4
T_DIRAC = 5
T_CUR = 6
T_DCT = 7
T_AT = 8
T_CMBLET = 9

TabTrans=intarr(NbrTotalTrans)
NbrTrans=0

T_LINEAR=1
T_MOM=2
T_MAD=3
T_EXP=4
T_FIT=5

TypeThreshold=T_LINEAR
NbrIterMax=10
if not keyword_set(niter) then niter=NbrIterMax
if keyword_set(mad) then TypeThreshold=T_MAD
if keyword_set(mom) then TypeThreshold=T_MOM
if keyword_set(expo) then TypeThreshold=T_EXP
if keyword_set(fit) then TypeThreshold=T_FIT


if not keyword_set(LastThreshold) then LastThreshold=0
if not keyword_set(SigmaNoise) then SigmaNoise=1.
if not keyword_set(NbrScale) then NbrScale=5
if not keyword_set(SelectTrans) then  SelectTrans = [2,4]   ; ALM and Wavelets are seclected by default
vs =size(SelectTrans)
NbrTrans = vs[1]
TabTrans = SelectTrans
for k=0,NbrTrans-1 do begin
  if TabTrans[k] LT 1 or TabTrans[k] GT NbrTotalTrans then begin
      print, 'Error: selected transforms must be between 1 and ', NbrTotalTrans
      goto, DONE 
  end
end
 
; DCTblocksize=nside/2
if not keyword_set(DCTblocksize) then  DCTblocksize=nside

if type_code(Data) EQ 8 then begin
   TabC = replicate(data, NbrTrans) 
   for k=0,NbrTrans-1 do tabC[k].t_sky = 0.
end else begin
   vs=size(Data)
   Npix=vs[1]
   TabC = fltarr(Npix, NbrTrans)
   TabC[*,*]=0.
end
if keyword_set(FirstGuess) then TabC = FirstGuess
; Compute the residual
Residual = mca_residual(Data, TabC, mask=mask)

print, 'NbrTrans = ', NbrTrans
TabMaxTrans=dblarr(NbrTrans)
if not keyword_set(FirstThreshold) then BEGIN
   for k=0,NbrTrans-1 do begin
       if TabTrans[k] EQ T_Wave then TabMaxTrans[k] = maxwt(Residual, NbrScale=NbrScale)
       if TabTrans[k] EQ T_PWave then TabMaxTrans[k] = maxpwt(Residual, NbrScale=NbrScale)
       if TabTrans[k] EQ T_OWT then TabMaxTrans[k] = maxowt(Residual, NbrScale=NbrScale)
       if TabTrans[k] EQ T_ALM then TabMaxTrans[k] = maxalm(Residual)
       if TabTrans[k] EQ T_CMBLET then TabMaxTrans[k] = maxcmblet(Residual,lmax=lmax, wtp=wtp)
       if TabTrans[k] EQ T_DIRAC then TabMaxTrans[k] = mrs_absmax(Residual)
       if TabTrans[k] EQ T_CUR then TabMaxTrans[k] = maxcur(Residual,nbrScale=nbrScale)
       if TabTrans[k] EQ T_DCT then TabMaxTrans[k] = maxdct(Residual,blocksize=DCTblocksize)
       if TabTrans[k] EQ T_AT then TabMaxTrans[k] = maxawt(Residual,NbrScale=NbrScale)
       print, 'Transform ',  tabNameTrans[ TabTrans[k] ], ' Max = ', TabMaxTrans[k]
   end
   FirstThreshold = double(mrs_get_mom(NbrTrans, TabMaxTrans))
   FirstThreshold = FirstThreshold / SigmaNoise
END

DeltaThreshold = FirstThreshold - LastThreshold
StepL = DeltaThreshold / float(NIter -1)
Lambda = FirstThreshold


for i=0,NbrTrans-1 do print, 'Selected Transform ', i+1,  ' ==> ', tabNameTrans[TabTrans[i]]
print, 'Nbr Iter = ', niter
if TypeThreshold EQ T_MAD then print, ' TypeThreshold = MAD Threshold'
if TypeThreshold EQ T_MOM then print, ' TypeThreshold = MOM Threshold'
if TypeThreshold EQ T_EXP then print, ' TypeThreshold = Exponential Decreasing Threshold'
if TypeThreshold EQ T_FIT then begin
     print, ' TypeThreshold = Fit alm threshold'
     threshold_fit = get_lambda_alm( data, lmax=lmax,niter=niter)
     print,threshold_fit
endif

print, 'FirstThreshold = ', FirstThreshold
print, 'LastThreshold = ', LastThreshold
print, 'StepLambda = ', StepL
print, 'FirstIterConstr = ', FirstIterConstr
 

; If IterResi is equal to 1 then the MCA iteration for the component k will be
;   x_k^{n+1} = x_k^{n} +  T_k^{-1}  HardThresh ( T_k Residual) 
;  If IterResi is equal to0 then the MCA iteration for the component k will be
;   x_k^{n+1} = T_k^{-1} [ HardThresh ( T_k  [ x_k^{n} + Residual] ) ]
if not keyword_set(IterResi) then IterResi=0
for i=0,niter-1 do BEGIN
   ; LambdaCst goes from 1 to zero. used in mca_proj_alm when /cmb option is set
   LambdaCst = float(niter-1-i) / float(niter-1)
   if i GT 0 then BEGIN 
      if TypeThreshold EQ T_LINEAR then Lambda = Lambda - StepL $
      else if TypeThreshold EQ T_MOM then begin
            ; In case of MOM and IterResi EQ 0, the MON threshold must be recalculated
            if IterResi EQ 0 then BEGIN
               for k=0,NbrTrans-1 do begin
                  if TabTrans[k] EQ T_Wave then TabMaxTrans[k] = maxwt(Residual, NbrScale=NbrScale)
                  if TabTrans[k] EQ T_PWave then TabMaxTrans[k] = maxpwt(Residual, NbrScale=NbrScale)
		  if TabTrans[k] EQ T_OWT then TabMaxTrans[k] = maxowt(Residual, NbrScale=NbrScale)
                  if TabTrans[k] EQ T_ALM then TabMaxTrans[k] = maxalm(Residual)
                  if TabTrans[k] EQ T_CMBLET then TabMaxTrans[k] = maxcmblet(Residual,lmax=lmax, wtp=wtp)
                  if TabTrans[k] EQ T_DIRAC then TabMaxTrans[k] = mrs_absmax(Residual)
                  if TabTrans[k] EQ T_Cur then TabMaxTrans[k] = maxcur(Residual,nbrscale=nbrscale)
		  if TabTrans[k] EQ T_DCT then TabMaxTrans[k] = maxdct(Residual,blocksize=DCTblocksize)
		  if TabTrans[k] EQ T_AT then TabMaxTrans[k] = maxawt(Residual, NbrScale=NbrScale)
                  print, 'Transform ',  tabNameTrans[ TabTrans[k] ], ' Max = ', TabMaxTrans[k]
              end ; for k
            END ; if IterResi      
            FirstThreshold = mrs_get_mom(NbrTrans, TabMaxTrans) / SigmaNoise
	    DeltaThreshold = FirstThreshold - LastThreshold
            StepL = DeltaThreshold / float(NIter-i)
	    if StepL LT SigmaNoise/2. then StepL = SigmaNoise/2.
            Lambda = MIN( [FirstThreshold, Lambda - StepL])
      end else if TypeThreshold EQ T_EXP or keyword_set(CMB) then begin   ; End if TypeThreshold EQ T_MOM
             Lambda =  LastThreshold  + DeltaThreshold  * (1.-erf(2.8*float(i) / float(niter)))
          end else if TypeThreshold EQ T_FIT then begin ; End if TypeThreshold EQ T_EXP 
               Lambda =threshold_fit(i)
        endif ; End if TypeThreshold EQ T_FIT 

   end ; END if i GT 0 
   
   if Lambda LT LastThreshold or i EQ niter-1 then Lambda = LastThreshold 
   num = 'ITER ' + STRCOMPRESS(string(i+1), /REMOVE_ALL) 
   if not keyword_set(CMB) then print,  num, ' Lambda = ', Lambda, ' MinResi = ', mrs_min(Residual), ' MaxResi = ', mrs_max(Residual), ' SigmaResi = ', mrs_sigma(Residual) $
   else print,  num, ' Lambda = ', LambdaCst, ' MinResi = ', mrs_min(Residual), ' MaxResi = ', mrs_max(Residual), ' SigmaResi = ', mrs_sigma(Residual)
   
   ; Iteration on the number of components
   for k=0,NbrTrans-1 do BEGIN
      ; print, '   ==> projection with ', tabNameTrans[TabTrans[k]]
      if GLESP EQ 1 then Ima = TabC[k] else Ima = TabC[*, k]
      if IterResi EQ 0 then Ima = mrs_add(Ima, Residual) ELSE Ima = Residual
       
      print,'########## SEUIL = ',Lambda
       
      if TabTrans[k] EQ T_Wave   then  if Lambda GT 0 then mca_proj_wave, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, SigmaNoise=SigmaNoise, soft=soft, FirstWTDetectScale=FirstWTDetectScale
      if TabTrans[k] EQ T_PWave  then  if Lambda GT 0 then mca_proj_wave, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, lmax=lmaxALM, SigmaNoise=SigmaNoise, soft=soft, /pyr
      if TabTrans[k] EQ T_OWT    then  if Lambda GT 0 then mca_proj_owt, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, SigmaNoise=SigmaNoise, soft=soft
      if TabTrans[k] EQ T_ALM    then  if Lambda GT 0 then mca_proj_alm, Ima, Lambda, S_ALM=S_ALM, Mad=Mad, MaxCoef=MaxCoef,  SigmaNoise=SigmaNoise , soft=soft, cmb=cmb, LambdaCst=LambdaCst
      if TabTrans[k] EQ T_CMBLET then  if Lambda GT 0 then mca_proj_cmblet, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef,  lmax=lmax, SigmaNoise=SigmaNoise, optwt=optwt, wtp=wtp
      if TabTrans[k] EQ T_DIRAC  then  if Lambda GT 0 then mca_proj_dirac, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, SigmaNoise=SigmaNoise 
      if TabTrans[k] EQ T_Cur    then  if Lambda GT 0 then mca_proj_cur, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, SigmaNoise=SigmaNoise , soft=soft
      if TabTrans[k] EQ T_DCT    then  if Lambda GT 0 then mca_proj_dct, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef,  SigmaNoise=SigmaNoise , soft=soft, blocksize=DCTblocksize, l2=l2
      if TabTrans[k] EQ T_AT    then  if Lambda GT 0 then mca_proj_at, Ima, Lambda, Mad=Mad, MaxCoef=MaxCoef, NbrScale=NbrScale, SigmaNoise=SigmaNoise, soft=soft, KillLastScale=KillLastScale

      TabMaxTrans[k] = MaxCoef
      if IterResi EQ 1 then BEGIN
        if GLESP EQ 1 then TabC[k] = mrs_add(TabC[k], Ima) else TabC[*,k] = TabC[*,k] + Ima
      END else  if GLESP EQ 1 then TabC[k] = Ima else TabC[*,k] = Ima
      
       
      ; Positivity of the solution
      if defined(Positivity) then begin
        ; print, 'POS', Positivity[k]
        if Positivity[k] EQ 1 then begin
        if GLESP EQ 1 then begin 
	     Ind = where (TabC[k].t_sky LT 0, c)
	     if c GT 0 then TabC[k].t_sky[Ind] = 0
        end else begin
	      Map = TabC[*, k]
	      Ind = where (Map LT 0, c)
	      if c GT 0 then Map[Ind] = 0
	      TabC[*, k] = Map
        end
      end 
      end; end positivity
      
      if keyword_set(TabSupport) then begin
        sup = TabSupport[*,k]
        ind = where(sup EQ 0, c)
	if c GT 0 then begin
	     if GLESP EQ 1 then TabC[k].t_sky[Ind] = 0 $
             else begin
 	      Map = TabC[*, k]
	      Map[Ind] = 0
	      TabC[*, k] = Map
             end
        end 
      end  ; end Tabsupport
      
      if keyword_set(Bounded) then begin
         if GLESP EQ 1 then begin 
	     Ind = where (TabC[k].t_sky GT BoundMax, c)
	     if c GT 0 then TabC[k].t_sky[Ind] = BoundMax
	     Ind = where (TabC[k].t_sky LT BoundMin, c)
	     if c GT 0 then TabC[k].t_sky[Ind] = BoundMin
         end else begin
 	      Map = TabC[*, k]
	      Ind = where (Map GT BoundMax, c)
	      if c GT 0 then Map[Ind] = BoundMax
	      Ind = where (Map LT BoundMin, c)
	      if c GT 0 then Map[Ind] = BoundMin
	      TabC[*, k] = Map
         end
      end ; End Bounded constraint

    ; tvs, TabC[*,0]

      
      ; Constraint inpainting ALM on the standard deviation; MUST be used only whith one transform
      ; print, '   Sigma IN = ', SigMaskOK, ' Sigma OUT = ', SigInp
      ; print,'ok0'
      if keyword_set(Mask) and keyword_set(CstSigma) and i GE  FirstIterConstr then begin
       print,'cst sigma',i
       if not glesp then begin
        lmax = 3.*npix2nside((size(tabc))(1)) 
	if lmax GT 3000 then lmax =3000
       end  else lmax =  tabc(0).nx /2 -2
       lmax = long(lmax)
       ; print,lmax
        
 	 ; mrs_mbtrans,reform(TabC[*,k]), wt, final, l_win=win;, ALM=S_ALM   ; NbrScale = NbrScale  
 	 
	; get_wp_meyer_filter, final,interval = cl_fine, lmax=lmax, nside=nside, hfilter=hfilter
        ; NbrScale = (size(final))(2)
        if GLESP EQ 0 then mrs_wptrans, reform(TabC[*,k]), wt,  NbrScale=NbrScale  $ ; , Wave=final,  l_win=win, lmax=lmax
        else mrs_wttrans, TabC[*,k], wt,  NbrScale=NbrScale
        

        
      ;  print,'ok1',nbrscale
        
        for j=0,NbrScale-1 do begin
	   Scale = mrs_wtget(wt, j)
	
            if GLESP EQ 1 then begin
              SigInp = sigma( Scale.t_sky[indMaskZero] )
	      SigMaskOK = sigma( Scale.t_sky[indMaskOK] )
 	      if siginp ne 0 and sigmaskok ne 0 then  Scale.t_sky[indMaskZero] = Scale.t_sky[indMaskZero] / SigInp * SigMaskOK
              siginp2 = sigma( Scale.t_sky[indMaskZero] )
              print, "Constraint: " , i,j,siginp,sigmaskok,SigMaskOK/ SigInp, SigMaskOK/ SigInp2 ,siginp2
           end else begin
              SigInp = sigma( Scale[indMaskZero] )
              MeanInp = mean( Scale[indMaskZero] )
	          SigMaskOK = sigma( Scale[indMaskOK] )
	          MeanMaskOK = mean( Scale[indMaskOK] )
	      ;  if siginp ne 0 and sigmaskok ne 0 then Scale[indMaskZero] = (Scale[indMaskZero] - MeanInp ) / SigInp * SigMaskOK +  MeanInp
	         if siginp ne 0 and sigmaskok ne 0 then Scale[indMaskZero] = Scale[indMaskZero] / SigInp * SigMaskOK 
	          moment=0
	          if keyword_set(moment) then begin
                Xin = Scale[indMaskOK]
                 K2 = sigma(Xin)
	            K1 = mean(Xin)
	            K3 = skewness(Xin)
 	            K4 = kurtosis(Xin)
 	            Xout = Scale[indMaskZero]
 	            print, j+1, ": K1 = ", K1, " K2 = ", K2, " K3 = ", K3, " K4 = ", K4
 	            Scale[indMaskZero] = MomentConst(Xout,K1,K2,K3,K4,20)
 	          end 
              if (SigInp GT 0) then print, "Constraint: ", i,j,siginp,sigmaskok,SigMaskOK/ SigInp else  print, "Constraint: ", i,j,siginp,sigmaskok
            end
	   mrs_wtput, wt, Scale , j
       end
       ;print,'ok2'
     mrs_wtrec,wt,recons
     
     if keyword_set(nomean) then recons = recons - mean(recons)
     ;help,tabc
     ;hs,tabc
     ;hs,recons
     ;tvs,mrs_diff(tabc(*,k),recons)
     ;tvs,tabc(*,k),title = strcompress('tabc '+string(i))
     ;tvs,recons(*,k),title = strcompress('recons '+string(i))
     tabc(*,k) = recons
     ;info,recons.t_sky
     ;info,tabc(*,k).t_sky
     end
 
    if keyword_set(film) then begin
         help,reform(tabc[*,0])
         if GLESP eq 1 then	mrs_write,strcompress('film_'+string(i)+'.fits',/remove_all),tabc[0] $
	  else mrs_write,strcompress('film_'+string(i)+'.fits',/remove_all),reform(tabc[*,0])
          endif
     
      ; print, '   Sigma IN = ', SigMaskOK, ' Sigma OUT = ', SigInp

      ; Compute the new residual
      Residual = mca_residual(Data, TabC, mask=mask)

      Tit='        ==> Rec ' + tabNameTrans[TabTrans[k]]
      if GLESP EQ 1 then mrs_info, TabC[k], mes=Tit $
      else mrs_info, TabC[*, k], mes=Tit
      
      END ; for k
    ; if i EQ 1 then mrs_tv, tabC[*,0], title=tabNameTrans[TabTrans[k]]
          
    if i GT 0 and i mod 1000 eq 0 then begin
       if GLESP EQ 1 then for k=0,NbrTrans-1 do  mrs_tv, tabC[k], title=tabNameTrans[TabTrans[k]], /healpix $
       else for k=0,NbrTrans-1 do  mrs_tv, tabC[*,k], title=tabNameTrans[TabTrans[k]]
    end 
    if Lambda EQ LastThreshold then i=niter

if keyword_set(disp) and i mod 5 EQ 0 then begin
  if GLESP EQ 1 then begin
     for k=0,NbrTrans-1 do  mrs_tv, tabC[k], title=tabNameTrans[TabTrans[k]]+string(i), /healpix, log=log 
     mrs_tv, residual, title='Residual '+string(i), /healpix
  end else  begin 
    for k=0,NbrTrans-1 do  mrs_tv, tabC[*,k], title=tabNameTrans[TabTrans[k]]+string(i), log=log
     mrs_tv, residual, title='Residual '+string(i)
  end
end

   END ; for i

   if GLESP EQ 0 and NbrTrans EQ 1 and Lambda EQ 0  then  TabC = TabC + residual

DONE:

end

;=====================================================

pro mrs_cmb_inpainting, Imag, Mask=Mask, ImagOut, niter=niter, to_glesp=to_glesp, FNImag=FNImag, FNMask=FNMask, FNOut=FNOut, FirstIterConstr=FirstIterConstr
COMMON C_PLANCK

DEF_ALM_FAST = 0
DEF_ALM_NITER = 0

if not keyword_set(niter) then niter=40

if keyword_set(cxx) then BEGIN 
   a = 0
END ELSE  BEGIN
  if keyword_set(FNMask) then Mask = mrs_read(FNMask)
  if keyword_set(FNImag) then Imag = mrs_read(FNImag)

  if not keyword_set(Mask)  then BEGIN
      Mask = Imag
      if  type_code(Imag) EQ 8 then begin
          Mask.t_sky[*] = 0
          ind = where(Imag.t_sky NE 0, c)
         if c GT 0 then Mask.t_sky[ind] = 1
      end else begin
         ind = where(Imag NE 0, c)
         if c GT 0 then Mask[ind] = 1
     end
    END
END

 ; We remove the mean value from the available input data
 if  type_code(Imag) NE 8 then BEGIN
  ind = where(mask NE 1, c)
  if c GT 0 then mask[ind] = 0
  indMask = where( mask NE 0, cMask)
  Map = Imag
  MeanVal = mean(map[indMask])
  if cMask NE 0 then map[indMask] = map[indMask] -  MeanVal
  map = map * mask
END

if keyword_set(to_glesp) then begin
  GImag = healpix2glesp(Map, /alm)
  GMask = healpix2glesp(mask, optnx=GImag.nx, optnp=GImag.np)
  mrs_mca, GImag, ImagOut, mask=GMask,  niter=niter , residual=residual, selecttrans=[4],/expo,/cstsigma;/fit, , FirstIterConstr=FirstIterConstr    ;,/cstsigma;,/disp;/expo;,/disp;,/cstsigma,/disp
end else mrs_mca, Map, ImagOut, mask=Mask,  niter=niter, residual=residual, selecttrans=[4],/expo ,/cstsigma , /nomean, FirstIterConstr=FirstIterConstr

if  type_code(Imag) NE 8 then BEGIN
if cMask NE 0 then ImagOut[indMask] = ImagOut[indMask] +  MeanVal
END

if keyword_set(FNOut) then mrs_write, FNOut, ImagOut
end

;=====================================================

pro run_may2012,  tabC, filter=filter
DIR= '/mnt/iscsi2-cosmo/sparseastro/JB/Planck/FFP4/'

TabFN = ['CMB_FFP4_Band0_sameres_inverted_inversion.fits', $   ; 33 min
'CMB_FFP4_Band1_sameres_inverted_inversion.fits', $  ; 24 min
'CMB_FFP4_Band2_sameres_inverted_inversion.fits', $  ; 14 min
'CMB_FFP4_Band3_sameres_inverted_inversion.fits', $  ; 10 min
'CMB_FFP4_Band4_sameres_inverted_inversion.fits']  ; 5 min
vs = size(TabFN)
Nf = vs[1]

; i = rim('CMB_FFP4_Band4_sameres_inverted_inversion.fits')*1e6
m = rims('../mask74pc.fits')
PSMask = rims('MaskPS_HFI_LFI_3sigma.fits')
COMask = 1. - rims('mask_co_v3.fits')
MGal = rims('/dsm/cosmo02/planck/starck/PL_DATA/DX8/NG/MASK/mask_gal_DX7_857GHz_5pc.fits')
Mask = PSMask * COMask * MGal
m512 = mrs_resize(Mask, nside=512)
ind = where(m512 ne 1, c)
if c GT 0 then m512[ind] = 0

m1024 = mrs_resize(Mask, nside=1024)
ind = where(m1024 ne 1, c)
if c GT 0 then m1024[ind] = 0

for f=1,Nf-1 do begin
   i = rim(DIR + TabFN[f]) * 1e6
   nside = gnside(i)
   if nside EQ 512 then m = m512 $
   else if nside EQ 1024 then m = m1024 $
   else m = Mask
   
   N = N_ELEMENTS(i)
   TabSupport = fltarr(N,2)
   TabSupport[*,0] = 1. - m
   TabSupport[*,1] = 1
   FirstGuess = TabSupport * 0
   NbrScale=5

   ; if not keyword_set(filter) then 
   killcmb, i, filter, nsigma=5, nbrscale=nbrscale, fdr=fdr, Use_FdrAll=Use_FdrAll
  ;tvs, filter, min=-500, max=500, tit='Filter'
  ; tvs, i-filter, min=-500, max=500, tit='Data - Filter'
   FirstGuess[*,0] = filter
  ; save, filename='xx2.xdr', filter
   mrs_mca, i, tabC, SelectTrans=[1,4], niter=40, /expo, LastThreshold=0, SigmaNoise=1, NbrScale=NbrScale,  tabsupport=tabsupport, FirstGuess= FirstGuess
   ; save, filename='xx_wmap-nilc.xdr', i, filter, tabC
   FN = 'xx_mca_' + strc(f+1) + '.xdr'
  save, filename=FN, i, filter, tabC
  FN = 'MCA_' +TabFN[f]  
  mrs_write, FN, tabC[*,1]
end

end

 ;=====================================================

pro killcmb2, ima, filter, nsigma=nsigma, nbrscale=nbrscale, fdr=fdr, Use_FdrAll=Use_FdrAll
if not keyword_set(nsigma) then nsigma=5
print,'NSigma',nsigma
if not keyword_set(nbrscale) then nbrscale=5

mrs_wtfilter, ima, filter, /kill, /mad, nsigma=nsigma, nbrscale=nbrscale, fdr=fdr, Use_FdrAll=Use_FdrAll, niter=5, /pos
z1 = filter-ima
mrs_wtfilter, z1, nfilter, /kill, /mad, nsigma=nsigma, nbrscale=nbrscale, fdr=fdr, Use_FdrAll=Use_FdrAll, niter=5, /pos
filter = filter - nfilter
tvs, filter, min=-500, max=500
end
;=====================================================


pro run_test_jan09_a
l = rims('lensed_map.fits')
m = rims('../mask_apcch2.fits')
d = m * l
; mrs_cmb_inpainting, d, i1, Mask=m, niter=40, FNOut='inp_idlv1.fits'
; mrs_cmb_inpainting, d, i2, Mask=m, niter=40, FNOut='inp_idlv2.fits', FirstIterConstr=20
i1 = rims('inp_idlv1.fits')
i2 = rims('inp_idlv2.fits')

pl = mrs_powspec(l)
pd = mrs_powspec(d)
pi1 = mrs_powspec(i1)
pi2 = mrs_powspec(i2)
save, filename='res_jan09.xdr', i1, i2, m, l, pl, pd, pi1, pi2
end
;=====================================================

pro run_test_jan09
l = rims('lensed_map_Lobb_hfi6chan-25001.fits')
m = rims('../mask_apcch2.fits')
d = m * l
mrs_cmb_inpainting, d, i1, Mask=m, niter=40, FNOut='inp_idlv1_lensed_map_Lobb_hfi6chan-25001.fits' 
pl = mrs_powspec(l)
pd = mrs_powspec(d)
pi = mrs_powspec(i1)
save, filename='res_jan09b.xdr',  pl, pd, pi
end

;=====================================================

pro old_mrs_alm_inpainting, Imag, Mask=Mask, ImagOut, niter=niter, to_glesp=to_glesp
COMMON C_PLANCK

DEF_ALM_FAST = 0
DEF_ALM_NITER = 0

if not keyword_set(niter) then niter=100
if not keyword_set(Mask) then BEGIN
  Mask = Imag
  if  type_code(Imag) EQ 8 then begin
     Mask.t_sky[*] = 0
     ind = where(Imag.t_sky NE 0, c)
     if c GT 0 then Mask.t_sky[ind] = 1
  end else begin
     ind = where(Imag NE 0, c)
     if c GT 0 then Mask[ind] = 1
  end
END

if keyword_set(to_glesp) then begin
  GImag = healpix2glesp(Imag, /alm)
  GMask = healpix2glesp(mask, optnx=GImag.nx, optnp=GImag.np)
  mrs_mca, GImag, ImagOut, mask=GMask,  niter=niter , residual=residual, selecttrans=[4,1], /expo ;/fit;,/cstsigma;,/disp;/expo;,/disp;,/cstsigma,/disp
end else mrs_mca, Imag, ImagOut, mask=Mask,  niter=niter, residual=residual, selecttrans=[4],/expo 

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro test_feb07
nside=64
cmb = getcmb(Cl=Cl, nside=nside, Fwhm=Fwhm, SigmaNoise=SigmaNoise, Cmb=Cmb)
mask = mrs_read('$MRS/WMAP3/mask_kp00_n0256.fits')
mask = mrs_resize(mask, nside=nside)
ind = where(mask NE 1, c)
if c GT 0 then mask[ind] = 0
c1 = cmb * mask
mrs_cmb_inpainting, C1, Mask=Mask, Inp, niter=30
 
 
gcmb = healpix2glesp(cmb, /alm)
gmask = healpix2glesp(mask, optnx=gcmb.nx, optnp=gcmb.np)

p = mrs_powspec(cmb, /nonorm)
h = healpix_window(nside)
plotcl, p
oplotcl, info

gdata = gcmb
gdata.t_sky = gdata.t_sky * gmask.t_sky

pglesp = mrs_powspec(gcmb, /nonorm)
pdata = mrs_powspec(gdata, /nonorm)

nbrPix = N_ELEMENTS(gdata.t_sky)
Ind = where(gmask.t_sky EQ 1, NbrPixOK)
Bias = float(nbrPix) / float(NbrPixOK)
pbiais = pdata*Bias
oplotcl, pdata
oplotcl, pdata*Bias, line=3


mrs_mca, gdata, gres, selecttrans=[4], niter=30, mask=gmask, t=t, /expo, /CstSigma, /disp
  
pres = mrs_powspec(gres, /nonorm)

end


;=====================================================


pro old_cmb_inp,  residual=residual, data, Rec,  niter=niter,  disp=disp,  Mask=Mask 

CstSigma=1
SelectTrans=[4]
expo=1
IterResi=1
if not keyword_set(disp) then disp=0

nside=0
if type_code(Data) EQ 8 then  GLESP=1 else begin
   GLESP=0
   npix = (size(Data))[1]
   nside = npix2nside(npix)
end

if not keyword_set(Mask) then begin
     if GLESP EQ 1 then begin
 	indMaskZero =  where (data.t_sky EQ 0)
	mask=data
	mask.t_sky[*]=1
	mask.t_sky[indMaskZero]=0
      end else begin
 	indMaskZero =  where (Mask EQ 0)
	mask=data
	mask[*]=1
	mask[indMaskZero]=0
      end
end
     if GLESP EQ 1 then begin
        indMaskOK = where (Mask.t_sky NE 0)
	SigMaskOK = sigma( data.t_sky[indMaskOK] )
	indMaskZero =  where (Mask.t_sky EQ 0)
      end else begin
        indMaskOK = where (Mask NE 0)
	SigMaskOK = sigma( data[indMaskOK] )
	indMaskZero =  where (Mask EQ 0)
      end
      
      
 Rec = Data
 mrs_set, Rec, 0
 FirstThreshold = maxalm(Data)
 LastThreshold = 0.
 DeltaThreshold = FirstThreshold - LastThreshold
  SigmaNoise=1.
for i=0, niter-1 do begin
   Lambda =  LastThreshold  + DeltaThreshold  * (1.-erf(2.8*float(i) / float(niter)))
   if i NE 0 then  mrs_almrec, ALMSol, Rec
   Resi = mrs_mult(mrs_diff(Data,Rec), Mask)
   mrs_almtrans, Resi, S_ALM, /norm, lmax=lmax
   if i NE 0 then  Scale = S_ALM.alm[1:*,*] + ALMSol.alm[1:*,*] $
   else begin
      ALMSol = S_ALM
      ALMSol.alm[*] = 0
      Scale = S_ALM.alm[1:*,*]
   end
   TScale = Scale
   NormVal = 1. / sqrt(2.)
   ThresholdLevel = Lambda*SigmaNoise*NormVal 
   TScale = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
   ALMSol.alm[1:*,*] = TScale
   mrs_almrec, ALMSol, Rec 

    if GLESP EQ 1 then begin
         SigInp = sigma( Rec.t_sky[indMaskZero] )
        end else begin
        SigInp = sigma( Rec[indMaskZero] )
        end
    print, i+1, ' Lambda = ' , Lambda, '  R(Sigma OUT/IN) % = ',  SigInp / SigMaskOK * 100.  
      
    if keyword_set(disp) and (i+1) mod disp eq 0 then tvs, rec  
    
end
 

end


;==========================================
pro mca_test_gauss_alm, ca, cw, d, r, mad=mad, mom=mom
n = randomn(seed, 64L^2*12)
mrs_almtrans, n, a
a.alm(*)=0
a.alm[50,0]=10.
mrs_almrec, a, ca

f = h2f(n)
f(*)=0
f[32,32,8] = mygauss(32,32,1.5) * 10.
cw = f2h(f) 
gca = healpix2glesp(ca, /alm)
gcw = healpix2glesp(cw, /alm)

gd = mrs_add(gca, gcw)
d = ca + cw
Colt=15

mrs_tv, d, title='Data', gif='test_gauss_alm.gif', colt=Colt
; mrs_mca, d, tabC, SelectTrans=[1,4], niter=50, /mom, tabNameTrans=tabNameTrans, LastThreshold=LastThreshold, SigmaNoise=1, NbrScale=NbrScale, lmax=lmax, Mask=Mask, /disp 

mrs_mca, d, tabC, SelectTrans=[1,4], niter=20, /mad, tabNameTrans=tabNameTrans, LastThreshold=LastThreshold, SigmaNoise=1, NbrScale=NbrScale, lmax=lmax, Mask=Mask, /disp, tabsupport=tabsupport



mrs_tv, r, title='MCA Residual', gif='test_gauss_alm_mca_resi.gif', colt=Colt
; mrs_tv, tabC[*,0], title='MCA Rec WT', gif='test_gauss_alm_mca_wt.gif', colt=Colt
; mrs_tv, tabC[*,1], title='MCA Rec ALM', gif='test_gauss_alm_mca_alm.gif', colt=Colt

 ; mrs_mca, s1, TabC2, selectTrans=[3,1],  /mom, niter=100,                                                                                                  

;tvs, tabC[*,0]
;tvs, tabC[*,1]

end

 

pro mca_test_cmbmask
; tabNameTrans = ['X', 'Wavelet',  'Pyramidal Wavelet', 'CMBLET', 'ALM', 'DIRAC']

; restore, /verb, 'd64.xdr'
c64 = rims('real_0_256.fits')
m64 = rims('mask256.fits')
d64 = c64 * m64
invm = 1. - m64

; ALM MCA
GlespTest = 1 

if GlespTest EQ 0 then mrs_mca, d64, tc, selecttrans=[3], niter=50, mask=m64, t=t, /wtp, /expo

; LOOK at the phase and module
InfoPhase = 0
if INfoPhase EQ 1 then begin
  mrs_almtrans, d64, a, /psp, /tab
  a1=a
  a1.alm[*,*,0]=1.
  mrs_almrec, a1, recphase  
  a2=a
  a2.alm[*,*,1]=0.
  mrs_almrec, a2, recpowspec
  mrs_tv, d64, title='CMB + Mask'
  mrs_tv, recphase, title='Phase only'
  mrs_tv, recpowspec, title='Power spectrum only'
end

; GLESP
if GlespTest EQ 1 then begin
  gc64 = healpix2glesp(c64, /alm)
  gm64 = healpix2glesp(m64, optnx=gc64.nx, optnp=gc64.np)
  gd64 = gc64
  gd64.t_sky = gd64.t_sky * gm64.t_sky
  mrs_mca, gd64, gtc, selecttrans=[3], niter=50, mask=gm64, t=t, /expo, /wtp, /disp
  mrs_almtrans, gc64, ap
  mrs_almtrans, gd64, apd
  mrs_almtrans, gtc, api
  p = mrs_alm2spec(ap)
  pd = mrs_alm2spec(apd)
  pi = mrs_alm2spec(api)
  mrs_pswtfil, p, fp, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul
  mrs_pswtfil, pd, fpd, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul 
  mrs_pswtfil, pi, fpi, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul 
  save, filename='ResMask.xdr', d64, m64, c64, gtc, ap, apd, api, p , pd, pi, fp, fpd, fpi
end

end

pro wres, ct=ct
restore, /verb, 'd64.xdr'
inp = readfits('cmb64_inp_alm_200iter.fits')
d64 = c64 * m64
mrs_tv, d64, png='fig_cmb64_mask.png', title='CMB nside=64 + mask', colt=ct
mrs_tv, inp, png='fig_cmb64_inp_alm_200iter.png', title='ALM inpainting CMB nside=64, 200 iter', colt=ct
mrs_tv, c64, png='fig_cmb64.png', title='CMB nside=64' , colt=ct

mrs_almtrans, c64, ap
mrs_almtrans, d64, apd
mrs_almtrans, inp, api
p = mrs_alm2spec(ap)
pd = mrs_alm2spec(apd)
pi = mrs_alm2spec(api)
 
mrs_pswtfil, p, fp, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul
mrs_pswtfil, pd, fpd, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul 
mrs_pswtfil, pi, fpi, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul 
mrs_pswtfil, px, fpx, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=5, NbrIter=10, mask=mask, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, /L1Regul 

plotcl, p
oplotcl, pd
oplotcl, pi

plotcl, fp
oplotcl, fpd
oplotcl, fpi

gc64 = healpix2glesp(c64, /alm)
nx=gc64.nx
np=gc64.np
gm64 = healpix2glesp(m64, optnx=nx, optnp=np)
gd64 = gc64
gd64.t_sky = gc64.t_sky * gm64.t_sky 
end       


pro mca_test_gauss_cmb, ca, cw, d, r, tabc, mad=mad, mom=mom
 restore, /verb, 'd64.xdr'
ca = c64
NbrScale=5
f = h2f(ca)
f(*)=0
Sig=1.
nside = 64
f[*,*,10] = mygauss(nside,nside,Sig) * 500.
f[*,*,9] = mygauss(nside,nside,Sig) * 400.
f[*,*,8] = mygauss(nside,nside,Sig) * 300.
f[*,*,7] = mygauss(nside,nside,Sig) * 200.
f[*,*,6] = mygauss(nside,nside,Sig) * 100.
f[*,*,5] = mygauss(nside,nside,Sig) * 80.
f[*,*,4] = mygauss(nside,nside,Sig) * 50
f[*,*,3] = mygauss(nside,nside,Sig) * 40
f[*,*,2] = mygauss(nside,nside,Sig) * 30
f[*,*,1] = mygauss(nside,nside,Sig) * 20
f[*,*,0] = mygauss(nside,nside,Sig) * 10
f = f * 10
cw = f2h(f) 

gca = healpix2glesp(ca, /alm)
gcw = healpix2glesp(cw, /alm)
gd = mrs_add(gca, gcw)
d = ca + cw
mrs_tv, gd, /log, /healpi

Colt=15
mrs_tv, d, title='Data', gif='test_gauss_cmb.gif', colt=Colt
mrs_mca, residual=r, gd, tabC, SelectTrans=[4,1], niter=50, /expo, /mad, cmb=2,  LastThreshold=LastThreshold, SigmaNoise=1, NbrScale=NbrScale
mrs_tv, r, title='MCA Residual', gif='test_gauss_cmb_mca_resi.gif', colt=Colt
mrs_tv, tabC[*,1], title='MCA Rec WT', gif='test_gauss_cmb_mca_wt.gif', colt=Colt
mrs_tv, tabC[*,0], title='MCA Rec CMB', gif='test_gauss_cmb_mca_alm.gif', colt=Colt

mrs_mca, residual=r, d, tabC, SelectTrans=[4,3], niter=50, /expo, /mad, cmb=2,  LastThreshold=LastThreshold, SigmaNoise=1, NbrScale=5

 ; mrs_mca, s1, TabC2, selectTrans=[3,1],  /mom, niter=100,                                                                                                  

;tvs, tabC[*,0]
;tvs, tabC[*,1]

end

pro tt1, gd, tabc, tabc1, r
sig=160.
mrs_mca, residual=r, gd, tabC, SelectTrans=[4,1], niter=10, /mad, SigmaNoise=sig, NbrScale=5, cmb=1, pos=[0,1]
mrs_mca, residual=r, gd, tabC, SelectTrans=[1,4], niter=10, /mad, SigmaNoise=sig, NbrScale=5, cmb=1, /soft, /disp
mrs_mca, residual=r, gd, tabC, SelectTrans=[1,4], niter=10, /mad, SigmaNoise=sig, NbrScale=5, cmb=1, /disp
mrs_mca, residual=r, gd, tabC, SelectTrans=[1,4], niter=50, /expo, SigmaNoise=sig, NbrScale=5,  /disp

mrs_mca, gd, SelectTrans=[4], mask=mask2, nbrscale=6, /expo, niter=100,  res1

end
 
