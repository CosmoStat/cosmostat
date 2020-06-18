
;; WTNscale : number of scales for the wavelet decomposition, used in
;; the sparisty connstraint
;; lmax : number of l used to recover cl
;; Nscale : number of scales used in the needlet decomposition used to
;; estimate the angular power spectrum


;;============================================== ntousi ===================================================

function cl_recons_test, Map, SigmaNoise, angSpec, filters=bkj, Mask=Mask, Nscale=Nscale, Firstl=Firstl, lmax=lmax, tousil=tousil, WTNscale=WTNscale, Niter=Niter, SDCT=SCDT, SWT=SWT, OptWT=OptWT, BestBasis=BestBasis, Pos=Pos, TabInterval=TabInterval, Cl=Cl, Proj_Resi=Proj_Resi, True_Resi=True_Resi, part_nMSE=part_nMSE, Mcl=Mcl, uresult=uresult

COMMON C_PLANCK

npix = N_elements(Map)
nside = npix2nside(npix)

;;---------------------------------------- Set default parameters ----------------------------------------
if not keyword_set(Niter) then Niter = fix(10*exp(float(Nscale)/10)) ; Niter = 500 ;

if not keyword_set(lmax) then begin
   lmax = long( nside )  * 3l
   if lmax GT P_Lmax then  lmax = P_Lmax  
endif
;; if we don't have many l, we don't want to use the DCT, and we consider more wavelet scales
if lmax lt 1000l then OWT = 1

if not keyword_set(WTNscale) then begin
   if not keyword_set(OWT) then  WTNscale = fix(alog(float(lmax)))  - 1 $
   else  WTNscale = fix(alog(float(lmax)))  + 1
endif

if not keyword_set(Firstl) then Firstl=0
if not keyword_set(Nscale) then Nscale = round(4*alog(lmax))
if not keyword_set(tousil) then tousil = min([max([100,2*Firstl]), lmax/2.])

if arg_present(BestBasis) then BestBasis = 1
if not arg_present(gen2) then gen2 = 1
if not arg_present(SWT) then SWT = 1
if not arg_present(HWT) then HWT = 1
if not arg_present(SDCT) then SDCT = 0
if not arg_present(HDCT) then HDCT = 0

;;------------------------------------------ Print Toucan options ------------------------------------------------------

print,  "SWT = ", SWT
print, "HWT = ", HWT

Bmat = compute_Bmat(filters=bkj)

;;---------------- Normalizing the input --------------------
sigmaSpec = sigma(angSpec)
angSpec = angSpec/sigmaSpec
if keyword_set(Firstl) then Bmat = Bmat[Firstl:*,*]
Phi = Bmat ## transpose(Bmat)
tau = maxeigval(Phi) ; step size for the gradient
sigl = lmax+1 - Firstl ; length of the Cl to recover
;;------------------------------------------------------------

;;------------------ Choose threshold level ------------------
Dirac = dblarr(sigl)
dirac[sigl/2] = 1
w  = uwt1d(dirac, Nscale=WTNscale, Gen2=Gen2, info=info, OptWT=OptWT)
TabPSWTNorm = dblarr(WTNscale)
for j=0, WTNscale-1 do TabPSWTNorm[j] = sqrt( total(w[*,j]^2) )

LevelDetect = fltarr(WTNscale)
Nsigma = 7 ;total(sigmaNoise)/N_elements(sigmaNoise) ; better way to determine the threhsold at each scale ?
for j = 0, WTNscale-1 do BEGIN
    if j NE 0 then LevelDetect(j) = NSigma * TabPSWTNorm(j)  $
    else LevelDetect(j) = (NSigma+1) * TabPSWTNorm(j)
 endfor
print, 'threshold chosen'
;;------------------------------------------------------------

;;------------ Set the masks for the HWT and HDCT ----------------
if not keyword_set(Mcl) then begin
   print, 'Compute estimated Power Spectrum with Master to estimate HWT and/or HDCT masks...'
   if not keyword_set(mask) then Mcl = mrs_powspec(map, lmax=lmax) else Mcl = mrs_master_powspec(map,mask,lmax=lmax)
endif
Signal =double(Mcl[Firstl:*])
SigmaSignal = sigma(Signal)
Signal =  Signal  / SigmaSignal
if N_elements(sigmaNoise) gt 1 then pn = mrs_powspec(randomn(seed, npix)*sigmaNoise) ; *****
powspecInit = mrs_tousi(Mcl[0:2*tousil], NoisePs=pn) ; ***** 

l = findgen(sigl) + firstl & lfactor = l * (l+1) / 2. / !dpi
if firstl EQ 0 then l[0]=1
if firstl EQ 0 then lfactor[0]=1
LastThreshold=0
if  keyword_set(SDCT)  then begin 
  if not keyword_set(BestBasis) then Coeff = transdico(Signal* lfactor)  $
  else Coeff = bb_dct(Signal* lfactor, LDCT=LDCT, UseFFT=UseFFT,TabInterval=TabInterval)
  FirstThreshold = max(Coeff[1:*])
endif else FirstThreshold=  max( sqrt(Signal) * lfactor)
DeltaThreshold = FirstThreshold - LastThreshold
DCTDeltaThreshold =  FirstThreshold

indl = 2. * l + 1.
if firstl EQ 0 then indl[0]=1
NormL = (1.+1/indl)
Fa =  -0.25839012       
Fb = 1.1820104
x = findgen(501) / 500.
yy = Fa*x + Fb
ind = where(yy lt 1, c)
if c GT 0 then yy[ind] = 1.
NormL[0:min([sigl-1,500])] = NormL[0:Min([sigl-1,500])] * yy

if keyword_set (HWT) then begin
   WTMask = fltarr(sigl, WTNscale) 
   WTMask[*, WTNscale-1] = 1.
   
   MR_Trans = uwt1d(Signal, gen2=gen2, Nscale=WTNscale, info=info, OptWT=OptWT) ;***
   for j = 0, WTNscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
      Scale     = MR_Trans[*,j]   
      
      ScaleMask = WTMask[*,j]
;;      Level = NormL * LevelDetect[j]
;;       Level = max([LevelDetect,0.1])
      Level = LevelDetect[j]
      
      ind       = where(abs(Scale) GT Level , count) ; Sig should be 1 in the GuassData space 
      if count GT 0 then ScaleMask[ind] = 1  
      WTMask[*,j] = ScaleMask  
   END 
   print, 'HWT mask computed'
endif


if  keyword_set(HDCT)  then begin
  GaussData = mrs_variance_stabilization(Signal, mu=mu, psi1=psi1,firstl= firstl)
  
  if not keyword_set(BestBasis) then Coeff = transdico(GaussData)  $
  else Coeff = bb_dct(GaussData, LDCT=LDCT, UseFFT=UseFFT,TabInterval=TabInterval)
  
  DCT_HThresh = NSigma
  ind = where( ABS(Coeff) GT  DCT_HThresh, c)
  DCTMask = Coeff * 0
  if c GT 0 then DCTMask[ind] = 1
  print, 'HDCT mask computed'
endif   
;;---------------------------------------------------------------

;;------------------ Initialization --------------------------
Resi =  (transpose(Bmat) ## angSpec)/(sqrt(tau)) ; Residual for the gradient step
Resi = Resi[*]
SigmaResi_old = sigma(Resi)
Result = fltarr(sigl) ; Initial guess for the power spectrum
step = 1
StopDecrease = 0
proj_resi = fltarr(Niter)
true_resi = fltarr(Niter)
part_nMSE = fltarr(Niter)
init = min([tousil, lmax-100])

;;-------------------------------------------------------------

for iter=1, Niter do begin

   if StopDecrease eq 0 then Lambda = ( 1.-erf(2.5*float(iter) / min([Niter,100])) )

   Resi_old = Resi
   Result_old = Result
   step_old = step
   step = (1. + sqrt(1.+4.*step^2))/2.
   Z = Result + (step_old - 1.)*(Result-Result_old)/step
   Result_old = Result
   Resi_old = Resi
   Resi =  (transpose(Bmat) ## (angSpec - Bmat ## Z) )/tau
   Resi = Resi[*]  
   
   ;;========= WT - HARD Thresholding

   if keyword_set (HWT)  then begin
      MR_Trans = uwt1d(Resi, gen2=gen2, Nscale=WTNscale, info=info, OptWT=OptWT)
      if iter EQ 1 then MaxHWT = max(ABS(MR_Trans[*,0: WTNscale-2]))

      for j = 0, WTNscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
         Scale     = MR_Trans[*,j]   
         ScaleMask = WTMask[*,j]
         MR_Trans[*,j] = Scale * ScaleMask
      END 
      Result =  Result  + iuwt1d(MR_Trans, gen2=gen2, info=info)
   endif
   
 
   ;;========= WT - SOFT Thresholding

   if keyword_set(SWT) then begin  
      MR_Trans = star1d(Result*lfactor, gen2=gen2, Nscale=WTNscale)
      if iter EQ 1 then MaxWT = max(ABS(MR_Trans[*,0: WTNscale-2]))

      for j = 0, WTNscale-2 do begin 
         SoftThreshold =   MaxWT  *  Lambda *  TabPSWTNorm[j] 
         Scale  = reform(MR_Trans[*,j]) 
         softthreshold, Scale, SoftThreshold
         MR_Trans[*,j] = Scale 
      endfor
      Result = istar1d(MR_Trans, gen2=gen2)/lfactor
   endif
   
   ;;========= DCT - HARD  Thresholding
   ;;   window, 1 & plot, Resi[2:*], title='resi'
   ;;   wait, 1
   if  keyword_set(HDCT)  then begin 
      step_old = step
      step = (1. + sqrt(1.+4.*step^2))/2.
      Z = Result + (step_old - 1.)*(Result-Result_old)/step
      Result_old = Result
      Resi =  (transpose(Bmat) ## (angSpec - Bmat ## Z) )/tau
      Resi = Resi[*]  
      
      PS_Estim = Resi
      
      if not keyword_set(BestBasis) then Coeff = transdico(PS_Estim)  $
      else Coeff = bb_dct(PS_Estim, LDCT=LDCT, UseFFT=UseFFT,TabInterval=TabInterval)
      
      CoeffT = Coeff * DCTMask
      
      if not keyword_set(BestBasis) then  PS_Estim    = (transdico(CoeffT, /inverse))   $ 
      else   PS_Estim = bb_dct(CoeffT, LDCT=LDCT, /inverse, UseFFT=UseFFT)
      
      Result = Result + PS_Estim
   endif   
   
   ;;========= DCT - SOFT Thresholding
    
   if  keyword_set(SDCT)  then begin 
      
      if not keyword_set(BestBasis) then Coeff = transdico(Result*lfactor)  $
      else  Coeff = bb_dct(Result*lfactor, LDCT=LDCT, UseFFT=UseFFT)
      
      ;Level =   DeltaThreshold  * Lambda
      Level = SigmaNoise * lambda
      ; bb_softthreshold, Coeff, Level,  LDCT=LDCT
      softthreshold, Coeff, Level
      
      if not keyword_set(BestBasis) then  Result  = (transdico(Coeff, /inverse))/lfactor   $ 
      else   Result = bb_dct(Coeff, LDCT=LDCT, /inverse, UseFFT=UseFFT)/lfactor
      
   endif   
  
   ;;========= Positivity constrain
   if keyword_set(Pos) then begin
      ind = where(Result LT 0, c)
      if c GT 0 then Result[ind] = 0 
   endif
   
   sigmaResi = sigma(Resi)
   proj_resi(iter-1) = sigmaResi
;;    true_resi(iter-1) = sigma(result[50-firstl :lmax-firstl]*sigmaSpec - cl[50 :lmax])
   true_resi(iter-1) = sigma(result[0 :lmax-firstl]*sigmaSpec - cl[firstl :lmax])
;;    part_nMSE(iter-1) = total(((cl[50:2*nside]-result[50-firstl :2*nside-firstl]*sigmaSpec)/cl[50:2*nside])^2)
   part_nMSE(iter-1) = total(((cl[init :lmax]-result[init-Firstl :lmax-firstl]*sigmaSpec)/cl[init :lmax])^2)
  
;;    print, Iter, '/', Niter, ':  Lambda = ', Lambda, ' | SigmaResi = ', sigma(Resi), ' | SigmaResult = ', sigma(Result)    
;;    print, Iter, ':  Linf(ResiDiff)= ', max(abs(Resi-Resi_old))

   if sigmaResi_old lt sigmaResi-(1.e-8) then begin
      print, 'Increasing residual -> iteration interrupted ==> try fewer iterations'
;;       interupted = 1
;;       goto, DONE
   endif
   sigmaResi_old = sigmaResi

endfor

DONE:
result = [fltarr(Firstl),result*sigmaSpec]
uresult = result
weight = 1+(exp(5.)-1)*(tousil-findgen(tousil+1))/tousil; **** change
weight = alog(weight)
weight2 = (1.-weight/5.)
result[0:tousil] = (weight2*result[0:tousil] + weight*powspecInit[0:tousil])/(weight+weight2) 

return, Result

end
