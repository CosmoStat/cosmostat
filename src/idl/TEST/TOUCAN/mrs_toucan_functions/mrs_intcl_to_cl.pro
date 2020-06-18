;+
; NAME:
;        mrs_intcl_to_cl
;
; PURPOSE:
;  Computes the CMB power spectrum Cl from the Cl integrated over
;  wavelet filters.
;
; CALLING:
;     toucanCl = mrs_intcl_to_cl( Map, intCl, bkj, Mask=Mask,
;     Firstl=Firstl, Bmat=Bmat, Niter=Niter, initl=initl,
;     WTNscale=WTNscale, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos,
;     Mcl=Mcl, Cl=Cl, Proj_Resi=Proj_Resi, True_Resi=True_Resi,
;     part_nMSE=part_nMSE)
;
; INPUTS:
;     Map        -- Healpix map = Input CMB image to estimate the CMB
;                   power spectrum from.
;     intCl      -- fltarr(Nscale) = Power Spectrum Cl integrated over the wavelet
;                   filters in bkj.
;     bkj        -- fltarr(lmax+1, Nscale) = wavelet filters in the
;                   multipole domain used to compute the integrated
;                   Cl.
;    
; OUTPUTS:
;     toucanCl   -- fltarr(lmax+1)  = Estimated Cl.
;
; INPUT KEYWORD: 
;     Mask        : Healpix map = mask applied to the input Map.
;     Firstl      : int = initial multipole value to be reconstructed
;                   in Cl.
;     Bmat        : IDL dblarr(lmax+1,Nscale) = measurement matrix
;                   constructed from the wavelet filters.
;     Niter       : int = number of iterations for the reconstruction.
;     initl       : int = multipole value to compute the partial nMSE
;                   from.
;     WTNscale    : int = number of wavelet scales used to impose
;                   sparsity during reconstruction.
;     SWT         : bool = if set, soft-thresholding of the wavelet
;                   coefficients of the power spectrum is performed
;                   during reconstruction.
;     HWT         : bool = if set, hard-thresholding of the wavelet
;                   coefficients of the power spectrum is performed
;                   during reconstruction.
;     OptWT       : options to be used in the wavelet decomposition of
;                   the power spectrum during reconstruction.
;     Pos         : bool = if set, a positivity constraint is applied
;                   to the power spectrum during reconstruction.
;     Mcl         : fltarr(lmax+1) = Cl estimated with Master
;     Cl          : fltarr(lmax+1) = true power Spectrum Cl.
;
; OUTPUT KEYWORD: 
;     Proj_Resi   : fltarr(Niter) = variance of the residual using
;                   projected onto Bmat.
;     True_Resi   : fltarr(Niter) = variance of the true residual if
;                   Cl is present in input.
;     part_nMSE   : fltarr(Niter) = true partial normalized mean
;                   square error if Cl is present in input.
;
; EXAMPLE:
;      Estimate the CMB power spectrum from the integrated power
;      spectrum intcl:
;      toucanCl = mrs_intcl_to_cl(Map, intCl, bkj,
;      Firstl=2, Niter=200)
;
;      Estimate the CMB power spectrum from the integrated power
;      spectrum intcl with mask:
;      toucanCl = mrs_intcl_to_cl(Map, intCl, bkj,
;      Mask=Mask. Firstl=2, Niter=300)
;         
; HISTORY:
;       Written:  Aurele Balavoine, 2012
;----------------------------------------------------------------------------------------------------------

 ;==============================================

function uwt1d, Signal, Nscale=Nscale, gen2=gen2, info=info, OptWT=OptWT
if not keyword_set(OptWT) then W = star1d(Signal, gen2=gen2, Nscale=Nscale) $
else begin
 mr1d_trans, Signal, info, OPT=OptWT, /nodel
 W=info.coef
end
return, w
end

 ;==============================================
 
function iuwt1d, Trans,  gen2=gen2, info=info
if not keyword_set(info) then  Rec = istar1d(Trans, gen2=gen2) $
else begin
  info.coef = Trans
  mr1d_recons, info, Rec
end
return, Rec
end

;;============================================== mrs_intcl_to_cl ===================================================

function mrs_intcl_to_cl, Map, intCl, bkj, Mask=Mask, Firstl=Firstl, Bmat=Bmat, Niter=Niter, initl=initl, WTNscale=WTNscale, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos, Mcl=Mcl, Cl=Cl, Proj_Resi=Proj_Resi, True_Resi=True_Resi, part_nMSE=part_nMSE

COMMON C_PLANCK

Result = -1
if N_PARAMS() LT 3  then begin 
   print, 'CALLING SEQUENCE: Res = mrs_intcl_to_cl(Map, intCl, bkj, Mask=Mask, Firstl=Firstl, Bmat=Bmat, Niter=Niter, initl=initl, WTNscale=WTNscale, SWT=SWT, HWT=HWT, OptWT=OptWT, Pos=Pos, Mcl=Mcl, Cl=Cl, Proj_Resi=Proj_Resi, True_Resi=True_Resi, part_nMSE=part_nMSE)'
   goto, DONE
endif

npix = N_elements(Map)
nside = npix2nside(npix)
Nscale = (size(bkj))[2]
lmax = (size(bkj))[1]-1

;;---------------------------------------- Set default parameters ----------------------------------------
if not keyword_set(Niter) then Niter = fix(10*exp(float(Nscale)/10)) ; Niter = 500 ;

;; if we don't have many l, we don't want to use the DCT, and we consider more wavelet scales
if lmax lt 1000l then OWT = 1

if not keyword_set(WTNscale) then begin
   if not keyword_set(OWT) then  WTNscale = fix(alog(float(lmax)))  - 1 $
   else  WTNscale = fix(alog(float(lmax)))  + 1
endif

if not keyword_set(Firstl) then Firstl = 0
if not keyword_set(initl) then initl = min([max([100,2*Firstl]), lmax/2., lmax-100])

if not arg_present(gen2) then gen2 = 1
if not arg_present(SWT) then SWT = 1
if not arg_present(HWT) then HWT = 1

;;------------------------------------------ Print Toucan options ------------------------------------------------------

print,  "SWT = ", SWT
print, "HWT = ", HWT

if not keyword_set(Bmat) then Bmat = mrs_compute_Bmat(bkj)

;;---------------- Normalizing the input --------------------
sigmaSpec = sigma(intCl)
intCl = intCl/sigmaSpec
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
Nsigma = 7 
for j = 0, WTNscale-1 do BEGIN
    if j NE 0 then LevelDetect(j) = NSigma * TabPSWTNorm(j)  $
    else LevelDetect(j) = (NSigma+1) * TabPSWTNorm(j)
 endfor
print, 'threshold chosen'
;;------------------------------------------------------------

;;------------ Set the mask for the HWT ----------------
if not keyword_set(Mcl) then begin
   print, 'Compute estimated Power Spectrum with Master to estimate HWT mask...'
   if not keyword_set(mask) then Mcl = mrs_powspec(map, lmax=lmax) else Mcl = mrs_master_powspec(map,mask,lmax=lmax)
endif
Signal =double(Mcl[Firstl:*])
SigmaSignal = sigma(Signal)
Signal =  Signal  / SigmaSignal

l = findgen(sigl) + firstl & lfactor = l * (l+1) / 2. / !dpi
if firstl EQ 0 then l[0]=1
if firstl EQ 0 then lfactor[0]=1
LastThreshold = 0
FirstThreshold =  max( sqrt(Signal) * lfactor)
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
      Level     = LevelDetect[j]
      
      ind       = where(abs(Scale) GT Level , count) ; Sig should be 1 in the GuassData space 
      if count GT 0 then ScaleMask[ind] = 1  
      WTMask[*,j] = ScaleMask  
   END 
   print, 'HWT mask computed'
endif

;;---------------------------------------------------------------

;;------------------ Initialization --------------------------
Resi =  (transpose(Bmat) ## intCl)/(sqrt(tau)) ; Residual for the gradient step
Resi = Resi[*]
SigmaResi_old = sigma(Resi)
Result = fltarr(sigl) ; Initial guess for the power spectrum
step = 1
StopDecrease = 0
proj_resi = fltarr(Niter)
true_resi = fltarr(Niter)
part_nMSE = fltarr(Niter)

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
   Resi =  (transpose(Bmat) ## (intCl - Bmat ## Z) )/tau
   Resi = Resi[*]  
   
   ;;========= WT - HARD Thresholding

   if keyword_set (HWT)  then begin
      MR_Trans = uwt1d(Resi, gen2=gen2, Nscale=WTNscale, info=info, OptWT=OptWT)
      if iter EQ 1 then MaxHWT = max(ABS(MR_Trans[*,0: WTNscale-2]))

      for j = 0, WTNscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
         Scale         = MR_Trans[*,j]   
         ScaleMask     = WTMask[*,j]
         MR_Trans[*,j] = Scale * ScaleMask
      END 
      Result =  Result  + iuwt1d(MR_Trans, gen2=gen2, info=info)
   endif
   
 
   ;;========= WT - SOFT Thresholding

   if keyword_set(SWT) then begin  
      MR_Trans = star1d(Result*lfactor, gen2=gen2, Nscale=WTNscale)
      if iter EQ 1 then MaxWT = max(ABS(MR_Trans[*,0: WTNscale-2]))

      for j = 0, WTNscale-2 do begin 
         SoftThreshold = MaxWT  *  Lambda *  TabPSWTNorm[j] 
         Scale         = reform(MR_Trans[*,j]) 
         softthreshold, Scale, SoftThreshold
         MR_Trans[*,j] = Scale 
      endfor
      Result = istar1d(MR_Trans, gen2=gen2)/lfactor
   endif
  
   ;;========= Positivity constrain
   if keyword_set(Pos) then begin
      ind = where(Result LT 0, c)
      if c GT 0 then Result[ind] = 0 
   endif
   
   sigmaResi = sigma(Resi)
   proj_resi(iter-1) = sigmaResi
   if keyword_set(cl) then begin
      true_resi(iter-1) = sigma(result[0 :lmax-firstl]*sigmaSpec - cl[firstl :lmax])
      part_nMSE(iter-1) = total(((cl[initl :lmax]-result[initl-Firstl :lmax-firstl]*sigmaSpec)/cl[initl :lmax])^2)
   endif
  
;;    print, Iter, '/', Niter, ':  Lambda = ', Lambda, ' | SigmaResi = ', sigma(Resi), ' | SigmaResult = ', sigma(Result)    
;;    print, Iter, ':  Linf(ResiDiff)= ', max(abs(Resi-Resi_old))

   if sigmaResi_old lt sigmaResi-(1.e-8) then begin
      print, 'Increasing residual -> iteration interrupted ==> try fewer iterations'
;;       interupted = 1
;;       goto, DONE
   endif
   sigmaResi_old = sigmaResi

endfor

result = [fltarr(Firstl),result*sigmaSpec]

DONE:
return, Result

end
