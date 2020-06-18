;+
; NAME:
;        mrs_tousi
;
; PURPOSE:
;  Apply a sparse denoising on a spectrum using wavelets, DCT and variance stabilisation.
;  If a noise power spectrum is given, then the output denoised power spectrum will take it into account.
;  By default, the Cl starts at l=0. If the input keyword NoMonoDip is set, then the input Cl[0] corresponds to l=2.
;
; CALLING:
;     DenCl = mrs_tousi(ClSpec, [ClCrosSPec, ] Niter=Niter, NScale=NScale, Nsigma= Nsigma, NoisePs=NoisePs, DenNoisePS=DenNoisePS, NoMonoDip=NoMonoDip)
;
; INPUTS:
;     DenCl -- IDL 1D array = Cl array to denoise  of a healpix map: Input image to be eroded 
;    
; OUTPUTS:
;     DenCl -- IDL 1D array  = Denoised Cl array
;
; INPUT KEYWORD: 
;     Niter: int  = Number of iterations, default is 50.
;     Nsigma: float = Detection level. Default is 7.
;     NoisePS: IDL 1D array  = input noise power spectrum.
;
; OUTPUT KEYWORD: 
;     DenNoisePS: IDL 1D array  = output denoised noise power spectrum
;
; EXAMPLE:
;      Remove the noise from a CMB power spectrum.
;      Denoised_Cl =  mrs_get_cl_theo_powspec(Cl)
;
;      Remove the noise from a CMB power spectrum, taking into account Cosmic Variance + Instrumental noise
;      DenNoisePS will contains the denoised spectrum of the instrumental noise.
;      Denoised_Cl =  mrs_get_cl_theo_powspec(Cl, NoisePS=NoisePS, DenNoisePS=DenNoisePS)
;         
; HISTORY:
;       Written: Paniez  Paykari & Jean-Luc Starck, 2011
;--------------------------------------------------------------------------------------------------------
  
 ;===========================
  
PRO tousi_denoise_cross, InputPowSpec, Result, Nscale=Nscale, NSigma=NSigma, FirstL=FirstL, BestBasis = BestBasis,  $
                  Niter_Den=Niter_Den, HWT=HWT, SWT=SWT,  gen2=gen2, Mask= Mask,  HDCT=HDCT, SDCT=SDCT, NoisePS=NoisePS, OptWT=OptWT, WTMask=WTMask, Pos=Pos, TabInterval=TabInterval, NoDebias=NoDebias, DCTMask=DCTMask
                      
IF N_PARAMS() LT 2 THEN BEGIN        
PRINT,'2 inputs needed.'  
GOTO, CLOSING 
ENDIF


if not keyword_set(NSigma)    then NSigma=7
if not keyword_set(Niter_Den) then Niter_Den=50
if not keyword_set(NoisePS)   then NoisePS=0.
if not keyword_set(FirstL)   then FirstL =0
vs = size(InputPowSpec)
Nl = vs[1]
if not keyword_set(Nscale) then Nscale = fix(  alog(float(Nl)) / alog(2.))  - 2
Signal =double( InputPowSpec)
SigmaSignal = sigma(Signal)
Signal =  Signal  / SigmaSignal
 
; TabPSWTNorm = [0.813558, 0.363782, 0.224453,0.137378,0.100910,0.0842879,0.114327,0.0699925, 0.0275717,0.0139204]
Dirac = Signal * 0.
vs = size(dirac)
dirac[Nl/2] = 1
w  = uwt1d(dirac, Nscale=Nscale, Gen2=Gen2, info=info, OptWT=OptWT)
TabPSWTNorm = dblarr(Nscale)
for j=0, Nscale-1 do TabPSWTNorm[j] = sqrt( total(w[*,j]^2) )

DataResult = Signal * 0
 Np = (size(Signal))(1)
l = findgen(Np) + firstl & lfactor = l * (l+1) / 2. / !dpi 
if firstl EQ 0 then l[0]=1
if firstl EQ 0 then lfactor[0]=1

;------------------Threshold Level---------------------------
LevelDetect = fltarr(Nscale)
for j = 0, Nscale-1 do BEGIN
    if j NE 0 then LevelDetect(j) = NSigma * TabPSWTNorm(j)  $
    else LevelDetect(j) = (NSigma+1) * TabPSWTNorm(j) 
END

;------------------------------------------------------------

GaussData = Signal
  
; plot, GaussData
;atwt1d, GaussData, MR_Data, Nscale=Nscale; performs WT - MRS/idl/STAT/at1dwt.pro

;------------------Denoising by Iteration--------------------
Resi   = GaussData ; Residual = Signal - Reconstructed Signal (for first iteration: Reconstructed Signal = 0)
Result = dblarr(Np)
 
LastThreshold=0
if  keyword_set(SDCT)  then begin 
         if not keyword_set(BestBasis) then Coeff = transdico(Signal* lfactor)  $
         else Coeff = bb_dct(Signal* lfactor, LDCT=LDCT, UseFFT=UseFFT,TabInterval=TabInterval)
         FirstThreshold = max(Coeff[1:*])
 end else  FirstThreshold=  max(sqrt( Signal) * lfactor)
 DeltaThreshold = FirstThreshold - LastThreshold
 DCTDeltaThreshold =  FirstThreshold  

indl = 2. * l + 1.
if firstl EQ 0 then indl[0]=1

MR_SignalTrans = uwt1d(Signal, gen2=gen2, Nscale=Nscale, info=info, OptWT=OptWT)

StopDecrease=0
for iter = 1, Niter_Den do BEGIN
   
    if StopDecrease EQ 0 then  Lambda =   (1.-erf(2.8*float(iter) / float(Niter_Den)))
       if iter EQ Niter_Den then Lambda = 0
    ;========= WT - HARD Thresholding
    
    Resi = GaussData - Result   
    
    if keyword_set (HWT)  then begin
    print, "HWT"
        MR_Trans = uwt1d(Resi, gen2=gen2, Nscale=Nscale, info=info, OptWT=OptWT)
        for j = 0, Nscale-2 do  MR_Trans[*,j] = MR_Trans[*,j] * WTMask[*,j]
        Result =  Result  + iuwt1d(MR_Trans, gen2=gen2, info=info)
    endif  
  
  ;========= WT - SOFT Thresholding

    if keyword_set(SWT) then begin  
    print, "SWT"
         MR_Trans = star1d(Result, gen2=gen2, Nscale=Nscale)
         if Iter EQ 1 then MaxWT = max(ABS(MR_Trans[*,0: Nscale-2]))  
         for j = 0, Nscale-2 do BEGIN
             SoftThreshold =   MaxWT  *  Lambda *  TabPSWTNorm[j]
             Scale  = reform(MR_Trans[*,j])
             softthreshold, Scale, SoftThreshold
             MR_Trans[*,j] = Scale
       END
       Result = istar1d(MR_Trans, gen2=gen2)
    END
    
          ;========= DCT - HARD  Thresholding
          
    Resi = GaussData - Result   
    if  keyword_set(HDCT)  then begin 
    print, "HDCT"
        Resi = Resi  * lfactor
       if not keyword_set(BestBasis) then Coeff = transdico(Resi)  $
       else Coeff = bb_dct(Resi, LDCT=LDCT, UseFFT=UseFFT,TabInterval=TabInterval)
       CoeffT = Coeff * DCTMask
       if not keyword_set(BestBasis) then  Resi    = (transdico(CoeffT, /inverse))   $ 
       else   Resi = bb_dct(CoeffT, LDCT=LDCT, /inverse, UseFFT=UseFFT)
       Result = Result + Resi /  lfactor
       Resi = GaussData - Result   
   endif   
   
    ;========= DCT - SOFT Thresholding
    
    if  keyword_set(SDCT)  then begin 
    print, "SDCT"
          Result = Result  * lfactor
          if not keyword_set(BestBasis) then Coeff = transdico(Result)  $
         else Coeff = bb_dct(Result, LDCT=LDCT, UseFFT=UseFFT)
          Level =   DeltaThreshold  * Lambda
         softthreshold, Coeff, Level
         if not keyword_set(BestBasis) then  Result    = (transdico(Coeff, /inverse))   $ 
        else   Result = bb_dct(Coeff, LDCT=LDCT, /inverse, UseFFT=UseFFT)
        Result = Result /  lfactor
      endif   
  
   ;===========================
   ; Positivity constrain
   if keyword_set(Pos) then begin
        ind = where(Result LT 0, c)
       if c GT 0 then Result[ind] = 0 
   end
   
   Resi = GaussData - Result   
   print, Iter, ':  Lambda = ', Lambda, ', SigmaResi = ', sigma(Resi), ', SigmaResult = ', sigma(Result)
   plotcl, Result
    ;  print, SigmaResi
    ; if SigmaResi LT 1 then StopDecrease = 1  ; GOTO, ENDLOOP
    ; if SigmaResi LT 1 then  GOTO, ENDLOOP
END

ENDLOOP:

 Result =  Result  * SigmaSignal
 
CLOSING: 
RETURN 
END

 ;===========================

function tousi_cross, InputPowSpec, NoMonoDip=NoMonoDip,  Niter=Niter, NScale=NScale, Nsigma=Nsigma,  SDCT= SDCT,  OptWT=OptWT,  SWT=SWT, WTMask=WTMask, BestBasis = BestBasis,  OWT=OWT, ODCT=ODCT, TabInterval=TabInterval, DCTMask=DCTMask

if not keyword_set(Niter) then Niter = 50
if not keyword_set(Nsigma) then Nsigma = 7
D = InputPowSpec
MonoDipRem = 0
FirstL = 2
ResNoise = 0

if keyword_set(NoMonoDip)  then begin
; In this case the input data does NOT contain l=0 and l=1 values.
  FirstUseIndL = 1
  DN=0
  D = D[FirstUseIndL:*]
  if keyword_set(NoisePS) then DN = NoisePs[FirstUseIndL:*]
end else begin
   FirstUseIndL = FirstL
   D = D[FirstUseIndL:*]
   DN = 0
   if keyword_set(NoisePS) then DN = NoisePs[FirstUseIndL:*]
end

; if we don't have many l, we don't want to use the DCT, and we consider more wavelet scales
vs = size(D)
Nl = vs[1]
if Nl LT 1000L then  OWT = 1
 
if not keyword_set(Nscale) then begin
   if not keyword_set(OWT) then  Nscale = fix(  alog(float(Nl)))  - 1 $
   else  Nscale = fix(  alog(float(Nl)))  + 1
end
  
if not keyword_set(ODCT) then HWT=1
if not keyword_set(OWT) then HDCT=1
if keyword_set(BestBasis) then BestBasis = 1

tousi_denoise_cross, D, Res, Niter_Den= Niter, FirstL= FirstL , Nsigma=Nsigma, Nscale= Nscale,  /gen2,  $
             HWT=HWT,  SWT=SWT, HDCT=HDCT, SDCT=SDCT, NoisePS= ResNoise, $
             BestBasis=BestBasis, OptWT=OptWT, WTMask=WTMask, /Pos, DCTMask=DCTMask

if FirstUseIndL GE 1 then Res =[InputPowSpec[0:FirstUseIndL-1], Res] 

return, res
end
 
  ;===========================
   ;===========================

 function mrs_tousi, InputPowSpec, InputCrossSpec, NoMonoDip=NoMonoDip,  Fit=Fit, Niter=Niter, NScale=NScale, Nsigma=Nsigma, NoisePs=NoisePs, SolFit=SolFit,  SDCT= SDCT, DenNoisePS=DenNoisePS, OptWT=OptWT,  SWT=SWT, WTMask=WTMask, BestBasis = BestBasis, Pos=Pos, OWT=OWT, ODCT=ODCT, TabInterval=TabInterval, DCTMask=DCTMask 
  
if not keyword_set(NFit) then NFIT = 5
if not keyword_set(Niter) then Niter = 50
if not keyword_set(Nsigma) then Nsigma = 7
Res=-1

  if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: Res =  mrs_tousi(InputPowSpec, [ InputCrossSpec, ] NoMonoDip=NoMonoDip,  Fit=Fit, Niter=Niter, NScale=NScale, Nsigma=Nsigma, NoisePs=NoisePs, SolFit=SolFit,  SDCT= SDCT, DenNoisePS=DenNoisePS, OptWT=OptWT,  SWT=SWT, WTMask=WTMask, BestBasis = BestBasis, Pos=Pos, OWT=OWT, ODC=ODCT, TabInterval=TabInterval, DCTMask=DCTMask)'
        goto, DONE
        end

Res =  mrs_get_cl_theo_powspec(InputPowSpec, NoMonoDip=NoMonoDip,  Fit=Fit, Niter=Niter, NScale=NScale, Nsigma=Nsigma, NoisePs=NoisePs, SolFit=SolFit,  SDCT= SDCT, DenNoisePS=DenNoisePS, OptWT=OptWT,  SWT=SWT, WTMask= WTMask, BestBasis = BestBasis, Pos=Pos, OWT=OWT, ODC=ODCT, TabInterval=TabInterval, DCTMask=DCTMask)
  
  ; If we want to denoise the cross-spec, then we use the mask WTMask and DCTMask computed on the spectrum to estimate the cross-spectrum.
 if N_PARAMS() GT 1 then begin 
      Res = tousi_cross(InputCrossSpec, NoMonoDip=NoMonoDip,  Niter=Niter, NScale=NScale, Nsigma=Nsigma,  /SDCT,  OptWT=OptWT, /SWT, WTMask=WTMask, BestBasis = BestBasis,   OWT=OWT, ODC=ODCT, TabInterval=TabInterval, DCTMask=DCTMask)
 end

DONE:

return, Res
end
 
  ;===========================
 ; compile first :   .r   mrs_get_cl_theo_powspec.pro
 
 pro run_ex_tousi
 
c   = getcmb(nside=512, cl=cl)
n1 = randomn(seed, 512L^2*12) * 50.
n2 = randomn(seed, 512L^2*12) * 50.

; make two CMB observations of the same realization
d1 = c + n1
d2 = c + n2

; Full data set
d = 0.5*(d1 + d2)
; Jack Knife map
jk = 0.5 *(d1-d2)

; TOUSI Estimation
pd = mrs_powspec(d)
pjk = mrs_powspec(jk)
tpd = mrs_tousi(pd, NoisePs=pjk)

; TOUSI Cross spectrum estimation
 pd1d2 = mrs_crosspowspec(d1,d2)
 tpd1d2 = mrs_tousi(pd, pd1d2)
  
 end
  
 ;===========================

