PRO mrs_pswtfil, Signal, Result, bias=bias, Nscale=Nscale, SigmaNoise=SigmaNoise, $
NSigma=NSigma, NbrIter=NbrIter, mask=mask, firstscale=firstscale, beam=beam, $
IdealBeam=IdealBeam, L1Regul1=L1Regul1, L1RegulDCT=L1RegulDCT, Verb=Verb, PSFBeam=PSFBeam, Eps=Eps, pos=pos, LastThreshold= LastThreshold, WTDen= WTDen

;+ 
; NAME: 
;        MRS_PSWTFIL
;
; PURPOSE: 
;   Filter a 1D power spectrum by the multiresolution support method and with a variance stabilisation.
;
; CALLING SEQUENCE: 
;   mrs_pswtfil, Signal, Result, Nscale=Nscale, SigmaNoise=SigmaNoise, NSigma=NSigma, NbrIter=NbrIter,
;   mask=mask, nobiais=nobiais, firstscale=firstscale, beam=beam, IdealBeam=IdealBeam, L1Regul=L1Regul, Verb=Verb 
;
; INPUTS: 
;   Signal -- 1D IDL array: power spectrum to filter. MUST be positive. Cl not l(l+1)Cl.
;    
; KEYED INPUTS: 
;   Nscale -- scalar: number of scales (default is 5)
;   SigmaNoise -- float: noise standard deviation in case of Gaussian noise.
;                      The default value is 1. This is on the Signal not map.
;   NSigma -- float: detection at NSigma at each scale (default is 3)
;   NbrIter -- int: number of iterations (default is 6)
;   firstscale -- int: first detection scale in the wavelet space
;   L1Regul -- scalar: if set, a l1 regularization is applied on a solution
;   bias -- scalar: if set, noise bias correction is applied - bias in the Cl space
;   beam -- 1D IDL float array: PSF corrupting the input power spectrum. Default is 1.0 (no PSF)
;   IdealBeam -- 1D IDL float array: ideal beam applied to the power spectrum. Default is 1.0 (no idealbeam)
;   Verb -- scalar: if set, print results on the screen.
;
; OUTPUTS: 
;   Result -- 1D IDL array: filtered power spectrum
;
; KEYED OUTPUTS: 
;   mask  -- 2D IDL array: position where signal is found at each scale.
;
; MODIFICATION HISTORY: 
;    8-Jan-1996  written by JL Starck with template_gen 
;-
     
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------

IF N_PARAMS() LT 2 THEN BEGIN
PRINT,'2 inputs needed.' 
GOTO, CLOSING
ENDIF

IF keyword_set(L1Regul1) and keyword_set(L1RegulDCT) THEN BEGIN
PRINT,'set either L1Regul1 or L1RegulDCT'
GOTO, CLOSING
ENDIF

TabPSWTNorm = [1.02228,0.407707,0.254645,0.168092,0.106421,0.0807738,0.0515747,0.0393627,1.20018]                
TabLevel = TabPSWTNorm

if not keyword_set(Nscale) then Nscale=5
if not keyword_set(NSigma) then NSigma=3;3sigma detection
if not keyword_set(NbrIter) then NbrIter=6
if not keyword_set(SigmaNoise) then SigmaNoise=1.
if not keyword_set(firstscale) then firstscale=1.
if not keyword_set(PSFBeam) then  PSFBeam=1.
if not keyword_set(IdealBeam) then IdealBeam=1.
if not keyword_set(Eps) then Eps =1.

FS=firstscale-1
 DataResult = Signal
 DataResult[*] = 0
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------

Np = (size(Signal))(1)
Mask = bytarr(Nscale, Np)
Scale = fltarr(Np)

Mask[Nscale-1, *] = 1.

;------------------Threshold Level---------------------------
LevelDetect = fltarr(Nscale)
for j = 0, Nscale-1 do BEGIN
   if j NE 0 then LevelDetect(j) = NSigma * TabLevel(j)  $
   else LevelDetect(j) = (NSigma+1) * TabLevel(j) 
END
;------------------------------------------------------------

ll = findgen(Np) + 2
lfactor = ll * (ll+1) / 2. / !pi

;------------------Take to Gauss Space-----------------------
NoisePowerSpectrum = SigmaNoise
    
   P = Signal
   
   GaussData = mrs_variance_stabilization(P, mu=mu, psi1=psi1)
;------------------------------------------------------------

atwt1d, GaussData, MR_Data, Nscale=Nscale; performs WT - MRS/idl/STAT/at1dwt.pro

;------------------Denoising by Iteration--------------------
Resi = GaussData ; Residual = Signal - Reconstructed Signal (for first iteration: Reconstructed Signal = 0)
Result = fltarr(Np)
for Iter = 1, NbrIter do BEGIN
    if keyword_set(WTDen) then begin
       atwt1d, Resi, MR_Trans, Nscale=Nscale
       for j = FS, Nscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
           Scale = MR_Trans[j,*]
           ScaleMask = Mask[j,*]
           ind = where(abs(Scale) GT LevelDetect(j), count)
          if count GT 0 then ScaleMask(ind) = 1
          Mask[j,*] = ScaleMask
      END ; for j = FS, Nscale-2
       MR_Trans = MR_Trans * Mask
       ; Reconstruction by summation over j - page 22/23 of notes
       rec = reform(MR_Trans[0,*])
       for j = 1, Nscale-1 do rec = rec + reform(MR_Trans[j,*])
    ; mrs_info, rec, mes=' REC '
      Result = Result + rec 
    end else Result = Result + Resi 
 
    ; L1 norm regularization on the solution
    if keyword_set(L1Regul1) then begin
       SigmaSignal = sigma(Result);-160.
       ; help, result, MR_Trans
       atwt1d, Result, MR_Trans, Nscale=Nscale
       Result = reform(MR_Trans[Nscale-1,*])
       for j = 0, Nscale-2 do BEGIN
           Level = TabLevel(j) * SigmaSignal * (NbrIter-Iter) / NbrIter
           Scale = reform(MR_Trans[j,*])
           softthreshold, Scale, Level
           Result = Result + Scale
       END ; for j
    end   ; end if keyword_set(L1Regul1)
    ;plot, GaussData
    ;oplot, Result, line=2
    ; mrs_info, Result, mes=' Result '
    ; if iter EQ 1 then save, filename='pb.sav', Rec, result, mu, psi
    PS_Estim = mrs_inv_variance_stabilization(Result, mu, psi1) 
    
    
    ;plotcl, P
    ;oplotcl, PS_Estim
    ;nn = p
    ;nn(*) = 1.
    ;oplotcl, nn
    ; mrs_info, PS_Estim, mes=' PS_Estim '
    ;print,PS_Estim
    
    ;++++++++++++++L1RegulDCT++++++++++++++++++++
    ; DCT L1 norm regularization on the solution
    ; but in l(l+1)/2pi * Cl
    ll = findgen(Np) + 2
    lfactor = ll * (ll+1) / 2. / !pi
    if keyword_set(bias) then PS_Estim =  PS_Estim - bias
    if keyword_set(Beam) then PS_Estim = PS_Estim / PSFBeam^2.
 
    if keyword_set(L1RegulDCT) then begin
       SigmaSignal = sigma(PS_Estim)
       Coeff = dct(PS_Estim * lfactor, T=C)
       Level = SigmaSignal * (NbrIter-Iter) / NbrIter 
       if keyword_set(LastThreshold) and (Level LT LastThreshold) then Level = LastThreshold
       softthreshold, Coeff, Level
       PS_Estim = dct(Coeff, C, /inverse) / lfactor
    end
    
    if keyword_set(Pos) then begin
       ind = where(PS_Estim LT 0, c)
       if c GT 0 then PS_Estim[ind] = 0
    end
    ;++++++++++++++++++++++++++++++++++++++++++++   
    
    if keyword_set(Beam) then  PS_Estim = PS_Estim* PSFBeam^2.

    if keyword_set(bias) then begin  
       DataResult = (1.-Eps)*DataResult  + Eps * (PS_Estim  + bias)
    endif else begin 
       DataResult = (1.-Eps)*DataResult  + Eps * PS_Estim  
       ;DataResult = PS_Estim + NoisePowerSpectrum
    endelse
    
    ;Oplot,Dataresult,color=120
    ;print,dataresult
    ; mrs_info, DataResult, mes=' DataResult '
    Result = mrs_variance_stabilization(DataResult , mu=mu, psi1=psi1) * PSFBeam^2. / IdealBeam^2.
    resi = GaussData - Result
    ;mrs_info, resi, mes=' RESI '
    ;if keyword_set(print) then $
    ;print, "Iter ", Iter, ",  Xi^2 / Np = ", (sigma(resi))^2/SigmaNoise^2
    if keyword_set(Verb) then print, "Iter ", Iter, ",  SigmaResi = ", sigma(resi)
    ;spec, resi
    ;wait,3
END ; for Iter = 1, NbrIter
    ;plotcl,ps_estim
;???????????? correct for PSFbeam and IdealBeam   ????????????
if keyword_set(bias) then begin  
   Result = DataResult - bias
endif else begin 
   Result = DataResult 
endelse
 if keyword_set(Beam) then Result = Result / PSFBeam^2.
    if keyword_set(Pos) then begin
       ind = where(Result LT 0, c)
      if c GT 0 then Result[ind] = 0
   end
; info, Result
plotcl, Result
;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
CLOSING:
RETURN 
END
