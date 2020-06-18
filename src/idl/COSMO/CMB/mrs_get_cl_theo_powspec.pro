;+
; NAME:
;        mrs_get_cl_theo_powspec
;
; PURPOSE:
;  Apply a sparse denoising on a spectrum using wavelets, DCT and variance stabilisation.
;  If a noise power spectrum is given, then the output denoised power spectrum will take it into account.
;  By default, the Cl starts at l=0. If the input keyword NoMonoDip is set, then the input Cl[0] corresponds to l=2.
;
; CALLING:
;     DenCl = mrs_get_cl_theo_powspec(Cl,  Niter=Niter, NScale=NScale, Nsigma= Nsigma, NoisePs=NoisePs, DenNoisePS=DenNoisePS, NoMonoDip=NoMonoDip)
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
 
 ;================================

function transdico, t, inverse=inverse, UseFFT=UseFFT, Dirac=Dirac
vs = size(t)
np = vs[1]

if keyword_set(Dirac) then Ret = t $
else begin
if keyword_set(UseFFT) then begin
     if not keyword_set(inverse) then begin
        Ret = dft(t)
        RealPart = double(Ret)
        ImaPart = double (imaginary(Ret))
        Ret=[RealPart, ImaPart]
        Ret = Ret / sqrt( double(Np)) * sqrt(2.)
     end else begin
        Np = Np / 2
        T = T * sqrt( double(Np)) / sqrt(2.)
        RealPart = t[0:np-1]
        ImaPart = t[np:*]
        cf = complex(RealPart, ImaPart, /double)
        z = dfti(cf)
        Ret = double(z)
end
end else Ret = dct(t, inverse=inverse)
END

return, Ret
end
    
;================================
 ; Compute the best DCT basis using top botton approach and a l1 criterion
pro get_best_basis, cl, lmin, last, nresol

if nresol GT 0 then begin
mid = (lmin+last)/2
Cl1 = cl[lmin:mid-1]
Cl2 = cl[mid:last]
Cl3 = cl[lmin:last]
D1 = transdico(cl1)
D2 = transdico(cl2)
D = transdico(Cl3)
Norm1 = total(abs(D1))
Norm2 = total(abs(D2))
NormD = total(abs(D))
if NormD LT Norm1+Norm2 then print, "BASIS: ", lmin, ' to ', last $
else begin
get_best_basis, cl, lmin, mid-1, nresol-1
get_best_basis, cl, mid, last, nresol-1
end  
end else   begin
print, "BASIS: ", lmin, ' to ', last 
end
end
  
;================================

function get_ldct_filter, TabInterval=TabInterval,  lmax=lmax 
if not keyword_set(lmax) then lmax = 3200
; if not keyword_set(TabInterval) then TabInterval =[92, 186, 279, 1500, 2250, Lmax]
; if not keyword_set(TabInterval) then TabInterval =[1500, 2250, Lmax]
if not keyword_set(TabInterval) then TabInterval =[70, 512, 1024, 1500, 2000, 2500, 3000]

 ind = where (TabInterval LT lmax) 
 TabInterval = [TabInterval[ind], lmax]
 
; if not keyword_set(TabInterval) then TabInterval =[2048, 3072, 3584, Lmax]  pour lmax = 4096 (norm l1)
; if not keyword_set(TabInterval) then TabInterval =[1500, 2250, 2625, Lmax]  pour lmax = 3000 (norm l1)

Np = Lmax + 1
 get_wp_meyer_filter, interval=TabInterval,  Filter, win, lmax=Lmax
; help, filter
vs = size(TabInterval)
 nb = vs[1]
 TabFirst = intarr(Nb)
 TabLast = intarr(Nb)
 TabFilter = dblarr(Np, Nb)
 for i=0, nb-1 do begin
    b = i+1
    if i EQ 0 then F = Filter[*,0] + Filter[*,1] $
    else F = Filter[*,b]
    TabFilter[*,i] = F
    ind = where (F NE 0, c)
    First = ind[0] - 1
    if First LT 0 then First = 0
    vs = size(ind)
    Last = ind[ vs[1]-1] + 1
    if Last EQ Np then Last = Last -1 
    TabFirst[i] = First
    TabLast[i] = Last
 end
 LDCT = {Nb:Nb, Np:Np, TabFilter: TabFilter, Lmax: Lmax, TabFirst:TabFirst, TabLast:TabLast}
 return, LDCT
 end


;================================
; Best Basis Coefficient soft thresholding
pro bb_softthreshold, Data, Lambda,  LDCT=LDCT
if not keyword_set(LDCT) then begin
    vs = size(Data)
    lmax = vs[1]-1
    LDCT = get_ldct_filter(TabInterval=TabInterval,  lmax=lmax)
 end

 for b=0, LDCT.nb-1 do begin
    F = LDCT.TabFilter[*,b]
    First= LDCT.TabFirst[b]
    Last = LDCT.TabLast[b]
    Local_Cl = Data[First:Last]
    C0 = Local_Cl[0]
    softthreshold, Local_Cl, Lambda
    ; Local_Cl[0] = C0
    Data[First:Last] =  Local_Cl    
end 
end

;================================

 function bb_dct, Cl,  inverse=inverse, TabInterval=TabInterval,   LDCT=LDCT, UseFFT=UseFFT, Dirac=Dirac 
 
 if not keyword_set(LDCT) then begin
    vs = size(cl)
    lmax = vs[1]-1
    ; print, lmax
    LDCT = get_ldct_filter(TabInterval=TabInterval,  lmax=lmax)
 end

 
 IndDCT  =0L
 Step = 1
 if keyword_set(UseFFT) then Step=2
 if keyword_set(inverse) then Res = dblarr(LDCT.Np) 
 
 for b=0, LDCT.nb-1 do begin
    F = LDCT.TabFilter[*,b]
    First= LDCT.TabFirst[b]
    Last = LDCT.TabLast[b]
   ;  print, First , ' ,' , Last
    
    if not keyword_set(inverse) then begin
         Local_Cl = double(Cl) * sqrt(F)
         Local_Cl = Local_Cl[First:Last]
         DC  = transdico(Local_Cl, UseFFT=UseFFT, Dirac=Dirac)
        ;  help, Local_Cl
        if b EQ 0 then Res = DC $
         else  Res = [Res, DC]
     end else begin
        Np = Last - First + 1
        Local_Cl  = Cl[ IndDCT: IndDCT+ Step*Np-1]
       ; help, Local_Cl
       DC  = transdico(Local_Cl, /inverse, UseFFT=UseFFT, Dirac=Dirac)
        U = dblarr(LDCT.Np)
        U[First:Last] = DC
        Res = Res + U * sqrt(F)
        IndDCT = IndDCT + Np* Step
     end
 end
  
return, Res
end
 ;================================

function uwt1d, Signal, Nscale=Nscale, gen2=gen2, info=info, OptWT=OptWT
if not keyword_set(OptWT) then W = star1d(Signal, gen2=gen2, Nscale=Nscale) $
else begin
 mr1d_trans, Signal, info, OPT=OptWT, /nodel
 W=info.coef
end
return, w
end

 ;================================
 
function iuwt1d, Trans,  gen2=gen2, info=info
if not keyword_set(info) then  Rec = istar1d(Trans, gen2=gen2) $
else begin
  info.coef = Trans
  mr1d_recons, info, Rec
end
return, Rec
end

 ;================================

PRO cl_denoise, InputPowSpec, Result, Nscale=Nscale, NSigma=NSigma, FirstL=FirstL, BestBasis = BestBasis,  $
                  Niter_Den=Niter_Den, HWT=HWT, SWT=SWT,  gen2=gen2, Mask= Mask,  HDCT=HDCT, SDCT=SDCT, NoisePS=NoisePS, OptWT=OptWT, WTMask=WTMask, Pos=Pos, TabInterval=TabInterval, NoDebias=NoDebias, DCTMask=DCTMask
                      
IF N_PARAMS() LT 2 THEN BEGIN        
PRINT,'2 inputs needed.'  
GOTO, CLOSING 
ENDIF

   if  keyword_set(HDCT)  then print,  "HDCT"
   if  keyword_set(SDCT)  then print, "SDCT"

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
NoisePS = NoisePS / SigmaSignal

; TabPSWTNorm = [0.813558, 0.363782, 0.224453,0.137378,0.100910,0.0842879,0.114327,0.0699925, 0.0275717,0.0139204]
Dirac = Signal * 0.
vs = size(dirac)
dirac[Nl/2] = 1
w  = uwt1d(dirac, Nscale=Nscale, Gen2=Gen2, info=info, OptWT=OptWT)
TabPSWTNorm = dblarr(Nscale)
for j=0, Nscale-1 do TabPSWTNorm[j] = sqrt( total(w[*,j]^2) )

 ;The numbers above are specific to the type of WT used here. They are to get the Nsigma to the true Nsigma, as there are some
;normalisations in the WT which need to be corrected for.

Signal = Signal     
DataResult = Signal
DataResult[*] = 0
Np = (size(Signal))(1)
Mask = fltarr(Np, Nscale) 
Mask[*, Nscale-1] = 1. 
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

GaussData = mrs_variance_stabilization(Signal, mu=mu, psi1=psi1,firstl= firstl);
mu0 = mu * 0
if keyword_set(NoDebias) then mu0 = mu

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
NormL = (1.+1/indl)
 Fa =  -0.25839012       
 Fb = 1.1820104
 x = findgen(501) / 500.
 yy = Fa*x + Fb
 ind = where(yy lt 1, c)
 if c GT 0 then yy[ind] = 1.
 NormL[0:min([Nl-1,500])] = NormL[0:Min([Nl-1,500])] * yy
MR_SignalTrans = uwt1d(Signal, gen2=gen2, Nscale=Nscale, info=info, OptWT=OptWT)

StopDecrease=0
for iter = 1, Niter_Den do BEGIN
   
    if StopDecrease EQ 0 then  Lambda =   (1.-erf(2.8*float(iter) / float(Niter_Den)))
       if iter EQ Niter_Den then Lambda = 0
    ;========= WT - HARD Thresholding
    
    Resi = GaussData - Result   
    
    if keyword_set (HWT)  then begin
        MR_Trans = uwt1d(Resi, gen2=gen2, Nscale=Nscale, info=info, OptWT=OptWT)
     ;  HalfWin = 2L
        for j = 0, Nscale-2 do BEGIN ; last scale includes the mean (increasing j -> coarser scale) - hence ignored
           Scale     = MR_Trans[*,j]   
         ;  Win = HalfWin * 2 + 1
           
           ScaleMask = Mask[*,j]
           Level = NormL * LevelDetect[j]
  
           ; lrms = mr1d_get_local_rms(Scale, WindowSize =Win)
           ; ind       = where(lrms GT Level , count)
           ind       = where(abs(Scale) GT Level , count) ; Sig should be 1 in the GuassData space 
           if count GT 0 then ScaleMask[ind] = 1  
           Mask[*,j] = ScaleMask  
           MR_Trans[*,j] = Scale * ScaleMask
        ;   HalfWin = HalfWin * 2L
       END 
       Result =  Result  + iuwt1d(MR_Trans, gen2=gen2, info=info)
    endif  
  
  ;========= WT - SOFT Thresholding

    if keyword_set(SWT) then begin  
          WaveletInDirectSpace=1
       if keyword_set(WaveletInDirectSpace) then begin
              PS_Estim    = mrs_inv_variance_stabilization(Result, mu, psi1) 
              PS_Estim = PS_Estim * lfactor
              MR_Trans = star1d(PS_Estim, gen2=gen2, Nscale=Nscale)
      end else MR_Trans = star1d(Result, gen2=gen2, Nscale=Nscale)
       if Iter EQ 1 then MaxWT = max(ABS(MR_Trans[*,0: Nscale-2]))  
       for j = 0, Nscale-2 do BEGIN
           SoftThreshold =   MaxWT  *  Lambda *  TabPSWTNorm[j]
           Scale  = reform(MR_Trans[*,j])
           softthreshold, Scale, SoftThreshold
           MR_Trans[*,j] = Scale
       END
         if  keyword_set(WaveletInDirectSpace) then begin
              PS_Estim = istar1d(MR_Trans, gen2=gen2)
              PS_Estim = PS_Estim  /  lfactor
              Result = mrs_variance_stabilization(PS_Estim, mu=m, psi1=psi1,firstl= firstl) ;* PSFBeam^2. / IdealBeam^2.
         end else  Result = istar1d(MR_Trans, gen2=gen2)
    END
    
          ;========= DCT - HARD  Thresholding
          
    Resi = GaussData - Result   
   ;   plot, Resi[2:*], title='resi'
   ;   wait, 1
   if  keyword_set(HDCT)  then begin 
       PS_Estim = Resi
   
      if not keyword_set(BestBasis) then Coeff = transdico(PS_Estim)  $
      else Coeff = bb_dct(PS_Estim, LDCT=LDCT, UseFFT=UseFFT,TabInterval=TabInterval)
       
       DCT_HThresh = NSigma
       ind = where( ABS(Coeff) GT  DCT_HThresh, c)
       if Iter EQ 1 then  DCTMask = Coeff * 0
      if c GT 0 then DCTMask[ind] = 1
      
      CoeffT = Coeff * DCTMask
     ;  if  keyword_set(BestBasis) then  for b=0, LDCT.nb-1 do  CoeffT[LDCT.TabFirst[b]] = 0  $; Coeff[LDCT.TabFirst[b]]
      ; else CoeffT[0] = 0
      ; CoeffT[0] = 0
      ; if not keyword_set(BestBasis) then  PS_Estim    = dct(CoeffT, CC1,  /inverse)  $;  + l2ClR $
      ; else   PS_Estim = bb_dct(CoeffT, LDCT=LDCT, /inverse, UseFFT=UseFFT)
      if not keyword_set(BestBasis) then  PS_Estim    = (transdico(CoeffT, /inverse))   $ 
       else   PS_Estim = bb_dct(CoeffT, LDCT=LDCT, /inverse, UseFFT=UseFFT)

      Result = Result + PS_Estim
      Resi = GaussData - Result   
   endif   
   
    ;========= DCT - SOFT Thresholding
    
    if  keyword_set(SDCT)  then begin 
          PS_Estim    = mrs_inv_variance_stabilization(Result, mu0, psi1) 
          PS_Estim = (PS_Estim - NoisePS) * lfactor
        
          if not keyword_set(BestBasis) then Coeff = transdico(PS_Estim)  $
         else Coeff = bb_dct(PS_Estim, LDCT=LDCT, UseFFT=UseFFT)
        
        Level =   DeltaThreshold  * Lambda
        ; bb_softthreshold, Coeff, Level,  LDCT=LDCT
         softthreshold, Coeff, Level

        if not keyword_set(BestBasis) then  PS_Estim    = (transdico(Coeff, /inverse))   $ 
        else   PS_Estim = bb_dct(Coeff, LDCT=LDCT, /inverse, UseFFT=UseFFT)
  
       PS_Estim = PS_Estim /  lfactor
      Result      = mrs_variance_stabilization(PS_Estim+ NoisePS, mu= mu0, psi1=psi1,firstl= firstl) ;* PSFBeam^2. / IdealBeam^2.
     ; plotcl, PS_Estim
     endif   
  
   ;===========================
   ; Positivity constrain
   if keyword_set(Pos) then begin
     PS_Estim    = mrs_inv_variance_stabilization(Result, mu0, psi1) - NoisePS
     ind = where(PS_Estim LT 0, c)
     if c GT 0 then PS_Estim[ind] = 0 
     sResult      = mrs_variance_stabilization(PS_Estim+ NoisePS, mu= mu0, psi1=psi1,firstl= firstl) ;* PSFBeam^2. / IdealBeam^2.
   end
   
   Resi = GaussData - Result   
   print, Iter, ':  Lambda = ', Lambda, ', SigmaResi = ', sigma(Resi), ', SigmaResult = ', sigma(Result)
    ;  print, SigmaResi
    ; if SigmaResi LT 1 then StopDecrease = 1  ; GOTO, ENDLOOP
    ; if SigmaResi LT 1 then  GOTO, ENDLOOP
END

ENDLOOP:

Result = mrs_inv_variance_stabilization(Result, mu0, psi1)  - NoisePS
Result =  Result  * SigmaSignal
NoisePS = NoisePS * SigmaSignal
WTMask = Mask

CLOSING: 
RETURN 
END

 ;===========================
  
function mrs_get_cl_theo_powspec, InputPowSpec, NoMonoDip=NoMonoDip,  Fit=Fit, Niter=Niter, NScale=NScale, Nsigma=Nsigma, NoisePs=NoisePs, SolFit=SolFit,  SDCT= SDCT, DenNoisePS=DenNoisePS, OptWT=OptWT,  SWT=SWT, WTMask=WTMask, BestBasis = BestBasis, Pos=Pos, OWT=OWT, ODC=ODCT, TabInterval=TabInterval, DCTMask=DCTMask
if not keyword_set(NFit) then NFIT = 5
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
 
if keyword_set(NoisePs) then begin 
   cl_denoise, DN, ResNoise, Niter_Den= Niter, FirstL=FirstL , Nsigma=Nsigma, Nscale= Nscale,  /gen2,  /HWT, /SWT, /NoDebias
   if keyword_set(NoMonoDip) then  DenNoisePS = ResNoise $
   else DenNoisePS=[NoisePs[0:FirstL-1], ResNoise] 
end


if not keyword_set(ODCT) then HWT=1
if not keyword_set(OWT) then HDCT=1
if keyword_set(BestBasis) then BestBasis = 1

cl_denoise, D, Res, Niter_Den= Niter, FirstL= FirstL , Nsigma=Nsigma, Nscale= Nscale,  /gen2,  $
             HWT=HWT,  SWT=SWT, HDCT=HDCT, SDCT=SDCT, NoisePS= ResNoise, $
             BestBasis=BestBasis, OptWT=OptWT, WTMask=WTMask, Pos=Pos, DCTMask=DCTMask


if keyword_set(Fit) then begin 
   SolFit = Res
    FirstFit = 0
  ;  if FirstL LT 2 then FirstFit = 2
   Np = NFIT + 5
   DD = D - ResNoise
   DD = DD[0:Np-1] 
   l = findgen(Np) + FirstL & lfactor = l * (l+1) / 2. / !dpi 
   ; DD = DD * lfactor
   x = alog(l)
   DD = alog(DD)
   resfit = LINFIT( X, DD, meas=sqrt(DD))
   Fa = resfit(1)
   Fb =  resfit(0)
   yy = Fa*x + Fb
   ; yy = yy / lfactor
   yy = exp(yy)
   SolFit[0:Nfit-2] = yy[0:Nfit-2]
   SolFit[Nfit-1] = 0.5*(res[Nfit-1]+ yy[Nfit-1])
 
   Res = SolFit
   if FirstUseIndL GE 1 then SolFit =[InputPowSpec[0:FirstUseIndL-1], SolFit] 
end

if FirstUseIndL GE 1 then Res =[InputPowSpec[0:FirstUseIndL-1], Res] 

return, res
end
 
  ;===========================
