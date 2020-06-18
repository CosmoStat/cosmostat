;+
; NAME: 
;       MRS_ISO_INPAINTING
;
; PURPOSE:
;        Apply an inpainting to a spherical map using an isotropic constraint on the spherical harmonics coefficients.
;        If a mask is not provided, all pixels with a zero value are considered as missing pixels.
;        A c++ program ($MRS/cxx/mrs_alm_inpainting) is called.
;        If file names are given for the input image and mask, these two images are not loaded into IDL.
;
; CALLING:
;     InpaintMap = mrs_alm_inpainting(Imag, Mask=Mask, Niter=Niter, OutPowSpec=OutPowSpec, lmax=lmax, gauss=gauss)
;
; INPUTS:
;     Imag -- IDL 1D array: Input Healpix image to be inpainted 
;     
; OUTPUTS:
;     InpaintMap -- IDL 1D array: Output inpainted Healpix image   
;          
; INPUT KEYWORDS:
;      niter: int: number of iterations used in the reconstruction
;      Lmax      : Number of spherical harmonics computed in the decomposition
;                                       (default is 3*nside, should be between 2*nside and 4*nside)
;      FNin: Filename containing the input image. The input image won't be read.
;      FNMask: Filename containing the mask. The mask won't be read.
;      FNOut: Filename containing the results. By default, nothing is written on the disk.
;      FNPS: Filename containing the mask. The mask won't be read.
;
; OUTPUT KEYWORDS:
;     OutPowSpec: IDL 1D array: Cl of the inpainted map. 
;       
; EXTERNAL CALLS:
;       mrs_alm_inpainting (C++ program)
;
; EXAMPLE:
;      
; HISTORY:
;       Written : Jean-Luc Starck   2009.
;-
;-----------------------------------------------------------------


;==============================================================================================

;=====================================================

function xerf,   X
 
   ;  double Z,A0,A1,B1,B2,Val_Return;
    
    Z = X
    if (X GT 0.5) then  Z = 1 - X
    if (Z LT 1e-12) then  Val_Return  = 0. $
     else begin
          Z = sqrt (-2. * alog (Z))
         A0 = 2.30753
         A1 = 0.27061
         B1 = 0.99229
         B2 = 0.04481
         Val_Return = Z - (A0 + A1 * Z) / (1 + ((B1 + B2 * Z) * Z))
         if (X GT 0.5) then Val_Return = - Val_Return
    end
    return, (Val_Return)
end

;=====================================================

 function xerfc1, F
     Nu = 0.
      P = F
   if (P GT 1.) then P = 1. $
   else if (P LT  0) then P = 0.
   if (P GT 0.5) then Nu = sqrt(2.)*inverfc( (1.-P)/0.5) $
   else Nu = - sqrt(2.)*xerfc( P/0.5)
   return, Nu
end

;=====================================================

function xerf1, y 
return,  1. - xerf(y)
end
;=====================================================

function phi, x, mu=u, sigma=sigma
if not keyword_set(mu) then mu = 0.
if not keyword_set(sigma) then sigma =1.

return, 1. / (double(2)) * (1. + errorf( (x-mu) / (sqrt(2.)*sigma)))
end
;=====================================================

function invphi, alpha
return, xerf(1 - alpha)
end

;=====================================================

function eps2nsig, alpha
InvPhi = invphi(double(1) - alpha/2.)
return, InvPhi
end

;=====================================================

function nsig2eps, InvPhi
P = 2. * ( 1. - phi(InvPhi))
return, P
end
;=====================================================

function eps_gaussian_constraint, TabMperL, alpha, minus=minus
TabM = TabMperL * 2 - 1
vs = size(TabM)
Nm = vs[1]
InvPhi = invphi(double(1) - alpha/2.)
if not keyword_set(minus) then Eps = sqrt(double(3)/4.) - sqrt(1. - 1. / (4*TabM))  + sqrt( (1. + 1. /TabM)/ 8.)  * InvPhi $
else  Eps = sqrt(double(3)/4.) - sqrt(1. - 1. / (4*TabM))  - sqrt( (1. + 1. /TabM)/ 8.)  * InvPhi
return, eps
end

;=====================================================

pro softcf, Re, Ima, T, ReS, ImaS
if N_ELEMENTS(T) NE N_ELEMENTS(Ima) then begin
SoftThresold = dblarray(N_ELEMENTS(Ima))
SoftThresold[*] = T[0]
end else SoftThresold = T

ReS = Ima
ReS[*] = 0
ImaS = ReS

M = sqrt(Re^2 + Ima^2)
ind = where ( M EQ 0., c)

if (c GT 0) then begin
   ReS[ind] = 0
   ImaS[ind] = 0
END 

ind = where ( M GT 0., c)

if (c GT 0) then begin
  Z = 1. - SoftThresold[ind] / M[ind]
  indZ0 = where (Z lt 0, cz)
  if cz GT 0 then Z[indZ0] = 0.
  ReS[ind]  = Z * Re[ind] 
  ImaS[ind] = Z * Ima[ind] 
END
end


;==============================================================================================


pro deconv_ps, PDat, Mask, MatMask, PRes, niter=niter, TrueCl=TrueCLm, NoisePS=NoisePS

if not keyword_set(niter) then niter=20
if not keyword_set(NoisePS) then NoisePS=0.

Fsky = float( N_ELEMENTS(Mask) / float ( total(mask)) )

; Firt initialisation: Fsky correction and noise subtraction
PRes = PDat*Fsky - NoisePS

Tmat = transpose(MatMask)

for i=0,Niter-1 do begin
   resi = PDat - MatMask # ( PRes + NoisePS )
   PRes += Tmat # resi
   
   ; Positivity constraint
   ind = where(Pres LT 0, c)
   if c GT 0 then Pres[ind] = 0
   
   if keyword_set(TrueCl) then begin
      plotcl, TrueCl, thick=2
      oplotcl, Pres
   end
end
end

;================================================================

pro isotropy_cst, S_ALM,  MuPowSpec=MuPowSpec, OutPowSpec=OutPowSpec, alpha=alpha, Threshold=Threshold, pscst=pscst

Ro = 1.  ; entre 0 et 2 (strictement)
TM = S_ALM.TABNBRM
if not keyword_set(Alpha) then Alpha = 0.05
;Alpha = 0.05   ; seuil a 2sigma
;Alpha = 0.0001   ; seuil a 2sigma
              ; 
if not keyword_set(threshold) then begin
  EpsP = eps_gaussian_constraint(tm,Alpha)
  EpsM = eps_gaussian_constraint(tm, Alpha, /min)
end else begin
  EpsP = TM
  EpsM = EpsP
  EpsP[*] = Threshold
  EpsM[*] = Threshold
end
; print, 'Alpha = ', Alpha

MedianPS=0
if not keyword_set(MuPowSpec) then P = mrs_alm2spec(s_alm, StdPS=StdPS) $
else P = MuPowSpec
Cl = P

TabL = float( 2*s_alm.TABNBRM-1)
; StdPS = StdPS / sqrt(TabL)
NormL = sqrt(1.+1/TabL) 
StdPS = sqrt(P) * NormL 
  
LP = fltarr(S_ALM.lmax+1, S_ALM.lmax+1)
LM = LP

SA = s_alm
Zero = reform( S_ALM.alm[0,0,*] )

A = s_alm.alm

; Isotropy constraint  
P = mrs_alm2spec(s_alm, StdPS=StdPS1)
DeltaP = sqrt(P) - sqrt(Cl)
for l=0, S_ALM.LMAX do BEGIN
      Nm = s_alm.TABNBRM[l]
      Mu = sqrt( CL[l] )
      m2n = sqrt(A[l,*,0]^2 + A[l,*,1]^2)
      DeltaMu = M2n - Mu
      indP = where(DeltaMu GE 0, cp)
      indM  = where(DeltaMu LT 0, cm)
      Z1 = dblarr( N_ELEMENTS(m2n) )
      Z2 = Z1
      T = Z1
      if cp GT 0 then begin
          Z1 =  DeltaMu - EpsP[l] * StdPS[l]
          ind = where(Z1 LT 0, c)
          if c gt 0 then Z1[ind] = 0.
          T[indP] = Z1[indP]
      end
      if cm GT 0 then begin
          Z2 =  DeltaMu  + EpsM[l]  * StdPS[l]
          ind = where(Z2 GT 0, c)
          if c gt 0 then Z2[ind] = 0.
          T[indM] = Z2[indM]
      end
                  
      ValRe = s_alm.alm[l,*,0]
      ValIma = s_alm.alm[l,*,1]
      softcf, ValRe, ValIma, T, ReS, ImaS
      A[l,*,0] = ReS
      A[l,*,1] = ImaS
ENDFOR ; l
S_ALM.alm = A




; Cl constraint
if keyword_set(pscst) then BEGIN
P = mrs_alm2spec(s_alm, StdPS=StdPS1)
DeltaP = sqrt(P) - sqrt(Cl)
SigMu = StdPS1 / sqrt(TabL)
Z1 = DeltaP -  EpsP * SigMu
Z2 = DeltaP +  EpsP * SigMu

Ind = where (Z1 LE 0, c)
if c GT 0 then Z1[ind] = 0
Ind = where (Z2 GT 0, c)
if c GT 0 then Z2[ind] = 0
for l=0, S_ALM.LMAX do BEGIN
      if DeltaP[l] GT 0 then s_alm.alm[l,*,*] = s_alm.alm[l,*,*] * sqrt(Z1[l] / P[l])  $
      else  s_alm.alm[l,*,*] = s_alm.alm[l,*,*] * sqrt(Z2[l] / P[l])  
      ENDFOR ; l
S_ALM.alm = A
END

 
OutPowSpec = mrs_alm2spec(s_alm)
 
END

;=====================================================

function ihaar1d, Haar, NScale=Nscale

Low = Haar
vs = size(Low)
np = vs[1]
Signal = fltarr(np)
c = 1. / sqrt( double(2))
MaxScale =  long(alog(float(np)) / alog(2.)) + 1
if not keyword_set(Nscale) then Nscale = MaxScale

TabNp = lonarr(Nscale)
TabResol = TabNp
Nx = Np
for j=0,NScale-2 do begin
  TabResol[j] = Nx
  Nhi = nx  / 2
  TabNp[j] = Nhi
  Nlo = (nx + 1) / 2
  nx = Nlo
end
TabNp[NScale-1] = Nlo
TabResol[NScale-1] = Nlo
; print, TabResol
; print, TabNp

Nx = TabResol[NScale-1]
Low = Haar[0:Nx-1]
for j=NScale-2,0,-1 do begin
  LastLow = Low[nx-1]
  Nx = TabResol[j]
  Nhi = Nx  / 2
  Nlo = (Nx + 1) / 2
  ie=indgen(nhi)*2
  io = ie + 1
  Det = Haar[Nlo:Nlo+TabNp[j]-1]
  EvenCoef = c * (Low[0:nhi-1] + Det[0:nhi-1])
  OddCoef = c * (Low[0:nhi-1] - Det[0:nhi-1])
  Low = fltarr(Nx)
  Low[ie] = EvenCoef
  Low[io] = OddCoef
  if Nx mod 2 EQ 1 then begin
     ; Low[Nx-1] = 1. / c * LastLow  - Low[0] 
     Low[Nx-1] = c * (LastLow + Det[0])
  end
end   
; plot, Low

return, Low
end


;=====================================================

function haar1d, Signal, NScale=Nscale, TabNp=TabNp, TapPos=TapPos

High = Signal
vs = size(High)
np = vs[1]
Haar = fltarr(np)
c = 1. / sqrt( double(2))
MaxScale =  long(alog(float(np)) / alog(2.)) + 1
if not keyword_set(Nscale) then Nscale = MaxScale

; print, np, ' ' , Nscale
TabNp = lonarr(Nscale)
TabPos = TabNp

for j=0,NScale-2 do begin
  vs = size(High)
  np = vs[1]
  H2 = shift(High, -1)
  Nhi = np  / 2
  TabNp[j] = Nhi
  Nlo = (np + 1) / 2
  TabPos[j] = Nlo
  ie=indgen(nlo)*2
  Low = c*(High[ie] + H2[ie])
  Det = c*(High[ie] - H2[ie])
  Haar[0:Nlo-1] = Low[0:Nlo-1]
  Haar[Nlo:Np-1] =  Det[0:Nhi-1]
  High = Haar[0:Nlo-1]
end   
TabNp[NScale-1] = Nlo
TabPos[NScale-1]= 0

;print, TabNp
;print, TabPos
;plot, Haar

return, Haar
end

;==============================================================================================

pro tth,n , h, r
n = randomn(seed, 301)
; n = rim('$HOME/Main/img/s1d.fits')
;n[*] = 0
;n[17] = 2
plot, n 
h = haar1d(n)
r = ihaar1d(h)
plot, n-r
info, n-r
end

;==============================================================================================

pro mrs_haar1d_alm, ALM, Nscale=Nscale
A = ALM.alm
for l=5, ALM.LMAX do BEGIN
      Nm = alm.TABNBRM[l]
      Al_re = reform(A[l,0:Nm-1,0])
      Al_im = reform(A[l,0:Nm-1,1])
      h_re =  haar1d(Al_re, Nscale=Nscale)
      h_im = haar1d(Al_im, Nscale=Nscale)
      A[l,0:Nm-1,0] = h_re
      A[l,0:Nm-1,1] = h_im
end
for l=0, ALM.LMAX-5 do BEGIN
      Nm = alm.TABNBRM[l]
      Al_re = reform(A[l:*,l,0])
      Al_im = reform(A[l:*,l,1])
      h_re =  haar1d(Al_re, Nscale=Nscale)
      h_im = haar1d(Al_im, Nscale=Nscale)
      A[l:*,l,0] = h_re
      A[l:*,l,1] = h_im
end
ALM.alm = A
end

pro mrs_ihaar1d_alm, ALM, Nscale=Nscale
A = ALM.alm
for l=0, ALM.LMAX-5 do BEGIN
      Nm = alm.TABNBRM[l]
      Al_re = reform(A[l:*,l,0])
      Al_im = reform(A[l:*,l,1])
      h_re =  ihaar1d(Al_re, Nscale=Nscale)
      h_im = ihaar1d(Al_im, Nscale=Nscale)
      A[l:*,l,0] = h_re
      A[l:*,l,1] = h_im
end
for l=5, ALM.LMAX do BEGIN
      Nm = alm.TABNBRM[l]
      Al_re = reform(A[l,0:Nm-1,0])
      Al_im = reform(A[l,0:Nm-1,1])
      h_re =  ihaar1d(Al_re, Nscale=Nscale)
      h_im = ihaar1d(Al_im, Nscale=Nscale)
      A[l,0:Nm-1,0] = h_re
      A[l,0:Nm-1,1] = h_im
end
ALM.alm = A
end

;==============================================================================================

function normalm, Alm, pow=pow, map=map
if keyword_set(map) then mrs_almtrans, map, Alm, /tab

Cl = mrs_alm2spec(Alm)
Coef = sqrt(Cl[2:*])
ind = where(Coef EQ 0, c)
if c GT 0 then Coef[ind] = 1.
Coef = 1. / Coef
TabAlm = ALM.alm
for l=2, ALM.LMAX do TabAlm[l,*,*] = TabAlm[l,*,*] * Coef[l-2]

if keyword_set(pow) then ima =  TabAlm[*,*,0]^2 +  TabAlm[*,*,1]^2 $
else begin
  ima = tabalm[*,*,0]
  for l=1, ALM.LMAX do ima[l, l:*] = TabAlm[ALM.LMAX-l, 0:(ALM.LMAX-l), 1]
end
ima[0,0]=0

return, ima
end

;==============================================================================================

function idl_mrs_iso_inpainting, dat, Mask=Mask, MatMask=MatMask, niter=niter, FNImag=FNImag, FNMask=FNMask, FNOut=FNOut, Cl=Cl, OutPowSpec=OutPowSpec, plot=plot, lmax=lmax, GKsig= GKsig, gauss=gauss, verb=verb, Haar=Haar, PRea=Prea
COMMON C_PLANCK

if keyword_set(GKsig) then Gauss_Cst_KSigma = GKsig $
else Gauss_Cst_KSigma = 4.

if not keyword_set(lmax) then begin
   nside = gnside(dat)
   if nside NE 2048 then lmax = long(nside * 3.)  $
   else lmax = P_LMAX
end
print, 'Lmax = ', lmax
ll = findgen(lmax)
ll =2*ll + 1

if not keyword_set(Niter) then Niter=10

if not keyword_set(Mask) then begin
   Mask = Dat 
   ind = where( Dat EQ 0, c)
   Mask[*] = 1
   if c GT 0 then Mask[ind] = 0
end

print, "Percentage of valid pixels = ", float(total(Mask)) / float(N_ELEMENTS(Mask)) * 100.

indm = where(mask NE 0)
Rec = dat
mrs_almtrans, Rec, ALM, lmax=lmax, /tab

if keyword_set(plot) then winbs, win=1

; Compute the true Power spectrum  from the mask matrix and the data power spectrum
if not keyword_set(Cl) then begin
   pm = mrs_powspec(mask, /nonorm, lmax=lmax)
   if not keyword_set(MatMask) then  MatMask  = mrs_matmask(Mask, lmax=lmax)
   pdat = mrs_alm2spec(ALM)
   deconv_ps, PDat, Mask, MatMask, PRes, niter=niter, TrueCl=TrueCL
   Cl = Pres
end else Cl = Cl[0:lmax]

indm = where(mask NE 0)
indz = where(mask EQ 0)

if keyword_set(plot) then begin
   winbs, win=1
   ; tvs, dat, window=2, tit='dat'
end


Eps = 0.5

mm = 1. - mask

; This parameter defines the solution space in the wavelet packet domain
; and in the Spherical harmonic domain:
;    1) a coeff w must verify:  | w | < Lambda * sigma
;    2)  | Alm^2 - Cl | < Lambda * sqrt(Cl)
; the constraint is relaxed with the iteration, lambda = 0.5sigma -> 5sigma

LambdaFirst = 0.5
StepLambda = 4.5  / float(Niter)
Lambda = LambdaFirst

; Main loop
for i=0, niter-1  do begin 
   OldRec = Rec
   if i GT 0 then mrs_almtrans, rec, ALM, lmax=lmax, /tab
   NsH = 3
   if keyword_set(Haar)  then  mrs_haar1d_alm, ALM, Nscale=NsH
   
   P = mrs_alm2spec(ALM)
   Alpha = nsig2eps(Lambda)
   
   if keyword_set(plot) then begin
      tit = 'Iter ' + strc(i+1)
      tvs, rec, title=tit, window=2
      wset, 1    
      plotcl, cl, line=2
      oplotcl, P
   end

 ; z = reform( ALM.alm[517,0:517,*])
 ; wset, 0
 ; plot, z, title=' Bef const '
 ; wait, 5
  ; Isotropy constraint in the Alm domain
  ThresholdLevel=Lambda
  isotropy_cst, ALM, MuPowSpec=Cl, OutPowSpec=ops, threshold=ThresholdLevel, pscst=pscst
 ; z = reform( ALM.alm[517,0:517,*])
 ; wset, 0
 ; plot, z, title=' After const '
 ; wait, 5
   
   if keyword_set(Haar)  then  mrs_ihaar1d_alm, ALM, Nscale=NsH
 ;   z = reform( ALM.alm[517,0:517,*])
 ; wset, 0
 ; plot, z, title=' After Ihaar '
 ; wait, 5


  if keyword_set(Prea) then begin 
       PSol = mrs_alm2spec(ALM)
       ind = where(Psol EQ 0, c)
       if c GT 0 then Psol[ind] = 1
       ColeFilter = sqrt( Prea / Psol)
       ; info, ColeFilter
      ; stop
       ;ind = where(ColeFilter le 0, c)
       ;if c GT 0 then ColeFilter[ind] = 1
       ;ColeFilter[0] = 1
       if c GT 0 then ColeFilter[ind] = 1
       for l=0,alm.lmax do  alm.alm[l,*,*] = alm.alm[l,*,*] * ColeFilter[l]
   end
       
    
   
  mrs_almrec, ALM, Rec
    
   ; We keep only the reconstruction in the masked area      
   Rec = dat * Mask + Rec * mm
   
 
   
   ; Constraint in the Wavelet packet domain
   waveWP_Cst=0
 
   if keyword_set(waveWP_Cst) then begin
   
       mrs_wptrans, Rec, W, nbrscale=nbrscale
        TabPsi2 = fltarr(W.NbrScale)
      for j=0,W.NbrScale-2 do BEGIN
          Scale = W.coef[*,j]
          
          ; Compute the expected power in the band
          PClj = Cl * W.TabPsi[*,j]^2
          SigBand = sqrt(total(PClj*ll) /  (4. * !DPI))
          TabPsi2(j) = SigBand
          ; projection. All  coeff > Gauss_Cst_KSigma*SigBand are set to zero
          ind = where(ABS(Scale*mm) GT Gauss_Cst_KSigma * SigBand, c)
          if c GT 0 then Scale[ind] = 0
          
          ; Gaussianity constraint
          if keyword_set(Gauss) then begin
             SigMask = sigma(Scale[Indm])
             DeltaMu = ABS(Scale) - SigMask
             SoftT = DeltaMu - Lambda *  SigBand
             ind = where(SoftT LT 0, c)
             if c gt 0 then SoftT[ind] = 0.
             SoftT[Indm] = 0.
             softthreshold, Scale, SoftT
          end
          
          ; band renormalisation
          Sigm =  sigma(Scale*Mask)
          Sigmm = sigma(Scale*mm)
          CoefMask =  (SigBand^2. - Sigm^2) / Sigmm^2
          if CoefMask GT 0 then CoefMask = sqrt(CoefMask) else CoefMask = 1.
          if keyword_set(Verbose) then print, '    Band ', j+1, ':  CoefM = ', CoefMask
          Scale[indz] = Scale[indz] * CoefMask
          
          W.coef[*,j] = Scale
       END
      mrs_wtrec, W, Rec
      writefits, 'xx_psi2.fits', TabPsi2
  END
  
   ; We keep only the reconstruction in the masked area      
   rec = dat * Mask + Rec * mm
  
   Err = sigma(oldrec -rec)
   print, ' Iter ', i+1, ', Lambda = ',  Lambda, ', DeltaResi = ',  Err 

   ; relax the constraint
   Lambda = Lambda + StepLambda
end


 mrs_almtrans, rec, ALM, lmax=lmax, /tab
 OutPowSpec = mrs_alm2spec(ALM)
 if keyword_set(plot) then begin
      tit = 'Rec  '  
      tvs, rec, title=tit, window=2
      wset, 1    
      plotcl, cl, line=2
      oplotcl, OutPowSpec
  end
  
 return, Rec

end

;==============================================================================================
;==============================================================================================


function mrs_iso_inpainting,  Ima, Mask=Mask, Niter=Niter,  OutPowSpec=OutPowSpec,  lmax=lmax,  Opt=Opt,  FNMask=FNMask, FNOut=FNOut, FNIn=FNin, FNPS=FNPS, gauss=gauss, verb=verb, idl=idl, MatMask=MatMask, plot=plot, Haar=Haar, Cl=CL, Prea=Prea, Log=Log, PNiter=PNiter

if N_PARAMS() LT 1 then begin 
        print, 'CALL SEQUENCE: InpaintMap = mrs_alm_inpainting(Imag, Mask=Mask, Niter=Niter,  OutPowSpec=OutPowSpec,  lmax=lmax,  FNMask=FNMask, FNOut=FNOut, FNIn=FNin, FNPS=FNPS, gauss=gauss'
        goto, DONE
        end

if not keyword_set(Mask)  then BEGIN
     if keyword_set(FNIn) then Imag = mrs_read(FNIn)
      Mask = Imag
      Mask[*] = 0
      ind = where(Imag NE 0, c)
      if c GT 0 then Mask[ind] = 1
END

EpsLog = 1.
if  keyword_set(log) then Imag = alog(Ima+EpsLog) $
else Imag = Ima


if not keyword_set(PRea) then begin
      Pi = mrs_powspec(Imag)
      P = mrs_deconv_powspec( Pi,  Mask,  MatMask=MatMask,  Niter= PNiter,  NoisePS=NoisePS, lmax=lmax)
end else P = Prea
if not keyword_set(Cl) then  Cli = mrs_get_cl_theo_powspec(P) else Cli= Cl

if not keyword_set(Opt) then Opt = '  '
OptI = '-m4 ' + Opt
if keyword_set(log) then OptI = OptI + ' -C ' +  '  '

NameCl =gettmpfilename()
writefits, NameCl, double(Cli)
OptI = OptI + ' -S '  +  NameCl + '  '
 
; IDL = 1 
if keyword_set(IDL) then begin 
     InpaintMap = idl_mrs_iso_inpainting(Imag, Mask=Mask,  niter=niter, FNImag=FNin, FNMask=FNMask, FNOut=FNOut, OutPowSpec=OutPowSpec, MatMask=MatMask, gauss=gauss, plot=plot, Haar=Haar, verb=verb, Cl=CLe, Prea=P)
end else BEGIN
 
 ind = where(mask NE 1, c)
 if c GT 0 then mask[ind] = 0
 indMask = where( mask NE 0, cMask)
 MeanVal = mean(Imag[indMask])


if keyword_set(Niter) then OptI = OptI + ' -i ' + strc(Niter) +  '  '
if keyword_set(lmax) then OptI = OptI + ' -l ' + strc(lmax) +  '  '
if keyword_set(gauss) then OptI = OptI + ' -G '  +  '  '
if keyword_set(Verb) then OptI = OptI + ' -v '  +  '  '

if not keyword_set(FNin)  then NameIma =gettmpfilename()   else NameIma = FNin
if not keyword_set(FNMask)  then NameMask =gettmpfilename()  else NameMask = FNMask
if not keyword_set(FNOut)  then NameResult = gettmpfilename()  else NameResult = FNOut
if not keyword_set(FNPS)  then NameResultPS = gettmpfilename()   else NameResultPS = FNPS

if not keyword_set(FNin)  then mrs_write,  NameIma,  Imag
if not keyword_set(FNMask)  then mrs_write,  NameMask,  Mask

com = '$MRS/cxx/mrs_alm_inpainting '+' '+ OptI + ' '+ NameIma  + '  ' +  NameMask  + '  ' +NameResult + '  ' +  NameResultPS   
print, com
spawn, com

InpaintMap = mrs_read(NameResult)
OutPowSpec = mrdfits(NameResultPS, /silent)

if not keyword_set(FNin)  then delete, NameIma
if not keyword_set(FNMask)  then delete, NameMask
if not keyword_set(FNOut)  then delete, NameResult
if not keyword_set(FNPS)  then delete, NameResultPS

end

if keyword_set(FNOut) then mrs_write, FNOut, InpaintMap

if  keyword_set(log) then InpaintMap = exp(InpaintMap) - EpsLog 

return, InpaintMap

 DONE:
 
end

;==============================================================================================

pro run_oct11
m = rims("$ISAP"+'/data/ISW_DATA/nvss/mask_cmb_n0032_ring.fits')
cmb = getcmb(cl=cl, nside=32)
i = cmb * m
Pi = mrs_powspec(i)
P = mrs_deconv_powspec( Pi,  M,  MatMask=MatMask,  Niter=Niter,  NoisePS=NoisePS, lmax=lmax)
Cle = mrs_get_cl_theo_powspec(P)

z = idl_mrs_iso_inpainting(i, Mask=m, MatMask=MatMask, niter=10, PRea=P, /plot)

 z = mrs_iso_inpainting(i, Mask=m, Niter=Niter,  Cl=CL1, Prea=Prea1)

end


;==============================================================================================

pro runt, dat, Mask=Mask, MatMask=MatMask, z1, z2, pz1, pz2, z3, z4, pz3, pz4
z1 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pz1,/idl, /plot, /verb)
z2 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pz2,/idl, /plot, /verb, /gauss)
z3 = mrs_alm_inpainting(dat, Mask=Mask, OutPowSpec=pz3)
z4 = mrs_alm_inpainting(dat, Mask=Mask, OutPowSpec=pz4, /gauss)
z5 = mrs_alm_inpainting(dat, Mask=Mask, OutPowSpec=pz5, /gauss)

r1 =  mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pz1,/idl, /plot, /verb, /haar)
save, filename='xxz.xdr', z1, z2, z3, z4, pz1, pz2, pz3, pz4
end

;==============================================================================================

pro run_inp_may09

mask = mrs_read('../mask_apcch2.fits')
;MatMask = mrs_matmask(mask)
;writefits, 'matmask.fits', MatMask
MatMask = readfits('matmask.fits')


OK=0
if OK EQ 1 then begin

FN = 'lensed_map-11187_Lobb_hfi6chan-29532.fits'
dat = mrs_read(FN)
dat = dat * mask

FN_Out='V4_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l1 = mrs_alm_inpainting(Dat, Mask=Mask, OutPowSpec=pl1, FNOut=FN_Out, FNPS=FNPS, /gauss)
wl1 = mrs_whitening(l1)
save, filename='l1.xdr', l1, pl1, wl1

FN_Out='V5_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l2 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pl2, FNOut=FN_Out, /idl, niter=5, /gauss)
wl2 = mrs_whitening(l2)
save, filename='l2.xdr', l2, pl2, wl2

FN_Out='V6_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l3 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pl3,FNOut=FN_Out, /idl, niter=10, /gauss)
wl3 = mrs_whitening(l3)
save, filename='l3.xdr', l3, pl3, wl3

FN_Out='V7_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l4 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pl4,FNOut=FN_Out, /idl, niter=20, /gauss)
wl4 = mrs_whitening(l4)
save, filename='l4.xdr', l4, pl4, wl4

FN_Out='V8_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l5 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pl5,FNOut=FN_Out, /idl, niter=40, /gauss)
wl5 = mrs_whitening(l5)
save, filename='l5.xdr', l5, pl5, wl5

FN = 'unlensed_map-11187_Lobb_hfi6chan-29532.fits'
dat = mrs_read(FN)
dat = dat * mask

FN_Out='V4_unlensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u1 = mrs_alm_inpainting(Dat, Mask=Mask, OutPowSpec=pu1, FNOut=FN_Out, FNPS=FNPS, /gauss)
wu1 = mrs_whitening(u1)
save, filename='u1.xdr', u1, pu1, wu1

FN_Out='V5_unlensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u2 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pu2, FNOut=FN_Out, /idl, niter=5, /gauss)
wu2 = mrs_whitening(u2)
save, filename='u2.xdr', u2, pu2, wu2

FN_Out='V6_unlensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u3 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pu3,FNOut=FN_Out, /idl, niter=10, /gauss)
wu3 = mrs_whitening(u3)
save, filename='u3.xdr', u3, pu3, wu3

FN = 'unlensed_map-11187_Lobb_hfi6chan-29532.fits'
dat = mrs_read(FN)
dat = dat * mask

print, "OK"
FN_Out='V7_unlensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u4 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pu4,FNOut=FN_Out, /idl, niter=20, /gauss)
wu4 = mrs_whitening(u4)
save, filename='u4.xdr', u4, pu4, wu4

FN_Out='V8_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u5 = mrs_iso_inpainting(dat, Mask=Mask, MatMask=MatMask, OutPowSpec=pu5,FNOut=FN_Out, /idl, niter=40, /gauss)
wu5 = mrs_whitening(u5)
save, filename='u5.xdr', u5, pu5, wu5



FN = 'lensed_map-11187_Lobb_hfi6chan-29532.fits'
dat = mrs_read(FN)
dat = dat * mask

FN_Out='V10_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l10 = mrs_alm_inpainting(Dat, Mask=Mask, OutPowSpec=pl10, FNOut=FN_Out, FNPS=FNPS, opt='  -v -W ')
wl10 = mrs_whitening(l10)
save, filename='l10.xdr', l10, pl10, wl10

FN_Out='V11_lensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
l11 = mrs_alm_inpainting(Dat, Mask=Mask, OutPowSpec=pl11, FNOut=FN_Out, FNPS=FNPS, opt='  -v  -W -B  ')
wl11 = mrs_whitening(l11)
save, filename='l11.xdr', l11, pl11, wl11


FN = 'unlensed_map-11187_Lobb_hfi6chan-29532.fits'
dat = mrs_read(FN)
dat = dat * mask

FN_Out='V10_unlensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u10 = mrs_alm_inpainting(Dat, Mask=Mask, OutPowSpec=pu10, FNOut=FN_Out, FNPS=FNPS, opt=' -v  -W ')
wu10 = mrs_whitening(u10)
save, filename='u10.xdr', u10, pu10, wu10

FN_Out='V11_unlensed_map-11187_Lobb_hfi6chan-29532_inp.fits'
u11 = mrs_alm_inpainting(Dat, Mask=Mask, OutPowSpec=pu11, FNOut=FN_Out, FNPS=FNPS, opt='  -v  -W -B ')
wu11 = mrs_whitening(u11)
save, filename='u11.xdr', u11, pu11, wu11
end 

;else begin
;end

end



pro run_oct09, ps, ips, t10e, t10,t30, t30e, tdct, m 
ps = rims('PSOnCmb_IR_milliKthermo_from_MJy_sr_LP_UnionMask.fits')  * 1.e3
; pps = mrs_powspec(ps)
ind = where(ps EQ 0)
m = ps
m[*] = 1
m[ind] = 0
ips = mrs_alm_inpainting(ps, Mask=m, OutPowSpec=pips, FNOut=FN_Out, FNPS=FNPS, opt='  -v  -W -B ')

dips = mrs_dct_inpainting(ps, m) 

mrs_mca, ps, Mask=m, t10e, niter=10, selecttrans=[2,4], /CstSigma, /expo
mrs_mca, ps, Mask=m, t10, niter=10, selecttrans=[2,4], /CstSigma

mrs_mca, ps, Mask=m, t30e, niter=30, selecttrans=[2,4], /CstSigma, /expo
mrs_mca, ps, Mask=m, t30, niter=30, selecttrans=[2,4], /CstSigma


mrs_mca, ps, Mask=m, tdct, niter=40, selecttrans=[7], /CstSigma, /expo

end



