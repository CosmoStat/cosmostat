;+
; NAME:  
;       CMB_LOWL_ALM_INPAINTING
;
; PURPOSE:
;        Inpaint a map and return its Alm coeff.
;        Several regularization methods can be applied: Energy, sparsity or isotropy.
;        This routine is optimized for low l (l < 50) inpainting. It has been test with map of size nside=32.
;        All tests have been done using CMB maps in muK units.
;
; CALLING:
;     InpaintAlm= cmb_lowl_alm_inpaiting(Imag, Mask, Niter=Niter, lmax=lmax, Energy=Energy, 
;                                                  Sparsity=Sparsity, Isotropy=Isotropy, InpMap=InpMap, Prea=Prea, Cl=Cl) 
; INPUTS:
;     Imag -- IDL 1D array: Input Healpix image to be inpainted 
;     
; OUTPUTS:
;     InpaintAlm -- IDL 1D array: Alm coeff of the inpainting map   
;          
; INPUT KEYWORDS:
;      Niter: int: number of iterations used in the reconstruction
;      Lmax      : Number of spherical harmonics computed in the decomposition (default is 50)
;      Energy    : Use a L2 regularization
;      Sparsity   : Use a sparse regularization
;      Isotropy   : Use an isotropy regularization
;      Prea         : Estimated power spectrum of the unmasked data (used by Isotropy)
;      Cl            : Theoretical Cl  (used by Energy constraint)
;
; OUTPUT KEYWORDS:
;     InpMap: IDL Healpix map: Inpainted map. 
;    
;   EXTERNAL CALLS:
;       mrs_alm_inpainting (C++ program) for isoptropie and Iterative Hard Thresholding.
;
; EXAMPLE:
;     InpAlm =  cmb_lowl_alm_inpaiting(Map, Mask, /sparsity, InpMap=InpMap)
;     mrs_tv, InpMap
;
; HISTORY:
;       Written : Jean-Luc Starck   2012.
;-
;-----------------------------------------------------------------
;================================================================

;=====================================================

pro softcf, Re, Ima, T, ReS, ImaS
 
if N_ELEMENTS(T) NE N_ELEMENTS(Ima) then begin
SoftThresold = dblarr(N_ELEMENTS(Ima))
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

;=====================================================

function mrs_l2_inp, Dat, Mask, Clin, Niter=Niter, lmax=lmax, cld=cld, FB=FB, beta=beta, AlmRes= AlmRes, plot=plot, verb=verb, l1=l1, BandLimitFilter=BandLimitFilter, monodipol=monodipol, eps=eps, log=log, firstguess=firstguess
; FB = Forward-Backward algorithm
; Cl = Prior on the power spectrum



if not keyword_set(niter) then Niter =10
if not keyword_set(lmax) then lmax =64
if not keyword_set(eps) then eps = 1e-10

RegulParam=0
if not keyword_set(l1) then Cl = Clin[0:lmax]

Rec = Dat * 0
Res = Rec
 if keyword_set(firstguess) then Res= firstguess
Alpha=1.
if not keyword_set(Beta) then Beta = 50.
  
 mrs_almtrans, Dat, AlmDat, /tab, lmax=lmax 
 mrs_almrec, AlmDat, D1
Nv = AlmDat.NORMVAL
if not keyword_set(l1) then RegulParam = Cl *Nv^2/ (Cl *Nv^2+ Beta)
; RegulParam [*] = 1. / 1e3

AlmRes = AlmDat
AlmProj = AlmDat
Res = D1
OldAlm = AlmProj.alm*0
indMask = where(Mask Eq 1)
SigData = sigma( Dat[IndMask])
OldSigResi = 1.
OldResi = Dat * Mask
;===ITER

i=0  
repeat begin

; for i=0, Niter-1  do begin  
   Resi = (Dat - Res) * Mask
    Z = (1 - mask) * Res + Dat
   ;  print, i+1, ': SigmaResi = ', sigma(Resi)
   if keyword_set(plot) then tvs, Z, window=1, title=strc(i+1)

   mrs_almtrans, Z, AlmProj, /tab, lmax=lmax 
   SigResi = sigma(Resi[indMask]) / SigData
  DeltaResi = ABS(OldSigResi - SigResi) / SigData
   Dresi = (OldResi - Resi) / SigData
   SigDresi = sigma( Dresi[indMask])
   if keyword_set(Verb) then  print, i+1, ': SigmaResi = ', SigResi, ', SigmaDiffResi = ', SigDresi
   OldSigResi = SigResi
   OldResi = Resi
   ; AlmResi.alm = AlmDat.alm   - AlmResi.alm
  ;  AlmProj.alm = AlmRes.alm + AlmResi.alm
  ;   Res = Res + Resi
  ; mrs_almtrans, Res, AlmProj, /tab, lmax=lmax
    if i EQ 0 then begin
       AlmProx = AlmProj
       AlmRes = AlmProx
       AlmRes.alm[*]=0 
   end
  AlmRes0 = AlmRes
  
    if keyword_set(FB) then begin
          RegulParam = Cl / (Cl + Beta)
          for l=0, AlmRes.LMAX do  AlmRes.alm[l,*,*] = AlmRes.alm[l,*,*] * RegulParam[l]
          AlmRes.alm = AlmRes.alm + alpha * AlmResi.alm
    end else begin
         AlmProx.alm = 2. * AlmProj.alm - AlmRes.alm   
        ; Prox Constraint
        if not keyword_set(l1) then  for l=0, AlmProx.LMAX do  AlmProx.alm[l,*,*] = AlmProx.alm[l,*,*] * RegulParam[l]  $
        else begin
            Re = AlmProx.alm[*,*,0] 
            Ima = AlmProx.alm[*,*,1]
            softcf, Re, Ima, Beta/NV, ReS, ImaS
            AlmProx.alm[*,*,0] = ReS
            AlmProx.alm[*,*,1] = ImaS
         end ; L1 solution
        ; Update the solution
         AlmRes.alm = AlmRes.alm + alpha * (AlmProx.alm - AlmProj.alm)
      end  ; if ForwardBackward solution
     if not keyword_set(monodipol) then AlmRes.alm[0:1,*,*]  = 0
   
     if keyword_set(Cld)  and i GT 1 then begin
         pp = mrs_alm2spec(AlmRes0)
        W = pp
        W[0:1] = 0
        W[2:*] = sqrt(Cld[2:*]) / sqrt(pp[2:*])
        for l=10, AlmRes.LMAX do  $
              for l=0, AlmProj.LMAX do  AlmRes0.alm[l,*,*] = AlmRes0.alm[l,*,*] * W[l]
       AlmRes.alm = 0.5 * (AlmRes.alm + AlmRes0.alm)
    end
    
    if keyword_set(BandLimitFilter) then begin
         for l=0, AlmRes.LMAX do  $
              for l=0, AlmProj.LMAX do  AlmRes.alm[l,*,*] = AlmRes.alm[l,*,*] * BandLimitFilter[l]
    end
    
    mrs_almrec, AlmRes, Res
    
      
;  end
i = i + 1
 
endrep until i EQ Niter ;  or SigDresi LT eps

return, Res
end

;================================================================
;  c = [ 1.00071, 1.00539, 1.00129, 0.99860, 1.00000, 0.99867, 0.99109, 1.21025, 1.00000]
function cmb_lowl_alm_inpainting, Ima, Mask, lmax=lmax, niter=niter, InpMap=InpMap, Prea=Prea, Sparsity=Sparsity, Isotropy= Isotropy, Galaxies=Galaxies, Energy=Energy, IHT=IHT, pniter=pniter, noinp= noinp, MatMask=MatMask, DeconvPSMatMask=DeconvPSMatMask, ImatMask=IMatmask, EstimPowSpecCMB=EstimPowSpecCMB,  ScalePowSpec=ScalePowSpec, BandLimitFilter=BandLimitFilter, Cl=Cl, GLS=GLS, KeepData=KeepData, log=log, eps=eps, nomonodipol=nomonodipol, beta=beta, NormInputMap=NormInputMap, silent=silent
 
if not keyword_set(lmax) then lmax=50
if not keyword_set(niter) then niter=150
if not keyword_set(eps) then eps = 1e-6
if not keyword_set(beta) then beta =1.
if not keyword_set(silent) then verb =1.
if not keyword_set(pniter) then pniter =500L

monodipol=1
if keyword_set(nomonodipol) then monodipol = 0

alm=-1
IF N_PARAMS() LT 2 THEN BEGIN        
PRINT,'Error: 2 inputs needed. ...'  
GOTO, DONE 
ENDIF
  
n1 = gnside(Ima)
n2 = gnside(Mask)
if n1 NE n2 then begin
   print, 'Error: input data and input mask must have the same size ... '
   GOTO, DONE 
ENDIF

; function mrs_deconv_powspec,  PowSpecMaskedData,  Mask,  MatMask=MatMask, PRes, Niter=Niter, TrueCl=TrueCL, NoisePS=NoisePS, lmax=lmax, mu=mu, zeromonodip=zeromonodip, verb=verb, fsky=fsky, inv=inv, ImatMask=IMatmask

if keyword_set(PBGalaxies) then begin
   Deconvlmax=0
   Pd = mrs_powspec(Ima*Mask)
   ; MasterPowSpec = mrs_deconv_powspec(pd,  mask, MatMask=DeconvPSMatMask, Niter=pniter)
   MasterPowSpec = mrs_deconv_powspec(Ima,  mask, MatMask=DeconvPSMatMask, Niter=pniter, /map)
   print, "PS Quad and Oct  (LOWL) = ", MasterPowSpec[2], MasterPowSpec[3]
end

EpsLog = 1. + eps
if  keyword_set(log) then begin
    Map = alog(double(Ima)+EpsLog) 
end else Map = Ima
ind = where(Mask  GE 1.)
MeanMap = mean( Map[ind])
if not keyword_set(Galaxies) then MeanMap=0.
; MeanMap=0.
Map[ind] = Map[ind] - MeanMap

D = double(Map*Mask)
NormMap = 100. / sigma(D[ind])
if not keyword_set(NormInputMap) then NormMap=1.
D = D * NormMap

if not keyword_set(Isotropy)  and not keyword_set(Galaxies)  and not keyword_set(Energy) and not keyword_set(IHT)  and not keyword_set(noinp)  and not keyword_set(GLS)  then Sparsity=1

if  keyword_set(Prea) then EstimPowSpecCMB=Prea $
else if  keyword_set(Isotropy)  or  keyword_set(Galaxies)  or  keyword_set(Energy) or keyword_set(Sparsity) or keyword_set(PowSpecCstr) or  (keyword_set(GLS) and not  keyword_set(Cl))    then begin
      Deconvlmax=0
      ; pd = mrs_powspec(d)
      ; EstimPowSpecCMB = mrs_deconv_powspec(pd,  mask, MatMask=DeconvPSMatMask, Niter=pniter)
      EstimPowSpecCMB = mrs_deconv_powspec(d,  mask, MatMask=DeconvPSMatMask, Niter=pniter, /map)
      ; print, "PS Quad and Oct = ", EstimPowSpecCMB[2], EstimPowSpecCMB[3]
end

if keyword_set(Sparsity) then begin
      ; firstguess = mrs_alm_inpainting(d, mask=mask, /set, outA=alm, lmax=lmax, niter=niter, pniter=pniter, verb=verb)
      InpMap = mrs_l2_inp(d, Mask, EstimPowSpecCMB, nit=niter, lmax=lmax, AlmRes=alm, verb=verb, /l1, beta=beta, BandLimitFilter=BandLimitFilter, eps=eps, log=log, monodipol=monodipol, firstguess=firstguess)
      PowSpecCstr=1
     if keyword_set(PowSpecCstr) then begin
        Palm = mrs_alm2spec(alm)
        for l=2,lmax do alm.alm[l,*,*] = alm.alm[l,*,*]  / sqrt(Palm[l]) * sqrt(EstimPowSpecCMB[l])
        mrs_almrec, alm, InpMap
    end 
end


 if keyword_set(IHT) then begin
      InpMap = mrs_alm_inpainting(d, niter= niter, opt='-m1 -W  ', lmax=lmax, outA=alm, mask=mask, verb=verb)
end

if keyword_set(Isotropy) then begin
      InpMap = mrs_alm_inpainting(d, niter= niter, opt='-m4 ', lmax=lmax, outA=alm, mask=mask, /set, prea=EstimPowSpecCMB, verb=verb)
end

if keyword_set(Galaxies) then begin
       ; InpMap = mrs_alm_inpainting(d, mask=mask, /set, outA=alm, lmax=lmax, niter=niter, pniter=pniter, verb=verb, PRea=EstimPowSpecCMB)
       InpMap = mrs_alm_inpainting(d, mask=mask, /set, outA=alm, lmax=lmax, niter=niter, pniter=pniter, verb=verb)
       ; mrs_almrec, alm, InpMap
end

if keyword_set(Energy) then begin
    if keyword_set(Cl) then InpMap = mrs_l2_inp(d, Mask, Cl, nit=niter, lmax=lmax, AlmRes=alm, verb=verb) $
    else InpMap = mrs_l2_inp(d, Mask, EstimPowSpecCMB, nit=niter, lmax=lmax, AlmRes=alm, eps=eps, verb=verb)
    PowSpecCstr=1
    if keyword_set(PowSpecCstr) then begin
        Palm = mrs_alm2spec(alm)
        for l=2,lmax do alm.alm[l,*,*] = alm.alm[l,*,*]  / sqrt(Palm[l]) * sqrt(EstimPowSpecCMB[l])
        mrs_almrec, alm, InpMap
    end 
end

if keyword_set(GLS) then begin
    PowSpecCstr=0
    if keyword_set(Cl) then  InpMap = mrs_gls_inpainting(d, Mask,  nit=niter, lmax=lmax, AlmRes=alm, verb=verb, cl=cl) $
   else InpMap = mrs_gls_inpainting(d, Mask,  nit=niter, lmax=lmax, AlmRes=alm, verb=verb, cl=EstimPowSpecCMB)
end

InpMap = InpMap /  NormMap
alm.alm = alm.alm / NormMap
InpMap = InpMap + MeanMap

if  keyword_set(log) then begin
    InpMap = exp(InpMap) - EpsLog 
end
if keyword_set(noInp) then InpMap = d*Mask

if keyword_set(KeepData) then  InpMap = InpMap*(1-mask) + mask*d
 
mrs_almtrans, InpMap, alm, lmax=lmax, /tab

if keyword_set(PBGalaxies) then begin
     PowSpecCstr=1
     EstimPowSpecCMB = MasterPowSpec
     if keyword_set(PowSpecCstr) then begin
        Palm = mrs_alm2spec(alm)
        for l=2,lmax do alm.alm[l,*,*] = alm.alm[l,*,*]  / sqrt(Palm[l]) * sqrt(EstimPowSpecCMB[l])
        mrs_almrec, alm, InpMap
    end 
end


if keyword_set(noinp) then begin
     fsky = float(N_ELEMENTS(mask)) / float(N_ELEMENTS(d))
     for l=0,lmax do alm.alm[l,*,*] = alm.alm[l,*,*]  / sqrt(fsky)
     mrs_almrec, alm, InpMap
end

DONE:

return, alm
end
;================================================================
 
