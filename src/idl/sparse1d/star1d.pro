;+
;
; NAME:
;        STAR1D
;
; PURPOSE:
;       Computes the starlet transform of a signal. 
;       The output is a 2D IDL array. Reconstruction can be done by simple 
;        co-addition of all frames: Rec = total(TransSignal, 2)
;       
; CALLING:
;
;      TransSignal = STAR1D(Signal, Nscale=Nscale) 
;       
; INPUTS:
;     Signal -- 1D IDL array: signal we want to transform
;
; OUTPUTS:
;     TransSignal -- 2D IDL array: Wavelet Transform
;                                 TransSignal(i, *) = ith band of the 
;                                                    wavelet transform
;  
; KEYWORDS:
;      Nscale -- int: Number of scales. Default is 4.
;                     TransSignal is an array of Nscale signals: 
;                                      TransSignal(0: Nscale-1, *)
;			There is no test of the validity of the number of scales with respect
;			to the size of the input images. 
; 
; EXAMPLE:
;
;       Compute the multiresolution of the signal I with default options
;       (i.e. a starlet transform with 4 scales).  
;               Output = stard1d_trans(I)
;
; REFERENCE:
;    J.L. Starck and F. Murtagh, 
;    "Image Restoration with Noise Suppression 
;    Using the Wavelet Transform",
;    Astronomy and Astrophysics, 288, pp-343-348, 1994.
;
; AUTHOR:
;    J.L. Starck
;    Service d'Astrophysique, Centre d'Etudes de SACLAY,
;    Orme des Merisiers, 91191 GIF-Sur-YVETTE CEDEX, France 
;    Email: jstarck@cea.fr        Tel:  +33 (0) 169 085 764
;    http://jstarck.free.fr       Fax:  +33 (0) 169 086 577
;-
;-------------------------------------------------------------------------------

function test_ind_1D, ind, N
	;this function implements mirror like limit conditions on the edges of the input signal
	;ATTENTION : the output may still be out of range ie not in [0, N-1]
	;refinements may be necessary, although not meaningful in practice
	
	ret = ind
	if ind LT 0 then ret = -ind $
	else if ind GE N then ret = 2*N-ind-2
	
	return, ret
	
end

;-------------------------------------------------------------------------------

pro b3spline_star1d, Im_in, Im_Out,  Step
; /edge_wra
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   KSize = 4*Step+1
   KS2 = KSize/2
   Kernel = fltarr(KSize)
   Kernel[0] = C1
   Kernel[KSize-1] = C1
   Kernel[KS2+Step] = C2
   Kernel[KS2-Step] = C2
   Kernel[KS2]=C3
   Im_Out = convol(Im_in,Kernel, /EDGE_TRUNCATE) 
 end


pro b3spline_1D, Sig_in, Sig_out,  Step
; help, Sig_in
   N = double( (size(Sig_in))(1) )
   
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Sig_out  = Sig_in

   for i = 0., N-1. do begin
      im = test_ind_1D(i-Step, N)
      ip = test_ind_1D(i+Step, N)
      im2 = test_ind_1D(i-2*Step, N)
      ip2 = test_ind_1D(i+2*Step, N)
      Sig_out(i) = C3 * Sig_in(i) + C2 * (Sig_in(im) + Sig_in(ip)) + C1 * (Sig_in(im2) + Sig_in(ip2)) 
   endfor
end

;-------------------------------------------------------------------------------
pro softth, Data, Lambda
ind = where ( Data GT 0, c)
if c GT 0 then begin
   Pos = Data[ind] - Lambda
   indN = where (Pos LT 0, c)
   if c GT 0 then Pos[indN] = 0
   Data[ind] = Pos
end 
ind = where ( Data LT 0, c)
if c GT 0 then begin
   Neg = Data[ind] + Lambda
   indP = where (Neg GT 0, c)
   if c GT 0 then Neg[indP] = 0
   Data[ind] = Neg
end 
end
;-------------------------------------------------------------------------------

function star1d, Signal, Nscale=Nscale, Gen2=Gen2

        TransSignal = 1
	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
          print, 'CALLING SEQUENCE: TransSignal = star1d(Signal, Nscale=Nscale, Gen2=Gen2)'
        goto, DONE
	end

	if not keyword_set(Nscale) then Nscale =4

	if Nscale LT 2 then  Nscale =4
	
	;-------------------------------------------------------------------------------
	;recovering a few parameter values, initializing the transform loop ...
	N = (size(Signal))(1)

	NStep=Nscale-1
	TransSignal = fltarr(N, Nscale)
	Signal_in = Signal
	Signal_Out = 0.* Signal
	Step_trou = 1

	for i=0., NStep-1 do begin
		;B3 Spline smoothing
		b3spline_1d, Signal_in, Signal_Out,  Step_trou
        ; b3spline_star1d, Signal_in, Signal_Out,  Step_trou
		; Wavelet coefficient calculations
		if keyword_set(Gen2) then begin
		   Im_Aux =  Signal_Out ; just to Im_Aux to the right dimension
		   b3spline_1d,  Signal_Out, Im_Aux, Step_trou 
           ; b3spline_star1d, Signal_Out, Im_Aux, Step_trou
		   TransSignal[*,i] = Signal_in - Im_Aux
		end else  TransSignal[*,i] = Signal_in - Signal_Out
		
		Signal_in = Signal_Out
   
		; New distance between two pixels
		if i NE NStep-1 then Step_trou = Step_trou * 2.
        endfor


	; Smooth array
 	TransSignal[*,NStep] = Signal_Out 

DONE:
   return, TransSignal
end

;-------------------------------------------------------------------------------
; *********************************************************************/
pro b3spline_istar1d, Sig_in, Sig_out,  Step
   N = double( (size(Sig_in))(1) )
   
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Sig_out  = Sig_in

   for i = 0., N-1. do begin
      im = test_ind_1D(i-Step, N)
      ip = test_ind_1D(i+Step, N)
      im2 = test_ind_1D(i-2*Step, N)
      ip2 = test_ind_1D(i+2*Step, N)
      Sig_out(i) = C3 * Sig_in(i) + C2 * (Sig_in(im) + Sig_in(ip)) + C1 * (Sig_in(im2) + Sig_in(ip2)) 
   endfor
end


function istar1d,  WT, Gen2=Gen2
RecIma = -1

        ;-------------------------------------------------------------------------------
        ;a few verifications
        if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: Rec = istar1d(WT, Gen2=Gen2)'
        goto, DONE
        end

        vsize = size(WT)
        if vsize[0] NE 2 then begin
        print, 'Error: First parameter is not a 2D array ...'
        print, 'CALLING SEQUENCE:  WT = istar1d(WT, Gen2=Gen2)'
        goto, DONE
        end
        ;-------------------------------------------------------------------------------
                
vs = size(WT)
Nscale = vs[2]
NStep=Nscale-1
Step_trou = 1
for i=0, NStep-2 do Step_trou = Step_trou * 2
        
RecIma = WT[*,Nscale-1]
if keyword_set(Gen2) then Im_Out = fltarr(vs[1])

for j=Nscale-2,0,-1 do begin
  ; print, j, Step_trou
  if keyword_set(Gen2) then  begin
      b3spline_istar1d, RecIma, Im_Out, Step_trou
      ; b3spline_istar1d, RecIma, Im_Out, Step_trou
      RecIma = Im_Out + WT[*,j]
  end else RecIma = RecIma + WT[*,j]
  Step_trou = Step_trou / 2
end

DONE:
    return,  RecIma
end

;-------------------------------------------------------------------------------





function get_fdr_pvalue, TabPVal,  Alpha=Alpha,   Correl=Correl
vs = size(TabPVal)
N = vs[1]
if not keyword_set(Alpha) then Alpha = 0.05

    PVal = dblarr(N)
    p_cutoff=0d
   
    TabS = sort(TabPVal)
   TabS = TabPVal[TabS]
   j_alpha = Alpha * dindgen(N) / double(N)
   for  i = 0, N-1 do begin
       if  TabS[i] LE j_alpha[i]  then p_cutoff = TabS[i]
   end
  print, "p_cutoff =" , p_cutoff 
  for i = 0, N-1  do   if TabPVal[i] LT p_cutoff  then TabPVal[i] = 1 else TabPVal[i] = 0
  
return,  p_cutoff
end

; *********************************************************************/

function get_fdr_detect_level,  Coef,  Alpha=Alpha,  Correl=Correl, SigmaNoise=SigmaNoise, Nsig=Nsig, Prob=Prob
if not keyword_set(Alpha) then Alpha = 0.05
if not keyword_set(SigmaNoise) then SigmaNoise = 1.
if keyword_set(Nsig) then Alpha =1d - errorf(double(Nsig) / sqrt(double(2.)))

if keyword_set(Prob) then PVal = Coef $
else PVal = 1d - errorf(double(abs(Coef) / SigmaNoise) / sqrt(double(2.)))
Prob = get_fdr_pvalue(PVal,  Alpha=Alpha,   Correl=Correl)
Nsigma = ABS(inverf(1-2* Prob))*sqrt(2.)
return,  Nsigma*SigmaNoise
end

; *********************************************************************/

function mad, data
  n = n_elements(data)
  d = fltarr(n)
  d[*] = abs(data-median(data))
  ind = sort(d)
  d = d[ind]
  return, d[n/2] / 0.6747
end


FUNCTION GET_NOISE, Data, Niter=Niter 
IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=GET_NOISE(Data)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------

if not keyword_set(Niter) then Niter = 3.
vsize = size(Data)
dim = vsize[0]
sigma=-1
case dim of 
  3: BEGIN
        nco=vsize[1]
        nli=vsize[2]
        npz=vsize[3]
        indices=indgen(npz-2)+1
        D_cube=fltarr(nco,nli,npz)
        c1=-1./sqrt(6.)
        c2=2./sqrt(6.)
        D_cube[*,*,1:npz-2]=c1*(Data[*,*,indices-1] + Data[*,*,indices+1]) + $
                            c2*(Data[*,*,indices])
        D_cube[*,*,0]=c2*(Data[*,*,0]-Data[*,*,1])
        D_cube[*,*,npz-1]=c2*(Data[*,*,npz-1]-Data[*,*,npz-2])
        sigma=sigma_clip(D_cube,Niter=Niter)
     END
  2: BEGIN
        im_smooth, Data, ima_med, winsize=3, method='median'
        sigma = sigma_clip(Data-ima_med,Niter=Niter) / 0.969684
     END
  1: BEGIN
        ; Sigma = sigma_clip(Data - median(Data,3), Niter=Niter) / 0.893421
        w = star1d(Data, nscale=2) 
        Sigma = MAD( w[*,0] ) / 0.727707
     END
  else: BEGIN
        print, 'Error: bad dimension of the input parameter'
        END
 END

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, Sigma
 
 END

; ===============================

function pstar1d, Ima, Nscale=Nscale, Niter=Niter, Rec=Rec, soft=soft, SigmaNoise=SigmaNoise, Nsigma=Nsigma, wpos=wpos, killlastscale=killlastscale, fdr=fdr, pos=pos
Gen2=1 
Resi  = Ima

if not keyword_set(SigmaNoise) then SigmaNoise = get_noise(Ima)
if not keyword_set(Niter) then Niter = 10
if not keyword_set(Nsigma) then Nsigma = 3
AlphaFDR = 1d - errorf(double(Nsigma) / sqrt(double(2.)))

print, "SigmaNoise = ", SigmaNoise


dirac = Ima
dirac [*] = 0
vs = size(dirac)
Nx = vs[1]
dirac[Nx/2] = 1
w  = star1d(dirac, Nscale=Nscale, Gen2=Gen2)
TabNorm = fltarr(Nscale)
for j=0, Nscale-1 do TabNorm[j] = sqrt( total(w[*,j]^2) )

w  = star1d(Resi, Nscale=Nscale, Gen2=Gen2)
MW = max(w) / 2.
TabDetectLevel = fltarr(Nscale)
for j=0, Nscale-1 do begin
    if not keyword_set(fdr) then  begin
        Nsig= Nsigma
        if j EQ 0 then Nsig= Nsigma+1
        TabDetectLevel[j] = SigmaNoise*TabNorm[j]*Nsig 
    end else begin
      Coef = W[*,j]
      print, AlphaFDR
      info, Coef / (SigmaNoise*TabNorm[j])
      TabDetectLevel[j] = get_fdr_detect_level(Coef, Alpha=AlphaFDR, SigmaNoise=SigmaNoise*TabNorm[j])
      print, 'Scale ', j+1, ': FDR Detect Level = ', TabDetectLevel[j], ', Nsig = ', TabDetectLevel[j] / (SigmaNoise*TabNorm[j])
    end
end


for i=0,Niter+10 do begin
   Lambda = MW * ( 1. - (i+1.) / Niter)
   if Lambda LT 0 then Lambda = 0
   w  = star1d(Resi, Nscale=Nscale, Gen2=Gen2)
   
    if i EQ 0 then WT = w $ 
   else WT = WT + w

   for j=0, Nscale-2 do begin
     Scale = WT[*,j]
     ind = where( ABS(Scale)  LT TabDetectLevel[j], c)
     if c GT 0 then Scale[ind] = 0
     WT[*,j] = Scale
   end
   if keyword_set(killlastscale) then WT[*,Nscale-1] = 0
   
    
   if keyword_set(soft) then begin
      for j=0, Nscale-2 do begin
          Scale = WT[*,j]
          softth, Scale,  Lambda* TabNorm[j]
          WT[*,j] = Scale
      end
  end
  
  if keyword_set(wpos) then begin
     ind = where(WT LT 0, c)
     if c GT 0 then WT[ind] = 0
  end
  
   Rec = istar1d(WT, Gen2=Gen2)
   
   if keyword_set(Pos) then begin
      ind = where(Rec LT 0, c)
      if c GT 0 then Rec[ind] = 0
   end
   
   Resi = Ima - Rec
   print, "Iter " , i+1, Lambda , ", Resi =  " ,  sigma(Resi)
end

DONE:

 return, WT
 
end


;============================================================

function star1d_filter, Ima, Nscale=Nscale, Niter=Niter, Rec=Rec, soft=soft, SigmaNoise=SigmaNoise, Nsigma=Nsigma, wpos=wpos, killlastscale=killlastscale, fdr=fdr, pos=pos
Gen2=0
UpdateSupport=0
Resi  = Ima

if not keyword_set(SigmaNoise) then SigmaNoise = get_noise(Ima)
if not keyword_set(Niter) then Niter = 10
if not keyword_set(Nsigma) then Nsigma = 3
AlphaFDR = 1d - errorf(double(Nsigma) / sqrt(double(2.)))

print, "SigmaNoise = ", SigmaNoise


dirac = Ima
dirac [*] = 0
vs = size(dirac)
Nx = vs[1]
dirac[Nx/2] = 1
w  = star1d(dirac, Nscale=Nscale, Gen2=Gen2)
TabNorm = fltarr(Nscale)
for j=0, Nscale-1 do TabNorm[j] = sqrt( total(w[*,j]^2) )

w  = star1d(Resi, Nscale=Nscale, Gen2=Gen2)
MW = max(w) / 2.
TabDetectLevel = fltarr(Nscale)
TabSupport = fltarr(Nx, Nscale) + 1

for j=0, Nscale-1 do begin
    Coef = W[*,j]
    if not keyword_set(fdr) then  begin
        Nsig= Nsigma
        if j EQ 0 then Nsig= Nsigma+1
        TabDetectLevel[j] = SigmaNoise*TabNorm[j]*Nsig 
    end else begin
      TabDetectLevel[j] = get_fdr_detect_level(Coef, Alpha=AlphaFDR, SigmaNoise=SigmaNoise*TabNorm[j])
      print, 'Scale ', j+1, ': FDR Detect Level = ', TabDetectLevel[j], ', Nsig = ', TabDetectLevel[j] / (SigmaNoise*TabNorm[j])
    end
    if not keyword_set(wpos) then  Coef = ABS(Coef) 
    ind = where(Coef  LT TabDetectLevel[j], c)
    Coef[*]=1
    if c GT 1 then Coef[ind] = 0
    TabSupport[*,j] = Coef
end
if  keyword_set(killlastscale) then TabSupport[*,Nscale-1] = 0

Rec = Ima * 0.
for i=0,Niter+10 do begin
   Lambda = MW * ( 1. - (i+1.) / Niter)
   if Lambda LT 0 then Lambda = 0
   w  = star1d(Resi, Nscale=Nscale, Gen2=Gen2)
   
   if keyword_set(UpdateSupport) then begin
      for j=0, Nscale-2 do begin
          if keyword_set(wpos) then  Scale = W[*,j] $
         else Scale = ABS(W[*,j])
         ind = where(Scale  GE TabDetectLevel[j], c)
         Scale = TabSupport[*,j] 
         if c GT 0 then Scale[ind] = 1
         TabSupport[*,j] =  Scale
    end
   end
   
   w = w * TabSupport
   ResiRec = istar1d(W, Gen2=Gen2)
   Rec = Rec + ResiRec

   if keyword_set(soft) then begin
      WT  = star1d(Rec, Nscale=Nscale, Gen2=Gen2)
      for j=0, Nscale-2 do begin
          Scale = WT[*,j]
          softth, Scale,  Lambda* TabNorm[j]
          WT[*,j] = Scale
      end
      Rec = istar1d(W, Gen2=Gen2)
  end
  
   if keyword_set(Pos) then begin
      ind = where(Rec LT 0, c)
      if c GT 0 then Rec[ind] = 0
   end
   
   Resi = Ima - Rec
   print, "Iter " , i+1, Lambda , ", Resi =  " ,  sigma(Resi)
end

DONE:

 return, Rec
 
end
 
;============================================================

function spectrum_decomp, Spectrum, Nscale=Nscale, Niter=Niter, Rec=Rec, soft=soft, SigmaNoise=SigmaNoise, Nsigma=Nsigma, wpos=wpos,   fdr=fdr, pos=pos, v2=v2

Niter = 50
killlastscale=1
if not keyword_set(v2) then Rec = star1d_filter(Spectrum, Nscale=Nscale, Niter=Niter, Rec=Rec, soft=soft, SigmaNoise=SigmaNoise, Nsigma=Nsigma, wpos=wpos, killlastscale=killlastscale, fdr=fdr, pos=pos)  $
else z = pstar1d(Spectrum, niter=Niter, nscale=Nscale, rec=Rec, nsig= Nsigma, killlastscale=killlastscale, pos=pos, wpos=wpos, SigmaNoise=SigmaNoise, fdr=fdr, /soft)

w  = star1d(Spectrum-Rec, Nscale=Nscale, Gen2=Gen2)
Continuum = w[*, Nscale-1]
Resi = Spectrum - Rec - Continuum

if not keyword_set(v2) then z = star1d(Rec, Nscale=Nscale, Gen2=Gen2)

Res = {Data: Spectrum, FilterData: Rec, Resi: Resi, Continuum: Continuum, Coef: z, SigmaNoise:SigmaNoise}
plot,  Res.data - Res.FILTERDATA - Res.CONTINUUM

info, Resi
return, Res
end

;*********************************************************************/

pro test
d = rim('specnoiseeq2.fits')
s = rim('purespec.fits')
p = rim('purespec.fits')
SigmaNoise = 0.131827
SigmaNoise = 0.12139959

opt='-k -f3 -K -n7 -i20 -C2 -p -v -s3 -g0.12'
mr1d_filter, d, f, opt=opt

z1 = spectrum_decomp(D, Nscale=7, Niter=10,   Nsigma=3, /fdr, /pos, /wpos)
plot, z1.filterdata
z2 = spectrum_decomp(D, Nscale=7, Niter=10,   Nsigma=3, /fdr, /pos, /wpos, /v2)
c = z2.Continuum
plot, z2.filterdata


end

;*********************************************************************/




