;+
; NAME:
;        mrs_cbfilter
;
; PURPOSE:
;	Combined filtering using Wavelet and Curvelet of an image on the sphere (Healpix pixel representation).
;       By default Gaussian noise is considered. If the keyword SigmaNoise is not 
;       set, then the noise standard deviation is automatically estimated.
;       If the keyword MAD is set, then a correlated Gaussian noise is considered,
;       and the noise level at each scale is derived from the Median Absolution Deviation (MAD)
;       method. If the keyword KillLastScale is set, the coarsest resolution is set to zero.
;       If the "undec" keyword is used, then a undecimated decomposition is used instead of the
;       pyramidal WT.
;
; CALLING:
;
;     mrs_cbfilter, Image, Filter, NbrScale=NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, 
;                   mad=mad, KillLastScale=KillLastScale, Undec=Undec, FirstBlockSize=FirstBlockSize,
;                   niter=niter, pos=pos, FirstScale=FirstScale       
;    
; INPUT:
;     Imag -- IDL array of healpix map: Input image be filtered 
;
; OUTPUT:
;     Filter -- IDL array of healpix map: reconstructed image from the thresholded wavelet coefficients   
;
; INPUT KEYWORDS:
;      NbrScale : int = Number of scales (default is 4)
;      NSigma: float = Level of thresholding (default is 3)
;      SigmaNoise: float = Noise standard deviation. Default is automatically estimated
;      MAD: int: if set, then the noise level is derive at each scale using the MAD of the
;                        wavelet coefficient. MAD = median ( ABS( WaveletScale) ) / 0.6745
;      KillLastScale: if set, the last scale is set to zero
;      FirstScale: int: Consider only scales larger than FirstScale. Default is 1 (i.e. all scales are used).
;      niter: int = number of iterations. Default is 10.
;      pos: scalar = if set, a positivity constraint is added
;
;
; OUTPUT KEYWORDS:
;
; EXTERNAL CALLS:
;       mrs_curtrans
;   	mrs_wcurrec
;       mrs_curget  
;   	mrs_curput
;
; EXAMPLE:
;       Filter an image with 5 scales. The result is stored in Filter 
;               mrs_cbfilter, Data, Filter, NbrScale=5 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;-

pro mrs_cbfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale,  Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos, FirstScale=FirstScale

 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_cbfilter, Data, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, Mad=Mad, KillLastScale=KillLastScale, Undec=Undec'
        goto, DONE
        end
;if keyword_set(Mask) then begin
;   vs = size(Mask)
;   Nm = vs[1]
;   if Nm NE npix then begin
;       print, 'Error: the mask must have the same dimension as the input map'
;       goto , DONE
;   end
;end

if not keyword_set(niter) then niter=10
if not keyword_set(FirstScale) then FirstScale=1
if not keyword_set(NbrScale) then NbrScale=4
if not keyword_set(NSigma) then NSigma=3.
Pyr=1   
if keyword_set(Undec) then Pyr=0

if Pyr EQ 0 then mrs_wttrans, Imag, TransWT, NbrScale=NbrScale $
else mrs_pwttrans, Imag, TransWT, NbrScale=NbrScale

TabWTThreshold = fltarr(NbrScale)
TabWTNoise = fltarr(NbrScale)

for j =0,NbrScale-2 do begin
    Scale = mrs_wtget(TransWT,j,NormVal=NormVal)
    
       ; Noise standard deviation estimation from the first scale
    if j EQ 0 and not keyword_set(MAD) and not keyword_set(SigmaNoise) then begin
          SigmaNoise=get_noise(Scale) / NormVal
          print, 'Noise standard deviation estimation = ', SigmaNoise
       end
       
    if j LT FirstScale-1 then Scale[*] = 0 $
    else BEGIN
        ; Threshold levele at scale j
       if j eq 0 then NS=NSigma+1 else NS=NSigma
       if not keyword_set(MAD) then ThresholdLevel = NS*SigmaNoise*NormVal $
       else ThresholdLevel = NS * median( abs(Scale)) / 0.6745
       TabWTThreshold[j] = ThresholdLevel
       TabWTNoise[j] = ThresholdLevel / NS
    
       index = where ( ABS(Scale) LT ThresholdLevel, count)
       if count GT 0 then Scale[index] = 0
    END
    mrs_wtput, TransWT, Scale, j
endfor
LastWTScale  = mrs_wtget(TransWT,NbrScale-1)
if keyword_set(KillLastScale) then LastWTScale[*]=0

; CURVELET PART
mrs_curtrans, Imag, TransCur, NbrScale=NbrScale, Undec=Undec, FirstBlockSize=FirstBlockSize, /overlap

vs = size(TransCur.TabNorm)
NxTNorm = vs[1]
NyTNorm = vs[2]

TabCurThreshold = fltarr(NxTNorm, NyTNorm)
TabCurNoise = fltarr(NxTNorm, NyTNorm)

for s2d = 0,TransCur.NbrScale-2 do begin
  for s1d = 0,TransCur.TabNbrScaleRid[s2d]-1 do begin
    Band = mrs_curget(TransCur, s2d, s1d)
    if s2d LT FirstScale-1 then Band[*] = 0 $
    else BEGIN
       s2d1 = s2d
       s1d1 = s1d
       if s2d1 GE NyTNorm then s2d1 = NyTNorm - 1
       if s1d1 GE NxTNorm then s1d1 = NxTNorm - 1
       if s1d1 EQ TransCur.TabNbrScaleRid[s2d]-1 and TransCur.TabNbrScaleRid[s2d] lt NxTNorm then s1d1 = s1d1-1
       NormVal = TransCur.TabNorm[s1d1,s2d1] * sqrt(TransCur.TabBlockSize[s2d])
       ; Threshold level  
       if s2d EQ 0 and s1d EQ 0 then NS=NSigma+1 else NS=NSigma
       if not keyword_set(MAD) then ThresholdLevel = NS*SigmaNoise*NormVal $
       else ThresholdLevel = NS * median( abs(Band)) / 0.6745
       TabCurThreshold[s1d1,s2d1] = ThresholdLevel
       TabCurNoise[s1d1,s2d1] = ThresholdLevel / NS
       index = where ( ABS(Band) LT ThresholdLevel, count)
       print, s2d, s1d, sigma(Band), ' Threshold = ', ThresholdLevel, ' Count = ',  count
       if count GT 0 then Band[index] = 0
    END
    mrs_curput, TransCur, Band, s2d, s1d
  end
end

if keyword_set(KillLastScale) then TransCur.LastScale[*] = 0
mrs_currec, TransCur, Filter 
if keyword_set(pos) then begin
  ind = where ( Filter LT 0, c)
  if c GT 0 then  Filter[ind] = 0
end

if keyword_set(niter) then begin
   Lambda = 1.
   delta = Lambda / float(niter-1)
   TabCurEps = TabCurNoise / 2
   LastCurScale = TransCur.LastScale
   TabWTEps = TabWTNoise / 2
   TabCurEps = TabCurNoise / 2

   for i=1,niter do begin
     ; WT PROJECTION
      if Pyr EQ 0 then mrs_wttrans, Filter, TransWTFil, NbrScale=NbrScale $
      else mrs_pwttrans,  Filter, TransWTFil, NbrScale=NbrScale 

      for j =0,NbrScale-2 do begin
         Scale = mrs_wtget(TransWT,j)
         ScaleFil = mrs_wtget(TransWTFil,j)
	 ScaleResi = Scale - ScaleFil
	 ind = where (Scale NE 0, c)
	 if c GT 0 then begin
	   ResiSignif = ScaleResi[ind]
	   indr = where (ABS(ResiSignif) LT TabWTEps[j], c)
	   if c GT 0 then ResiSignif[indr] = 0
	   ScaleFil[ind] = ScaleFil[ind] + ResiSignif
	 end
	   
	 ; Now soft thresholding
	 if Lambda GT 0 then softthreshold, ScaleFil, Lambda*TabWTNoise[j]
         mrs_wtput, TransWTFil, ScaleFil, j
      end
      mrs_wtput, TransWTFil, LastWTScale, NbrScale-1
      
      if Pyr EQ 0 then mrs_wtrec,  TransWTFil, Filter $
      else mrs_pwtrec, TransWTFil, Filter      
      if keyword_set(pos) then begin
          ind = where ( Filter LT 0, c)
         if c GT 0 then  Filter[ind] = 0
      end
   
   
     ; CURVELET PROJECTION
     mrs_curtrans, Filter, TransCurFil, NbrScale=NbrScale, Undec=Undec, FirstBlockSize=FirstBlockSize, /overlap, /silent
     for s2d = 0,TransCur.NbrScale-2 do begin
     for s1d = 0,TransCur.TabNbrScaleRid[s2d]-1 do begin
         Scale = mrs_curget(TransCur, s2d, s1d)
         ScaleFil = mrs_curget(TransCurFil,s2d, s1d)
	 ScaleResi = Scale - ScaleFil
	 ind = where (Scale NE 0, c)
	 if c GT 0 then begin
	   ResiSignif = ScaleResi[ind]
	   indr = where (ABS(ResiSignif) LT TabCurEps[s1d,s2d], c)
	   if c GT 0 then ResiSignif[indr] = 0
	   ScaleFil[ind] = ScaleFil[ind] + ResiSignif
	 end
 	 ; Now soft thresholding
	 if Lambda GT 0 then softthreshold, ScaleFil, Lambda*TabCurNoise[s1d,s2d]
         mrs_curput, TransCurFil, ScaleFil, s2d, s1d 
     end
     end
     
     TransCurFil.LastScale = LastCurScale
     mrs_currec, TransCurFil, Filter 
     if keyword_set(pos) then begin
         ind = where (Filter LT 0, c)
         if c GT 0 then  Filter[ind] = 0
      end
      resi = Imag - Filter 
      ; tvs, resi 
      print, 'Iter ', i, '   ==>  SigmaResi = ', sigma(Resi), '  Lambda = ', Lambda
      Lambda = Lambda - delta 
    end
end


DONE:

END

