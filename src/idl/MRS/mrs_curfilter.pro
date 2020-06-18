;+
; NAME:
;        mrs_curfilter
;
; PURPOSE:
;	Curvelet denoising of an image on the sphere (Healpix pixel representation).
;       By default Gaussian noise is considered. If the keyword SigmaNoise is not 
;       set, then the noise standard deviation is automatically estimated.
;       If the keyword MAD is set, then a correlated Gaussian noise is considered,
;       and the noise level at each scale is derived from the Median Absolution Deviation (MAD)
;       method. If the keyword KillLastScale the coarsest resolution is set to zero.
;       If the UNDEC keyword is used, then a undecimated decomposition is used instead of the
;       pyramidal WT. If the keyword CYCLE is set, the denoising is performed three times,
;       by shifting the data by PI/4 and -PI/4, denoising the shifted version, and averaging
;       the unshifted denoising maps. This procedure also us to remove the block effect
;       which may appear on the border of the Healpix faces.
;       The threshold curvelet coefficient can be obtained using the keyword Trans.
;  
;       If the input keyword NITER is set, then an iterative algorithm is applied and
;       if the POS keyword is also set, then a positivity constraint is added.
;   
;
; CALLING:
;
;     mrs_curfilter, Image, Filter, NbrScale=NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, 
;                   mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, 
;                   FirstBlockSize=FirstBlockSize, niter=niter, pos=pos, cycle=cycle, FirstScale=FirstScale  
;    
; INPUT:
;     Imag -- IDL array of healpix map: Input image be filtered 
;
; OUTPUT:
;     Filter -- IDL array of healpix map: reconstructed image from the thresholded wavelet coefficients   
;
; INPUT KEYWORDS:
;		NbrScale : int = Number of scales (default is 4)
;		NSigma: float = Level of thresholding (default is 3)
;		SigmaNoise: float = Noise standard deviation. Default is automatically estimated
;		MAD: int: if set, then the noise level is derive at each scale using the MAD of the
;                        wavelet coefficient. MAD = median ( ABS( WaveletScale) ) / 0.6745
;		KillLastScale: if set, the last scale is set to zero
;		Undec: if set, an undecimated WT is used instead of the the pyramidal WT
;		niter: int: number of iterations. By default, there is no iterations
;		pos: int: if set and if niter, then a positivity constaint is added.
;		cycle: int: if set, then a cycle spanning is applied.
;		FirstScale: int: Consider only scales larger than FirstScale. Default is 1 (i.e. all scales are used).
;		FirstBlockSize: int, Block size in the ridgelet transform at the finest scale (default is 16)
;
; OUTPUT KEYWORDS:
;		Trans -- IDL structure: Thresholded curvelet decomposition of the input image
;
; EXTERNAL CALLS:
;       mrs_curtrans
;   	mrs_currec
;       mrs_curget  
;   	mrs_curput
;
; EXAMPLE:
;       Filter an image with 5 scales. The result is stored in Filter 
;               mrs_curfilter, Data, Filter, NbrScale=5 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	September, 2005 File creation
;-

pro mrs_curfilter_one_cycle, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos, FirstScale=FirstScale

if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_curfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, niter=niter, pos=pos, FirstScale=FirstScale'
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

if not keyword_set(NbrScale) then NbrScale=4
if not keyword_set(NSigma) then NSigma=3.
if not keyword_set(FirstScale) then FirstScale=1

mrs_curtrans, Imag, Trans, NbrScale=NbrScale, Undec=Undec, FirstBlockSize=FirstBlockSize, /overlap, opt=' -O '

; Noise estimation
if not  keyword_set(MAD)  and not keyword_set(SigmaNoise) then begin
   Scale = mrs_wtget(Trans.WT, 0,NormVal=NormVal)
   SigmaNoise=get_noise(Scale) / NormVal
   print, 'Estimated SigmaNoise  = ', SigmaNoise
end

vs = size(Trans.TabNorm)
NxTNorm = vs[1]
NyTNorm = vs[2]

TabThreshold = fltarr(NxTNorm, NyTNorm)
TabNoise = fltarr(NxTNorm, NyTNorm)

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
      if s2d EQ 0 and s1d EQ 0 then NS=NSigma+1 else NS=NSigma
      if not keyword_set(MAD) then ThresholdLevel = NS*SigmaNoise*NormVal $
      else ThresholdLevel = NS * median( abs(Band)) / 0.6745
      TabThreshold[s1d1,s2d1] = ThresholdLevel
      TabNoise[s1d1,s2d1] = ThresholdLevel / NS
      index = where ( ABS(Band) LT ThresholdLevel, count)
      print, s2d, s1d, sigma(Band), ' Threshold = ', ThresholdLevel, ' Count = ',  count
      if count GT 0 then Band[index] = 0
      END
    mrs_curput, Trans, Band, s2d, s1d
  end
end

if keyword_set(KillLastScale) then Trans.LastScale[*] = 0

print, 'Reconstruction '
mrs_currec, Trans, Filter 
if keyword_set(pos) then begin
  ind = where ( Filter LT 0, c)
  if c GT 0 then  Filter[ind] = 0
end

if keyword_set(niter) then begin
   Lambda = 1.
   delta = Lambda / float(niter-1)
   TabEps = TabNoise / 2
   LastScale = Trans.LastScale
   
   for i=1,niter do begin
     mrs_curtrans, Filter, TransFil, NbrScale=NbrScale, Undec=Undec, FirstBlockSize=FirstBlockSize, /overlap, opt=' -O ', /silent
     for s2d = 0,Trans.NbrScale-2 do begin
     for s1d = 0,Trans.TabNbrScaleRid[s2d]-1 do begin
         Scale = mrs_curget(Trans, s2d, s1d)
         ScaleFil = mrs_curget(TransFil,s2d, s1d)
	 ScaleResi = Scale - ScaleFil
	 ind = where (Scale NE 0, c)
	 if c GT 0 then begin
	   ResiSignif = ScaleResi[ind]
	   indr = where (ABS(ResiSignif) LT TabEps[s1d,s2d], c)
	   if c GT 0 then ResiSignif[indr] = 0
	   ScaleFil[ind] = ScaleFil[ind] + ResiSignif
	 end
 	 ; Now soft thresholding
	 if Lambda GT 0 then softthreshold, ScaleFil, Lambda*TabNoise[s1d,s2d]
         mrs_curput, TransFil, ScaleFil, s2d, s1d 
     end
     end
     
     TransFil.LastScale = LastScale
     mrs_currec, TransFil, Filter 
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

;==================================================================

pro mrs_curfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos, cycle=cycle, FirstScale=FirstScale

if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_curfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, 	  niter=niter, pos=pos, cycle=cycle, FirstScale=FirstScale'
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

if not keyword_set(FirstScale) then FirstScale = 1

if not keyword_set(cycle) then begin
    mrs_curfilter_one_cycle, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos 
end else begin

   npixel = (size(Imag))[1]
   nside = npix2nside(npixel)
   mak_map,nside,interpole,t_interpol = 1

   rotate_map_nest,Imag,!dpi/4,0,0,Imag2
   rotate_map_nest,interpole,!dpi/4,0,0,interpole2
   rotate_map_nest,Imag,0,!dpi/2,0,Imag3
   rotate_map_nest,interpole,0,!dpi/2,0,interpole3

   mrs_curfilter_one_cycle, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos, FirstScale=FirstScale

   mrs_curfilter_one_cycle, Imag2, Filter2, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos, FirstScale=FirstScale

    mrs_curfilter_one_cycle, Imag3, Filter3, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, FirstBlockSize=FirstBlockSize, $
		  niter=niter, pos=pos, FirstScale=FirstScale

   rotate_map_nest,Filter2,-!dpi/4,0,0,F2
   rotate_map_nest,Filter3,0,-!dpi/2,0,F3

   FilterFinal = (  (Filter * interpole) + (Filter2 * interpole2)+ (Filter3 * interpole3) ) / (interpole+interpole2+interpole3 )
   Filter = FilterFinal
   end

DONE:


END
