;+
; NAME:
;        mrsp_wtfilter
;
; PURPOSE:
;	Wavelet denoising of a POLARIZED image on the sphere (Healpix NESTED pixel representation).
;       By default Gaussian noise is considered. If the keyword SigmaNoise is not 
;       set, then the noise standard deviation is automatically estimated.
;       If the keyword MAD is set, then a correlated Gaussian noise is considered,
;       and the noise level at each scale is derived from the Median Absolution Deviation (MAD)
;       method. If the keyword KillLastScale  is set, the coarsest resolution is set to zero.
;       The thresholded wavelet coefficients can be obtained using the keyword Trans.
;       If the input keyword NITER is set, then an iterative algorithm is applied and
;       if the POS keyword is also set, then a positivity constraint is added.
;
; CALLING:
;
;		mrsp_wtfilter, Imag, Filter, NbrScale=NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, mad=mad, KillLastScale=KillLastScale,
;						Trans=Trans, niter=niter, pos=pos, FirstScale=FirstScale, soft=soft, fdr=fdr, Use_FdrAll=Use_FdrAll,
;						lmax=lmax, FilterLast=FilterLast, mask=mask
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
;      niter: int: number of iterations used in the reconstruction
;      pos: if set, the solution is assumed to be positive
;      FirstScale: int: Consider only scales larger than FirstScale. Default is 1 (i.e. all scales are used).
;		soft: if set, use soft thresholding instead of hard thresholding
;		fdr: float between 0 (default) and 1 (max, if greater or equal to 1, set to 0.05),
;						used to estimate a threshold level instead of a NSigma threshold,
;						threshold is applied from scale j=FirstScale to the last.
;		Use_FdrAll: same as fdr but applied to all scales.
;		FilterLast: if set, the last scale is filtered.
;		mask: IDL array of healpix map, input mask applied.
;		pos: if set, the solution is assumed to be positive
;
; INPUT/OUTPUT:
;		lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;
; OUTPUT KEYWORDS:
;      Trans -- IDL structure: Threshold wavelet decomposition of the input image
;
; EXTERNAL CALLS:
;       mrsp_wttrans
;   	mrsp_wtrec
;       mrs_pwttrans
;   	mrs_pwtrec
;       mrs_wtget
;   	mrs_wtput
;
; EXAMPLE:
;       Filter an image with 5 scales. The result is stored in Filter 
;               mrsp_wtfilter, Data, Filter, NbrScale=5 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;-

pro mrsp_wtfilter, Imag, Filter, NbrScale=NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, KillLastScale=KillLastScale, Trans=Trans, Undec=Undec, niter=niter, $
		  pos=pos, FirstScale=FirstScale, soft=soft, fdr=fdr, Use_FdrAll=Use_FdrAll, $
		  lmax=lmax, FilterLast=FilterLast, mask=mask, atrou=atrou ; , NSigMask=NSigMask

 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
print, 'CALLING SEQUENCE: mrsp_wtfilter, Imag, Filter, NbrScale=NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, mad=mad, KillLastScale=KillLastScale, Trans=Trans, niter=niter, pos=pos, FirstScale=FirstScale, soft=soft, fdr=fdr, Use_FdrAll=Use_FdrAll, lmax=lmax, FilterLast=FilterLast, mask=mask'
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

nside=0
Pyr=1   
if keyword_set(Undec) then Pyr=0
if keyword_set(atrou) then Pyr=0 

if not keyword_set(NbrScale) then NbrScale=4
if not keyword_set(NSigma) then NSigma=3.
if not keyword_set(NSigMask) then NSigMask=NSigma

if not keyword_set(FirstScale) then FirstScale=1

mrsp_wttrans, Imag, Trans, NbrScale=NbrScale, lmax=lmax

for comp =0,2 do begin

TabThreshold = fltarr(NbrScale)
TabNoise = fltarr(NbrScale)

if keyword_set(FilterLast) then LS = NbrScale-1 $
else  LS = NbrScale-2

for j =0,LS do begin
    Scale = mrsp_wtget(Trans,j,comp,NormVal=NormVal)
    
    if j LT FirstScale-1 then Scale[*] = 0 $
    else BEGIN
      ; Noise standard deviation estimation from the first scale
      if j EQ 0 and not keyword_set(MAD) and not keyword_set(SigmaNoise) then begin
         SigmaNoise=get_noise(Scale) / NormVal
       print, 'Noise standard deviation estimation = ', SigmaNoise
      end
    
      ; Threshold level at scale j
      if j eq 0 then NS=NSigma+1 else NS=NSigma
      if j eq 0 then NSM=NSigMask+1 else NSM=NSigMask
      if not keyword_set(MAD) then ThresholdLevel = NS*SigmaNoise*NormVal $
      else begin
         if keyword_set(mask) then begin
            Ind = where( Mask EQ 1)
            ThresholdLevel = NS * mad(Scale[Ind]) 
          end else ThresholdLevel = NS * mad(Scale)
      end
      TabThreshold[j] = ThresholdLevel
      TabNoise[j] = ThresholdLevel / NS
      ThresholdLevelMask = TabNoise[j] * NSM 
      
      if not keyword_set(fdr) then begin
         ScaleIn = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
	 if keyword_set(mask) and  (NSigMask NE NSigma) then begin
	    ScaleOut = mrs_absthreshold(Scale, ThresholdLevelMask, soft=soft)
	    Ind = where( Mask EQ 0)
	    ScaleIn[Ind] = ScaleOut[Ind]
	 end   
	 Scale=ScaleIn
       end else if not keyword_set(Use_FdrAll) then begin
         coef  = mrs_getpix(Scale) / TabNoise[j]
         vs = size(coef)
         nx = vs[1]
         ny = vs[2]
         nel = N_elements(coef)
         Coef1D = reform(coef, nel)
         if fdr GE 1 then fdr = 0.05
         opt = ' -v -g 1. -a   ' + STRCOMPRESS(string(fdr), /REMOVE_ALL) 
         mr_prog, 'im_fdr' , Coef1D, Res, opt=opt
         ind = where (Res EQ 0, c)
         if c gt 0 then  coef[ind] = 0
	 mrs_putpix, Scale, Coef
      end
      
      ;index = where ( ABS(Scale) LT ThresholdLevel, count)
      ;if count GT 0 then Scale[index] = 0
    END
    mrsp_wtput, Trans, Scale, j,comp
endfor



if keyword_set(KillLastScale) then begin
    Scale = mrsp_wtget(Trans,NbrScale-1,comp)
    mrs_set, Scale, 0.
    mrsp_wtput, Trans, Scale, NbrScale-1,comp
end

endfor ; end component
 mrsp_wtrec,  Trans, Filter 




if keyword_set(pos) then mrs_pos, Filter



DONE:

END

