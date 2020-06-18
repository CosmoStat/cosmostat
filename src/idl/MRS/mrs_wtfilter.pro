;+
; NAME:
;        mrs_wtfilter
;
; PURPOSE:
;	Wavelet denoising of an image on the sphere (Healpix pixel NESTED representation) or Glesp Data representation.
;       By default Gaussian noise is considered. If the keyword SigmaNoise is not 
;       set, then the noise standard deviation is automatically estimated.
;       If the keyword MAD is set, then a correlated Gaussian noise is considered,
;       and the noise level at each scale is derived from the Median Absolution Deviation (MAD)
;       method. If the keyword KillLastScale is set, the coarsest resolution is set to zero.
;       If the "Pyr" keyword is used, then the pyramidal WT is used instead of the undecimated WT.
;		If the "atrou" keyword is used, then the "a trou" WT is used instead of the undecimated WT.
;		If the keyword CYCLE is set, the denoising is performed three times,
;       by shifting the data by PI/4 and -PI/4, denoising the shifted version, and averaging
;       the unshifted denoising maps. This procedure also us to remove the block effect
;       which may appear on the border of the Healpix faces.
;       The thresholded wavelet coefficients can be obtained using the keyword Trans.
;       If the input keyword NITER is set, then an iterative algorithm is applied and
;       if the POS keyword is also set, then a positivity constraint is added.
;
; CALLING:
;
;	mrs_wtfilter, Imag, Filter, NbrScale=NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, lmax=lmax, TabNSigma=TabNSigma,
;				  mad=mad, localmad=localmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale, Trans=Trans, Pyr=Pyr,
;				  niter=niter, pos=pos, cycle=cycle, FirstScale=FirstScale, Soft=Soft, fdr=fdr, Use_FdrAll=Use_FdrAll,
;				  FilterLast=FilterLast, mask=mask, atrou=atrou, OutMask=OutMask
;    
; INPUT:
;     Imag -- IDL array of healpix map or Glesp image IDL structure: Input image be filtered 
;
; OUTPUT:
;     Filter -- IDL array of healpix map or Glesp image IDL structure: reconstructed image from the thresholded wavelet coefficients   
;
; INPUT KEYWORDS:
;		NbrScale: int = Number of scales (default is 4)
;		NSigma: float = Level of thresholding (default is 3)
;		SigmaNoise: float = Noise standard deviation. Default is automatically estimated
;		TabNSigma: float array = Level of thresholding at each scale
;		MAD: int if set, then the noise level is derive at each scale using the MAD of the
;                        wavelet coefficient. MAD = median ( ABS( WaveletScale) ) / 0.6745
;		localmad: int if set, similar to keyword MAD but with one value for each patch of the
;						image and each instead of one value for the full image at each scale
;						DO NOT WORK WITH GLESP
;
;		WinMinSize: int minimal size of patches (default is 8)
;		KillLastScale: if set, the last scale is set to zero
;		niter: int: number of iterations used in the reconstruction
;		pos: if set, the solution is assumed to be positive
;		Pyr: if set, a pyramidal WT is used instead of the the undecimated WT
;		cycle: int if set, then a cycle spanning is applied.
;						DO NOT WORK WITH GLESP
;
;		FirstScale: int Consider only scales larger than FirstScale. Default is 1 (i.e. all scales are used).
;		Soft: if set, use soft thresholding instead of hard thresholding
;		fdr: float between 0 (default) and 1 (max, if greater or equal to 1, set to 0.05),
;						used to estimate a threshold level instead of a NSigma threshold,
;						threshold is applied from scale j=FirstScale to the last.
;		Use_FdrAll: same as fdr but applied to all scales.
;		FilterLast: if set, the last scale is filtered.
;		mask: IDL array of healpix map, input mask applied.
;		atrou: if set, a "a trou" WT is used instead of the the undecimated WT
;						DO NOT WORK WITH GLESP
;
; INPUT/OUTPUT:
;		lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;
; OUTPUT KEYWORDS:
;		Trans -- IDL structure: Thresholded wavelet decomposition of the input image
;		OutMask: IDL array of healpix map, part of Imag that were set to 0 via the filtering, including keyword mask if used.
;
; EXTERNAL CALLS:
;       mrs_wttrans
;   	mrs_wtrec
;       mrs_pwttrans
;   	mrs_pwtrec
;       mrs_wtget
;   	mrs_wtput
;
; EXAMPLE:
;       Filter an image with 5 scales. The result is stored in Filter 
;               mrs_wtfilter, Data, Filter, NbrScale=5 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;
;---------------------------------------------------------------------------------------------------------------------------------------------------------

pro mrs_wtfilter_one_cycle, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, $
                  mad=mad, localmad=localmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale,Trans=Trans, Undec=Undec, niter=niter, $
		  pos=pos, FirstScale=FirstScale, soft=soft, fdr=fdr, Use_FdrAll=Use_FdrAll, $
		  lmax=lmax, FilterLast=FilterLast, mask=mask, atrou=atrou, TabNSigma=TabNSigma, UseTabThreshold=UseTabThreshold  ; , NSigMask=NSigMask

 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_wtfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, lmax=lmax, TabNSigma=TabNSigma, mad=mad, localmad=localmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale, Trans=Trans, Pyr=Pyr, niter=niter, pos=pos, cycle=cycle, FirstScale=FirstScale, Soft=Soft, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, atrou=atrou, OutMask=OutMask'
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
if type_code(Imag) EQ 8 then begin
	GLESP=1
end else begin
	GLESP=0
	npix = (size(Imag))[1]
	nside = npix2nside(npix)
end

Pyr=1   
if keyword_set(Undec) then Pyr=0
if keyword_set(atrou) then Pyr=0 

if not keyword_set(NbrScale) then NbrScale=4
if not keyword_set(NSigma) then NSigma=3.
if not keyword_set(NSigMask) then NSigMask=NSigma
if not keyword_set(WinMinSize ) then  WinMinSize = 8
MinBinSize = WinMinSize
bin_size = MinBinSize

if not keyword_set(Fdr) then AlphaFDR = 0 else AlphaFDR = fdr


if not keyword_set(FirstScale) then FirstScale=1

if Pyr EQ 1 then mrs_pwttrans, Imag, Trans, NbrScale=NbrScale, lmax=lmax $
else begin
   if keyword_set(atrou) then mrs_attrans, Imag, Trans, /healpix, NbrScale=NbrScale $
   else mrs_wttrans, Imag, Trans, NbrScale=NbrScale, lmax=lmax
end

TabThreshold = fltarr(NbrScale)
TabNoise = fltarr(NbrScale)

if keyword_set(FilterLast) then LS = NbrScale-1 $
else  LS = NbrScale-2

for j =0,LS do begin
    Scale = mrs_wtget(Trans,j,NormVal=NormVal)
    if GLESP EQ 1 then begin
       GlData = Scale
       Scale = Scale.T_sky
    end
     
    if j LT FirstScale-1 then Scale[*] = 0 $
    else BEGIN
      if keyword_set(localmad) then BEGIN
          ScaleData = H2F(Scale)
          nside =gnside(Scale)
          if not keyword_set(pyr) then bin_size = MinBinSize*2.^double(j) $
          else bin_size = MinBinSize
          if bin_size GT nside then bin_size = nside
          nb_bin = double(nside)/bin_size
          BlockNoise = fltarr(nb_bin, nb_bin, 12)     
          for f=0, 11 do $
          for ii=0, nb_bin-1 do $
          for jj=0, nb_bin-1 do $    
              BlockNoise[ii,jj, f] = mad(ScaleData[ii*bin_size:(ii+1)*bin_size-1, jj*bin_size:(jj+1)*bin_size-1, f])
          ScaleNoise = F2H(BlockNoise)
          ScaleNoise = mrs_resize(ScaleNoise, nside=nside)
          ; tvs, scalenoise
          ; help, nside, ScaleNoise
       END  ; if Loc_mad EQ 1
      if keyword_set(localmad) then print, "LOCAL MAD ", bin_size, nb_bin
      ; Noise standard deviation estimation from the first scale
      if j EQ 0 and not keyword_set(MAD) and not keyword_set(UseTabThreshold) and not keyword_set(SigmaNoise) and not keyword_set(localmad) then begin
         SigmaNoise=get_noise(Scale) / NormVal
       print, ' ==> Noise standard deviation estimation = ', SigmaNoise
      end
    
      ; Threshold level at scale j
      ThresholdLevel = 0.
      NS = 0.
      NSM = 0.
      if j eq 0 then NS=NSigma+1 else NS=NSigma
      if j eq 0 then NSM=NSigMask+1 else NSM=NSigMask
      if keyword_set(TabNSigma) then NS=TabNSigma[j]
      if not keyword_set(MAD) and not keyword_set(localmad) and not keyword_set(UseTabThreshold) then ThresholdLevel = NS*SigmaNoise*NormVal $
      else if not keyword_set(localmad) then begin
         if keyword_set(mask) then begin
            Ind = where( Mask EQ 1)
            ThresholdLevel = NS * mad(Scale[Ind]) 
          end else ThresholdLevel = NS * mad(Scale)
      end
      if keyword_set(UseTabThreshold) then ThresholdLevel = UseTabThreshold[j]
      TabThreshold[j] = ThresholdLevel
      TabNoise[j] = ThresholdLevel / NS
      ThresholdLevelMask = TabNoise[j] * NSM 
      
      if not keyword_set(fdr) then begin
         if not keyword_set(localmad) then begin
           ScaleIn = mrs_absthreshold(Scale, ThresholdLevel, soft=soft)
	       if keyword_set(mask) and (NSigMask NE NSigma) then begin
	           ScaleOut = mrs_absthreshold(Scale, ThresholdLevelMask, soft=soft)
	           Ind = where( Mask EQ 0)
	           ScaleIn[Ind] = ScaleOut[Ind]
 	       end   
	       Scale=ScaleIn
	     end else begin
	       if not keyword_set(pos) then Ind = where ( ABS(Scale) LT NS*ScaleNoise,c) $
	       else Ind = where (Scale LT NS*ScaleNoise,c)
	       print, sigma(Scale), sigma(ScaleNoise)
	       if c GT 0 then Scale[Ind] = 0
	     end
       end else if not keyword_set(Use_FdrAll) then begin
         coef = mrs_getpix(Scale)
         if not keyword_set(localmad) then coefnorm  = mrs_getpix(Scale) / TabNoise[j] $
         else coefnorm  = mrs_getpix(Scale) / scalenoise
         vs = size(coef)
         nx = vs[1]
         ny = vs[2]
         nel = N_elements(coef)
         Coef1D = reform(coefnorm, nel)
         if Alphafdr GE 1 then Alphafdr  = 0.05
         opt = ' -v  -a   ' + STRCOMPRESS(string(Alphafdr), /REMOVE_ALL) 
         print, "FDR : Alpha = ", Alphafdr
         mr_prog, 'im_fdr' , Coef1D, Res, opt=opt
         ind = where (Res EQ 0, c)
         if c gt 0 then  coef[ind] = 0
         if keyword_set(pos) then begin
             Ind = where ( coef LT 0,c)
             if c gt 0 then  coef[ind] = 0
         end
	     mrs_putpix, Scale, Coef
      end ; if not keyword_set(Use_FdrAll) 
    END
    if GLESP EQ 1 then begin
       GlData.T_sky = Scale
       Scale = GlData
    end
    if not keyword_set(Use_FdrAll) then mrs_wtput, Trans, Scale, j
    
 ;   AlphaFDR = AlphaFDR * 2.
 ;   if (AlphaFDR GT 0.1) then AlphaFDR = 0.1
endfor


if keyword_set(Use_FdrAll) then begin
   print, "FDR"
   coef = trans.coef
   for j =0,NbrScale-2 do  coef[*,j] = coef[*,j] / TabNoise[j]
   vs = size(coef)
   nx = vs[1]
   ny = vs[2]
   nel = N_elements(coef)
   Coef1D = reform(coef, nel)
   if fdr GE 1 then fdr = 0.05
   opt = ' -v  -a   ' + STRCOMPRESS(string(fdr), /REMOVE_ALL) 
   mr_prog, 'im_fdr' , Coef1D, Res, opt=opt
   ind = where (Res EQ 0, c)
   if c gt 0 then trans.coef[ind] = 0
end


if keyword_set(KillLastScale) then begin
    Scale = mrs_wtget(Trans,NbrScale-1)
    mrs_set, Scale, 0.
    mrs_wtput, Trans, Scale, NbrScale-1
end

if Pyr EQ 1 then mrs_pwtrec, Trans, Filter $
else begin
   if keyword_set(atrou) then mrs_atrec, Trans, Filter  $
   else mrs_wtrec,  Trans, Filter 
end


if keyword_set(pos) then mrs_pos, Filter

if keyword_set(niter) then begin
   LastScale = mrs_wtget(Trans,NbrScale-1)
   
   Lambda = 10.
   delta = Lambda / float(niter-1)
   TabEps = TabNoise / 2
   for i=1,niter do begin
      if Pyr EQ 1 then mrs_pwttrans, Filter, TransFil, NbrScale=NbrScale, lmax=lmax $
      else begin
        if keyword_set(atrou) then mrs_attrans, Filter, TransFil, /healpix, NbrScale=NbrScale $
        else mrs_wttrans, Filter, TransFil, NbrScale=NbrScale, lmax=lmax
      end

      for j =0,LS do begin
         Scale = mrs_wtget(Trans,j)
         ScaleFil = mrs_wtget(TransFil,j)
	 ScaleResi = mrs_diff(Scale,ScaleFil)
         if GLESP EQ 1 then begin
            GScale = Scale
            Scale = Scale.T_sky
	    GScaleFil = ScaleFil
	    ScaleFil= ScaleFil.T_sky
	    GScaleResi =  ScaleResi
	    ScaleResi = ScaleResi.T_sky
         end	 
	 
	 ind = where (Scale NE 0, c)
	 if c GT 0 then begin
	   ResiSignif = ScaleResi[ind]
	   indr = where (ABS(ResiSignif) LT TabEps[j], c)
	   if c GT 0 then ResiSignif[indr] = 0
	   ScaleFil[ind] = ScaleFil[ind] + ResiSignif
	 end
	   
	 ; Now soft thresholding
	 if Lambda GT 0 then softthreshold, ScaleFil, Lambda*TabNoise[j]
	 if GLESP EQ 1 then begin
	    GScaleFil.T_sky = ScaleFil
	    ScaleFil = GScaleFil
	 end
         mrs_wtput, TransFil, ScaleFil, j
      end
      mrs_wtput, TransFil, LastScale, NbrScale-1
      
      if Pyr EQ 1 then mrs_pwtrec, TransFil, Filter  $
      else begin
         if keyword_set(atrou) then mrs_atrec, TransFil, Filter  $
      else mrs_wtrec,  TransFil, Filter 
      end
 
      if keyword_set(pos) then mrs_pos, Filter
      
      resi = mrs_diff(Imag, Filter)
      ; tvs, resi 
      if GLESP EQ 1 then print, 'Iter ', i, '   ==>  SigmaResi = ', sigma(Resi.t_sky), '  Lambda = ', Lambda $
      else  print, 'Iter ', i, '   ==>  SigmaResi = ', sigma(Resi), '  Lambda = ', Lambda
        
      Lambda = Lambda - delta 
    end
end


DONE:

END


;==================================================================


pro mrs_wtfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, lmax=lmax, TabNSigma=TabNSigma, $
                  mad=mad, localmad=localmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale, Trans=Trans, Pyr=Pyr, niter=niter, UseTabThreshold=UseTabThreshold, $
		  pos=pos, cycle=cycle, FirstScale=FirstScale, Soft=Soft, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, atrou=atrou, OutMask=OutMask ; , NSigMask=NSigMask

 
if N_PARAMS() LT 2 or N_PARAMS() GE 3 then begin 
        print, 'CALLING SEQUENCE: mrs_wtfilter, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, lmax=lmax, TabNSigma=TabNSigma, mad=mad, localmad=localmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale, Trans=Trans, Pyr=Pyr, niter=niter, pos=pos, cycle=cycle, FirstScale=FirstScale, Soft=Soft, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, atrou=atrou, OutMask=OutMask'
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
if type_code(Imag) EQ 8 then begin 
	GLESP=1
	npix = (size(Imag.t_sky))[1]
	nx = Imag.nx
	np = Imag.np
	x_sky = Imag.x_sky
	y_sky = Imag.y_sky
end else begin
	GLESP=0
	npixel = (size(Imag))[1]
	nside = npix2nside(npixel)
end

if keyword_set(atrou) AND GLESP EQ 1 then begin
	print, 'Error: GlESP format is not supported by the "a trou" algorithm'
	goto , DONE
end

if keyword_set(localmad) AND GLESP EQ 1 then begin
	print, 'Error: Localmad keyword cannot be used with GlESP format images'
	goto , DONE
end

if keyword_set(cycle) AND GLESP EQ 1 then begin
	print, 'Error: Cycle keyword cannot be used with GlESP format images'
	goto , DONE
end

if not keyword_set(Pyr) then Undec=1 else Undec=0
if not keyword_set(FirstScale) then FirstScale=1
if not keyword_set(cycle) then begin
   mrs_wtfilter_one_cycle, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, TabNSigma=TabNSigma, $
                  mad=mad, localmad=localmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale,Trans=Trans, Undec=Undec, niter=niter, pos=pos, soft=soft, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, lmax=lmax, atrou=atrou, UseTabThreshold=UseTabThreshold 
end else begin

   mak_map,nside,interpole,t_interpol = 1
   rotate_map_nest,Imag,!dpi/4,0,0,Imag2
   rotate_map_nest,interpole,!dpi/4,0,0,interpole2
   rotate_map_nest,Imag,0,!dpi/2,0,Imag3
   rotate_map_nest,interpole,0,!dpi/2,0,interpole3

   mrs_wtfilter_one_cycle, Imag, Filter, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, TabNSigma=TabNSigma, $
                  mad=mad, localmad=lmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale,Trans=Trans, Undec=Undec, niter=niter, pos=pos, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, lmax=lmax, atrou=atrou, NSigMask=NSigMask, UseTabThreshold=UseTabThreshold

   mrs_wtfilter_one_cycle, Imag2, Filter2, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, TabNSigma=TabNSigma, $
                  mad=mad, localmad=lmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale,Trans=Trans, Undec=Undec, niter=niter, pos=pos, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, lmax=lmax, atrou=atrou, NSigMask=NSigMask, UseTabThreshold=UseTabThreshold

   mrs_wtfilter_one_cycle, Imag3, Filter3, NbrScale= NbrScale, NSigma=NSigma, SigmaNoise=SigmaNoise, TabNSigma=TabNSigma, $
                  mad=mad, localmad=lmad, WinMinSize=WinMinSize, KillLastScale=KillLastScale,Trans=Trans, Undec=Undec, niter=niter, pos=pos, fdr=fdr, Use_FdrAll=Use_FdrAll, FilterLast=FilterLast, mask=mask, lmax=lmax, atrou=atrou, NSigMask=NSigMask, UseTabThreshold=UseTabThreshold

   rotate_map_nest,Filter2,-!dpi/4,0,0,F2
   rotate_map_nest,Filter3,0,-!dpi/2,0,F3

   FilterFinal = (  (Filter * interpole) + (Filter2 * interpole2)+ (Filter3 * interpole3) ) / (interpole+interpole2+interpole3 )
   Filter = FilterFinal
end

if GLESP EQ 1 then ind = where( abs(Filter.t_sky) NE 0, c) else ind = where( abs(Filter) NE 0, c)

;if keyword_set (mask) then ind = where( abs(Filter) + (1-mask) NE 0, c)

if GLESP EQ 1 then OutMask = Imag.t_sky else OutMask = Imag

OutMask[*] = 1

if c GT 0 then OutMask[ind]=0

if keyword_set (mask) then OutMask = OutMask*mask

;resi = mrs_diff(Imag, Filter)
;OutMask  = mrs_getpix(resi)
;ind = where( abs(OutMask) NE 0, c)
;OutMask[*]=1
;if c GT 0 then OutMask[ind]=0

DONE:

END


pro killcmb, ima, filter, nsigma=nsigma, nbrscale=nbrscale, fdr=fdr, Use_FdrAll=Use_FdrAll
if not keyword_set(nsigma) then nsigma=4
print,'NSigma',nsigma
if not keyword_set(nbrscale) then nbrscale=6

mrs_wtfilter, ima, filter, /kill, /mad, nsigma=nsigma, nbrscale=nbrscale, fdr=fdr, Use_FdrAll=Use_FdrAll
tvs, filter
end


pro ttf, i, s, fs, fi
s = rims('SyncMap.fits')
i  = rims('FreeFreeMap.fits')
mrs_wtfilter, s, fs, niter=10, /mad, nbrscale=8, nsigma=4, firstscale=2, trans=trans, /pos
save, filename='xx_fs.xdr', fs
mrs_wtfilter, i, fi, niter=10, /mad, nbrscale=8, nsigma=4, firstscale=2, trans=trans, /pos
save, filename='xx_fi.xdr', fi
end





