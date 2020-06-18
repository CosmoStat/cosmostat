;+
; NAME:
;        df_get_spectra_components
;
; PURPOSE:
;   Decompose a spectrum or an array of spectra into several components (i.e. emission lines, absorption lines, baseline, noise)
;
; CALLING:
;     DecSpectra =  df_get_spectra_components(TabSpectra, RMS=RMS, nscale=nscale, nsigma=nsigma)
;
; INPUTS:
;     TabSpectra -- IDL 1D or 2D array :  TabSpectra[*, 0:N-1]   N input spectra
;    
; OUTPUTS:
;     DecSpectra -- IDL structure : decomposition of all specrta
;                                   BASELINE   - IDL Array[*, 0:N-1]   :   baseline or continuum
;                                   POSLINE    - IDLArray[*, 0:N-1]     :   emission lines
;                                   NEGLINE   - IDLArray[*, 0:N-1]     : absorption lines
;                                   NOISE   - IDLArray[*, 0:N-1]         : noise[*,i]  = TabSpectra[*,i] - (BASELINE[*,i]+ POSLINE[*,i] + NEGLINE[*,i])
;                                   DATA_NOBASELINE   - IDLArray[*, 0:N-1] :   data_nobaseline[*,i]  = TabSpectra[*,i] - (BASELINE[*,i]
;                                   DATA_CLEAN   - IDLArray[*, 0:N-1]  :  data_clean[*,i] = POSLINE[*,i] + NEGLINE[*,i]
;
; INPUT KEYWORDS:
;  RMS -- IDL 1D array:  noise rms per pixel
;  nscale -- scalar : number of scales used in the decomposition. Default is 6.
;  nsigma -- double:  detection level in the FDR. nsigma=2 by default, corresponding to alpha_FDR = 0.05, i.e. 5% of detections may be false detections.
;
; EXAMPLE:
;       Decompose all spectra in  TabSpectra, assuming the noise per pixel is given by RMS[0:Npix]
;               DecSpectra =  df_get_spectra_components(TabSpectra, RMS=RMS) 
;         
; HISTORY:
;	Written: Daniel Machado & Jean-Luc Starck, 2013
;	Sept, 2013 File creation

;--------------------------------------------------------------------------------------------------------
 
function  df_get_spectra_components, Spectra, RMS=RMS, nscale=nscale, nsigma=nsigma, RegulParam=RegulParam, Gen2Starlet=Gen2Starlet, $
                                     RegulFilter=RegulFilter

if N_PARAMS() LT 1  then begin 
        Dec =-1
        print, 'CALLING SEQUENCE: DecSpectra =  df_get_spectra_components(TabSpectra, RMS=RMS, nscale=nscale, nsigma=nsigma)'
        goto, DONE
        end

input_data =   Spectra
vs = size(input_data)
if vs[0] EQ 1 then nspecs=1 $
else nspecs =vs[2]
Nx = vs[1]

errorcurve=0
if keyword_set(RMS) then errorcurve=RMS
if not keyword_set(nsigma) then nsigma =2.   ; detection level in the FDR; nsigma=2 ==> alpha_FDR = 0.05, i.e. 5% of detections may be false detections.
if not keyword_set(RegulParam) then RegulParam = 0.05
if not keyword_set(nscale) then nscale=fix(alog(Nx)) - 1

;; First remove the continuum from the spectra
Contsubdata  = df_del_baseline(Spectra,  nscale=nscale , rms=rms,  baseline=continua)

optstring2 =  ' -n ' +  STRC(nscale) + ' -K -k -C2 -i20  -s'  + STRC(nsigma) + ' '

if keyword_set(Gen2Starlet) then optstring2 = optstring2 + ' -t9 ' $
else optstring2 = optstring2 + ' -t3 '
if not keyword_set(RegulFilter) then  optstring2 = optstring2 + ' -f3 ' $
else optstring2 = optstring2 + ' -f4 -G'  +  STRC(RegulParam) + ' '

optstring2 = optstring2 + ' -p ' 
;; Now denoise the spectra using a positivity constraint to obtain
;; emission lines
mr1d_filter, contsubdata, g1, rmsnoise=errorcurve, OPT= optstring2
;; And now denoise the spectra using a negativity constraint to obtain
;; absorption lines
mr1d_filter, -contsubdata, g2, rmsnoise=errorcurve, OPT= optstring2
g2 = -g2
;; The denoised spectrum is the sum of these two components
g = g1 + g2

Noise = input_data - continua - g
Dec = {baseline: continua, posline: g1, negline: g2, noise: Noise, data_nobaseline: contsubdata, data_clean: g}

;stop

DONE:

return, Dec

END



