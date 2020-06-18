;+
; NAME:
;        df_peaks
;
; PURPOSE:
;   Count the number of emission and absoption lines in a spectrum or a set of spectra
;
; CALLING:
;     TabPeak = df_peaks(DecSpectra,    noise=noise, nsigma=nsigma)
;
; INPUTS:
;     DecSpectra -- IDL structure: decomposed spetra derived from df_get_spectra_components
;    
; OUTPUTS:
;     TabPeak --   IDL array of structures : Tab[i] is a structure relative to the ith spectrum, with the following fields
;                                        abspeaks: number of absoption lines
;                                         empeaks:  number of emission lines
;                                         allpeaks: number of peaks (i.e. abspeaks + empeaks)
;
; INPUT KEYWORDS:
;  noise -- IDL array : RMS per pixel
;  nsigma  -- double :  peaks with a maximum amplitude small than nsigma*noise are not taken into account. Default value is 0.1.
;
; EXAMPLE:
;       Remove the baseline of all spectra
;               TabPeak = df_peaks(DecSpectra,  noise=RMSnoisePerPixe, nsigma=0.1) 
;         
; HISTORY:
;	Written: Daniel Machado & Jean-Luc Starck, 2013
;	Sept, 2013 File creation

;--------------------------------------------------------------------------------------------------------
 
function df_get_peaks,  NegLine, PosLine,  noise=noise, nsigma=nsigma, AmplitudeEmPeak=AmplitudeEmPeak, AmplitudeAbsPeak=AmplitudeAbsPeak, loc=loc
 
;create arrays to count peaks
gabs = ABS(NegLine)
gem = PosLine
InfoPeaks= {abspeaks: 0L, empeaks: 0L, allpeaks: 0L}

if not keyword_set(noise) then begin
   threshcutoffp=(5d/100d)                  ; user-defined cutoff for defining whether a tiny peak really is a peak or just a denoising artifact,
   threshcutoffn=(5d/100d)                  ;this should remain small. Here p and n stand for positive and negative domains of the spectrum.
   absthreshlevel = max(gabs)*threshcutoffn ;calculate an absorption theshold for a peak
   emthreshlevel = max(gem)*threshcutoffp   ;alculate an emission theshold for a peak
end else begin
   absthreshlevel = noise* nsigma
   emthreshlevel= noise* nsigma
end

allpeaks=0 
abspeaks=0
empeaks=0

if total(gabs+ gem) NE 0  and finite(total(gabs+gem)) NE 0  then  begin
   absd0 = gabs - shift(gabs,1)    ; peak-counting procedure, based on changes of slope,  doublets will count twice
   absd1 = gabs - shift(gabs,-1)   ; if they have two separate peaks but are otherwise blended.
   abspk = where(absd0 gt 0 and absd1 gt 0, absnpk)
  ; stop
   if absnpk  LE 0 then AmplitudeAbsPeak = 0 $
   else AmplitudeAbsPeak = gabs[abspk]
   emd0 = gem - shift(gem,1)
   emd1 = gem - shift(gem,-1)
   empk = where(emd0 gt 0 and emd1 gt 0, emnpk)
   if emnpk  LE 0 then AmplitudeEmPeak = 0 $
   else AmplitudeEmPeak = gem[empk]
   
   if absnpk gt 0 then begin
      absbigind = where(gabs[abspk] gt absthreshlevel[abspk], absnpk) ;count number of peaks higher than threshold.  
      if absnpk gt 0 then loc=[abspk[absbigind]]
   endif
   if emnpk gt 0 then begin
      embigind = where(gem[empk] gt emthreshlevel[empk], emnpk)         
      if emnpk gt 0 then if absnpk gt 0 then loc=[loc,empk[embigind]] else loc=[empk[embigind]]
   endif
     
   abspeaks= long(absnpk)            ;record the absorption peak count
   empeaks= long(emnpk)              ;record the emission peak count
   allpeaks = long(absnpk+emnpk)     ;combined count.
   InfoPeaks = {abspeaks: absnpk, empeaks: emnpk, allpeaks:  absnpk+emnpk}
   peakloc = loc[sort(loc)]
end
return, InfoPeaks
end  

;====================================================================================

; TabOk is interesting to study the peak distribution on simulations for gal where the redshift estimation is not correct
function df_peaks, DecSpectra, TabOk=TabOk, noise=noise, nsigma=nsigma, peakloc=peakloc
TabPeak=-1
  
if N_PARAMS() LT 1  then begin 
   Dec =-1
   print, 'CALLING SEQUENCE: TabPeak = df_peaks(DecSpectra, noise=noise, nsigma=nsigma)'
   goto, DONE
end

if not keyword_set(nsigma) then nsigma = 0.01

vs= size(DecSpectra.POSLINE)
NSpectra = vs[2]

for i=0,NSpectra-1 do begin
                                ; if keyword_set(TabOK)  then print, 'Gal ', i, ' ==> ', TabOk[i]
   GalPeaks = df_get_peaks(DecSpectra.NegLine[*,i],  DecSpectra.PosLine[*,i],   noise=noise,  nsigma=nsigma)
   
   if i EQ 0 then TabPeak = replicate(GalPeaks, NSpectra) $
   else TabPeak[i] = GalPeaks
   if keyword_set(TabOK) then if TabOK[i] EQ 0   then begin
      allpeaks = InfoPeask.allpeaks
      empeaks  = InfoPeask.empeaks
      abspeaks = InfoPeask.abspeaks
      print, 'Gal ', i+1, ' ==> ', STRC(TabOk[i]), ", #absPeak = ", STRC(abspeaks), ", #emPeak = ",STRC( empeaks), ", #Peak = ", STRC(allpeaks)
   endif
endfor

DONE: 

return,    TabPeak
end

;====================================================================================



