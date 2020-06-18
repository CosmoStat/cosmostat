;+
; NAME:
;        darth_fader.pro
;
; PURPOSE:
;    Run the Darth Fader pipeline
;
; INPUTS: 
;    Input options are specified in df_input_params.pro
;
; OUTPUTS: 
;    
;    Eigentemplates computed from a training set
;    Components of test spectra (continua and line features)
;    Redshift estimates for test galaxies
;    Catalogue cleaned according to feature counting criterion
;    Table of peak counts for each galaxy as S_TabPeaks.allpeaks
;    S_TabPeaks.empeaks and S_TabPeaks.abspeaks. 
;    Input data
;    Input training set
;    Input training set after continuum-subtraction
; 
; HISTORY: 
;       Written: Adrienne Leonard, Oct 2013
;-

pro darth_fader, DecSpectra, template, estred, clean_catalogue, verbose=verbose, S_TabPeaks=S_TabPeaks, Data=Data, Training=Training, TnoBaseline=TnoBaseline, iterz=iterz, flag = flag

;; Read in the data
if keyword_set(verbose) then print, 'Reading data...'
input = df_input_params()

if strmatch(input.incat,'*.fits') or strmatch(input.incat,'*.fits.gz') then $
   data = readfits(input.indir + input.incat) else message,'Input catalogue not specified!'
if strmatch(input.trainingcat,'*.fits')  or  strmatch(input.trainingcat,'*.fits.gz') then $
   training = readfits(input.indir + input.trainingcat) else $
if strmatch(input.templatecat,'*.fits') eq 0 and strmatch(input.templatecat,'*.fits.gz') eq 0 then $
   message,'No templates or training set specified!'

if strmatch(input.templatecat,'*.fits') or  strmatch(input.templatecat,'*.fits.gz') then $
   template = readfits(input.indir + input.templatecat)
if strmatch(input.rmscurve,'*.fits') or  strmatch(input.rmscurve,'*.fits.gz')  then $
   rms = readfits(input.indir + input.rmscurve) else begin
   rms = 0.
   print,'RMS curve not set; assuming white Gaussian noise'
   rms = get_noise(data)
endelse
if keyword_set(verbose) then print,'done!'
;stop
;; Compute eigentemplates if required 
if input.ComputeEigenTemplates eq 1 then begin
   if keyword_set(verbose) then print, 'Computing eigentemplates from training data...'
   Training_noBaseline = df_del_baseline(Training, nscale = input.nscale)
   TnoBaseline = Training_noBaseline
   Template = df_get_eigenvector(Training_noBaseline, AllEigen = AllEigen, EnergPercent = input.EnergPercent, $
                                 ntemp = input.Ntemplates)
   if input.OutputTemplates eq 1 then writefits, input.outdir + input.TempFileName, Template
   if input.OutputAllEigen eq 1 then writefits, input.outdir + input.AllEigenName, AllEigen
   if keyword_set(verbose) then print, 'done!'
endif else if strmatch(input.templatecat,'*.fits') eq 0 and strmatch(input.templatecat,'*.fits.gz') eq 0 $
then message, 'No templates specified!'

;; Decompose test spectra into components
if keyword_set(verbose) then print, 'Getting spectra components...'
if input.ComputeComponents eq 0 then begin
   dum = readfits(input.indir+input.componentcat)
   DecSpectra={baseline:dum[*,*,0], posline:dum[*,*,1], negline:dum[*,*,2], noise:dum[*,*,3], data_nobaseline:dum[*,*,4], data_clean:dum[*,*,5]}
endif else begin
   DecSpectra = df_get_spectra_components(data, rms=rms, nsigma=input.fdr, gen2starlet = input.gen2,regulfilter=input.regul,regulparam=input.regulparam)
endelse
if keyword_set(verbose) then print, 'done!'

if input.OutputComponents eq 1 then begin
   outs = [[[DecSpectra.baseline]],[[DecSpectra.posline]],[[DecSpectra.negline]],[[DecSpectra.noise]],[[DecSpectra.data_nobaseline]], [[DecSpectra.data_clean]]]
   writefits, input.outdir + input.CompFileName, outs
endif 

;; Compute redshifts by cross-correlation
if keyword_set(verbose) then print, 'Computing redshifts...'
if input.estred eq 0 then estred = readfits(input.indir+input.zcat) else begin   
   if input.EstimRedshiftFromNoisyData eq 1 then EstRed = df_get_redshift(DecSpectra.data_nobaseline, Template, $
                                                                          lstep = input.lstep, shiftpar = input.shiftpar) $
   else Estred = df_get_redshift(DecSpectra.data_clean, Template, lstep = input.lstep, shiftpar = input.shiftpar)
endelse
if keyword_set(verbose) then print, 'done!'
   
if input.OutputRedshifts eq 1 then writefits, input.outdir + input.RedshiftFileName, Estred

;; Clean catalogue, if required
if input.CleanCatalogue eq 1 then begin
   if keyword_set(verbose) then print, 'Cleaning catalogue...'
   if input.Linematch eq 0 then begin
      S_TabPeaks = df_peaks(DecSpectra,  noise=rms, Nsigma=input.NsigmaPeak)
      TabPeak = S_TabPeaks.allpeaks
      
      incl = where(TabPeak ge input.NpeakCrit, count)
      if input.OutputCleanCatalogue eq 1 then begin
         clean_catalogue = [[incl],[Estred[incl]]]
         writefits, input.outdir+input.CleanFileName, clean_catalogue
      endif 
      clean_catalogue={indices:incl, redshifts:estred[incl]}
   endif else begin
      clean_catalogue = df_match_lines(DecSpectra, estred, rms, input)
   endelse
 ; help,clean_catalogue,/str
 ;  stop
   if keyword_set(verbose) then print, 'done!'
endif

return
end
