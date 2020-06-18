
pro test_darth_fader, first=first
;EstimRedshiftFromNoiseData = 0
;lstep=2.17e-4

rms=0
print,'Getting inputs..'
ins=df_input_params()
print,'Done!'

;DIR = '~/Dropbox/DanielMachado/Code/'
if keyword_set(first) then begin

   print,'Reading data..'
   data = readfits(ins.indir+ins.incat+'.gz')
   vs = size(data)
   Nx = vs[1]
   Training = readfits(ins.indir+ins.trainingcat+'.gz')
   ; Tred = readfits(ins.indir+ins.truecat)
   Tred =readfits('~/Dropbox/DanielMachado/Code/trueredshifts.fits')
   ngals = n_elements(Tred)

   if strmatch(ins.noisetype,'lambda') then RMS  = readfits(ins.indir+ins.rmscurve+'.gz') else RMS=0.
   print,'Done!'
   
   print,'Getting spectra components..'
   DecSpectra = df_get_spectra_components(Data, RMS=RMS)
   print,'Done!'
   fdrvals = [2.0,0.5,0.8,1.0,1.2,1.5,1.8,2.2,2.5,2.8,3.0]
   b=10
   fdr = fdrvals[10]
   print, 'FDR = ', fdr   ; ==> FDR = 3 and not 2 !!!!
   namesave2 = string(double(fdr),Format='(F8.1)')
   optstring1= '-M -f3 -n6 -t11 -k -K -P -i20 -s4'
   optstring2='-M -f3 -n6 -k -K -i20 -C2 ' + strcompress('-s' + namesave2,/remove_all)
    ; mrfiltermod, data, filtered_zs, noise=rms, OPT =optstring1
    mr1d_filter, data, filtered_zs, rmsnoise= rms, OPT=optstring1
    ;  fd = readfits('~/Dropbox/DanielMachado/Code/filtered_test.fits'.gz)
    ;    mrfiltermod, testing, filtered_zs, OPT = optstring1    
    ; ==> OK, fd =    filtered_zs
    istherespec = total(abs(filtered_zs),1,/double)
    zsnd = data-filtered_zs        
    continuum=0d * zsnd
    nscales=7
    for p=0, ngals-1 do continuum[*,p] = (star1d(zSNd[*,p],nscale=nscales))[*,(nscales-1)]
   ; cd  = readfits('~/Dropbox/DanielMachado/Code/continuum.fits.gz')
    zs = data-continuum   ; continnum is ALWAYS equal to zerp.
    mr1d_filter,zs, cleanspec, rmsnoise =errorcurve, OPT = strcompress(optstring2 + ' -P')
    cs = readfits('~/Dropbox/DanielMachado/Code/cleanspec.fits.gz')
    mr1d_filter,zs, emcleanspec, rmsnoise =errorcurve, OPT = strcompress(optstring2 + ' -p')
    em = readfits('~/Dropbox/DanielMachado/Code/emspec.fits.gz')
    mr1d_filter, -zs, abscleanspec, rmsnoise =errorcurve, OPT = strcompress(optstring2 + ' -p')
    as = readfits(~/Dropbox/DanielMachado/Code/absspec.fits.gz')
   g = cleanspec;/noisescaler
   gabs = abscleanspec;/noisescaler
   gem = emcleanspec;/noisescaler
  

  
   print,'Removing baseline from training set..'
   nscale=6
   Training_noBaseline = del_baseline(Training, nscale=nscale)
   
   optstring1= '-M -f3 -n6 -t11 -k -K -P -i20 -s4'
   ; mrfiltermod, training, S, OPT= optstring1    
   mr1d_filter, training, Sn, OPT= optstring1
   ;  td = readfits('~/Dropbox/DanielMachado/Code/training_filtered.fits.gz')   ; ==> OK td and S, and Sn are the same
   
   TemplateS = df_get_eigenvector(Sn,  AllEigen=AllEigen)
   Template = TemplateS[*, 0:19]
   ;  ed = readfits('~/Dropbox/DanielMachado/Code/eigentemplates.fits')

   print,'Done!'
   print,'Running PCA to compute eigenvectors..'
   ; Template = df_get_eigenvector(Training_noBaseline,  AllEigen=AllEigen)
   ntrain = (size(training))[2]  
   trains = ntrain-1
   tbins = (size(training))[1]


;   NbrTraining = (size(Training_noBaseline))[2]
;   T = Training_noBaseline
;   for i = 0, NbrTemp-1 do begin T[*,i] = ( T[*,i] - mean(T[*,i] ) ) / total( T[*,i]^2 )
;   mr_prog, "mr1d_gmca", T, GMCA_template, opt='-S10'   ; to get 20 templates using GMCA
;   print,'Done!'
   
   save, filename='xxjlstest.xdr', DecSpectra, data, training, Template, Tred, rms, ins
end else restore, 'xxjlstest.xdr', /verb

;stop
NbrTest = (size(data))[2]
NbrTemp = (size(Template))[2]

print,'Estimating redshifts..'
; if ins.EstimRedshiftFromNoiseData EQ 1 then EstRed = df_get_redshift(data,  Template,lstep = ins.lstep, diffpar = ins.diffpar) $
; else EstRed = df_get_redshift(DecSpectra.DATA_CLEAN[*,0: NbrTest-1],  Template,lstep = ins.lstep, diffpar = ins.diffpar)
;EstRed = dblarr(NbrTest)
;for i=0, NbrTest-1 do begin
;    if ins.EstimRedshiftFromNoiseData EQ 1 then EstRed[i]  = df_get_redshift1(data[*,i],  Template,lstep = ins.lstep, diffpar = ins.diffpar) $
;    else EstRed[i] = df_get_redshift1(DecSpectra.DATA_CLEAN[*,i],  Template,lstep = ins.lstep, diffpar = ins.diffpar)
;end

;if ins.EstimRedshiftFromNoiseData EQ 1 then GMCAEstRed = df_get_redshift(data,  GMCA_template,lstep = ins.lstep, diffpar = ins.diffpar) $
;else GMCAEstRed = df_get_redshift(DecSpectra.DATA_CLEAN[*,0: NbrTest-1],  GMCA_template,lstep = ins.lstep, diffpar = ins.diffpar)
;T = Template
;for i = 0, NbrTemp-1 do begin T[*,i] = ( T[*,i] - mean(T[*,i] ) ) / total( T[*,i]^2 )
;D = DecSpectra.DATA_CLEAN
;for i = 0, NbrTemp-1 do begin D[*,i] = ( D[*,i] - mean(D[*,i] ) ) / RMS^2
;EstRed = df_get_redshift(D,  T, lstep = ins.lstep, diffpar = ins.diffpar) 

print,'Done!'
; DecSpectra1 = df_get_spectra_components(Data, RMS=RMS, RegulParam=0.1)
; EstRed = df_get_redshift(DecSpectra.DATA_CLEAN[*,0: NbrTest-1],  Template,lstep = ins.lstep, diffpar = ins.diffpar)
EstRed = df_get_redshift(DecSpectra.DATA_CLEAN[*,0: NbrTest-1],  Template,lstep = ins.lstep, shiftpar = ins.SHIFTPAR)



diff =  Tred[0: NbrTest-1] - EstRed

; i=17
diffpar = ins.SHIFTPAR
lstep = ins. lstep
pixtol=4
delz = double((10d ^ (pixtol * Lstep)) - 1d )
ind = where(abs(diff) GT delz, c)
print, "Percentage of outliers = ", c / double(NbrTest) * 100.


; a =  df_sdss_get_redshift(DecSpectra.DATA_CLEAN[*,i],  Template, lstep= lstep,  diffpar=diffpar, rms=rms, tred=tred[i])

goto, DONE

; plot, Tred[0: NbrTest-1] ,  EstRed, psym=1, xr=[0,2], yr=[0,2]
; diff =  Tred[0: NbrTest-1] - EstRed
; ind = where(abs(diff) GT ins.delz, c)
;indOK =  where(abs(diff) LE ins.delz, c)
;TabOK = intarr(NbrTest) + 1
;TabOK[ind] = 0
;df_peaks_allgal, DecSpectra1, TabOk=TabOk, TabPeak, TabPosPeak, TabNegPeak, MaxPeak_BadGal= MaxPeak_BadGal, noise= rms, nsigma=nsigma

;TabPeakPB = TabPeak[ind]
;TabPeakOK = TabPeak[indOK]
;TabMaxChiSquare = TabMaxChiSquare  / stddev(TabMaxChiSquare) * 100.

;TabMPB = TabMaxChiSquare[ind]
;info, TabMPB
;TabMOK = TabMaxChiSquare[indOK]
;info, TabMOK

;print, "GMCA Percentage of outliers = ", c / double(NbrTest) * 100.
;plot, Tred[0: NbrTest-1] ,  GMCAEstRed, psym=1, xr=[0,2], yr=[0,2]




;ind1 = where(abs(diff[0: NbrTest-1] ) GT 0.005*(1+ Tred[0: NbrTest-1] ) , c1)
;print, "Percentage of outliers ( Errz > 0.005*(1+z)) = ", c1 / double(NbrTest) * 100.

;save, filename='xxjlsest.xdr', NbrTest, EstRed, ind, ind1, /verb
;restore, /verb, 'xxjlsest.xdr'
stop

DONE:
end


pro ttx


df_get_redshift1, g,  Training, lstep= lstep, chisquare=chisquare


 gi = 2002
 gd = DecSpectra.DATA_NOBASELINE[*,gi]
 print, Tred[gi], '  ', EstRed[gi]
 ;  g = data[*, gi]
 g = DecSpectra.DATA_CLEAN[*, gi]
 d1 = g
Dec = df_get_spectra_components(D1, RMS=RMS)
print, 'Gal ', gi, ' ==> ', TabOk[gi]
df_peaks,  Dec, allpeaks=allpeaks, empeaks=empeaks, absnpk=absnpk

gg= [[g],[g]]

CleanSpectrum = DecSpectra.DATA_CLEAN[*,gi]
z1 = df_get_redshift(g,  Template, lstep = ins.lstep, diffpar = ins.diffpar, chisquare=c)
z = df_get_redshift(gg,  Template, lstep = ins.lstep, diffpar = ins.diffpar)
print,  Tred[gi], '  ', EstRed[gi], ' ' , z[0],  z1

EstRed1 = dblarr(NbrTest)
for i=0, NbrTest-1 do  EstRed1[i] = df_get_redshift(DecSpectra.DATA_CLEAN[*,i],  Template,lstep = ins.lstep, diffpar = ins.diffpar)
plot, Tred[0: NbrTest-1] ,  EstRed1, psym=1, xr=[0,2], yr=[0,2]
diff1 =  Tred[0: NbrTest-1] - EstRed1
indx = where(abs(diff1) GT ins.delz, c)
print, "Percentage of outliers = ", c / double(NbrTest) * 100.
 print, Tred[gi], '  ', EstRed[gi], EstRed1[gi]

  
EstRed = df_get_redshift(Data[*, 0: NbrTest-1],  Template)

print, EstRed
plot, Tred[0: NbrTest] - EstRed, psym=1
plot, Tred[0: NbrTest] ,  EstRed, psym=1, xr=[0,2], yr=[0,2]
diff =  Tred[0: NbrTest] - EstRed
ind = where(abs(diff) GT 1, c)
print, "Percentail of outloiers = ", c / double(NbrTest) * 100.

end


