;+
; NAME:
;   extract_object
;
; PURPOSE:
;
;   Performs all object extraction tasks
;      0) Locate bright fibers, and test image background
;      1) 4 Step Optimal extraction
;      2) Tweak to sky lines
;      3) Sky subtraction
;      4) Flux calibration
;      5) Telluric correction
;
; CALLING SEQUENCE:
;   extract_object, outname, objhdr, image, invvar, plugsort, wset, $
;    xarc, lambda, xtrace, fflat, fibermask, proftype=, color=, $
;    [ widthset=, dispset=, skylinefile=, plottitle=, superflatset=, $
;    /do_telluric ]
;
; INPUTS:
;   outname    - Name of outputs FITS file
;   objhdr     - Header cards from object image
;   image      - Object image [nx,ny]
;   invvar     - Inverse Variance of object image [nx,ny]
;   plugsort   - Plugmap structure for [ntrace] spectra
;   wset       - Wavelength solution from arc line spectra
;   xarc       - centroids measured in arc line spectra
;   lambda     - air wavelengths corresponding to xarc
;   xtrace     - spatial traces from flat field
;   fflat      - 1d flat field vectors
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   proftype   - Which type of profile should we use, (default=1 Gaussian)
;   superflatset- If present, then divide by median superflat! ???
;   color      - ???
;   widthset   - ???
;   dispset    - ???
;   skylinefile- ???
;
; REQUIRED KEYWORDS:
;   color      - camera color (red or blue)
;
; OPTIONAL KEYWORDS:
;   plottitle  - Prefix for titles in QA plots.
;   do_telluric- If set, then perform telluric-corrections for the red CCDs;
;                the default is to no longer do this, because the v5 code
;                does this in the later flux-calibration steps
;
; OUTPUTS:
;   A fits file is output in outname, which contains
;      FLOAT flambda [NX,NTRACE]
;      FLOAT flambda_invvar [NX,NTRACE]
;      LONG finalmask [NX,NTRACE]
;      STRUCT vacuum wavelengths
;      STRUCT wavelength dispersion
;      STRUCT plugmap [NTRACE]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   calcscatimage()
;   divideflat
;   djs_median()
;   djs_oplot
;   djs_plot
;   doscatter()
;   extract_boxcar()
;   extract_image
;   fibermask_bits()
;   fitsn()
;   fitvacset()
;   get_tai
;   heliocentric()
;   locateskylines
;   mwrfits
;   pixelmask_bits()
;   qaplot_scatlight
;   qaplot_skydev
;   qaplot_skyline
;   qaplot_skyshift
;   qaplot_skysub
;   skyline_dispersion()
;   skysubtract
;   spadd_guiderinfo
;   splog
;   sxaddpar
;   sxpar()
;   telluric_corr
;   traceset2xy
;   find_whopping()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-Jan-2000  Written by S. Burles, Chicago
;   26-Jul-2001  Also pass proftype and superflatset
;-
;------------------------------------------------------------------------------
pro extract_object, outname, objhdr, image, invvar, plugsort, wset, $
 xarc, lambda, xtrace, fflat, fibermask, color=color, proftype=proftype, $
 widthset=widthset, dispset=dispset, skylinefile=skylinefile, $
 plottitle=plottitle, superflatset=superflatset, do_telluric=do_telluric

   objname = strtrim(sxpar(objhdr,'OBJFILE'),2) 
   flavor  = strtrim(sxpar(objhdr,'FLAVOR'),2) 
   camera  = strtrim(sxpar(objhdr,'CAMERAS'),2) 
 
   ;------------------
   ; Identify very bright objects
   ; Do a boxcar extraction, and look for fibers where the median
   ; counts are 10000 ADU per row.

   fextract = extract_boxcar(image*(invvar GT 0), xtrace)
   scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
   whopping = find_whopping(scrunch, 10000.0, whopct)
   scrunch_sort = sort(scrunch)
   i5 = n_elements(scrunch)/20
   i95 = i5 * 19

   splog, 'Whopping fibers: ', whopping
   splog, 'Median counts in all fibers = ', djs_median(scrunch)
   splog, 'Number of bright fibers = ', whopct

   ; Assume that all fiber-mask bits are fatal for selecting sky fibers???
   iskies = where(strtrim(plugsort.objtype,2) EQ 'SKY' $
    AND plugsort.fiberid GT 0 AND (fibermask EQ 0), nskies)

   if (nskies LT 2) then begin
      splog, 'ABORT: Only '+ string(nskies) + ' sky fibers found' 
      return
   endif 

   skymedian = djs_median(scrunch[iskies])
   splog, 'Sky fiber median '+string(skymedian)
;   if (skymedian GT 2000) then begin
;      splog, 'ABORT: Median sky flux is brighter than 2000 e-'
;      return
;   endif

   splog, '5% and 95% count levels ', scrunch[scrunch_sort[i5]], $
                                      scrunch[scrunch_sort[i95]]

   if (whopct GT 20) then begin
      splog, 'WARNING: Disable whopping terms ' + objname
      whopping = -1
      whopct = 0
   endif

   ;------------------------------------------------------------
   ;  Check for bad pixels within 3 pixels of trace

   badcheck = extract_boxcar((invvar LE 0), xtrace, radius=2.5)
   badplace = where(badcheck GT 0)

   nx = (size(fextract,/dim))[0] 
   ny = (size(fextract,/dim))[1] 
   pixelmask = lonarr(nx,ny)

   badcolumns = where(total(badcheck GT 0,1) GT 0.1 * nx)

   if (badplace[0] NE -1) then pixelmask[badplace] = $
                pixelmask[badplace] OR pixelmask_bits('NEARBADPIXEL')

   if (badcolumns[0] NE -1) then fibermask[badcolumns] = $
                fibermask[badcolumns] OR pixelmask_bits('MANYBADCOLUMNS')

   if (whopping[0] NE -1) then begin
      ; Set the mask bit for whopping fibers themselves
      fibermask[whopping] = fibermask[whopping] OR pixelmask_bits('WHOPPER')

      ; Set the mask bit for fibers near whopping fibers, excluding the
      ; whopping fibers themselves.  Note that a fiber could still have both
      ; WHOPPER and NEARWHOPPER set if it is both bright and near another
      ; bright fiber.
      wp = [whopping - 2 , whopping -1, whopping+1 , whopping+2]
      wp = wp[ where(wp GE 0 AND wp LT ny) ]
      fibermask[wp] = fibermask[wp] OR pixelmask_bits('NEARWHOPPER')
   endif

   ;-----------------------------------------------------------------------
   ;  This is a kludge to fix first and last column ???
   ;-----------------------------------------------------------------------
   image[0,*] = image[0,*]*0.7
   image[2047,*] = image[2047,*]*0.7

   ;
   ;  First we should attempt to shift trace to object flexure
   ;

   xnow = match_trace(image, invvar, xtrace)
   bestlag = median(xnow-xtrace)

   splog, 'Shifting traces by match_trace ', bestlag

   if (abs(bestlag) GT 1.0) then begin
      splog, 'WARNING: pixel shift is large!'
   endif

   highrej = 10  ; just for first extraction steps
   lowrej = 10   ; just for first extraction steps
                 ; We need to check npoly with new scattered light backgrounds
   npoly = 16 ; maybe more structure, lots of structure
   nrow = (size(image))[2]
   yrow = lindgen(nrow) 
   nfirst = n_elements(yrow)

   splog, 'Extracting frame '+objname+' with 4 step process'

   traceset2xy, widthset, xx, sigma2
   ntrace = (size(sigma2,/dimens))[1]
   wfixed = [1,1] ; Fit gaussian height + width (fixed center position)
   nterms = n_elements(wfixed)

   splog, 'Step 1: Initial Object extraction'

   extract_image, image, invvar, xnow, sigma2, flux0, flux0ivar, $
    proftype=proftype, wfixed=wfixed, yrow=yrow, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
    ansimage=ansimage, chisq=firstchisq, ymodel=ymodel0, /relative

   ; (1a) Calculate scattered light, just for curiosity
   splog, 'Step 2: Find scattered light image'
   scatfit0 = calcscatimage(ansimage[ntrace*nterms:*,*], yrow, nscatbkpts=npoly)

   qaplot_scatlight, scatfit0, yrow, $
    wset=wset, xcen=xtrace, fibermask=fibermask, $
    title=plottitle+'Scattered Light before halo removal on '+objname

   ; (3) Calculate per-camera modelled scattering surface
   splog, 'Step 3: Calculate halo image'
   scatter = doscatter(camera, ymodel0, invvar, wset, sigma=0.9)

   ; (4) Re-extract with scattering model
   ; subtracted in order to measure
   ;     any remaining scattered light
   splog, 'Step 4: Extracting with halos removed'
   image_sub = image - scatter

   extract_image, image_sub, invvar, xnow, sigma2, tempflux, tempfluxivar, $
    proftype=proftype, wfixed=wfixed, yrow=yrow, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
    ansimage=ansimage2, chisq=secondchisq, ymodel=tempymodel, /relative

   ; (4) Calculate remaining scattered light
   splog, 'Step 4: Find scattered light image'
   scatfit = calcscatimage(ansimage2[ntrace*nterms:*,*], yrow, nscatbkpts=npoly)

   qaplot_scatlight, scatfit, yrow, $
    wset=wset, xcen=xtrace, fibermask=fibermask, $
    title=plottitle+'Scattered Light on '+objname

   image_sub2 = image_sub - scatfit

   ;-----------------------------------------------------------------------
   ;  Now, subtract halo image and do final extraction with all rows
   ;-----------------------------------------------------------------------
   ; (6) Final extraction
   splog, 'Step 6: Final Object extraction'

   highrej = 8
   lowrej = 5
;   wfixed = [1,1]
;   reject = [0.2,0.6,0.6]
   wfixed = [1] ; Fit to height only (fixed width + center)
   nterms = n_elements(wfixed)
   reject = [0.2,0.2,0.2]
   npoly = 0

   extract_image, image_sub2, invvar, xnow, sigma2, flux, fluxivar,$
    proftype=proftype, wfixed=wfixed, ansimage=ansimage3, $
    highrej=highrej, lowrej=lowrej, npoly=npoly, whopping=whopping, $
    chisq=chisq, ymodel=ymodel, pixelmask=pixelmask, reject=reject, /relative

   ;----------------------------------------------------------------------
   ; Can we find cosmic rays by looking for outlandish ansimage ratios???
   ;
   ; a = where(ansimage[lindgen(ntrace)*nterms, *] LT $
   ;           (-2*ansimage[lindgen(ntrace)*nterms+1, *])

   ;------------------
   ; QA chisq plot for fit calculated in extract image (make QAPLOT ???)

   xaxis = lindgen(n_elements(chisq)) + 1
   djs_plot, xaxis, chisq, $
    xrange=[0,N_elements(chisq)], xstyle=1, $
    xtitle='Row number',  ytitle = '\chi^2', $
    title=plottitle+'Extraction chi^2 for '+objname

   djs_oplot, !x.crange, [1,1]
   djs_oplot, yrow, secondchisq[yrow], color='blue'
   djs_oplot, yrow, firstchisq[yrow], color='green'

   xyouts, 100, 0.05*!y.crange[0]+0.95*!y.crange[1], $
            'BLACK = Final chisq extraction'
   xyouts, 100, 0.08*!y.crange[0]+0.92*!y.crange[1], $
            'BLUE = Initial-scattered chisq extraction'
   xyouts, 100, 0.08*!y.crange[0]+0.89*!y.crange[1], $
            'GREEN = Initial chisq extraction'

   ;------------------
   ; Flat-field the extracted object fibers with the global flat

   divideflat, flux, invvar=fluxivar, fflat
 
   pixelmask = pixelmask OR ((fflat LT 0.5) * pixelmask_bits('LOWFLAT'))

   ;------------------
   ; Look for pixels where scattered light is dominating

   scatteredpix = where(extract_boxcar(scatfit, xnow, radius=2.0) GT 2.0 * flux)
   if (scatteredpix[0] NE -1) then pixelmask[scatteredpix] = $
                 pixelmask[scatteredpix] + pixelmask_bits('SCATTEREDLIGHT')

   ;------------------
   ; Tweak up the wavelength solution to agree with the sky lines.
   ; xshet contains polynomial coefficients to shift arc to sky line frame.

   locateskylines, skylinefile, flux, fluxivar, wset, $
    xarc, arcshift=arcshift, $
    xsky=xsky, skywaves=skywaves, skyshift=skyshift

   qaplot_skyshift, wset, xsky, skywaves, skyshift, $
    title=plottitle+'Sky Line Deviations for '+objname

   if (NOT keyword_set(arcshift)) then $
    splog, 'WARNING: Cannot shift to sky lines'

   ;------------------
   ; Fit for the widths of the sky lines (relative to those of the arcs)

   ; The values in XSKY are noisy measurements, so replace them with
   ; the predicted positions based upon the arc solution.
   ; We should also apply the arcshift here, but I haven't yet ???

   xsky = transpose(traceset2pix(wset, alog10(skywaves)))
   skydispset = skyline_dispersion(flux, fluxivar, xsky, iskies, dispset)
   splog, 'Not applying skyline adjusted line widths'

   ;------------------
   ; Apply heliocentric correction
   ; Correct LAMBDA, which is used to shift to vacuum wavelengths.

   helio=0.0
   get_platecenter, plugsort, ra, dec

   ;--------------------------------------------------------
   ; Call standard proc to determine time-stamps

   get_tai, objhdr, tai_beg, tai_mid, tai_end

   ; Set TAI equal to the time half-way through the exposure
   ; If all these keywords are present in the header, they will be either
   ; type FLOAT or DOUBLE.  Note that SDSS will put NaN in the header for
   ; these values if they are unknown.
   if (size(tai_mid, /tname) NE 'INT' $
       AND finite(ra) AND finite(dec) AND finite(tai_mid) ) then begin
      helio = heliocentric(ra, dec, tai=tai_mid)
      splog, 'Heliocentric correction = ', helio, ' km/s'
      sxaddpar, objhdr, 'HELIO_RV', helio, $
       ' Heliocentric correction (added to velocities)'
   endif else begin
      splog, 'WARNING: Header info not present to compute heliocentric correction'
   endelse
   if (size(tai_mid, /tname) EQ 'INT' OR finite(tai_mid) EQ 0) then begin
      splog, 'WARNING: Header info not present to compute airmass correction to sky level'
      tai_mid = 0
   endif

   ;------------------
   ; Shift to skylines and fit to vacuum wavelengths

   vacset = fitvacset(xarc, lambda, wset, arcshift, helio=helio, airset=airset)
; No longer make the following QA plot ???
;   qaplot_skydev, flux, fluxivar, vacset, plugsort, color, $
;    title=plottitle+objname
   sxaddpar, objhdr, 'VACUUM', 'T', ' Wavelengths are in vacuum'

   ;------------------
   ;  If present, reconstruct superflat and normalize

   sxaddpar, objhdr, 'SFLATTEN', 'F', ' Superflat has not been applied'
   superfit = 0

   if keyword_set(superflatset) AND keyword_set(airset) then begin
     superfit = float(smooth_superflat(superflatset, airset, $
      plottitle=plottitle+'Smooth superflat for '+objname))
     if keyword_set(superfit) then begin
       divideflat, flux, invvar=fluxivar, superfit 
       sxaddpar, objhdr, 'SFLATTEN', 'T', ' Superflat has been applied'
     endif
   endif  

   ;----------
   ; Sky-subtract

   nbkpt = color EQ 'blue' ? 3*nx/4 : nx
   skystruct = skysubtract(flux, fluxivar, plugsort, vacset, $
    skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
    fibermask=fibermask, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt)
   if (NOT keyword_set(skystruct)) then return

   ;----------
   ; If any of the sky-fibers are bad, then re-do sky-subtraction.

   ibadfib = where(djs_median(skysub[*,iskies]^2 * $
    skysubivar[*,iskies], 1) GT 2.0)               
   if (ibadfib[0] NE -1) then begin
      fibermask[iskies[ibadfib]] = fibermask[iskies[ibadfib]] OR $
       fibermask_bits('BADSKYFIBER')

      splog, 'Calling skysubtract again; masked skyfibers',$
       string(iskies[ibadfib])
      skystruct = skysubtract(flux, fluxivar, plugsort, vacset, $
       skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
       fibermask=fibermask, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt)
      if (NOT keyword_set(skystruct)) then return
   endif

   ;----------
   ; QA plots for chi^2 from 1D sky-subtraction.

   qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
    vacset, iskies, title=plottitle+objname+' 1D Sky-subtraction'

   ;----------
   ; Sky-subtract one final time, this time with dispset (PSF subtraction)
   ; (rejected sky fibers from above remain rejected).
   ; Modify pixelmask in this call.

   nskypoly = 3L
   skystruct = skysubtract(flux, fluxivar, plugsort, vacset, $
    skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
    fibermask=fibermask, upper=10.0, lower=10.0, tai=tai_mid, $
    dispset=dispset, npoly=nskypoly, nbkpt=nbkpt, $
    relchi2set=relchi2set, newmask=newmask)
   pixelmask = newmask
   if (NOT keyword_set(skystruct)) then return

   ;----------
   ; QA plots for chi^2 from 2D sky-subtraction.

   qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
    vacset, iskies, title=plottitle+objname+' 2D Sky-subtraction'

   ;----------
   ; QA for 2 skylines in the blue (specify vacuum wavelengths below)

   if (color EQ 'blue') then begin
      qaplot_skyline, 4359.5, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=4.0, $
       tai=tai_mid, title=plottitle+objname+' Skyline Flux at 4359.5'
      qaplot_skyline, 5578.9, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=5.0, $
       tai=tai_mid, title=plottitle+objname+' Skyline Flux at 5578.9'
   endif

   ;----------
   ; QA for 2 skylines in the red (specify vacuum wavelengths below)

   if (color EQ 'red') then begin
      qaplot_skyline, 7343.0, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=7.0, $
       tai=tai_mid, title=plottitle+objname+' Skyline Flux at 7343,0'
      qaplot_skyline, 8888.3, flux, fluxivar, skysub, skysubivar, $
       plugsort, vacset, iskies, fibermask=fibermask, dwave=7.0, $
       tai=tai_mid, title=plottitle+objname+' Skyline Flux at 8888.3'
   endif

   ;------------------------------------------
   ; Save the sky-subtracted flux values as is, and now modify flambda.

   flambda = skysub
   flambdaivar = skysubivar
   skyimg = flux - flambda
   sxaddpar, objhdr, 'PSFSKY', nskypoly, ' Order of PSF skysubtraction'
;   if (keyword_set(relchi2struct)) then skychi2 = mean(relchi2struct.chi2) $
;    else skychi2 = 0.0
;   sxaddpar, objhdr, 'SKYCHI2', skychi2, ' Mean chi^2 of sky-subtraction'

   ;------------------------------------------
   ; Telluric correction called for 'red' side

   if (keyword_set(do_telluric) AND color EQ 'red')  then begin
      telluricfactor = telluric_corr(flambda, flambdaivar, vacset, plugsort, $
       fibermask=fibermask, $
       plottitle=plottitle+'Telluric correction for '+objname)

      divideflat, flambda, invvar=flambdaivar, telluricfactor, minval=0.1
   endif

   ;----------
   ; Interpolate over masked pixels, just for aesthetic purposes

   flambda = djs_maskinterp(flambda, flambdaivar EQ 0, /const, iaxis=0 )

   ;----------
   ; Combine FIBERMASK and PIXELMASK to FINALMASK

   finalmask = pixelmask
   for itrace=0, ntrace-1 do $
    finalmask[*,itrace] = finalmask[*,itrace] OR fibermask[itrace]

   ;----------
   ; Get an estimate of the relative chi^2 at each pixel.
   ; Do this with a simple linear interpolation.

   if (keyword_set(relchi2set)) then begin
      xx = 0
      traceset2xy, vacset, xx, loglam
;      rchi2img = interpol(relchi2struct.chi2, relchi2struct.wave, loglam)
      rchi2img = bspline_valu(loglam, relchi2set)
      skychi2 = mean(rchi2img)
   endif else begin
      rchi2img = 0 * flambda + 1.
      skychi2 = 0.
   endelse
   sxaddpar, objhdr, 'SKYCHI2', skychi2, ' Mean chi^2 of sky-subtraction'

   ;----------
   ; Add keywords to object header

   sxaddpar, objhdr, 'VERS2D', idlspec2d_version(), $
    ' Version of idlspec2d for 2D reduction', after='VERSREAD'
   if (keyword_set(osigma)) then $
    sxaddpar, objhdr, 'OSIGMA',  sigma, $
     ' Original guess at spatial sigma in pix'
   sxaddpar, objhdr, 'PREJECT', reject[1], ' Profile area rejection threshold'
   sxaddpar, objhdr, 'LOWREJ', lowrej, ' Extraction: low rejection'
   sxaddpar, objhdr, 'HIGHREJ', highrej, ' Extraction: high rejection'
   sxaddpar, objhdr, 'SCATPOLY', npoly, ' Extraction: Order of scattered light polynomial'
   sxaddpar, objhdr, 'PROFTYPE', proftype, ' Extraction profile: 1=Gaussian'
   sxaddpar, objhdr, 'NFITPOLY', nterms, ' Extraction: Number of parameters in each profile'
   sxaddpar, objhdr, 'XCHI2', mean(chisq), ' Extraction: Mean chi^2'

   snvec = djs_median(flambda*sqrt(flambdaivar),1)
   if color EQ 'blue' then begin
        mag=plugsort.mag[1] 
        snmag = 20.2
   endif else begin
        mag=plugsort.mag[3]
        snmag = 19.9
   endelse

   fitmag = snmag + [-2.0,-0.5]
   sncoeff = fitsn(mag, snvec, fitmag=fitmag, sigma=sigma)
   sn2 = 10^(2.0 * poly([snmag],sncoeff))

   sxaddpar,objhdr,'FRAMESN2', sn2[0]
   sxaddpar,objhdr,'EQUINOX',2000.0,after='DEC'
   sxaddpar,objhdr,'RADECSYS', 'FK5', after='EQUINOX'
   sxaddpar,objhdr,'AIRMASS',$
    float(tai2airmass(ra, dec, tai=tai_mid)) $
            * (tai_mid GT 0), after='ALT'

   spadd_guiderinfo, objhdr

   sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'

   ;----------
   ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk

   mwrfits, flambda, outname, objhdr, /create
   mwrfits, flambdaivar, outname
   mwrfits, finalmask, outname
   mwrfits, vacset, outname
   mwrfits, dispset, outname
   mwrfits, plugsort, outname
   mwrfits, skyimg, outname
   mwrfits, superfit, outname
   mwrfits, skystruct, outname
   mwrfits, scatter, outname
   if (keyword_set(do_telluric)) then $
    mwrfits, telluricfactor, outname ; This array only exists for red frames.
   spawn, ['gzip', '-f', outname], /noshell

   heap_gc

   return
end
;------------------------------------------------------------------------------
