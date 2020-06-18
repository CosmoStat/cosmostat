;+
; NAME:  
;       GETGAL
;
; PURPOSE:
;        Calculate a overdensity map with a given power spectrum (defauit is 2mass galaxy survey power spectrum), 
;         then calculate a density field by ;        
;            DensityField = MeanGalPerPixel * ( 1.0d + OverDensityFIeld )
;        By default, MeanGalPerPixel is equal to 1.6  at nside=256 (this is the corresponding 2mass value)
;        Then return the density field or poissson realasition is the keyword /noise is set
;
; CALLING:
;     GalMap = getgal(powerspec=powerspec,  nolog=nolog,  noise=noise, cmean=cmean, poisson=poisson, lmax=lmax, nside=nside, overdens=overdens)
; 
; INPUTS:
;     
; OUTPUTS:
;     GalMap -- IDL 1D array: galaxy survey map   
;          
; INPUT KEYWORDS:
;      powerspec -- double array: power spectrum of the overdensity galaxy field. By default the 2mass theoretical power spectrum (in the PR1 Planck cosmology) is used.
;      Noise --  int : if set, a mean density is added (variable cmean), and Poisson noise is added.
;      cmean --  double : mean galaxy density. By default, the 2mass mean at nside=256 is used, i.e. cmean = 1.6  
;      nolog --  int : if set, a  Gaussian simulation is done instead of a lognormal one
;      lmax --  int :lmax   :  largest multipole.
;      nside -- int : nside of the output map
;
; OUTPUT KEYWORDS:
;      overdens -- IDL Healpix map:  overdensity map. If the keyword /noise is used, the noise-free overdensity map is return in overdens
;      DensityField -- IDL Healpix map:  density map. If the keyword /noise is not used, this is the sane map as the returned map.
;    
;   EXTERNAL CALLS:
 ;
; EXAMPLE:
;     Galmap  =  getgal(over=map, dens=d, /noise)
;     Calculate a noisy 2mass simulated map Galmap, its overdensity map and its density d  (i..e  Galmap = poisson(density)).
;
; HISTORY:
;       Written : F.X. Dupe and Jean-Luc Starck   2014.
;-
;-----------------------------------------------------------------

function getgal, powerspec=powerspec,  nolog=nolog,  noise=noise, cmean=cmean, poisson=poisson, lmax=lmax, nside=nside, overdens=overdens, verb=verb, DensityField= DensityField

; IF N_PARAMS() LT 1 THEN BEGIN        
;PRINT,'Error: 2 inputs needed. ...'  
;GOTO, DONE 
;ENDIF

if keyword_set(Verb) then begin


end

  if not keyword_set(nside) then nside = 512
  gmean = 1.6d ;; Mean number of galaxies per pixel in the 2MASS survey at NSIDE=256
  if  not  keyword_set(cmean) then cmean= gmean * 256.^2 / double(nside)^2

if not keyword_set(powerspec) then begin
   restore, /verb, '$ISAP/param/galaxies/1202_2mass_isw.sav'
   lmaxGal = long(max(cgg.l))
   galspec = dblarr(lmaxGal +1)
   galspec[1: lmaxGal] =  interpol(reform(CGG.cl), CGG.l, fix(findgen(lmaxGal)+1), /qua)
   if not keyword_set(lmax) then lmax = lmaxGal $
   else begin
       if lmax LT lmaxGal then galspec = galspec[0: lmax] $
       else lmax = lmaxGal
   end
   powerspec=galspec
end

  npix = nside2npix(nside)

  ;; First generate the overdensity
  if not keyword_set(nolog) then begin

     ;; First the normal field
     signal = dcomplexarr(lmax+1,lmax+1)
     ppln   = alog(1.0d + powerspec)
           
     for l=1ul,lmax do begin
           rv0  = randomn(seed,1,/normal,/double)*sqrt(ppln[l]) - dcomplex(ppln[l])
           rvar = randomn(seed,l,/normal,/double)*sqrt(ppln[l]/2.0d) - dcomplex(0.5d * ppln[l],0.5d * ppln[l])
           ivar = randomn(seed,l,/normal,/double)*sqrt(ppln[l]/2.0d) - dcomplex(0.5d * ppln[l],0.5d * ppln[l])
           signal[l,0]   = dcomplex(rv0,0.0d)
           signal[l,1ul:l] = dcomplex(rvar,ivar)
     endfor

     ;; Second convert to a log-normal field
     gaussf = dblarr(npix)
     mrs_almtrans,gaussf,galm,/tab,/complex,lmax=lmax
     galm.alm = signal
     mrs_almrec,galm,gaussf
     overdens = exp(gaussf) - 1.0d

  endif else begin

     signal = dcomplexarr(lmax+1,lmax+1)
     ppga   = powerspec
           
     for l=1ul,lmax do begin
              
        rv0  = randomn(seed,1,/normal,/double)*sqrt(ppga[l])
        rvar = randomn(seed,l,/normal,/double)*sqrt(ppga[l]/2.0d)
        ivar = randomn(seed,l,/normal,/double)*sqrt(ppga[l]/2.0d)

        signal[l,0]   = dcomplex(rv0,0.0d)
        signal[l,1ul:l] = dcomplex(rvar,ivar)

     endfor

     gaussf = dblarr(npix)
     mrs_almtrans,gaussf,galm,/tab,/complex,lmax=lmax
     galm.alm = signal
     mrs_almrec,galm,overdens

  endelse

  ;; Add a mean number of count before Poisson noise

  DensityField = cmean * ( 1.0d + overdens )
  
  ;; Add noise if needed
  if keyword_set(noise) then begin
     imagpois = DensityField
     for i=0ul,npix-1 do if DensityField[i] le 0 then imagpois[i] = 0.0d else imagpois[i] = randomn(seed,poisson=DensityField[i],/double)
     return, imagpois
  endif else $
     return,  DensityField

DONE:
     return, -1
end
