;+
; NAME: 
;       MR1D_DETECT
;
; PURPOSE: 
;       Detect the bands in a 1D signal by using a multiresolution transform
;
; CALLING:
;       MR_Detect, Signal, result, OPT=OPT, print=print, tabobj=tabobj, $
;              NbrObj=NbrObj, tabband=tabband, nodel=nodel, tabw=tabw 
;                   
;
; INPUT:
;
; KEYWORDS:
;      OPT: string which contains the differents options. Options are:
;   	     [-n number_of_scales]
;    	         number of scales used in the multiresolution transform 
;    	         for the filtering. 
;    	          default is 5
; 
;   	      [-s NSigma]
;   	           Detection level at NSigma * SigmaNoise
;   	           default is 3.000000
; 
;  	       [-a]
;  	            detection of Absorption lines. Default is no. 
; 
;   	      [-e]
;   	           detection of Emission lines. Default is no. 
; 
;    	      [-m type_of_noise]
;              1: Gaussian Noise 
;              2: Poisson Noise 
;              3: Poisson Noise + Gaussian Noise 
;              4: Multiplicative Noise 
;              5: Non uniform additive noise 
;              6: Non uniform multiplicative noise 
;              7: Undefined uniform Noise 
;              8: Undefined Noise 
;      	       Default is Gaussian noise.
;
;   	      [-g SigmaNoise]
;     	          gaussian noise standard deviation. 
; 
; 
;    	      [-c gain,sigma,mean]
;     	        gain = gain of the CCD
;     	        sigma = read-out noise standard deviation
;     	        mean = read-out noise mean
;      	         noise = poisson + readout noise. default is no (Gaussian)
;      	         
;             [-a]
;               detection of Absorption lines. Default is no. 
; 
;             [-e]
;               detection of Emission lines. Default is no. 
; 
;      	      [-f FirstScale]
;      	        first scale. Default is 1.
; 
;              [-l LastScale]
;               Last scale. Default is number_of_scales-2
; 
;              [-i IterNumber]
;          	   Number of iteration for the reconstruction. 
;           	  Default is 20
; 
;              [-M]
;           	  Use the multiresolution median transform 
;           	  instead of the a-trous algorithm. 
; 
;              [-A]
;              detect only negative  multiresolution coefficients. Default is no. 
; 
;              [-E]
;              detect only positive multiresolution coefficients. Default is no. 
; 
;              [-w ]
;              write other results:
;                tabadd.fits = sum of the reconstructed objects  
;                tabseg.fits = segmented wavelet transform  
;
;              [-v]
;              Verbose. Default is no.
;
;      print: if set, information about each detected object is printed.
; 
;      tabobj: IDL structure which contains the information about 
;              the objects. For each object, we have
;                 NumObj: Object number 
;                 Pos: Position in the signal (in wavelength unit).
;                 Sigma: standard deviation of the band (in wavelength unit).
;                 Fwhm:  Full width at half maximum (in wavelength unit).
;                 PosPix:  Position in the signal (in pixel unit).
;                 SigmaPix:standard deviation of the band(in pixel unit).  
;                 FwhmPix" Full width at half maximum  (in pixel unit).             
;                 ValMax: Maximum value of the band
;                 Flux: integrated flux of the band.
;
;      NbrObj: number of detected bands
;
;      nodel: if set, the created files (by the program mr1d\_detect)
;             are not deleted.
;
;      tabw=tabw: wavelength array
;
; OUTPUTS:
;           Result: result of the detection. The output image contains all
;                   detected  objects which are coadded.
;
; EXTERNAL CALLS
;           mr1d_detect (C++ program)
;
; EXAMPLE:
;           detection   with all default options 
;                mr1d_detect, Signal, Result
;
; HISTORY:
;       Written: Jean-Luc Starck 1997.
;       May, 1997 File creation
;-



pro mr1d_infoband, band, flux, sig, Fwhm, CenterPos, $
   MaxObj,ni=ni, xrange=xrange

if not keyword_set(ni) then ni = 0

MaxObj = max(band)
MinObj = min(band)
if abs(MinObj) GT abs(MaxObj) then sign = -1 else sign = 1
data = band*sign

flux = total(data)
vs = size(band)
Np = vs(1)



if ni GT 0 then BEGIN

x = float(indgen(Np))  + 1
C0 = total(x*data) / flux
sig0 = total(x*x*data) / flux - C0^2
print, 'sig0 = ', sig0
sig0 = sqrt( sig0 )
Max0 = max(data)

Cn = C0
Sign = sig0
Maxn = Max0

for iter = 1, ni do BEGIN
  if Sign GT 0.00001 then BEGIN
     y = exp( - (x - Cn)^2. / (2*Sign^2))*Maxn
     Maxt = max(y)
     flux_y = total(y)
     Ct = total(x*y) / flux_y
     Sigt = sqrt( total(x*x*y) / flux_y - Ct^2)
     Cn = Cn + C0 - Ct
     Sign = Sign + sig0 - Sigt
     Maxn = Maxn + Max0 - Maxt
  END ELSE BEGIN
     Cn = C0
     Sign = sig0
     Maxn = Max0
  END
  print, "Pos = ", Cn, "  Sigma = ", Sign, "   Max = ", Maxn
 END

CenterPos = Cn - 1
sig = Sign
MaxObj = Maxn
END ELSE BEGIN

if not keyword_set(xrange) then x = indgen(Np) else x= xrange
y = band
r = gaussfit(x,y,a)
CenterPos =  a(1)
sig = a(2)
MaxObj = a(0)
Fwhm = sig * (2.*sqrt(2.*alog(2.)))
 END
end

;-----------------------------------------------------------------

pro mr1d_detect, Signal, result, OPT=OPT, print=print, tabobj=tabobj, $
    NbrObj=NbrObj, tabband=tabband, nodel=nodel, tabw=tabw

if N_PARAMS() LT 2 then begin 
        spawn, 'mr1d_detect'
        print, 'CALL SEQUENCE: mr1d_detect, Signal, Signal_Out, OPT=Opt,'
        print, 'print=print, tabobj=tabobj, NbrObj=NbrObj, nodel=nodel, tabw=tabw'
        goto, DONE
        end
        
 if not keyword_set(Opt) then Opt = ' '  

       
Nx = (size(Sginal))(1)
NameSig = 'xx_signal.fits'
NameResult = 'xx_result.fits'
NameTabAdd =  'tabadd.fits'

writefits,  NameSig,  Signal
com = 'mr1d_detect ' + OPT + ' -w '+ NameSig  + ' ' +  NameResult
spawn, com
Result = readfits(NameTabAdd, /silent) ; result=1D array (co-added bands)
tabband = readfits(NameResult, /silent)   ; tabband=2d array (1 band per line)

vs = size(tabband)
if vs(0) eq 2 then NbrObj = vs(2) $
else if vs(0) eq 1 then  NbrObj = 1 else NbrObj = 0
if n_elements(tabband) EQ 1 then NbrObj = 0
 
print, "Number of detected bands = ", NbrObj
if NbrObj LT 1 then BEGIN
   tabobj = -1
   tabband = -1
   Result = -1
   goto, done
   END
info_obj = {   NumObj: 0, $
                Pos : 0., $
                Sigma: 0., $
                Fwhm: 0., $
                PosPix : 0., $
                SigmaPix: 0., $
                FwhmPix: 0., $             
                ValMax: 0., $
                Flux: 0.}
tabobj = replicate(info_obj, NbrObj)

for i=0,NbrObj-1 do begin
mr1d_infoband, tabband(*,i), flux, sig, Fwhm, CenterPos, MaxObj,xrange=tabw
tabobj(i).Sigma = sig
tabobj(i).Fwhm = Fwhm  
tabobj(i).Pos = CenterPos 
 
mr1d_infoband, tabband(*,i), flux, sig, Fwhm, CenterPos, MaxObj
tabobj(i).NumObj = i+1
tabobj(i).PosPix = CenterPos
tabobj(i).SigmaPix = sig
tabobj(i).FwhmPix = Fwhm
tabobj(i).Flux = flux
tabobj(i).ValMax = MaxObj
 
end

if not keyword_set(nodel) then delete, NameSig
if not keyword_set(nodel) then delete, NameResult
if not keyword_set(nodel) then delete, NameTabAdd
; delete, "tabseg.fits"
DONE:
end

; f = 'hh100_FINAL.fits'
; spec = readspectrum(f,nbrscale=7)
; s= spec.buff
; w= spec.wave
; mr1d_detect,s,ts
; mr1d_detect,s1,tb,opt='-s5 -n7 -e -m5 -f3',tabband=tb,tabobj=to,/nodel,tabw=w
