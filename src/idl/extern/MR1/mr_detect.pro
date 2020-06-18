;+
; NAME: 
;       MR_DETECT
;
; PURPOSE: 
;       Detect the sources in  an image by using a multiresolution transform. 
;       This routine is calling the C++ executable {mr\_detect}. 
;       The keyword "OPT" allows to pass to the executable all
;       options. By default, all files created by the executable
;       are deleted.
;
; CALLING:
;       MR_Detect, Imag, Result, OPT=Opt, print=print,
;                  tabobj=tabobj, NbrObj=NbrObj, nodel=nodel, RMS=RMS
;
; INPUT:
;       Imag: image 
;
; KEYWORDS:
;      print: if set, information about each detected object is printed.
;
;      tabobj: IDL structure which contains the information about 
;              the objects. For each object, we have
;                 ScaleObj: Scale where the object is detected.
;                 NumObj: Object number 
;                 PosX: X coordinate in the image.
;                 PosY: Y coordinate in the image.
;                 SigmaX: standard deviation in first main axis of the object.
;                 SigmaY: standard deviation in second main axis of the object.
;                 Angle: angle between the main axis and x-axis.
;                 ValPixMax: Value of the maximum of the object.
;                 Flux: integrated flux of the object.
;                 Magnitude: magnitude of the object.
;                 ErrorFlux: flux error.
;                 SNR_ValMaxCoef: signal to noise ration of 
;                                  the maximum of the wavelet coefficient.
;                 SNR_Obj: signal to noise ration of the object.
;                 PosCoefMaxX:  X coordinate of maximum wavelet coefficient
;                 PosCoefMaxY:  Y coordinate of maximum wavelet coefficient
;
;      NbrObj: number of detected objects
;
;      nodel: if set, the created files (by the program mr\_detect)
;             are not deleted.
;
;      RMS: RMS image related to the input image.
;
;      Opt: string which contains the differents options. Options are:
;
;
;         [-t type_of_multiresolution_transform]
;              1: bspline wavelet transform: a trous algorithm 
;              2: Half-pyramidal transform 
;              3: pyramidal bspline wavelet transform 
;              4: Mixed WT and PMT method (WT-PMT) 
;              5: Mixed Half-pyramidal WT and Median method (WT-HPMT) 
;             default is bspline wavelet transform: a trous algorithm
;
;          [-V Multiscale_Vision_Model]
;              1: No vision model 
;              2: Blinded Objects 
;              3: Rue-Bijaoui Vision Model for blinded + embedded Objects 
;              default is Rue-Bijaoui Vision Model for blinded + embedded Objects
;
;         [-n number_of_scales]
;             number of scales used in the multiresolution transform
;              default is 5
;
;         [-m type_of_noise]
;              1: Gaussian Noise 
;              2: Poisson Noise 
;              3: Poisson Noise + Gaussian Noise 
;              4: Multiplicative Noise 
;              5: Non uniform additive noise 
;              6: Non uniform multiplicative noise 
;              7: Undefined uniform Noise 
;              8: Undefined Noise 
;              9: Stationary correlated noise 
;              10: Poisson Noise with few events 
;             default is Gaussian noise
;
;         [-g sigma]
;             sigma = noise standard deviation
;              default is automatically estimated
;
;         [-c gain,sigma,mean]
;             gain = gain of the CCD
;             sigma = read-out noise standard deviation
;             mean = read-out noise mean
;               noise = poisson + readout noise. default is no (Gaussian)
;
;         [-s nsigma]
;              Thresholding at nsigma * SigmaNoise
;              default is 3.000000
;
;         [-E Epsilon]
;             Epsilon = precision for computing thresholds
;                       (only used in case of poisson noise with few events)
;             default is 1.000000e-03 
;
;         [-e minimum_of_events]
;             Minimum number of events for a detection.
;             default is 4
;
;         [-S SizeBlock]
;             Size of the  blocks used for local variance estimation.
;            default is 7 
;
;         [-N NiterSigmaClip]
;             iteration number used for local variance estimation.
;             default is 1 
;
;         [-F first_detection_scale]
;              first scale used for the detection 
;              default is 1
;
;         [-R RMS_Map_File_Name]
;               RMS Map (only used with -m 5 and -m 9 options). 
;
;         [-L last_detection_scale]
;              last scale used for the detection 
;
;         [-i number_of_iterations]
;              iteration number per object reconstruction
;              default is 10
;
;         [-u object_reconstruction_error]
;              default is: 0.000010
;
;         [-k]
;             keep isolated objects
;             default is no.
; 
;         [-K]
;             Keep objects at the border.
;             default is no.
; 
;         [-A FluxMult]
;              Flux in tex table are multiplied by FluxMul.
;              default is 1.
; 
;         [-w writing_parameter]
;              1: write each object separately in an image. 
;                 The image file name of the object will be: 
;                      ima_obj_xx_yy.fits 
;              2: simulated two images 
;                   xx_ellips.fits: an ellipse is drawn around each object 
;                   xx_simu.fits: image created only from the morphological ;parameters 
;              3: equivalent to 1 and 2 together 
;
;         [-U]
;              Sub Segmentation.
;
;         [-p]
;              Detect also negative structures 
;              default is no.
;
;         [-q]
;              Define the root of an object from the maximum position and its value 
;
;         [-d DistMax]
;              Maximum distance between two max positions
;              of the same object at two successive scales.
;              Default is 1.
;         [-o] 
;             Sub-object analysis. Default is no. 
;
;         [-D] 
;             Perform a deconvolution. Default is no. 
;
;         [-P PsfFileName] 
;             PSF file name.
;
;         [-f Fwhm] 
;             Full Width at Half Maximum.
;
;         [-O PSF_Sampling] 
;             PSF over-sampling value.
;
;         [-a BgrMethod] 
;             Aperture photometry.
;              1: Aperture photometry using a background image model 
;              2: Aperture photometry using an estimated background image 
;              3: Aperture photometry using a constant background 
;              default is no aperture photometry.
;
;         [-B BgrFileName] 
;             Background image file name.
;
;         [-G BgrValue] 
;             Constant background value.
;             default is: 0.000000.
;
;         [-b BGR_Size] 
;             Background image size for automatic background esitmation.
;             default is: 16.
;
;         [-l KSigmaAperture] 
;             Aperture photometry size parameter.
;             default is: 3.000000.
;	     
;         [-M object_reconstruction_method]
;              1: reconstruction from the fixed step gradient method 
;              2: reconstruction from the optimum step gradient method 
;              3: reconstruction from the conjugate gradient method 
;              4: reconstruction using the PSF 
;              default is: reconstruction from the conjugate gradient method
;
;         [-C RADEC_Table_Order]
;              1: tex table ordered by object number 
;              2: tex table ordered by the Right Ascension 
;              3: tex table ordered by object SNR 
;              default is: tex table ordered by object number.
; 
;         [-v] 
;              Verbose
;
; OUTPUTS:
;           Result: result of the detection. The output image contains all
;                   detected  objects which are coadded.
;
; EXTERNAL CALLS:
;           mr_detect (C++ program)
;
; EXAMPLE:
;           detection in an image with all default options 
;                mr_detect, Imag, Result
;
; HISTORY:
;	Written: Jean-Luc Starck 1997.
;	January, 1997 File creation
;       October, 1999 Header Update
;-

PRO mr_detect, imag, result, OPT=OPT, print=print, tabobj=tabobj, $
    NbrObj=NbrObj, nodel=nodel, RMS=RMS

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_detect'
        print, 'CALL SEQUENCE: mr_detect, Imag, Imag_Out, OPT=Opt, print=print, tabobj=tabobj, NbrObj=NbrObj, nodel=nodel'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
NameRMS = 'xx_rms.fits'

NameRin = 'xx_result'
NameResult = 'xx_result.fits'
TabDetect =  "xx_result.mes"
Tabmrg = 'xx_result.mrg'
Tabmro = 'xx_result.mro'
Tabtex = 'xx_result.tex'
Tabmes = 'xx_result.mes'
TabPS1 = 'xx_result_obj.ps'
TabPS2 = 'xx_result_subobj.ps'

writefits, NameImag, imag
if keyword_set(RMS) then begin
    writefits, NameRMS, RMS
    OPT = OPT + ' -R ' + NameRMS
    end

com = 'mr_detect ' + OPT + ' ' + NameImag + ' ' + NameRin
spawn, com
Result = readfits(NameResult, /silent)

if not keyword_set(nodel) then begin
delete, NameImag
delete, NameResult
delete, Tabmrg
delete, Tabmro
delete, Tabtex
delete, TabPS1
delete, TabPS2
if keyword_set(RMS) then delete, NameRMS
end

openr,1, TabDetect
NbrObj = 0
info_obj = {  ScaleObj: 0, $
                NumObj: 0, $
                PosX : 0., $
                PosY  : 0., $
                SigmaX: 0., $
                SigmaY: 0., $
                Angle: 0., $
		ValPixMax: 0., $
                Flux: 0., $
                Magnitude: 0., $
                ErrorFlux: 0., $
                SNR_ValMaxCoef: 0., $
                SNR_Obj: 0., $
		PosCoefMaxX: 0., $
		PosCoefMaxY: 0.}

tabobj=-1
while not eof (1) do begin
   readf,1, ScaleObj, NumObj
   readf,1, PosX, PosY, SigmaX, SigmaY
   readf,1, Angle, ValPixMax, Flux, Magnitude
   readf,1, ErrorFlux, SNR_ValMaxCoef, SNR_Obj
   readf,1, PCMX, PCMY
   
   NbrObj = NbrObj + 1
   NewObj = info_obj
   NewObj.NumObj = NumObj
   NewObj.ScaleObj = ScaleObj
   NewObj.PosX = PosX
   NewObj.PosY = PosY
   NewObj.SigmaX = SigmaX
   NewObj.SigmaY = SigmaY
   NewObj.Angle = Angle
   NewObj.ValPixMax = ValPixMax
   NewObj.Flux = Flux
   NewObj.Magnitude = Magnitude
   NewObj.ErrorFlux = ErrorFlux
   NewObj.SNR_ValMaxCoef = SNR_ValMaxCoef
   NewObj.SNR_Obj = SNR_Obj
   NewObj.PosCoefMaxX = PCMX
   NewObj.PosCoefMaxY = PCMY
   
   if NbrObj GT 1 then BEGIN
     tabobj = replicate(info_obj, NbrObj)
     tabobj[0:NbrObj-2] = TabOld
     tabobj[NbrObj-1] = NewObj
     TabOld = tabobj
   END ELSE BEGIN
       tabobj = NewObj
       TabOld = NewObj
   END

   if keyword_set(print) then BEGIN
      print, " ScaleObj = ", ScaleObj, " NumObj = ", NumObj
      print, " PosX = ", PosX, " PosY = ", PosY, " SigmaX = ", SigmaX, " SigmaY = ", SigmaY
      print, " Angle = ", Angle, " ValPixMax = ", ValPixMax, " Flux = ", Flux, " Magnitude = ", Magnitude
      print, " ErrorFlux = ", ErrorFlux, " SNR_ValMaxCoef = ", SNR_ValMaxCoef, " SNR_Obj = ", SNR_Obj
      print, " PosCoefMaxX = ", PCMX, " PosCoefMaxY = ",  PCMY
   END

endwhile
close, 1

if not keyword_set(nodel) then delete, Tabmes

DONE:
end
