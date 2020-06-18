;+
; NAME: 
;       MR1D_TRANS
;
; PURPOSE: 
;        One dimensional wavelet transform. 16 transforms are available, which are
;        grouped in 5 classes 
;                         Class 1: no decimation
;                                  transform 1 to 5 and 9 to 12
;                         Class 2: pyramidal transform
;                                  transform 6 to 8
;                         Class 3: orthogonal transform (13 and 14)
;                         Class 4: Wavelet packets (15)
;                         Class 5: Wavelet packets from the a-trous algorithm  (16)
;
;        Depending on the class, the transform does not contain the
;        same number of pixels, and the data representation differs.
;
; CALLING:
;       MR1D_Trans, Signal, result, OPT=OPT, BAND=BAND, NODEL=NODEL
;                   
;
; INPUT:
;      Signal: 1D array; input signal 
;
; KEYWORDS:
;      OPT: string which contains the differents options. Options are:
;         [-t type_of_multiresolution_transform]
;               1: linear wavelet transform: a trous algorithm 
;               2: b1spline wavelet transform: a trous algorithm 
;               3: b3spline wavelet transform: a trous algorithm 
;               4: Derivative of a b3spline: a trous algorithm 
;               5: undecimated Haar wavelet transform: a trous algorithm 
;               6: morphological median transform 
;               7: Undecimated (bi-) orthogonal wavelet transform 
;               8: pyramidal linear wavelet transform 
;               9: pyramidal b3spline wavelet transform 
;               10: pyramidal median transform 
;               11: Morlet's wavelet transform 
;               12: mexican hat wavelet transform 
;               13: french hat wavelet transform 
;               14: Gaussian Derivative wavelet transform 
;               15: bi-orthogonam transform   
;               16: bi-orthogonam transform via lifting sheme (CDF filters) 
;               17: Wavelet packets (CDF filters) 
;               18: Wavelet packets from lifting sheme 
;               19: Wavelet packets using the a-trous algorithm 
;               Default is 3.
;
;          [-r]
;               rebin all scales to the input size
;               (for pyramidal median transform only)
;
;          [-k]
;               set to 0 pixels contaminated by the border problem.     
; 
;   	   [-n number_of_scales]
;    	        number of scales used in the multiresolution transform 
;    	        Default value is automatically calculated.
;
;          [-T type_of_filters]
;              1: Antonini 7/9 filters 
;              2: Daubechies filter 4 
;              3: Biorthogonal 2/6 Haar filters 
;              4: Biorthogonal 2/10 Haar filters 
;              5: Odegard 7/9 filters 
;              6: User's filters 
;              default is Antonini 7/9 filters
;
;          [-L]
;              Use a L2 normalization. Default is L1.
;
;          [-l type_of_lifting_transform]
;              1: Lifting scheme: CDF WT 
;              2: Lifting scheme: median prediction 
;              3: Lifting scheme: integer Haar WT 
;              4: Lifting scheme: integer CDF WT 
;              5: Lifting scheme: integer (4,2) interpolating transform 
;              6: Lifting scheme: 7/9 WT 
;              7: Lifting scheme: integer 7/9 WT 
;             default is Lifting scheme: integer Haar WT
;
;         [-w InfoFileName]
;              write in xx_info.fits the size and the starting index of each band.
;
;
; 
;      BAND:    if set, a tag per band is created in the output structure
;
;      NODEL:   if set, the two created file are not deleted: 
;                xx_result.fits : wavelet coefficients file
;                xx_info.fits:    information about the transform
;          
; OUTPUTS:
;           Result: IDL structure which contains the wavelet transform
;            The structure contains the following tags:
;             N_BAND  : float  ; number of bands in the transfrom    
;                                                
;             INFO    : 2D float array (Array[2, NbrBand+3])
;                       info[0,0] = transform number
;                       info[1,0] = number of scales
;                       info[0,1] = transform class number (5 classes)
;                       info[1,1] = number of bands
;                                 it is not equal to the number of scales
;                                 for wavelet packets transform.
;                       info[0,2] = number of pixels
;                       info[1,2] = lifting scheme type
;                       info[0,3] = type of filter
;                       info[1,3] = type of normalization
;                       for i=4 to Number_of_bands + 3
;                       info[0,i] = number of pixels in the band i
;                       info[1,i] = position number of the pixel of the band
;                       If a user filter file is given (i.e. -T 6,filename), 
;                       with a filename of $L$ caracters, $L$ lines are added 
;                       to the array:
;                       info[1,Number_of_bands + 4] = number of caracters of 
;                                                    the filter file name
;                        for i=Number_of_bands+4 to  Number_of_bands+4+L-1
;                             info[0,i] = ascii number of the ith caracter.
;
;             FROM    : 1D array (NbrBand); position of the first pixel 
;             TO      : 1D array (NbrBand); position of the last pixel
;
;             COEF    : 1D or 2D FLOAT array; wavelet coefficients
;                       for non-redundant transform (13,14,15), it is a 1D array
;                       for other transform, it is 2D array
;                       class 1 and 5: coeff[*,i] = band i (i in [0..NbrBand-1])
;                       class 2: coeff[0:to[i],i] = band i
;                       class 3 and 4: coeff[from[i]:to[i]] = band i
; 
;         if BAND keyword is set, the array coef is also splitted into bands:
;            BAND1   : band 1  
;            BAND2   : band 2  
;                      ...
;            BANDi   : band i  
;             
;
; EXTERNAL CALLS
;           mr1d_trans (C++ program)
;
; EXAMPLE:
;   > mr1d_trans, Signal, Result
;   > plot, result.coef[*,1]
;       Wavelet transform using the a-trous algorithm, and plot the second band
;
;   > mr1d_trans, Signal, Result, /band
;   > plot, result.band2
;       Idem
;
;   > mr1d_trans, Signal, Result, OPT='-t 13 -n 4'
;   > B2First = result.from[1]
;   > B2End = result.to[1]
;   > plot, result.coef[B2First:B2End]
;       Orthogonal wavelet transform with 4 scales, and plot the second scales
;
;   > mr1d_trans, Signal, Result, OPT='-t 13 -n 4', /band
;   > plot, result.band2
;       Idem
;
;   > mr1d_trans, Signal, Result, OPT='-t 8 -n 4'
;   > B2First = 0
;   > B2End = result.to[1]
;   > plot, result.coef[B2First:B2End,1]
;       Pyramidal median transform with 4 scales, and plot the second scales
;
; HISTORY:
;       Written: Jean-Luc Starck 1998.
;       July, 1998 File creation
;-


;-----------------------------------------------------------------

PRO mr1d_struc, result, info, output, band=band
; print, "transform number = ", Info[0,0]
;print, "number of scales = ", Info[1,0]
;print, "transform class number = ", Info[0,1]
;print, "number of bands = ", Info[1,1]

  Np = info[0, 2]
  trans = info[0, 1]   ; transform class number
  N_Scale = info[1,1]  ; band number
  from = fltarr(N_Scale)
  to = fltarr(N_Scale)
  
  for i=1, N_Scale do begin
     first = info[1,i+3]
     if trans NE 3 and trans NE 4 then first = 0
      Size = info[0,i+3]
      last  = first + Size - 1
      from[i-1] = first
      to[i-1] = last
      ; print, i-1, first, last
 
     if keyword_set(band) then begin
        ma_commande = 'band'+string(i)
        IndexScale = i-1
        ; orthogonal transform  
        if trans EQ 3 or trans EQ 4 then  IndexScale = 0 
        ma_commande = ma_commande +'=result[first:last,IndexScale]'
        ma_commande = strcompress( ma_commande, /remove_all)
        ; print, 'cmd = ',  ma_commande
        ACK = EXECUTE( MA_COMMANDE)
     end
   endfor

  ma_commande = 'output = { N_Band: N_Scale,'
  ma_commande = ma_commande+'coef: result,'
  ma_commande = ma_commande+'info: info,'
  ma_commande = ma_commande+'from: from,'
  ma_commande = ma_commande+'to: to,'
  if keyword_set(band) then begin
    for i=1, N_Scale do begin
     ma_commande = ma_commande+'band'+string(i)+':band'+string(i)+','
    endfor
  end
  ma_commande = strcompress( ma_commande, /remove_all)
  ma_commande = strmid(ma_commande,0,strlen(ma_commande)-1)
  ma_commande =ma_commande+'}'
  ; print, ma_commande
  ACK = EXECUTE( MA_COMMANDE)
return
end

;-----------------------------------------------------------------
 
pro mr1d_trans, Signal, result, OPT=OPT, band=band, nodel=nodel

if N_PARAMS() LT 2 then begin 
        spawn, 'mr1d_trans'
        print, 'CALL SEQUENCE: mr1d_trans, Signal, Struct_Out, OPT=Opt, BAND=Band, NODEL=Nodel'
        goto, DONE
        end
        
Nx = (size(Signal))[1]
; print, 'npix = ', nx
NameSig = 'xx_signal.fits'
NameResult = 'xx_result.fits'
NameInfo = 'xx_info.fits'

writefits,  NameSig,  Signal
if not keyword_set(OPT) then OPT = " "
OPT1 = OPT + " -w " + NameInfo + " "

com = 'mr1d_trans ' + OPT1 + ' '+ NameSig  + ' ' +  NameResult
spawn, com
Result = readfits(NameResult, /silent); 
Info =  readfits(NameInfo, /silent);
 
delete, NameSig
if not keyword_set(NODEL) then begin
  delete, NameResult
  delete, NameInfo
END

mr1d_struc, result, info, output, band=band
result = output
 DONE:
end
