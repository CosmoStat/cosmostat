;+
; NAME:
;       CB_MCA
;
; PURPOSE:
;	      Inpainting by decomposition of an image on multiple bases
;
; CALLING:
;
;      CB_MCA, Imag, Struct, OPT=Opt
;       
;
; INPUTS:
;     Imag -- 2D IDL array: image we want to decompose 
;
; OUTPUTS:
;     Result -- Image inpainted
;     Struct -- Decomosition of Imag in each of the base
;
; KEYWORDS:
;
;      Opt -- string: string which contains the differents options. Options are:
;   where options =  
;         [-t TransformSelection]
;              1: A trous algorithm 
;              2: bi-orthogonal WT with 7/9 filters 
;              3: Ridgelet transform 
;              4: Curvelet transform 02 
;              5: Local Discrete Cosinus Transform 
;              6: Wavelet Packet basis 
;              7: Curvelet transform 05 (Pyramidal Construction) 
;              8: Curvelet transform 04 (Fast Transform) 
;
;         [-n number_of_scales]
;             Number of scales used in the WT, the a trous, the PMT and the curvelet transform.
;             default is 5.
;
;         [-b BlockSize]
;             Block Size in the ridgelet transform.
;             Default is image size. 
;
;         [-i NbrIter]
;             Number of iteration. Default is 30.
;
;         [-B DCT_BlockSize]
;             Local DCT block size.
;             By default, a global DCT is used. 
;
;         [-S FirstThresholdLevel]
;            First thresholding value.
;            Default is derived from the data.
;
;         [-s LastThresholdLevel]
;            Last thresholding value..
;            default is 3.000000.
;
;         [-N]
;             Minimize the L1 norm.
;             Default is L0 norm.
;
;         [-L]
;             Replacing the linear descent by a non linear descent.
;             (Should be used when one components is much larger than another one).
;
;         [-l]
;             Remove last scale. Default is no. 
;
;         [-g sigma]
;             sigma = noise standard deviation. 
;             Default is 0.5 (quantification noise).
;             if sigma is set to 0, noise standard deviation
;             is automatically estimated.
;
;         [-O]
;             Supress the block overlapping. Default is no. 
;
;         [-P]
;             Supress the positivity constraint. Default is no. 
;
;         [-G RegulVal[,NbrScale]]
;             Total Variation regularization term. Default is 0.
;             NbrScale = number of scales used in Haar TV regularizarion. 
;                        default is 2. 
;
;         [-H]
;             Data contained masked area (must have a zero value). Default is no. 
;
;         [-I]
;             Interpolate the data (super-resolution). Default is no. 
;
;         [-z]
;             Use virtual memory.
;                default limit size: 4
;                default directory: .
;
;         [-Z VMSize:VMDIR]
;             Use virtual memory.
;                VMSize = limit size (megabytes) 
;                VMDIR = directory name 
;
;         [-v]
;             Verbose. Default is no.
;
;      
; EXTERNAL CALLS:
;       cb_mca (C++ program)
;
; EXAMPLE:
;
;       Decomposition of an image I into 2 bases (Curvelet and DCT with block size of 16).
; 			The result is stored in the strcture Struct
;               CB_MCA, I, Struct, OPT='-t4 -t5 -B16'
;
; HISTORY:
;	Written: Sandrine Pires 2006.
;-

pro CB_MCA, Imag, Result, Struct, OPT=Opt

if N_PARAMS() LT 2 then begin 
        spawn, 'cb_mca'
        print, 'CALLING SEQUENCE: cb_mca, Imag,  Result, Struct,OPT=Opt'
        goto, DONE
        end

vsize = size(Imag)
if vsize(0) NE 2 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: cb_mca, Imag,  Result, Struct, OPT=Opt'
        goto, DONE
        end

deltmpmr=0
if not keyword_set(Opt) then Opt = ' '

filename = 'Result'
deltmpmr=1 

NameImag=strcompress('tmp'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all) 
writefits, NameImag, Imag

com = 'cb_mca ' + Opt + ' ' + NameImag + ' ' +  filename
spawn, com

sz = size(Imag)  
n1= sz[1]
n2= sz[2]

Result = readfits(filename+'.fits', header)

Imag1 = readfits(filename+'_atrou'+'.fits', /silent)
Imag2 = readfits(filename+'_wt'+'.fits', /silent)
Imag3 = readfits(filename+'_rid'+'.fits', /silent)
Imag4 = readfits(filename+'_cur'+'.fits', /silent)
Imag5 = readfits(filename+'_cos'+'.fits', /silent)
Imag6 = readfits(filename+'_wp'+'.fits', /silent)
Imag7 = readfits(filename+'_pcur'+'.fits', /silent)
Imag8 = readfits(filename+'_fcur'+'.fits', /silent)
Imag9 = readfits(filename+'_resi'+'.fits', /silent)
if total(Imag1) eq -1 then Imag1 =fltarr(n1, n2)
if total(Imag2) eq -1 then Imag2 =fltarr(n1, n2)
if total(Imag3) eq -1 then Imag3 =fltarr(n1, n2)
if total(Imag4) eq -1 then Imag4 =fltarr(n1, n2)
if total(Imag5) eq -1 then Imag5 =fltarr(n1, n2)
if total(Imag6) eq -1 then Imag6 =fltarr(n1, n2)
if total(Imag7) eq -1 then Imag7 =fltarr(n1, n2)
if total(Imag8) eq -1 then Imag8 =fltarr(n1, n2)



Struct={im1:Imag1,im2:Imag2,im3:Imag3,im4:Imag4,im5:Imag5,im6:Imag6,im7:Imag7,im8:Imag8,im9:Imag9, header:header}


if deltmpmr EQ 1 then begin
  delete, NameImag
  delete, filename+'_atrou'+'.fits'
  delete, filename+'_wt'+'.fits'
  delete, filename+'_rid'+'.fits'
  delete, filename+'_cur'+'.fits'
  delete, filename+'_cos'+'.fits'
  delete, filename+'_wp'+'.fits'
  delete, filename+'_pcur'+'.fits'
  delete, filename+'_fcur'+'.fits'
  delete, filename+'_resi'+'.fits'
endif

DONE:

END
