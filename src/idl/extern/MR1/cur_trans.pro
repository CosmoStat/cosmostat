;+
; NAME: 
;       CUR_TRANS
;
; PURPOSE:
;       Computes a curvelet transform of an image. A multiresolution file
;       has a .cur extension, if the parameter file name have not one, then
;       the extension is added to the file name. Result is store in DataTransf.
;       Dependiong on the options (see OPT), DataTransf can be a cube or an
;       image. If the keyword is not set, the created file is deleted.
;
; CALLING:
;     cur_trans, Imag, Trans, Opt=Opt 
;
; INPUTS:
;     Imag -- IDL 2D array: Input image be transformed 
;     
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;         NXRID     -- LONG: Number of pixels in the x-axis direction
;         NYRID     -- LONG: Number of pixels in the y-axis direction
;         COEF      -- 2D IDL array: Ridgelet coefficients
;                             Image which contains all ridgelet coefficients
;                        COEF[ TabDepX[j]:TabBandNx[j]-1, *] are the ridgelet 
;                        coefficients at the scale j
;         BSIZE     --  LONG: Block size used in the ridgelet transform
;         OVERLAP   -- LONG: is equal to 1 if blocks are overlapping
;         HEADTRANS -- IDL STRING Array: contains the Fits Header of the decomposition
;
; KEYWORDS:
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;
;          [-t type_of_ridgelet]
;              1: RectoPolar Ridgelet Transform using a standard bi-orthogonal WT 
;              2: RectoPolar Ridgelet Transform using a FFT based Pyramidal WT
;              3: Finite ridgelet transform
;              4: RectoPolar Ridgelet Transform using the Slant Stack Radon transform 
;              Default is RectoPolar Ridgelet Transform using a FFT based Pyramidal WT.
;   
;           [-n number_of_scales]
;                number of scales used in the wavelet transform
;                default is 4
;
;           [-N number_of_scales]
;                number of scales used in the Curgelet transform
;                default is automatically calculed
;
;           [-b BlockSize]
;                Block Size.
;                default is 16. 
;
;           [-r]
;                Inverse Curvelet transform.
;
;           [-i]
;                Print statistical information about.
;                each band. Default is no. 
;
;           [-O]
;                Use overlapping block. Default is no.  
;
;           [-x]
;                Write all bands separately as images with prefix 'band_j' (j being the band number)
;          
;           [-z]
;                Use virtual memory.
;                default limit size: 4
;                default directory: 
;
;           [-Z VMSize:VMDIR]
;                Use virtual memory.
;                VMSize = limit size (megabytes) 
;                VMDIR = directory name 
;
;           [-v]
;                Verbose. Default is no.
;
; EXTERNAL CALLS:
;       cur_trans (C++ program)
;
; EXAMPLE:
;
; HISTORY:
;       Written : Sandrine PIRES 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

pro cur_trans, Ima, CurTrans, OPT=OPT 

if N_PARAMS() LT 2 then begin 
        spawn, 'cur_trans'
        print, 'CALL SEQUENCE: cur_trans, Signal, Struct_Out, OPT=Opt'
        goto, DONE
        end

Nx = (size(Ima))[1]
Ny = (size(Ima))[2]

NameIma = 'xx_signal.fits'
NameResult = 'xx_result'
NameResultFits = 'xx_result.cur'

writefits,  NameIma,  Ima
if not keyword_set(OPT) then OPT = ' '
 
com = 'cur_trans'+' '+ OPT + ' '+ NameIma  + ' ' +  NameResult
spawn, com

Cur = readfits(NameResultFits, HeadTrans, /silent); 

Nxcur = (size(Cur))[1]
Nycur = (size(Cur))[2]
Nzcur = (size(Cur))[3]
Nl = FXPAR(HeadTrans, "NL")
NC = FXPAR(HeadTrans, "NC")
Bsize = FXPAR( HeadTrans, "BSIZE")
Overlap = FXPAR(HeadTrans, "OVERLAP")
TypeTrans = FXPAR(HeadTrans, "TYPE_TRA")                               

CurTrans = {Nxcur : Nxcur, Nycur:Nycur, Nzcur:Nzcur,Coef : Cur, Bsize : Bsize, $
            Overlap: Overlap, HeadTrans:HeadTrans}
DONE:
 
end
