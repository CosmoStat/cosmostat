;+
; NAME: 
;       MR1D_RECONS
;
; PURPOSE: 
;        Reconstruct a one dimensional signal from its wavelet transform.   
;
; CALLING:
;       MR1D_RECONS, WT_Struct, result
;                   
;
; INPUT:
;       WT_Struct : IDL structure; Wavelet transform structure (see MR1D_TRANS) 
;          
; OUTPUTS:
;      Result:  1D array
;             
;
; EXTERNAL CALLS
;           mr1d_recons (C++ program)
;
; EXAMPLE:
;
;   > mr1d_trans, Signal, WT, OPT='-t 13 -n 4'
;   > mr1d_recons, WT, Rec
;       Orthogonal wavelet transform with 4 scales, and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck 1998.
;       July, 1998 File creation
;-

 
;-----------------------------------------------------------------
 
pro mr1d_recons, WT, Rec

if N_PARAMS() LT 2 then begin 
        spawn, 'mr1d_recons'
        print, 'CALL SEQUENCE: mr1d_recons, WT, Rec_Signal'
        goto, DONE
        end
        
Nx = (size(Signal))[1]
NameSig = 'xx_signal.fits'
NameWT = 'xx_coef.fits'
NameInfo = 'xx_info.fits'

writefits,  NameWT,  WT.Coef
writefits,  NameInfo,  WT.Info
com = 'mr1d_recons ' + NameWT  + ' ' +  NameInfo + ' ' + NameSig
spawn, com
Rec = readfits(NameSig, /silent); 
  
delete, NameSig
delete, NameWT
delete, NameInfo
 
DONE:
end
