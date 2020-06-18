;+
; NAME: 
;       CUR_RECONS
;
; PURPOSE: 
;        Reconstruct an image from its ridgelet transform.   
;
; CALLING:
;       CUR_RECONS, Rid_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Ridgelet transform structure (see RID_TRANS) 
;          
; OUTPUTS:
;      Result:  2D array
;             
; EXTERNAL CALLS
;           cur_trans (C++ program)
;
; EXAMPLE:
;   cur_trans, Imag, Rid 
;   cur_recons, Rid, RecIma
;       Curvelet transform and reconstruction
;
; HISTORY:
;       Written: Sandrine PIRES 2005.
;       March, 2005 File creation
;-

 
;-----------------------------------------------------------------
 
pro cur_recons, Cur, Rec

if N_PARAMS() LT 2 then begin 
        spawn, 'cur_recons'
        print, 'CALL SEQUENCE: cur_recons, Cur_Struct, result'
        goto, DONE
        end
  
     
NameSig = 'xx_signal.fits'
NameCur = 'xx_coef.cur'
 
writefits,  NameCur,  Cur.Coef, Cur.HEADTRANS
com = 'cur_trans -r '+ ' ' + NameCur  + ' ' +  NameSig
spawn, com
Rec = readfits(NameSig, /silent); 
  
delete, NameSig
delete, NameCur
DONE:
end
