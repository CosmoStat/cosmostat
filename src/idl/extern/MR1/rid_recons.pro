;+
; NAME: 
;       RID_RECONS
;
; PURPOSE: 
;        Reconstruct an image from its ridgelet transform.   
;
; CALLING:
;       RID_RECONS, Rid_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Ridgelet transform structure (see RID_TRANS) 
;          
; OUTPUTS:
;      Result:  2D array
;             
; EXTERNAL CALLS
;           rid_trans (C++ program)
;
; EXAMPLE:
;   rid_trans, Imag, Rid 
;   rid_recons, Rid, RecIma
;       Ridgelet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       February, 2005 File creation
;-

 
;-----------------------------------------------------------------
 
pro rid_recons, Rid, Rec 

if N_PARAMS() LT 2 then begin 
        spawn, 'rid_recons'
        print, 'CALL SEQUENCE: rid_recons, Rid_Struct, result'
        goto, DONE
        end
        
NameSig = 'xx_signal.fits'
NameRid = 'xx_coef.rid'
 
writefits,  NameRid,  Rid.Coef, Rid.HEADTRANS
com = 'rid_trans -r ' + NameRid  + ' ' +  NameSig
spawn, com
Rec = readfits(NameSig, /silent); 
  
delete, NameSig
delete, NameRid 
DONE:
end
