;+
; NAME: 
;       RID3D_RECONS
;
; PURPOSE: 
;        Reconstruct a cube from its ridgelet transform.   
;
; CALLING:
;       RID3D_RECONS, Rid_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Ridgelet transform structure (see RID3D_TRANS) 
;          
; OUTPUTS:
;      Result:  3D array
;             
; EXTERNAL CALLS
;           rid3d_trans (C++ program)
;
; EXAMPLE:
;   rid3d_trans, Cube, Rid 
;   rid3d_recons, Rid, RecCube
;       Ridgelet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       September, 2005 File creation
;-

 
;-----------------------------------------------------------------
 
pro rid3d_rec, Rid, Rec 

if N_PARAMS() LT 2 then begin 
        spawn, 'rid3d_recons'
        print, 'CALL SEQUENCE: rid3d_recons, Rid_Struct, result'
        goto, DONE
        end
        
NameSig = 'xx_signal.fits'
NameRid = 'xx_coef.rid'
 
writefits,  NameRid,  Rid.Coef, Rid.HEADTRANS
com = 'rid3d_trans -r ' + NameRid  + ' ' +  NameSig
spawn, com
Rec = readfits(NameSig, /silent); 
  
delete, NameSig
delete, NameRid 
DONE:
end
