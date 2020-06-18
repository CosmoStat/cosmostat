;+
; NAME: 
;       BEAM3D_RECONS
;
; PURPOSE: 
;        Reconstruct a cube from its beamlet transform.   
;
; CALLING:
;       BEAM3D_RECONS, Beam_Struct, result
;
; INPUT:
;       WT_Struct : IDL structure; Beamlet transform structure (see BEAM3D_TRANS) 
;          
; OUTPUTS:
;      Result:  3D array
;             
; EXTERNAL CALLS
;           beam3d_trans (C++ program)
;
; EXAMPLE:
;   beam3d_trans, Cube, Beam 
;   beam3d_recons, Beam, RecCube
;       Beamlet transform and reconstruction
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       September, 2005 File creation
;-

 
;-----------------------------------------------------------------
 
pro beam3d_rec, Beam, Rec 

if N_PARAMS() LT 2 then begin 
        spawn, 'beam3d_recons'
        print, 'CALL SEQUENCE: beam3d_rec, Beam_Struct, result'
        goto, DONE
        end
        
NameSig = 'xx_signal.fits'
NameBeam = 'xx_coef.bet'
 
writefits,  NameBeam,  Beam.Coef, Beam.HEADTRANS
com = 'beam3d_trans -r ' + NameBeam  + ' ' +  NameSig
spawn, com
Rec = readfits(NameSig, /silent); 
  
delete, NameSig
delete, NameBeam 
DONE:
end
