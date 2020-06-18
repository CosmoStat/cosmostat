;+
; NAME:
;        mrs_wtget 
;
; PURPOSE:
;	Return a scale of the wavelet transform obtained by the command mrs_wttrans or by the command mrs_pwttrans.
;
; CALLING:
;
;     Scale = mrs_wtget( Trans, ScaleNumber, Face=Face, NormVal=NormVal )  
;       
; INPUTS:
;     Trans -- Trans: IDL structures containing the wavelet transform 
;     ScaleNumber -- int: Scale number, The scale number must be 
;                     between 0 and Trans.NbrScale-1
;
; KEYWORDS:
;      Face -- int : If set, the routines retuns a Cube[*,*,0:11] containg
;                    the twelve faces of the Healpix representation
;      NormVal -- float: Normalization value of the band. 
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;
;-------------------------------------------------------------------------------------------------------------------------------------

function mrs_wtget,  Trans, ScaleNumber, Face=Face, NormVal=NormVal

if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: Scale=mrs_wtget(Trans, ScaleNumber, Face=Face)'
        goto, DONE
        end

if ScaleNumber lt 0 or ScaleNumber ge Trans.NbrScale then begin 
                                 print,'Error: Number of scales should be between 0 and ', Trans.NbrScale-1
    	    	    	    	 goto, DONE
				 end
UseGlesp = Trans.UseGlesp
ScaleNumber = fix(ScaleNumber)
 
if Trans.pyrtrans EQ 0 then begin
         Scale = Trans.coef[*,ScaleNumber]
         if UseGlesp EQ  1 then  Scale = { T_SKY: Scale, nx:Trans.nx, np:Trans.np, x_sky:Trans.x_sky, y_sky:Trans.y_sky}
    end else begin
        my_command = 'Trans.scale'+strcompress(string(ScaleNumber+1), /remove_all)
        my_command = 'Scale = ' + my_command 
        my_command = strcompress( my_command, /remove_all)
        ACK = EXECUTE( my_command) 
    end 

TN = Trans.TabNorm
vs = size(TN)
Last = vs[1]
if ScaleNumber LT Last then NormVal = Trans.TabNorm[ScaleNumber] else NormVal = Trans.TabNorm[Last-1]

if keyword_set(Face) then begin
        get_all_faces, Scale, CubeFace
	return, CubeFace
end else return, Scale
    
    
DONE:
   return, -1 
end
