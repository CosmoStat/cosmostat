;+
; NAME:
;        MRSP_WTGET 
;
; PURPOSE:
;	Return a band of a transform for polarized data (wavelet, curvelet...) obtained by the command 
;       mrsp_trans.
;
; CALLING:
;
;     Band = MRSP_WTGET( Trans, Component, ScaleNumber, BandNumber=BandNumber, NormVal=NormVal )  
;       
; INPUTS:
;		Trans -- Trans: IDL structures containing the coefficients'transform (see mrsp_trans.pro)
;		ScaleNumber -- int: Scale number, The scale number must be 
;                     between 0 and Trans.NbrScale-1
;		Component -- int: choice of the component, 0 is for T, 1 for E and 2 for B
;     
;
; KEYWORDS:
;		NormVal -- float: Normalization value of the band.
;		BandNumber -- int:  Ridgelet band number (for curvelet transform)
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2008
;	May, 2008 File creation
;---------------------------------------------------------------------------------------------------------

function mrsp_wtget, Trans, Component, ScaleNumber, BandNumber=BandNumber, NormVal=NormVal

if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: Scale = mrsp_wtget( Trans, Component, ScaleNumber, BandNumber=BandNumber, NormVal=NormVal )'
        goto, DONE
        end

if ScaleNumber lt 0 or ScaleNumber ge Trans.NbrScale then begin 
                                 print,'Error: Number of scales should be between 0 and ', Trans.NbrScale-1
    	    	    	    	 goto, DONE
				 end

if component lt 0 or component GT 2 then begin 
                                 print,'Error: Component should be between 0 and 2' 
    	    	    	    	 goto, DONE
				 end

NormVal = 1.

FirstComponentStructField = Trans.FirstComponentStructField
TransChoice= Trans.TransChoice

if TransChoice EQ 'T_EBDEC' then Scale = Trans.( FirstComponentStructField +Component)
 
; BI-Orthogonal WT
if TransChoice EQ 'T_OWT' then BEGIN
   Scale = mrs_owtget(Trans.(FirstComponentStructField +Component), component, ScaleNumber, BandNumber)
END

; Pyramidal WT
if TransChoice EQ 'T_PyrWT'  then BEGIN
   scale = mrs_wtget(Trans.(FirstComponentStructField +Component), ScaleNumber, NormVal=NormVal)
END

; Undecimated WT
if TransChoice EQ 'T_UWT' then BEGIN
   scale = mrs_wtget(Trans.(FirstComponentStructField +Component), ScaleNumber, NormVal=NormVal)
END

; Module-Phase Decimated Transform
if TransChoice EQ 'T_MPDWT' then BEGIN
   Scale = 0
END

; Module-Phase Undecimated Transform
if TransChoice EQ 'T_MPUWT' then BEGIN
   Scale = 0
END

; Curvelet
if TransChoice EQ 'T_CUR' then BEGIN
   Scale = mrs_curget(Trans.(FirstComponentStructField +Component), ScaleNumber, BandNumber)
END

return, Scale
    
    
DONE:
    
end
