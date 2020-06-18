;+
; NAME:
;        mrsp_wtput
;
; PURPOSE:
;	Put a scale in a transform for polarized data (wavelet, curvelet...) obtained by the command 
;       mrsp_trans.
;
; CALLING:
;
;     mrsp_wtput, Trans, Scale, Component, ScaleNumber, BandNumber=BandNumber
;       
; INPUTS:
;		Trans -- IDL structure: IDL structures containing the coefficients'transform (see mrsp_trans.pro)
;		Scale -- IDL array: wavelet scale we want use in the decomposition. 
;		ScaleNumber -- int: Scale number, The scale number must be 
;                     between 0 and Trans.NbrScale-1
;		Component -- int: choice of the component, 0 is for T, 1 for E and 2 for B
;
; KEYWORDS:
;		BandNumber -- int:  Ridgelet band number (for curvelet transform)
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;------------------------------------------------------------------------------------

pro mrsp_wtput, Trans, Scale, Component, ScaleNumber, BandNumber=BandNumber

if N_PARAMS() LT 4 then begin 
        print, 'CALLING SEQUENCE: mrsp_wtput, Trans, Scale, Component, ScaleNumber, BandNumber=BandNumber'
        goto, DONE
        end

    
if ScaleNumber lt 0 or ScaleNumber ge Trans.NbrScale then begin 
                                 print,'Error: Number of scales should be between 0 and ', Trans.NbrScale-1
    	    	    	    	 goto, DONE
				 end
 
 
 
FirstComponentStructField = Trans.FirstComponentStructField
TransChoice= Trans.TransChoice

if TransChoice EQ 'T_EBDEC' then Trans.( FirstComponentStructField +Component) = Scale
 
; BI-Orthogonal WT
if TransChoice EQ 'T_OWT' then BEGIN
  t = Trans.(FirstComponentStructField +Component)
  mrs_owtput, t, Scale, component, ScaleNumber, BandNumber
  Trans.(FirstComponentStructField +Component) = t
END

; Pyramidal WT
if TransChoice EQ 'T_PyrWT'  then BEGIN
   t = Trans.(FirstComponentStructField +Component)
   mrs_wtput, t, Scale, ScaleNumber
   Trans.(FirstComponentStructField +Component) = t
END

; Undecimated WT
if TransChoice EQ 'T_UWT' then BEGIN
   t = Trans.(FirstComponentStructField +Component)
   mrs_wtput, t, Scale, ScaleNumber
   Trans.(FirstComponentStructField +Component) = t

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
   t = Trans.(FirstComponentStructField +Component)
   mrs_curput, t, Scale, ScaleNumber, BandNumber
   Trans.(FirstComponentStructField +Component) = t
END


DONE:
    
end
