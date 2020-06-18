;+
; NAME: 
;       MRS_CURGET
;
; PURPOSE: 
;        Extract a curvelet band from the curvelet transform (see mrs_curtrans).   
;        If the keyword NormMad is set, a normalization is applied to all coefficients.   
;
; CALLING:
;       result = MRS_CURGET( Cur_Struct, ScaleWT2D, ScaleRid, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad )
;
; INPUT:
;       Cur_Struct -- IDL structure: Curvelet transform structure (see MRS_CURTRANS) 
;       ScaleWT2D -- int:  Scale of the 2D WT
;       ScaleRid  -- int:  Ridgelet band number
;       NormMad -- scalar: if set, normalized the coefficients by the Median Absolution Deviation
;                          of all coefficients at a give position in the block.
;       ImaMad  -- 2D fltarr: Image containing the normalization parameter
;       ImaMean -- 2D fltarr: Image containing  the mean value for all coefficients at a given position.
;
; OUTPUTS:
;      Result -- 4D IDL float array[*,*,*,12]: Curvelet band
;             
; EXTERNAL CALLS
;           mrs_ridget
;
; EXAMPLE:
;   mrs_curtrans, Imag, Cur 
;   Band = mrs_curget(Cur, 1,2)
;       Extract the the third ridgelet scale inside the second 2D wavelet scale
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       September, 2005 File creation
;-
;-----------------------------------------------------------------

function mrs_curget, Cur, s2d, s1d, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad

  if N_PARAMS() NE 3  then begin 
        print, 'CALLING SEQUENCE: Band = mrs_curget(Cur, s2d, s1d, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad)'
        goto, DONE
        end
  if s2d LT 0 or s2d GE Cur.NBRSCALE then begin
     print, 'Bad first parameter: ', s2d, ', Min Val = 0', ' Max Val = ', Cur.NBRSCALE-1
     goto, DONE
        end

  if s2d LT Cur.NBRSCALE-1 then begin
    if s1d LT 0 or s1d GE Cur.TABNBRSCALERID[s2d]  then begin
     print, 'Bad second parameter: ', s1d, ', Min Val = 0', ' Max Val = ', Cur.TABNBRSCALERID[s2d]-1
     goto, DONE
     end
  
    my_command = 'Band = mrs_ridget(Cur.RidScale' + strcompress(string(s2d+1), /remove_all) + ', s1d, NormMad=NormMad, ImaMean=ImaMean, ImaMad=ImaMad)'
    my_command = strcompress( my_command, /remove_all)
    ACK = EXECUTE( my_command)   
  end else Band = Cur.LASTSCALE
  
  return, Band

DONE:
  
end
