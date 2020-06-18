;+
; NAME: 
;       MRS_RIDPUT
;
; PURPOSE: 
;        Insert a band in the ridgelet transform (see mrs_ridtrans).   
;
; CALLING:
;       MRS_RIDPUT, Rid_Struct, Band, ScaleRid
;
; INPUT:
;       Rid_Struct -- IDL structure: Ridgelet transform structure (see MRS_RIDTRANS)
;		Band -- 4D IDL float array[*,*,*,12]: Ridgelet coefficients band to insert.
;       ScaleRid  -- int: Ridgelet band number
;
; OUTPUTS:
;       Rid_Struct -- IDL structure: Ridgelet transform structure (see MRS_RIDTRANS) 
;             
; EXAMPLE:
;   mrs_ridtrans, Imag, Rid 
;   Band = mrs_ridget(Rid,1)
;   mrs_ridput, Rid, Band, 1
;       Extract and reinsert a band in the ridgelet transform  
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

pro mrs_ridput, Rid, Band, j

if N_PARAMS() NE 3  then begin 
        print, 'CALLING SEQUENCE: mrs_ridput, Rid, Band, j'
        goto, DONE
        end

 if j LT 0 or j GE Rid.NBRSCALE then begin
     print, 'Bad second parameter: ', j, ', Min Val = 0', ' Max Val = ', Rid.NBRSCALE-1
     goto, DONE
     end
     	
  if Rid.bin EQ 1 then Rid.coef[ Rid.TabDepX[j]:Rid.TabDepX[j]+Rid.TabBandNx[j]-1,*,*] = Band $
  else rid.coef.(11+j) = Band
  
  DONE:

end

