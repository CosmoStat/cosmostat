;+
; NAME: 
;       MRS_CURPUT
;
; PURPOSE: 
;        Insert a band in the curvelet transform (see mrs_curtrans).   
;
; CALLING:
;       MRS_CURPUT, Cur_Struct, Band, ScaleWT2D, ScaleRid
;
; INPUT:
;       Cur_Struct -- IDL structure: Curvelet transform structure (see MRS_CURTRANS)
;       Band      -- 4D IDL float array[*,*,*,12]: Curvelet band   
;       ScaleWT2D -- int:  Scale of the 2D WT
;       ScaleRid  -- int:  Ridgelet band number
;
; OUTPUTS:
;      Cur_Struct --  IDL structure: Curvelet transform structure (see MRS_CURTRANS)
;             
; EXAMPLE:
;   mrs_curtrans, Imag, Cur 
;   Band = mrs_curget(Cur,1,2)
;   mrs_curput, Cur, Band, 1,2
;       Extract and reinsert a band in the curvelet transform  
;
; HISTORY:
;       Written: Jean-Luc Starck 2005.
;       February, 2005 File creation
;-
;-----------------------------------------------------------------

pro mrs_curput, Cur, Band, s2d, s1d

if N_PARAMS() NE 4  then begin 
        print, 'CALLING SEQUENCE: mrs_curput, Cur, Band, s2d, s1d'
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
  
    if Cur.(4+s2d).bin EQ 1 then BEGIN
       my_command = 'Cur.RidScale'+strcompress(string(s2d+1), /remove_all) + '
       C1 = 'DepX = ' + my_command + '.TabDepX[s1d]'
       C2 = 'Nx = '   + my_command + '.TabBandNx[s1d]'
       C1 = strcompress( C1, /remove_all)
       C2 = strcompress( C2, /remove_all)
       ACK = EXECUTE( C1) 
       ACK = EXECUTE( C2) 
       my_command = my_command + '.coef[DepX:DepX+Nx-1,*,*] = Band'
       my_command = strcompress( my_command, /remove_all)
       ACK = EXECUTE( my_command) 
     END ELSE BEGIN
       Cur.(4+s2d).coef.(11+s1d) = Band
     END
  end else Cur.LASTSCALE = Band 
 
 
   DONE:

end

