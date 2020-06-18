;+
; NAME:
;        mrs_wtput
;
; PURPOSE:
;	Put a scale in the wavelet transform obtained by the command mrs_wttrans or by the command mrs_pwttrans.
;
; CALLING:
;
;     mrs_wtput, Trans, Scale, ScaleNumber, Face=Face   
;       
; INPUTS:
;     Trans -- IDL structure: IDL structures containing the wavelet transform 
;     Scale -- IDL array: wavelet scale we want use in the decomposition. 
;     ScaleNumber -- int: Scale number, The scale number must be 
;                     between 0 and Trans.NbrScale-1
;
; KEYWORDS:
;      Face -- int : If set, the Scale is a Cube[*,*,0:11] containg
;                    the twelve faces of the Healpix representation
;
; EXTERNAL CALLS:
;
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	February, 2005 File creation
;
;--------------------------------------------------------------------------------------------------------------------------------

pro mrs_wtput, Trans, Image, ScaleNumber, Face=Face  

if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: mrs_wtput, Trans, Scale, ScaleNumber, Face=Face'
        goto, DONE
        end


if  type_code(Image) EQ 8 then Scale = Image.T_Sky $
else begin
  if keyword_set(Face) then  put_all_faces, Image, Scale $
  else  Scale = Image 
end
    
if ScaleNumber lt 0 or ScaleNumber ge Trans.NbrScale then begin 
                                 print,'Error: Number of scales should be between 0 and ', Trans.NbrScale-1
    	    	    	    	 goto, DONE
				 end
 
 
if Trans.pyrtrans EQ 0 then Trans.coef[*,ScaleNumber] = Scale $
    else begin
        if Trans.useglesp eq 0 then begin
	  my_command = 'Trans.scale'+strcompress(string(ScaleNumber+1), /remove_all)
          my_command =  my_command + '= Scale '  
          my_command = strcompress( my_command, /remove_all)
          ; print, 'cmd = ',  my_command
          ACK = EXECUTE( my_command)
	endif else begin
	  my_command = strcompress('Trans.scale'+string(ScaleNumber+1)+'.t_sky', /remove_all)
          ;ack = execute ('help,'+my_command)
	  my_command =  my_command + '= Scale '  
          ;my_command = strcompress( my_command, /remove_all)
          ; print, 'cmd : ',  my_command
          
	  ACK = EXECUTE( my_command)
	  
	
	endelse
	 
    end 

DONE:
    
end
