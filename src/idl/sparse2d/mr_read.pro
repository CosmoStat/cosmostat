
PRO mr_make_struc, result, Nl, Nc, N_Scale, output
  scale1 = 0 & scale2 = 0 & scale3 = 0 & scale4 = 0
  scale5 = 0 & scale6 = 0 & scale7 = 0 & scale8 = 0
  first = [0,0]
  last  = [Nc-1, Nl-1]
  delta = [Nc, Nl]
  for i=1, N_Scale do begin
     ma_commande = 'scale'+string(i)
     ma_commande = ma_commande +'=result[first(0):last(0), first(1):last(1)]'
     ma_commande = strcompress( ma_commande, /remove_all)
     ;print, ma_commande
     ACK = EXECUTE( MA_COMMANDE)
     first = last + 1
     last = last + (delta+1)/2
     delta = (delta+1)/2
  endfor

  ma_commande = 'output = { N_Scale: N_Scale,'
  for i=1, N_Scale do begin
     ma_commande = ma_commande+'scale'+string(i)+':scale'+string(i)+','
  endfor
  ma_commande = strcompress( ma_commande, /remove_all)
  ma_commande = strmid(ma_commande,0,strlen(ma_commande)-1)
  ma_commande =ma_commande+'}'
  ;print, ma_commande
  ACK = EXECUTE( MA_COMMANDE)
return
end

PRO mr_make_struc1, result, Nli, Nci, N_Scale, output
  scale1 = 0 & scale2 = 0 & scale3 = 0 & scale4 = 0
  scale5 = 0 & scale6 = 0 & scale7 = 0 & scale8 = 0
  Nl = Nli
  Nc = Nci
  for i=1, N_Scale do begin
     ma_commande = 'scale'+string(i)
     ma_commande = ma_commande +'=result[0:Nc-1, 0:Nl-1,i-1]'
     ma_commande = strcompress( ma_commande, /remove_all)
     ;print, ma_commande
     ACK = EXECUTE( MA_COMMANDE)
     if i GT 1 then Nc = (Nc + 1)/ 2
     if i GT 1 then Nl = (Nl + 1)/ 2
    endfor

  ma_commande = 'output = { N_Scale: N_Scale,'
  for i=1, N_Scale do begin
     ma_commande = ma_commande+'scale'+string(i)+':scale'+string(i)+','
  endfor
  ma_commande = strcompress( ma_commande, /remove_all)
  ma_commande = strmid(ma_commande,0,strlen(ma_commande)-1)
  ma_commande =ma_commande+'}'
  print, ma_commande
  ACK = EXECUTE( MA_COMMANDE)
return
end
;+
; NAME:
;         MR_READ
;
; PURPOSE:
;    read a multiresolution file (extension ".mr")
;    If it is a pyramidal transform or half pyramidal transform, 
;    we have 3 different outputs:
;        1) an image  containing several sub-images ( flag raw set)
;        2) a cube  of interpolated or rebinned images ( flag interpol set)
;        3) a structure containing several sub-images  (default)
;
; CALLING:
;
;      output = MR_READ(filename, /interpol, /raw, /debug)
;
; INPUTS:
;     Imag -- string: file name
;
; KEYED INPUTS:
;     interpol -- integer: 0 -> the output will not be interpolated
;                          1 -> the output will be rebinned
;                          2 -> the output will be interpolated
;                  this for pyramidal
;     raw   -- flag : if set the output is the input of the fits file
;
; OUTPUTS:
;     output -- cube of float: the image(s)
;
; CALLED ROUTINES
;   READFITS
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;       Modification : 18-FEB-1997  R Gastaud for fits file

;-

function mr_read, filename, interpol=interpol, debug=debug, raw=raw
output=-1
if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: output=mr_read(filename,/interpol,/raw, /debug)'
        goto, DONE
end

taille = size(findfile(filename))
if ( taille(0) eq 0) then begin
    PRINT,'Error file '+ filename +' not found '
    GOTO, DONE
ENDIF

Result = READFITS( Filename, Header)
ABORT='MR_READ'
Nl               = FXPAR( Header, 'NL'  , ABORT  )
Nc               = FXPAR( Header, 'NC'  , ABORT  )
N_Scale          = FXPAR( Header, 'NBR_PLAN'  , ABORT  )
TypeTransform    = FXPAR( Header, 'Type_Tra'  , ABORT  )
SetTransform     =  FXPAR( Header, 'Set_Tran', ABORT  )
MedianWindowSize = FXPAR( Header, 'MEDIANWI', ABORT  )
FC               = FXPAR( Header, 'FC', ABORT  )
NBR_ITER         = FXPAR( Header, 'NBR_ITER', ABORT  )
EXACTPYR         = FXPAR( Header, 'EXACTPYR', ABORT  )
 
set_transform =['TRANSF_PAVE', 'TRANSF_PYR', 'TRANSF_SEMIPYR', $
                 'TRANSF_MALLAT', 'TRANSF_DIADIC_MALLAT', $
		 'TRANSF_FEAUVEAU']
		 
type_transform = ['TO_PAVE_LINEAR',$
                  'TO_PAVE_BSPLINE',$
                  'TO_PAVE_FFT', $
                  'TM_PAVE_MEDIAN',$
                  'TM_PAVE_MINMAX',$
                  'TO_PYR_LINEAR',$
                  'TO_PYR_BSPLINE',$
                  'TO_PYR_FFT_DIFF_RESOL',$
                  'TO_PYR_FFT_DIFF_SQUARE',$
                  'TM_PYR_MEDIAN',$
                  'TM_PYR_LAPLACIAN',$
                  'TM_PYR_MINMAX',$
                  'TM_PYR_SCALING_FUNCTION',$
                  'TO_MALLAT',$
                  'TO_FEAUVEAU',$
                  'TO_PAVE_FEAUVEAU',$
                  'TM_MIN_MAX',$  ;/* G transform */
                  'TO_HAAR',$
                  'TO_SEMI_PYR',$
                  'TM_TO_SEMI_PYR',$
                  'TO_DIADIC_MALLAT',$ 
                  'TM_TO_PYR',$
                  'TO_PAVE_HAAR',$
                  'TO_LIFTING']

IF KEYWORD_SET(debug) THEN BEGIN
print,format='("Nl = ",I4, "  Nc = ",I4,"  N_Scale = ",I2)', Nl, Nc,N_Scale
print, format='("TypeTransform= ",I2,A," SetTransform = ",I2,A)',$
        TypeTransform,Type_Transform(TypeTransform),SetTransform, $
                                                  set_transform(SetTransform)
print, format='("MedianWindowSize =",I2,"  FC   = ",F6.2," NBR_ITER  = ",I2)',$
       MedianWindowSize, FC, NBR_ITER

ENDIF

IF not KEYWORD_SET( interpol) then interpol = 0


IF KEYWORD_SET( raw) then begin
   output = result

ENDIF else begin

   case SetTransform of
      0: BEGIN ; transformation pave  cube Nl*Nc*N_Scale
         output = result
         END
      1: BEGIN ; transformation pyramidal image 2*Nl*2*Nc

         case interpol of
           0: begin
              mr_make_struc, result, Nl, Nc, N_Scale, output
              end
           1: begin
              output = fltarr(Nc, Nl, N_Scale)
              first = [0,0]
              last  = [Nc-1, Nl-1]
              delta = [Nc, Nl]
              for i=0, N_Scale-1 do begin
                image = result(first[0]:last[0], first[1]:last[1])
                output[*,*,i] = congrid(image, Nc, Nl)
                first = last + 1
                last = last + (delta+1)/2
                delta = (delta+1)/2
              endfor
              end

           2: begin
              output = fltarr(Nc, Nl, N_Scale)
              first = [0,0]
              last  = [Nc-1, Nl-1]
              delta = [Nc, Nl]
              for i=0, N_Scale-1 do begin
 		image = result(first[0]:last[0], first[1]:last[1])
		output[*,*,i] = congrid(image, Nc, Nl,/cubic)
                first = last + 1
                last = last + (delta+1)/2
                delta = (delta+1)/2
              endfor
              end
          endcase
         END
      2: BEGIN ; transformation pyramidal image 2*Nl*2*Nc

         case interpol of
           0: begin
              mr_make_struc1, result, Nl, Nc, N_Scale, output
              end
           1: begin
              output = fltarr(Nc, Nl, N_Scale)
              Ncs  = Nc
              Nls = Nl
              for i=0, N_Scale-1 do begin
                image = result[0:Ncs-1, 0:Nls-1, i]
                output[*,*,i] = congrid(image, Nc, Nl)
		if i GT 0 then BEGIN
		   Ncs = (Ncs+1)/2
		   Nls = (Nls+1)/2
 		END
              endfor
              end
           2: begin
              output = fltarr(Nc, Nl, N_Scale)
	      Ncs = Nc
              Nls = Nl
              for i=0, N_Scale-1 do begin
 		image = result[0:Ncs-1, 0:Nls-1, i]
		output[*,*,i] = congrid(image, Nc, Nl,/cubic)
		if i GT 0 then BEGIN
		   Ncs = (Ncs+1)/2
		   Nls = (Nls+1)/2
 		END
              endfor
              end
          endcase
         END
      3: BEGIN ; transformation mallat image  
         output = result
         END
      4: BEGIN ; transformation diadic mallat image  
         output = result
         END
      5: BEGIN ; transformation feauveau image  
         output = result
         END
      else: print, 'the file ', filename, ' has a bad SetTransform', SetTransform
   endcase 

ENDELSE 
    
;---------------------------------------------------
DONE:

return, output
END
