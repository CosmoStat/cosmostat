pro rdfits_struct, filename, struct, NODELETE=nodelete, SILENT = silent 
;+
; NAME:
;      RDFITS_STRUCT
; PURPOSE:
;      Read an entire FITS file (all extensions) into a single IDL structure. 
; EXPLANATION:
;      Each header, image or table array is placed in a separate structure 
;      tag.
;
; CALLING SEQUENCE:
;      RDFITS_STRUCT, filename, struct, [ /NODELETE, /SILENT ]
;
; INPUT:
;      FILENAME = Scalar string giving the name of the FITS file
;
; OPTIONAL KEYWORD: 
;      /NODELETE -  RDFITS_STRUCT creates a temporary file with the name 
;              temp_'fitsname'.pro which contains the IDL structure definition
;              Normally, this temporary file is deleted -- set the /NODELETE
;              keyword to keep it. 
;      /SILENT - Set this keyword to suppress informational displays at the
;               terminal.
;
; OUTPUT:
;      struct = structure into which FITS data is read.   The primary header
;             and image are placed into tag names HDR0 and IM0.   The ith
;             extension is placed into the tag names HDRi, TABi
;
; PROCEDURES USED:
;       FDECOMP, FITS_INFO, HEADFITS(), GETTOK(), READFITS(), STRN()
;
; METHOD:
;       The procedure FITS_INFO is used to determine whether a primary image
;       exists and the number of extensions.     The number and type of
;       structure tags required is written to a temporary file and assigned
;       to an appropiate HEADFITS or READFITS call.     The temporary file
;       is executed using CALL_PROCEDURE.   
;
; EXAMPLE:
;       Read the FITS file 'm33.fits' into an IDL structure, st
;
;       IDL> rdfits_struct, 'm33.fits', st
;       IDL> help, /str, st                   ;Display info about the structure
;
; RESTRICTIONS:
;       The current algorithm is not particularly efficient. 
;
;       Does not handle random groups
; MODIFICATION HISTORY:
;       Written K. Venkatakrishna, STX April 1992
;       Code cleaned up a bit  W. Landsman  STX  October 92
;       Modified for MacOS     I.  Freedman  HSTX April 1994
;       Work under Windows 95  W. Landsman   HSTX  January 1996
;       Use anonymous structures, skip extensions without data WBL April 1998
;       Converted to IDL V5.0, W. Landsman, April 1998
;       OS-independent deletion of temporary file  W. Landsman  Jan 1999
;-
 COMMON descriptor, file_descript

 if N_Params() LT 2 then begin 
        print,'Syntax - RDFITS_STRUCT, file, struct, [ /NODELETE, /SILENT ]'
        return
 endif

 fits_info, filename,/SILENT                ; Get the description of the file
 Nhdr_prim = gettok( file_descript, ' ')     ;Get header size & # of dimensions 

 Naxis_prim = gettok( file_descript, ' ')    ;of primary image, if any


 if Naxis_prim GT 0 then begin $
      for i = 0, Naxis_prim do Naxisn_im = gettok(file_descript, ' ')
 endif

 struct_im = ''

 if Naxis_prim GT 0 then begin 
          im_string = ", im0 : readfits( '"+filename + "' )" + " $"
 endif else im_string = ''
                
 fname = filename
 fdecomp, fname, disk, dir, fitsname

;   Form a string to read primary header

 hd_string = ["struct = {hdr0 : headfits( '" + fname + "' ) $" ]
                                                 
 N_ext = 1 

;   If there are extensions, form the strings defining the corresponding tags

 while strlen(file_descript) GT 0 do begin
   h_size = gettok(file_descript,' ')

    table_type = gettok( file_descript,' ')
  naxis = gettok( file_descript, ' ')

  if naxis GT 0 then im_string = [im_string,  ', tab' + strn(N_ext) + ' :' + $
        'readfits(' + '"' + fname + '", EXTEN_NO = '+ strn(N_ext) + ') $ ' ] 

   hd_string = [ hd_string, ", hdr"+ strn(N_ext)+ " : " + $
            'headfits(' +  '"' + fname + '", EXT = ' + strn(N_ext) + ') $ ' ] 

   if Naxis GT 0 then for i = 0, Naxis do $ 
                    Naxisn = gettok( file_descript,' ')

  N_ext = N_ext + 1
endwhile

;  Write the structure code into a file named temp_struct.pro 

 openw, unit, 'temp_' + fitsname+'.pro',/get_lun
 printf, unit, 'pro '+ ' temp_'+fitsname+' ,struct'

 for j = 0, N_elements( hd_string ) - 1 do $
        printf, unit, strtrim( hd_string[j] )

 for j = 0, N_elements( im_string ) - 1 do $
        printf, unit, strtrim( im_string[j] )  

 printf, unit, '        }'
 printf, unit, 'return'
 printf, unit, 'end'

 free_lun, unit
;                                      
 if not keyword_set( SILENT) then $ 
         message,' Creating a structure from the file ' + fname, /INF
 call_procedure, 'temp_' + fitsname, struct
 if not keyword_set( SILENT ) then help, /str, struct

 if keyword_set( NODELETE ) then begin 
       message,'File temp_' + fitsname+'.pro has been created',/INF  
       return 
  endif

; Delete the file temp_struct.pro in an OS-independent way

 openr,unit, 'temp_' + fitsname + '.pro',/delete
 close,unit
                             
 return
 end
