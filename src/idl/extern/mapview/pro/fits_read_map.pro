Pro FITS_READ_Map, File_Name, T, N_obs, EHeader, Pheader = pheader
;+
; NAME:
;      FITS_READ_MAP
; PURPOSE:
;      Read a WMAP temperature map FITS file
; EXPLANATION: 
;      The map is stored as a FITS binary table, which is read using FITS_READ.
;      The map resolution, number of pixels, pixel type, and other information 
;      is  given in the FITS header.
; 
; CALLING SEQUENCE:
;      FITS_READ_Map, File_Name, T, N_obs, [ EHeader, PHeader= ]
; INPUTS:
;       File_Name - Character name of archive map FITS file, scalar string
;                   A .fits extension is optional.    Compressed files (with a
;                   .gz extension) are allowed.
;                  
; OUTPUTS:
;       T  -  R*4 array with temperature of each pixel
;       N_obs   R*4 array with number of observations of each pixel
; OPTIONAL OUTPUT:
;       EHeader - FITS header with map information, string array.  This is
;	          the FITS header associated with the FITS binary table
;	          that contains the map.
; OPTIONAL KEYWORD OUTPUT:
;       Pheader - Primary FITS header, string array
; EXAMPLE:
;       Read the first year K band map 
;       IDL> fits_read_map,'map_k_imap_yr1_v1',t,N_obs,h
;
; PROCEDURES USED:
;       FITS_READ, TBINFO, TBGET()                ;IDL Astronomy Library
; REVISION HISTORY:
;       -- Al Kogut, RSTX       February 6, 1998
;        Use FITS_READ   W. Landsman    January 2003
;        Only read N_OBS if parameter supplied  W. Landsman February 2003
;        Added Pheader keyword, close file unit W. Landsman April 2003
;-----------------------------------------------------------------------------
on_error, 2
;
If N_params() LT 2 Then Begin
    print,'Syntax - FITS_READ_Map, File_Name, T [, N_obs, EHeader, Pheader=]'
    Return
EndIf
;
;   If we can't open the file name directly then try a .fits extension
;
fits_open,file_name,fcb, /no_abort,message = message
If (message NE '') Then $
       fits_open,file_name + '.fits',fcb,message = message 
If (arg_present(pheader)) Then fits_read,fcb,dummy,pheader,/header,exten=0     
fits_read,fcb, data, EHeader, exten_no=1
fits_close,fcb
tbinfo, EHeader, tb_str
T = tbget(tb_str,data,'TEMPERATURE')
If (N_params() GT 2) Then N_obs= tbget(tb_str,data,'N_OBS')
;
Return
End
