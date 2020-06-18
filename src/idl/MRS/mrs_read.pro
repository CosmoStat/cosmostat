;+
; NAME:
;        mrs_read
;
; PURPOSE:
;	
;   Read a spherical map, either in healpix or glesp format
;
;
; CALLING:
;
;     map = mrs_read(file)
;       
; INPUTS:
;    file : file to be read 
;    
; OUTPUTS:
;     
;    either healpix map or glesp structure. For a Healpix map, the map is setted to the NESTED format after reading.
;
; KEYWORDS:
;    none
;
; EXTERNAL CALLS:
;       g_read_fits_map (glesp)
;   	read_fits_map (healpix)
;
; EXAMPLE:
;
;    healmap = mrs_read( 'my_file_healpix.fits' )
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2006
;	February 2006
;
;
;		map = mrsp_read( FileName, noverb=noverb )	Read a HealPIX (no glesp) POLARIZED map, the map is setted to the NESTED format with the same calling as mrs_read.
;													The option /noverb prevent the printing on the screen of the format (RING or NESTED) of the read map and the number of pixels.
;
;
;-------------------------------------------------------------------------------------------

function rimsg, file
  g_read_fits_map,file,t_sky,x_sky,y_sky
  nx = (size(x_sky))[1]
  np = max(y_sky)
  data = {t_sky : t_sky , x_sky : x_sky , y_sky : y_sky , nx : nx , np : np }
  return, Data
end

function rims, FileName, noverb=noverb, glesp=glesp

read_fits_map, FileName, Imag,  hdr, xhdr, /silent
;tab = headfits(FileName, /exten)
;values = FTGET( xhdr, tab, "order")
TypeOrder = FXPAR( xhdr, "ordering")
if type_of(TypeOrder) EQ 'INTEGER' then TypeOrder='RING    '

Order=0
if TypeOrder  eq 'RING    '  or  TypeOrder  eq 'ring    ' then Order=1
if TypeOrder  eq 'NESTED  '  or  TypeOrder  eq 'nested  ' then Order=2
if Order EQ 0 then Order=1

if not keyword_set(noverb) then if Order EQ 1 then print, 'Ring order' else print, 'Nested order'
if Order ne 2 then imag =reorder(imag,in='ring',out='nested')

npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(noverb) then print, "Npix = ", npix, "  nside = " , nside  ; , " SQRT(Npix/12) = ", sqrt(npix/12.)

if keyword_set(glesp) then imag = healpix2glesp(imag)
return, imag
end

;=================================================================================

function mrs_read, file   
  hdr = headfits(file,ext=1)
  size_hdr = (size(hdr))[1]
  
  if size_hdr eq 0 then goto,done
  
  if total(strmatch(hdr,"*healpix*",/fold_case)) ge 1 then  map = rims(file, /noverb)
  if total(strmatch(hdr,"*glesp*",/fold_case)) ge 1 then map = rimsg(file)
  
  ;endfor

  return,map  
  DONE:

END

;=================================================================================

function mrsp_read, FileName, noverb=noverb
 
READ_TQU, FileName, TQU, HDR=Hdr, XHDR=Xhdr , NSIDE=nside, ORDERING=order, COORDSYS=coord

if not keyword_set(noverb) then print, order

if Order ne 'NESTED' then TQU = reorder(TQU,in='ring',out='nested')

return, TQU 
end

;=================================================================================

