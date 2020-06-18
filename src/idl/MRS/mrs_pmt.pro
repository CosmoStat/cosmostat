;+
; NAME:
;        mrs_pmt
;
; PURPOSE:
;	Computes the pyramidal median transform of an Healpix spherical image. 
;
; CALLING:
;
;     mrs_pmt, Imag, Trans, NbrScale=NbrScale
;
; INPUTS:
;     Imag -- IDL array of healpix image: Input image to be transformed
;
; INPUT/OUTPUT:
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;                  NbrScale : int = number of scales 
;                     nside : int = Healpix nside parameter (0 for a Glesp image)
;                      npix : long = Number of pixels of the input image (12*nside*nside)
;                      Scale1 : finest scale (highest frequencies). A IDL array of healpix map or IDL structure of a Glesp map
;                      Scale2 : Second scale    
;                      Scalej : j th scale
;                        ...
;                      ScaleJ : with J = NbrScale, coarsest resolution
;	                   Tab_nside[NbrScale] : int array Tab_nside[j] = nside parameter of the scale j+1, j=0..NbrScale-1 (Healpix input map)
;														 nx number of rings, Glesp parameter of the scale j+1, j=0..NbrScale-1 (Glesp input map)
;
; KEYWORDS:
;      NbrScale  : Number of scales (default is 4). If it is set to -1, then the number scales is:  log(lmax) / log(2)  - 2
;
; EXAMPLE:
;
;       Compute the pyramidal median transform of an image I  
;        The result is stored in Output
;               mrs_pmt, Imag, Output, NbrScale=5
;         
; HISTORY:
;	Written:   Jean-Luc Starck, 2011
;---------------------------------------------------------------------------------------------------------------------------------------------

pro mrs_pmt, Imag, out, NbrScale=NbrScale, WindowSize=WindowSize

COMMON C_PLANCK

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_pmt, Imag, out, NbrScale=NbrScale'
        goto, DONE
        end
	    
GLESP=0
nside=0
nx = 0
np = 0
x_sky = 0
y_sky = 0
Healpix_with_Glesp = 0
pyrtrans = 0

npix = (size(imag))[1]
nside = npix2nside(npix)
if not keyword_set(NbrScale) then NbrScale = 4
if NbrScale  EQ -1 then NbrScale  = ceil(alog(Lmax)/alog(2.)) - 2.

if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end

if not keyword_set(WindowSize) then Win = 3 else Win = WindowSize
BandName = strarr(NbrScale)
Tab_lmax = intarr(NbrScale)
Tab_nside = intarr(NbrScale)

Hscale = imag 							
for j=0,NbrScale-2 do begin
   print, ' Scale ', j+1
   Tab_nside[j] = nside
    LScale = mrs_median( Hscale,  windowsize=Win)
    LScale = mrs_resize( LScale, nside=nside/2)
    LRScale = mrs_resize( LScale, nside=nside)
	WScale = double(Hscale) - double(LRScale)
	Hscale = double(LScale)

    my_command = 'scale'+strcompress(string(j+1), /remove_all)
    BandName[j] = my_command 
    my_command = my_command +'=WScale'
    my_command = strcompress( my_command, /remove_all)
    ; print, 'cmd = ',  my_command
   ACK = EXECUTE( my_command) 
   
    Hscale = LScale
    nside = nside /2
endfor

j=NbrScale-1
Tab_nside[j] = nside

my_command = 'scale'+strcompress(string(j+1), /remove_all)
BandName[j] = my_command 
my_command = my_command +'=LScale'
my_command = strcompress( my_command, /remove_all)
; print, 'cmd = ',  my_command
ACK = EXECUTE( my_command) 


 nside = Tab_nside[0]
TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 
  pyrtrans = 1
  my_command = 'out = { NbrScale : NbrScale,'
  my_command = my_command+'UseGLESP: 0, '
  my_command = my_command+'nside : nside, '
  my_command = my_command+' npix:npix, '
  my_command = my_command+' lmax:0, '
  my_command = my_command+' Tab_lmax: fltarr(NbrScale), '
  my_command = my_command+' Tab_nside:Tab_nside, '
  my_command = my_command+' MeyerWave:0,'
  my_command = my_command+' pyrtrans :pyrtrans,'
  my_command = my_command+' TabNorm :TabNorm,'
  for j=0, NbrScale-1 do  my_command = my_command+BandName[j]+':'+BandName[j]+','
  my_command = my_command+' DifInSH:0'
  my_command = strcompress( my_command, /remove_all)
  my_command =my_command+'}'
  ; print, my_command
  ACK = EXECUTE( my_command)

DONE:

END

