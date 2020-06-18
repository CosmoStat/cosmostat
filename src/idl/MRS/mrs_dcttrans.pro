;+
; NAME:
;        mrs_dcttrans
;
; PURPOSE:
;	Compute the Discrete Cosite Transform (DCT) )on the sphere, 
;       using the healPix pixel representation (NESTED data representation). 
;       The DCT is applied successively on the 12 faces of the Healpix image.
;       The output is an IDL cube [0:nside-1, 0:Nside-1, 0:11].
;       If the keyword HEALPIX is set, then the ouput is an healpix image.
;
; CALLING:
;     mrs_dcttrans, Imag, TransDct, healpix=healpix
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image be transformed 
;    
; OUTPUTS:
;     TransDct -- IDL 3D array [0:nside-1, 0:Nside-1, 0:11]  --    TransDct[*,*,i] is the DCT the ith face of the input map.
;
; KEYWORDS:
;         healpix -- Scalar: if set, the output is a healpix image instead of a cube
;
; EXTERNAL CALLS:
;          
;
; EXAMPLE:
;       Compute the DCT transform of an image I.
;       The result is stored in DCT
;               mrs_dcttrans, Imag, DCT
;               tvscl, DCT[*,*,f] ; plot the fth face wavelet transform (f = 0..11) 
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2010
;-
;-----------------------------------------------------------------

pro mrs_dcttrans, data, TransDct, blockSize=BlockSize, healpix=healpix
COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE:  mrs_dcttrans,  Imag, TransDct,  blockSize=BlockSize, healpix=healpix'
        goto, DONE
        end

nside = gnside(data)
if not keyword_set(BlockSize) then BlockSize=nside

C=0
NbrFace=12
AllFace = h2f(data)
Opt='-b ' + STRCOMPRESS(string(BlockSize), /REMOVE_ALL)  

for f=0,11 do begin
	;if package mre is available
	if keyword_set(mr1ok) then im_dct, AllFace[*,*,f], Dct, opt=Opt  $
	else begin
 		  if f eq 0 then DctIma = dct(AllFace[*,*,f],  T=C) $
 		  else DctIma = dct(AllFace[*,*,f],  C) 
	end

   if f EQ 0 then begin
      vs = size(DctIma)
      nx= vs[1]
      ny= vs[2]
     TransDct = fltarr(nx,ny,12)
   end
   TransDct[*,*,f] = DctIma
end

if keyword_set(healpix) then TransDct = F2H(TransDct)

DONE:
end




