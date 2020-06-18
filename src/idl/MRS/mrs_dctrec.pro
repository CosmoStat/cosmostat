;+
; NAME:
;        mrs_dctrec
;
; PURPOSE:
;	Compute the inverse Discrete Cosite Transform (DCT) ) on the sphere, 
;       using the healPix pixel representation (NESTED data representation). 
;       The inverse DCT is applied successively on the 12 faces of the Healpix image.
;       The input is an IDL cube [0:nside-1, 0:Nside-1, 0:11] or an healpix image is the healpix keyword is set.
;       The output is an healpix image
;
; CALLING:
;     mrs_dctrec, TransDct,  RecIma, healppix= healppix
;
; INPUTS:
;     TransDct -- IDL 3D array [0:nside-1, 0:Nside-1, 0:11]  --    TransDct[*,*,i] is the DCT the ith face of the input map.
;
; OUTPUTS:
;     RecIma -- IDL array of healpix map: Input image be transformed 
;
; KEYWORDS:
;         healpix -- Scalar: if set, the input is a healpix image instead of a cube
;
; EXTERNAL CALLS:
;          im_dct.pro  -- MR1 IDL program
;
; EXAMPLE:
;       Compute the DCT transform of an image I.
;       The result is stored in DCT
;               mrs_dcttrans, Imag, DCT
;               tvscl, DCT[*,*,f] ; plot the fth face wavelet transform (f = 0..11) 
;      Reconstruction
;              mrs_dctrec, DCT, RecIma
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2010
;-
;-----------------------------------------------------------------
  
pro mrs_dctrec,  InputTransDct, Data, blockSize=BlockSize, healpix=healpix
COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE:  mrs_dctrec,  InputTransDct, RecIma,  blockSize=BlockSize, healpix=healpix'
        goto, DONE
        end

if keyword_set(healpix) then TransDct = H2F(InputTransDct) $
else TransDct = InputTransDct

vs = size(TransDct)
nside = vs[1]
if not keyword_set(BlockSize) then BlockSize=nside

NbrFace=12
Opt='-r -b ' + STRCOMPRESS(string(BlockSize), /REMOVE_ALL)  
for f=0,11 do begin
	if keyword_set(mr1ok) then  im_dct, TransDct[*,*,f], Rec, opt=Opt  $
	else begin
	      if f eq 0 then Rec = DCT(TransDct[*,*,f],  T=C, /inverse) $ 
 		  else Rec = DCT(TransDct[*,*,f],  C, /inverse) 
	end


if f EQ 0 then begin
   vs = size(Rec)
   nx= vs[1]
   ny= vs[2]
   AllFace = fltarr(nx,ny,12)
   end
AllFace[*,*,f] = Rec
end

data = f2h(AllFace)

DONE:

end


pro testdec
COMMON MR1ENV

n = getcmb()

mr1ok=0
mrs_dcttrans, n, d, /healpix
mrs_dctrec, d, r, /healpix
info, n-r

mr1ok=1
mrs_dcttrans, n, d, blocksize=16
mrs_dctrec, d, r, blocksize=16
info, n-r

end

