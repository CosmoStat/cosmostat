;+
; NAME:
;        mrsp_resize
;
; PURPOSE:
;   Resize polarized map in healpix representation 
;
; CALLING:
;     resize_map = mrsp_resize( map, nside=nside, ViaAlm=ViaAlm, teb=teb )
;
; INPUTS:
;     Imag -- 3D IDL array of healpix polarized map in TQU scheme: Input image to be transformed 
;    
; OUTPUTS:
;     resize_map -- 3D IDL array of healpix polarized map.
;
; KEYWORDS:
;      nside     : the nside of the healpix output map
;	   ViaAlm	 : use ALM transform
;	   teb		 : The input and output imag are in TEB scheme
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;       cl2map (glesp software)
;       reorder (healpix software)
;
; EXAMPLE:
;       resize an healpix map
;
;               map2 = mrsp_resize(map,nside = 256)
;
;         
; HISTORY:
;	Written:  Olivier Fourt, 2009
;	February, 2009 File creation
;--------------------------------------------------------------------------------------------------------


function mrsp_resize, Imag, nside=nside, ViaAlm = ViaAlm, teb = teb

if keyword_set(teb) then begin

	map_1 = Imag[*,0]
	map_2 = Imag[*,1]
	map_3 = Imag[*,2]
	
	map_1 = mrs_resize( map_1, nside=nside, ViaAlm = ViaAlm )
	map_2 = mrs_resize( map_2, nside=nside, ViaAlm = ViaAlm )
	map_3 = mrs_resize( map_3, nside=nside, ViaAlm = ViaAlm )
	
	NImag = [ [map_1], [map_2], [map_3] ]
	
end else begin

	if not keyword_set(ViaAlm) then ViaAlm = 0
	ImagNside = gnside(Imag[*,0])
	if not keyword_set(nside) then nside=ImagNside
	if ImagNside GE nside then ViaAlm=0
	
	if ViaAlm EQ 0 then begin
		
		if keyword_set(nside) then begin
        	ud_grade, Imag, NImag, nside_out=nside, order_in='nested',order_out='nested'
        end else begin
        	NImag = Imag
        end
		
	end else begin
	
		NImag = 0.*dblarr(nside^2.*12.,3)
		
		mrsp_almtrans, NImag, A, /tab
    	mrsp_almtrans, Imag, alm, /tab
    	A.alm[*]=0.
    
    	A.alm[0:alm.lmax, 0:alm.lmax, *, *] = alm.alm
    	mrsp_almrec, A, NImag
		
	end

end
	
return, NImag
        
end
   
    
