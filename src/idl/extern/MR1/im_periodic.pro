;function im_periodic, Im_in, Smooth=Smooth
;+ 
; NAME: 
;     im_periodic
;
; PURPOSE: 
;     Decompose an image into a two images, one with periodic border conditions  and the second one is the 
;     smooth compoenent
;        InputImage = PeriodicPart  + SmoothPart
;     Return the periodic part.
;
; CALLING SEQUENCE: 
;    PeriodIma = im_periodic(ImagIn, Smooth=Smooth)
;
; INPUTS: 
;   ImagIn -- 2D IDL array: image to transform into a periodic image
;
; OPTIONAL OUTPUT PARAMETER: 
;   Smooth-- 2D IDL array: 
;
; OUTPUTS: 
;   PeriodIma -- 2D IDL array: Periodic image 
; 
; OPTIONAL OUTPUT PARAMETERS: 
;   none
;
; MODIFICATION HISTORY: 
;    21-Dec-2009 Arnaud Woiselle 
;-

function repmat,M,x,y

	A=M
	MM=M

	for i=1,x-1 do begin
		A=[A,M]
		MM=[MM,M]
	endfor

	for i=1,y-1 do begin
		A=[[A],[MM]]
	endfor

	return,A
end


function get_smooth_part,o
	
	s=size(o)
	nx = s(1)
	ny = s(2)
	
	o=double(o)
	m=double(o*0.)
	
	m(0,*) += o(0,*) - o(nx-1,*)
	m(nx-1,*) += o(nx-1,*) - o(0,*)
	m(*,0) += o(*,0) - o(*,ny-1)
	m(*,ny-1) += o(*,ny-1) - o(*,0)
	
	x=repmat(indgen(s(1),1),1,s(2)) - s(1)/2
	y=repmat(indgen(1,s(2)),s(1),1) - s(2)/2
	
	kf = 4 - 2*cos(2*!pi*x/nx) - 2*cos(2*!pi*y/ny)
	kf(nx/2,ny/2) = 1
	
	sf = dft(m)
	sf(nx/2,ny/2) = 0
	s = real_part(dfti(sf/kf))
	
	return,s
end

function im_periodic, InMap, Smooth=Smooth

Smooth  = get_smooth_part(InMap)
P = InMap - Smooth
return, P

end
