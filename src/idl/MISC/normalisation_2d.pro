;-------------------------------------------------------------------------------
;+
; NAME:
;	normalisation_2d
;	
; PURPOSE:
;	Computes the wavelet transform of a dirac in a 2D map and from there computes
;   l2 normalisation coefficients.
;
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	normalisation_2d, example_map, nb_scales, l2norm_out, l1norm_out
;
; INPUTS:
;	example_map : 2D idl array, is used here to determine the size of the maps used in the current experiments
;						for which we need the normalizing coefficients
;
;   nb_scales:			number of scales for the wavelet transform of the dirac
;
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   l2norm_out : array giving the normalizing coeffs in the l2 norm sense
;   l1norm_out : array giving the normalizing coeffs in the l1 norm sense
;
; DEPENDENCIES:
;	
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		atwt2D : the 2d a trou wavelet transform
;
; EXAMPLE:
;	normalisation_2D, example_map, nb_scales, l2norm_out, l1norm_out
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   
pro normalisation_2d, example_map, nb_scales, l2norm_out, l1norm_out 

	;put all pixels to zero except one which is put equal to one and compute the wavelet transform
	map = 0. *double(example_map )
	tx = double( (size(map))(1) )
	ty = double( (size(map))(2) )
	map(floor(tx/2.), floor(ty/2.) ) = 1.
	
	atwt2d, map, temp_out, Nscale = nb_scales
	
	l2norm_out = dblarr(nb_scales)
	l1norm_out = dblarr(nb_scales)

	for num_scale= 0., nb_scales-1 do begin
		l1norm_out(num_scale) = total( abs( temp_out( *, *, num_scale) ), /double)  
		l2norm_out(num_scale) = sqrt( total( (  temp_out( *, *, num_scale) )^2. , /double)  ) 
	endfor
	

return
end
