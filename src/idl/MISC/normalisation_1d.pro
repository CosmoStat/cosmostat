;-------------------------------------------------------------------------------
;+
; NAME:
;	normalisation_1d
;	
; PURPOSE:
;	Computes the wavelet transform of a dirac in a 1D signal and from there computes
;   l2 normalisation coefficients.
;
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	normalisation_1d, example_map, nb_scales, l2norm_out, l1norm_out
;
; INPUTS:
;	example_map : 1D idl array, is used here to determine the size of the maps used in the current experiments
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
;		atwt1D : the 1D a trou wavelet transform
;
; EXAMPLE:
;	normalisation_1D, example_map, nb_scales, l2norm_out, l1norm_out
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   
pro normalisation_1d, example_map, nb_scales, l2norm_out, l1norm_out 

	;put all pixels to zero except one which is put equal to one and compute the wavelet transform
	map = 0. *double(example_map )
	T = double( (size(map))(1) )
	map( floor(T/2.) ) = 1.
	
	atwt1d, map, temp_out, Nscale = nb_scales
	
	l2norm_out = dblarr(nb_scales)
	l1norm_out = dblarr(nb_scales)

	for num_scale = 0., nb_scales-1 do begin
		l1norm_out(num_scale) = total( abs( temp_out( num_scale,*) ), /double )  
		l2norm_out(num_scale) = sqrt( total( (  temp_out( num_scale,*) )^2. , /double )  ) 
	endfor
	

return
end
