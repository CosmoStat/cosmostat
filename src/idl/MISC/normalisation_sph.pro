;-------------------------------------------------------------------------------
;+
; NAME:
;	normalisation_sph
;	
; PURPOSE:
;	Computes the wavelet transform of a dirac on the sphere and frome there computes the 
;   l2 normalisation coefficients.
;
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	normalisation_sph, example_file_name, nb_scales, nlmax, l2norm_out, l1norm_out
;
; INPUTS:
;	example_file_name : healpix file name, determines the size of the maps used in the current experiments
;						for which we need the normalizing coefficients
;
;   nb_scales:			number of scales for the wavelet transform of the dirac
;   nlmax :				integer : maximum multipole to be computed for the wavelet transform
;
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   l2norm_out : array giving the normalizaing coeffs in the l2 norm sense
;   l1norm_out : array givieng the normalizing coeffs in the l1 norm sense
;
; DEPENDENCIES:
;		THe Healpix package for handling spherical maps and computing forward and reverse 
;		spherical harmonics transforms is used.
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		
;
; EXAMPLE:
;	normalisation_sph, example_file_name, nb_scales, nlmax, l2norm_out, l1norm_out
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   
pro normalisation_sph, example_file_name, nb_scales, nlmax, l2norm_out, l1norm_out

	read_fits_map, example_file_name, map
		
	nlmax = ceil(nlmax)		;just in case, should not create any undesired effect ... 
			
	;put all pixels to zero except one which is put equal to one and compute the wavelet transform
	map = 0. *double(map )
	npix = double( (size(map))(1))
	nside = npix2nside(npix)
	map( round( npix/2. + 2.*nside) ) = 1.
	mrs_wttrans, map, wave_map, NbrScale=nb_scales, lmax = nlmax		;, /DifInSH
	
	l2norm_out = dblarr(nb_scales)
	l1norm_out = dblarr(nb_scales)

	for num_scale= 0., nb_scales-1 do begin
		l1norm_out(num_scale) = total( abs( wave_map.coef( *, num_scale) ), /double )  
		l2norm_out(num_scale) = sqrt( total( (  wave_map.coef( *, num_scale) )^2.  , /double )  ) 
		;mollview, reform( map_out.map(num_scale, *) )/normalisation_out(num_scale) , /online

	endfor
	

return
end
