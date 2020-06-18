;-------------------------------------------------------------------------------
;+
; NAME:
;	data2stats_sph_w
;	
; PURPOSE:
;	Computes mean covariance matrices in a wavelet representation. This implementation
;   computes one covariance on each scale and includes the possibility of handling correctly possible missing patches
;   specified by an optional input Mask. Variations are obviously possible.
;
;   The wavelet transform used here is the undecimated isotropic 'a trous' wavelet transform on the sphere.
;   Other algorithms could be used. And many extensions are therefore possible.
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	data2stats_sph_w, X, nb_scales, Stats, MASK = mask, NLMAX = nlmax, l2_norm = l2_norm, old_mask = old_mask, old_data = old_data 
;
; INPUTS:
;	X : is an array of strings containing the file names for the m observed mixtures in Healpix format
;		the maps should be centered.
;
;   nb_scales : integer specifying the number of wavelet scales (including smoothed array) to be computed.
;				We do not check that the specified number of scales is valid.
;
; OPTIONAL INPUTS:
;   MASK =  mask: a  filename, the corresponding file is a boolean Healpix map containing 0. 
;					to indicate invalid pixels and 1. to indicate valid pixels 
;
;			THE INPUT DATA X SHOULD HAVE BEEN MASKED PRIOR TO CALLING THIS PROCEDURE IF
;			NECESSARY (MULTIPLIED BY THE MASK PRIOR TO COMPUTING THE STATISTICS ).
;
;
;   NLMAX = nlmax : integer, maximum value of l for which the spherical harmonics are to be computed.
;			default is ceil (3.5*nside) where nside is the number of pixels on each side of the Healpix rhombi
;
;   l2_norm = a filename, the corresponding file is a fits file containing an array of l2 normalization coefficients 
;			 for the nromalisation of the wavelet coefficients. if the file is specified, then these coeffiecients are used
;			  otherwise they are computed from the wavelet transform of a dirac on the sphere.
;
; KEYWORD PARAMETERS:
;	old_mask = if not set, a mask on each scale is determined from the wavelet transform of the specified mask
;			   if set, the masks on each scale are read from the current directory, from files named 
;			    strcompress('wave_' + 'i' +'_'+ mask, /remove_all ) where mask is the filename of the specified mask.
;
;	old_data = if not set, the wavelet transforms of the specified data maps are computed and saved in the current directory
;			   if set, the wavelet coefficients of the data maps are read from the current directory, from files named 
;			   strcompress('wave_' + 'i' +'_'+ name, /remove_all ) where name stands for the filenames of the specified data files.
;
; OUTPUTS:
;	Stats : structure grouping spectral covariances computed from the data.
;			this is an output of the statistics computation step.
;			required fields for use in other procedures are:
;
;				Stats.covmat :  m*m*q 3D array of q m*m covariance matrices
;								q is the number of bands or scales, etc.
;								m is the number of data channels
;
;				Stats.weight :  statistical significance weights of the covariance matrices
;
;			other fields for graphical purposes, source estimation, etc :
;		
;				Stats.normalisation : 1*q array containing the l2 norm of the wavelet or smoothing filters
;										on each of the q scales
;									 USED TO PONDER THE COVARIANCE MATRICES
;	
;
; DEPENDENCIES:
;		The Healpix package for handling spherical maps and computing forward and reverse 
;		spherical harmonics transforms is used.
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;   mrs_wttrans
;
; EXAMPLE:
;	data2stats_sph_w,  X, nb_scales, Stats, MASK = mask, NLMAX = nlmax, l2_norm = l2_norm, old_mask = old_mask, old_data = old_data
;
; MODIFICATION HISTORY:
;		Written: Yassir Moudden, 2005
;-------------------------------------------------------------------------------   


pro data2stats_sph_w, X, nb_scales, Stats, MASK = mask, NLMAX = nlmax, l2_norm = l2_norm, old_mask = old_mask, old_data = old_data

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: data2stats_2d_w, X, nb_scales, Stats, MASK = mask, NLMAX = nlmax, l2_norm = l2_norm, old_mask = old_mask, old_data = old_data'
        goto, DONE
	end
	
	if keyword_set(MASK) eq 0 $
	then mask_on = 0 $
	else mask_on = 1
	
	read_fits_map, X(0), tmp
	nb_pixel = (size(tmp))(1)
	if not keyword_set(NLMAX) then begin
		nlmax =  ceil( 3.5 * npix2nside(nb_pixel) ) 
	endif else nlmax = ceil(nlmax)					;nlmax is expected to be an 
													;integer so this just a precaution
													;which should not affect the in/out variable
													; nlmax in any undesired way
													
	;------------------------------------------------------------------------
	
	
	;------------------------------------------------------------------------
	;recovering a few parameter values
   	q = nb_scales					;number of wavelet scales and of covariance matrices to be computed
	m = double((size(X))(1))		;number of channels

	;create anonymous structure in which to group the summary covariance statistics
	Stats = {covmat:dblarr(m,m,q), weight:dblarr(1,q), normalisation:dblarr(1,q)}



	;------------------------------------------------------------------------------
	;computing l2 normalisation coefficients 
	;this step may be suppressed if this computation has been previously performed
	;and saved to disk! use the l2norm keyword to specify a .fits file containing the right 
	;normalisation values
	;------------------------------------------------------------------------------

	if not keyword_set(l2_norm) then begin
	
		print, 'computing l2 normalisation coefficients'
		example_file_name = X(0)
		normalisation_sph, example_file_name, nb_scales, nlmax, normalisation, l1norm_out
		print, 'normalisation = ', normalisation
		fits_write, 'normalisation_sphere.fits', normalisation
		print, 'normalisation coefficients were saved to normalisation_sphere.fits'
	
	endif else begin
	
		print, 'reading l2 normalisation coefficients from the specified file: ', l2_norm		
		fits_read, l2_norm, normalisation
		print, 'normalisation = ', normalisation
		;verify that the specified array is the right size
		size_norm = size(normalisation)
		if (size_norm(1) LT q) then begin
				print, 'The specified array of normalisatoin coefficients is invalid.'
				goto, DONE
		endif
	endelse	

	;------------------------------------------------------------------------------
	;compute the wavelet transform of the mask in order to get the right
	;mask on each scale.
	;if the same mask is used several times, one might prefer to compute
	; the transform once and write it to disk in a separate file for each scale
	;then recall those when necessary. This is done using the keyword new_mask.
	;------------------------------------------------------------------------------


	if mask_on then begin
		if not keyword_set(old_mask) then begin

			print, 'BEGIN compute the wavelet transform of the mask' 
		
			read_fits_map, mask, masque
			;ones in mask indicate valid pixels, zeros indicate invalid pixels
		
			mrs_mask, masque, 'Sphere', 'atrous', q, wave_masque,  nlmax = nlmax
			;wave_masque is an array which gives the masques to be used on the different scales								
			;ie 1 indicates a valid pixel and 0 indicates an invalid pixel
		
			;it may be useful to save the results to disk
			for num_scale = 0., q-1. do begin
				write_fits_map, strcompress('wave_' + string(floor( num_scale) ) +'_'+ mask, /remove_all ), wave_masque(*, num_scale), /nested
			endfor	

			print, 'END compute the wavelet transform of the mask'
		
		
		endif else begin
			print, 'Will be using precomputed masks on each scale.'
			print, 'Mask on scale i will be read from file ', strcompress('wave_' + 'i' +'_'+ mask, /remove_all )
		endelse	  
	endif

	
	
	;------------------------------------------------------------------------------
	;compute the wavelet transform of the different observed maps
	;Clearly, the maps should be masked using the initial mask
	; (if a mask is used) before the wavelet transform
	;------------------------------------------------------------------------------

	if not keyword_set(old_data) then begin

		print, 'BEGIN computing the wavelet transform of the data maps' 
		for num_map = 0., m-1. do begin
		
			;recovering data map filename
			name = X(num_map)
			read_fits_map, name, tmp
			mrs_wttrans, tmp, wave_tmp, NbrScale = nb_scales, lmax = nlmax			;, /DifInSH
		
			;writing to disk the different wavelet maps for each data map
			for num_scale = 0., nb_scales - 1. do begin
				wave_tmp_name = strcompress('wave_' + string(floor( num_scale) ) +'_'+ name, /remove_all )
				write_fits_map, wave_tmp_name, wave_tmp.coef(*, num_scale), /nested
			endfor
	
		endfor
		print, 'END computing the wavelet transform of the data maps' 
		
	endif else begin
		print, 'Will be using precomputed wavelet transforms of the data maps.'
		print, 'Scale i of data map name will be read from', strcompress('wave_' + 'i' +'_'+ 'name', /remove_all )
	endelse	  

	

	;------------------------------------------------------------------------------
	;computing the covariance matrix at each scale
	
	print, 'BEGIN computing the covariance matrices at each scale' 

	for num_scale= 0., nb_scales-1. do begin     


		if mask_on then begin 
			mask_local = wave_masque(*, num_scale)
			iM_inter = where(  ( mask_local eq 1) , count )
		endif else begin 
			iM_inter = findgen( nb_pixel)
			count = 1.
		endelse
			
		if (count ne 0)  then begin 
			lsize = double(  (size(iM_inter))(1) )
			carte_ok = dblarr(m,lsize)
	
			for num_map = 0., m-1. do begin
			print, 'num_scale = ', num_scale, '  num_map = ', num_map
				name = X(num_map)
				wave_tmp_name = strcompress('wave_' + string(floor( num_scale)) +'_'+ name, /remove_all )
				read_fits_map, wave_tmp_name, carte
				carte_ok(num_map, *) = carte(iM_inter)
				carte_ok(num_map, *) = carte_ok(num_map, *) - mean(carte_ok(num_map, *))
			endfor
		
			Stats.covmat(*,*,num_scale) = real_part( carte_ok # transpose(  carte_ok ) ) / ( normalisation(num_scale)^2. * lsize ) 
			
			Stats.weight(num_scale) = (1.0/(4.0^num_scale))*lsize
			;this is not really right on the sphere, I should be counting the number of spherical harmonic modes in each wavelet band
			;lsize is a correction for sky coverage
			;need to evaluate more precisely the number of harmonic modes in each band
			;this is not very difficult since the wavelet used has compact support in fourier space.
			
			
		endif else begin
			print, 'this band was empty : ', num_scale
			goto, DONE
		endelse
	
	endfor		
	
	print, 'END computing the covariance matrices at each scale' 

	Stats.normalisation = normalisation
	Stats.weight(q-1) = Stats.weight(q-1)/3.0
	Stats.weight = Stats.weight/total(Stats.weight, /double)
	
	;the number of statistically significant coefficients is three times less in the smoothed image than in the
	;detail maps
	;this is not really right on the sphere, I should be counting the number of spherical harmonic modes in each wavelet band
	;minus the statistical significance of the missing pixels


DONE : 

return

end


