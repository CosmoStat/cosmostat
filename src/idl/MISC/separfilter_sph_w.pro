;-------------------------------------------------------------------------------
;+
; NAME:
;	separfilter_sph_w
;	
; PURPOSE:
;	Filters the data X using the optimized separating filter, in order to get estimates
;   of the component sources S as in X = AS + B
;
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	separfilter_sph_w,  X, Stats, Param, filter_type, S, MASK = mask
;
; INPUTS:
;	X : is an array of strings containing the file names for the m observed mixtures in Healpix format
;		the maps should be centered
;
;		Stats : structure grouping spectral covariances computed from the data.
;			this is an output of the statistics computation step.
;		NOT USED HERE BUT WE KEEP IT FORUNIFORMITY WITH OTHER separfilter PROCEDURES
;
;
;	Param : structure grouping the covariance matching model parameters
;			Param.mixmat : m*n matrix of mixing coefficients
;			Param.source : n*q matrix of source variance profiles
;			Param.noise  : m vector of white noise variance on each channel
;						   in the colored noise case, this is an m*q matrix of noise variance profiles
;			Param.type   : =1 in the white noise case
;						   =2 in the colored noise case 
;
;   filter_type :   string, either "wiener", "pinv", "pinvA", "passband" depending on desired 
;					technique for source signal reconstruction
;
; OPTIONAL INPUTS:
;
;   MASK = mask: filename, a boolean Healpix map containing 0. to indicate invalid pixels and 1. to indicate valid pixels 
;			
;			THE INPUT DATA X SHOULD HAVE BEEN MASKED PRIOR TO CALLING THIS PROCEDURE IF
;			NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO COMPUTING THE STATISTICS.
;			IT IS THEN NOT USEFUL TO HAVE mask AS AN INPUT HERE. HOWEVER WE KEEP IT TO HAVE A UNIFORMITY
;			BETWEEN ALL data2stats PROCEDURES.
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	S :  is an array of strings containing the file names for the n  reconstructed source maps in Healpix format
;
; DEPENDENCIES:
;		THe Healpix package for handling spherical maps and computing forward and reverse 
;		spherical harmonics transforms is used.
;
; RESTRICTIONS:
;   This procedure assumes the wavelet transform of the data maps computed while calling procedure data2stats_sph_w were saved to disk
;   with file names as in data2stats_sph_w.
;
; PROCEDURES USED:
;		makefilter
;
; EXAMPLE:
;	separfilter_sph_w,  X, Stats, Param, filter_type, S, MASK = mask
;
; MODIFICATION HISTORY:
;   Yassir Moudden , 2005
;-------------------------------------------------------------------------------   
pro separfilter_sph_w, X, Stats, Param, filter_type, S, MASK = mask

	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 5 $
	then begin 
        print, 'CALLING SEQUENCE: separfilter_2d_w, X, Stats, Param, filter_type, S, MASK = mask'
        goto, DONE
	end
	
	
	if keyword_set(mask) EQ 0 $
	then mask_on = 0 $
	else mask_on = 1
	;------------------------------------------------------------------------

	
	;------------------------------------------------------------------------
	;recovering a few model parameters
	A = Param.mixmat 
	Rp = Param.source 
	Rn = Param.noise  
	
	m = double( (size(A))(1) )
	n = double( (size(A))(2) )
	q = double( (size(Rp))(2) )

	if mask_on $
	then read_fits_map, mask, masque				;reading mask file for masking on exit
	
	;building the filter response for the q wavelet scales bands
	makefilter, Param, filter_type, F

	;------------------------------------------------------------------------
	;scalewise filtering of the data maps
	
	S = strarr(n)
	source_map = 'rec_sph_map_w_'

	name = X(0)
	wave_tmp_name = strcompress('wave_' + string(0) +'_'+ name, /remove_all )
	read_fits_map, wave_tmp_name, carte
	;we are using the wavelet coefficient maps computed while calling data2stats_sph_w so we do not need to compute again 
	; the wavelet transform of the data
	
	S_tmp = 0.* carte
	print, 'BEGIN scalewise filtering of the data maps'
	for num_source = 0, n-1 do begin		    ;looping over the components
		
		S_tmp(*) = 0.
		S(num_source) = strcompress( source_map + string(floor( num_source) ) + '.fits',  /remove_all)  
		
		for num_scale = 0, q-1 do begin			;looping over the scales
			
			for num_map = 0, m-1 do begin		;looping over the data maps
				name = X(num_map)
				wave_tmp_name = strcompress('wave_' + string(floor( num_scale ) ) +'_'+ name, /remove_all )
				read_fits_map, wave_tmp_name, carte
				S_tmp = S_tmp +  F(num_source, num_map, num_scale) * carte
			endfor
			
		endfor
		;since we are using the a trous algorithm, 
		;these last two loops perform the inverse transform
		
		if mask_on then begin
			S_tmp =  real_part( S_tmp ) * masque	;multiplication on exit by the initial mask
			write_fits_map, S(num_source), S_tmp, /nested
		endif else write_fits_map, S(num_source), real_part( S_tmp ) , /nested

	endfor 
	print, 'END scalewise filtering of the data maps'


DONE:

return
end



