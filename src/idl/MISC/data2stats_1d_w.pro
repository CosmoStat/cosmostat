;-------------------------------------------------------------------------------
; NAME:
;	data2stats_1d_w.pro
;	
; PURPOSE:
;	Computes mean covariance matrices in a wavelet representation. This implementation
;   computes one covariance on each scale and includes the possibility of handling correctly possible missing patches
;   specified by an optional input Mask. Variations are obviously possible.
;
;   The wavelet transform used here is the undecimated isotropic 'a trous' wavelet transform.
;   Other algorithms could be used. And many extensions are therefore possible.
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	data2stats_1d_w, X, nb_scales, Stats, MASK = mask
;
; INPUTS:
;	X : m*T matrix containing the observed signals of length T in m channels
;		the signals should be centered
;
;   nb_scales : integer specifying the number of wavelet scales (including smoothed array) to be computed.
;				We do not check that the specified number of scales is valid.
;
; OPTIONAL INPUTS:
;   MASK : length T boolean array with MASK(i) = 1 indicating valid samples in X(*,i)
;								  MASK(i) = 0 indicating invalid/masked samples in X(*,i)
;			
;			The pixels on the edges should be masked to leave 
;			edge wavelet coefficients out of the estimation of the covariance matrices.
;
;			IF A MASK IS SPECIFIED, THE INPUT DATA X SHOULD BE MASKED (MULTIPLIED BY THE MASK) 
;			PRIOR TO COMPUTING THE STATISTICS.
;
; KEYWORD PARAMETERS:
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
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		atwt1d : idl implementation of the 'a trous' wavelet transform for 1d data, in file 'at1dwt.pro'
;
; EXAMPLE:
;	data2stats_1d_w, X, nb_scales, Stats, MASK = mask
;
; MODIFICATION HISTORY:
;
; COMMENTS:
;
;		Written: Yassir Moudden 2005.
;-------------------------------------------------------------------------------   

pro data2stats_1d_w, X, nb_scales, Stats, MASK = mask

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: data2stats_1d_w, X, nb_scales, Stats, MASK = mask'
        goto, DONE
	end
	
	if keyword_set(MASK) EQ 0 $
	then mask_on = 0 $
	else mask_on = 1

	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	;recovering a few parameter values
   	q = nb_scales					;number of wavelet scales and of covariance matrices to be computed
	
	T = double((size(X))(2))		;signals are length T
	m = double((size(X))(1))		;number of channels
		
	;------------------------------------------------------------------------
	;create anonymous structure in which to group the summary covariance statistics
	Stats = {covmat:dblarr(m,m,q), weight:dblarr(1,q), normalisation:dblarr(1,q)}
	
	
	
	;------------------------------------------------------------------------------
	;computing a trous wavelet transform of the data and of the mask if necessary
	;SAVE THIS TO DISK?
	;------------------------------------------------------------------------------
	XW = dcomplexarr(q, T, m )
	for num_channel = 0. , m-1. do begin
		 atwt1d, reform(  X(num_channel, *) ), temp_out, Nscale = q
		 XW(*, *,num_channel) =  temp_out
	end	 

	if mask_on then begin
		mrs_mask, mask, '1D', 'atrous', q, MW
	endif
	
	
	
	
	;------------------------------------------------------------------------------
	;computing l2 normalisation coefficients 
	;this step may be suppressed by hand if this computation has been previously performed
	;and saved to disk!
	;------------------------------------------------------------------------------
	

	if 1 then begin
		print, 'computing l2 normalisation coefficients '
		example_map = reform(  X(0, *) )
		normalisation_1D, example_map, nb_scales, normalisation, l1norm_out
		print, 'normalisation = ', normalisation
		fits_write, 'normalisation_1D.fits', normalisation
	endif		

	
	;for instance
	if 0 then begin
		print, 'reading l2 normalisation coefficients '
		fits_read,'normalisation_1D.fits',normalisation
		print, 'normalisation = ', normalisation
	endif	

	;normalisation = [ 0.0226091,   0.00892032,   0.00556086,   0.00381947,   0.00268160,   0.00189282,   0.00133783,  0.000945888,  0.000717896, 0.000226322]
	;Stats.normalisation  = normalisation(0:q-1)

	Stats.normalisation  = normalisation

	
	;------------------------------------------------------------------------------			
	;computing covariance matrices
	;------------------------------------------------------------------------------
	for num_scale = 0., q-1. do begin
	
		Xq = transpose( XW( num_scale, *,*) )

		if mask_on then begin
			maskq = reform( MW(num_scale, *) )
			maskq_okay = where(   (maskq eq 1) )		
			weight_q = (size(maskq_okay))(1)
			if weight_q eq 0 $
			then begin
				print, 'Error : Some specified wavelet bands are empty.'
				goto, DONE
			endif
			Xq_okay = Xq( *, maskq_okay)
		endif else begin
			weight_q = T
			Xq_okay = Xq
		endelse

		Stats.weight(num_scale) = (1.0/(2.0^num_scale))*weight_q
		;this heuristic is because the a trous transform is redundant
		;the number of statistically significant pixels is divided by
		;four from one scale to the next scale (as in the orthogonal wavelet transform)
	
	
		for num_channel = 0, m-1 do Xq_okay(num_channel, *) = Xq_okay(num_channel, *) - mean( Xq_okay(num_channel, *) )
		Stats.covmat(*, *, num_scale) = real_part(  Xq_okay  # transpose(Xq_okay) ) / ( normalisation(num_scale)^2 * weight_q )
		;real part is not necessary since everything is real, but still...
		
	endfor

	Stats.weight(q-1) = Stats.weight(q-1)*2.0/2.0
	Stats.weight = Stats.weight/total(Stats.weight, /double)
	;the number of statistically significant coefficients is three times less in the smoothed image than in the
	;detail maps (heurisitcally speaking)


DONE:

return

end
		
		
		
	
		
		
		

