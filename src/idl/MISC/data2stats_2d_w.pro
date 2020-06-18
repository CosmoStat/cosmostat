;-------------------------------------------------------------------------------
; NAME:
;	data2stats_2d_w.pro
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
;	data2stats_2d_w, X, nb_scales, Stats, MASK = mask
;
; INPUTS:
;	X : tx*ty*m matrix containing the observed data maps of size tx*ty in m channels
;		the maps should be centered
;
;   nb_scales : integer specifying the number of wavelet scales (including smoothed array) to be computed.
;				We do not check that the specified number of scales is valid.
;
; OPTIONAL INPUTS:
;   MASK : tx*ty boolean array with MASK(i,j) = 1 indicating valid samples in X(i,j,*)
;								  MASK(i,j) = 0 indicating invalid/masked samples in X(i,j,*)
;			
;			The pixels on the edges should be masked to leave 
;			edge wavelet coefficients out of the estimation of the covariance matrices.
;
;			THE INPUT DATA X SHOULD BE MASKED IF NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO COMPUTING THE STATISTICS.
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
;									 IS THIS REALLY USEFUL?
;									 USED TO PONDER THE COVARIANCE MATRICES
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		atwt2d : idl implementation of the 'a trous' wavelet transform for 2d data, in file 'at2dwt.pro'
;
; EXAMPLE:
;	data2stats_2d_w, X, nb_scales, Stats, MASK = mask
;
; MODIFICATION HISTORY:
;		Written: Yassir Moudden 2005.
;
; COMMENTS:
;
;
;-------------------------------------------------------------------------------   

pro data2stats_2d_w, X, nb_scales, Stats, MASK = mask

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: data2stats_2d_w, X, nb_scales, Stats, MASK = mask'
        goto, DONE
	end
	
	if keyword_set(mask) eq 0 $
	then mask_on = 0 $
	else mask_on = 1

	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	;recovering a few parameter values
   	q = nb_scales					;number of wavelet scales and of covariance matrices to be computed
	tx = double((size(X))(1))		;data maps have size tx*ty
	ty = double((size(X))(2))
	m = double((size(X))(3))		;number of channels
		
	;------------------------------------------------------------------------
	;create anonymous structure in which to group the summary covariance statistics
	Stats = {covmat:dblarr(m,m,q), weight:dblarr(1,q), normalisation:dblarr(1,q)}
  
	;------------------------------------------------------------------------------  
	;computing a trous wavelet transform of the data and of the mask if necessary   
	;SAVE THIS TO DISK?
	;------------------------------------------------------------------------------

	XW = dcomplexarr(tx, ty, m, q)
	for num_channel = 0. , m-1. do begin
		print, 'computing 2d wavelet transform of input map number ', num_channel
		 atwt2d, X(*,*,num_channel), temp_out, Nscale = q
		 XW(*,*,num_channel, *) = temp_out
	end	 

	if mask_on then begin
		mrs_mask, mask, '2D', 'atrous', q, MW
	endif
	
	
	;------------------------------------------------------------------------------
	;computing l2 normalisation coefficients 
	;this step may be suppressed by hand if this computation has been previously performed
	;and saved to disk!
	;------------------------------------------------------------------------------
	
	
	;normalisation = [0.890796,   0.200664,  0.0855075,   0.0412174,  0.0204249,  0.0101896 , 0.00509159,0.00254474, 0.00132169, 0.0010303]
	;Stats.normalisation  = normalisation(0:q-1)


	if 1 then begin
		print, 'computing l2 normalisation coefficients '
		example_map = X(*,*, 0)
		normalisation_2D, example_map, nb_scales, normalisation, l1norm_out
		print, 'normalisation = ', normalisation
		fits_write, 'normalisation_2D.fits', normalisation
	endif		

	
	;for instance
	if 0 then begin
		print, 'reading l2 normalisation coefficients '
		fits_read,'normalisation_2D.fits',normalisation
		print, 'normalisation = ', normalisation
	endif	
	
	Stats.normalisation  = normalisation

		
	;------------------------------------------------------------------------------
	;computing covariance matrices
	;------------------------------------------------------------------------------
	for num_scale = 0., q-1. do begin
	
		Xq = reform( XW(*, *, *, num_scale), tx*ty, m )

		if mask_on then begin
			maskq = reform( MW(*,*, num_scale) ,  tx*ty, 1 )
			maskq_okay = where(   (maskq eq 1) )	
			weight_q = (size(maskq_okay))(1)
			if weight_q eq 0 $
			then begin
				print, 'Error : Some specified wavelet bands are empty.'
				goto, DONE
			endif
			Xq_okay = Xq( maskq_okay, *)
		endif else begin
			weight_q = tx*ty
			Xq_okay = Xq
		endelse

		Stats.weight(num_scale) = (1.0/(4.0^num_scale))*weight_q
		;this heuristic is because the a trou transform is redundant
		;the number of statistically significant pixels is divided by
		;four from one scale to the next scale (as in the orthogonal wavelet transform)
		
		for num_channel = 0., m-1. do Xq_okay(*,num_channel) = Xq_okay(*,num_channel) - mean( Xq_okay(*,num_channel) )
		Stats.covmat(*, *, num_scale) = real_part( transpose(Xq_okay) # Xq_okay ) / ( normalisation(num_scale)^2 * weight_q )
		;real part is not necessary since everything is real, but still...
		
	endfor

	Stats.weight(q-1) = Stats.weight(q-1)*4.0/3.0
	Stats.weight = Stats.weight/total(Stats.weight, /double)
	;the number of statistically significant coefficients is three times less in the smoothed image than in the
	;detail maps (heurisitcally speaking)


DONE:

return

end
		
		
		
	
		
		
		

