;-------------------------------------------------------------------------------
;
; NAME:
;	separfilter_1d_f.pro
;	
; PURPOSE:
;	Filters the data X using the optimized separating filter, in order to get estimates
;   of the component sources S as in X = AS + B
;
; EXPLANATION:
;   Further versions should allow to use differrent filters in different bands, especially in the case
;   where the set of frequency bins is not a partition of the frequency intervalle.
;
; CALLING SEQUENCE:
;	separfilter_1d_f,  X, Stats, Param, filter_type, S, MASK = mask
;
;
; INPUTS:
;	X : m*T matrix containing the observed data in m channels during "time" T
;
;	Stats : structure grouping spectral covariances computed from the data.
;			this is an output of the statistics computation step.
;			required fields used here are:
;				Stats.fbins :  2*q  array specifying frequency bins
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
;   MASK : 1*T boolean array with MASK(0, t) = 1 indicating valid samples in X(*, t)
;								  MASK(0, t) = 0 indicating invalid/masked samples in X(*, t)
;			
;			THE INPUT DATA X SHOULD BE MASKED IF NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO CALLING THE PROCEDURE.
;			ESTIMATED SOURCE SIGNALS ARE MULTIPLIED BY THE MASK BEFORE OUTPUT 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	S :  n*T matrix containing the estimated signals from n sources during "time" T
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		makefilter
;
; EXAMPLE:
;	separfilter_1d_f,  X, Stats, Param, filter_type, S, MASK = mask
;
; MODIFICATION HISTORY:
;   Yassir Moudden , 2005
;-------------------------------------------------------------------------------   
pro separfilter_1d_f, X, Stats, Param, filter_type, S, MASK = mask

	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 5 $
	then begin 
        print, 'CALLING SEQUENCE: separfilter_1d_f, X, Stats, Param, filter_type, S, MASK = mask'
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
	T = double(  (size(X))(2) )

	;building the filter response for the Q frequency bands
	makefilter, Param, filter_type, F
	
	
	;computing FFT of the data
	XF = dcomplexarr(m, T)																	;- mean( X(num_channel, * ) )
	for num_channel=0,  m-1 do XF(num_channel, * ) = T * fft(  X(num_channel, * ), -1, /double )		;is centering useful here? 
	;useless recomputation of XF since it was already computed once in data2stats_1d_f
	;we should try to keep it somewhere

	;computing a matrix containing the values of |f| in the 1d fft arrays i.e. a frequency index
	fx = dindgen(T , 1)/T
	fx = abs( fx - double(fx gt 0.5 ) ) 

	; bandwise filtering of the data maps
	tol = 0.						;possible extra overlapping between bins
	SF = 0.*dcomplexarr(n, T)		;fourier transform of the reconstructed maps
	poids = 0.* dblarr(n, T)		;weighting matrix so that reconstructed fourier modes at bin edges are averaged contributions from both sides.

	for num_q = 0, q-1 do begin
		fstart = double( 	Stats.fbins(0,num_q) )
		fstop  = double(    Stats.fbins(1,num_q) )			
		fx_okay   = WHERE(  (  fx GE (fstart-tol)) AND ( fx LE (fstop+tol)  )  )
		SF(*,fx_okay) = SF(*, fx_okay) + F(*,*,num_q) # XF(*,fx_okay)				
		poids(*,fx_okay) = poids(*,fx_okay) +1.
	endfor 
	
	
	;STILL TO DO
	;what if some fourier modes do not fall in the specified modes?
	;do something like:
	;fx_not_in_the_bands = 
	;SF(*,fx_not_in_the_bands) = SF(*, fx_not_in_the_bands) + invert ( transpose(A)# A, /double) # transpose(A) # XF(*,fx_not_in_the_bands) 
	;poids(*,fx_not_in_the_bands) = poids(*,fx_not_in_the_bands) +1.
	
	
	poids = poids + (poids LT 1.)
	SF = SF/poids
  
	;inverse fourier transform to get the signals in the time representation
	S = 0. * dblarr(n, T)
	if mask_on $
	then for num_source= 0 , n-1 do S(num_source, *) = mask * real_part(  fft( SF(num_source, * ) , 1 ,/double) )$   
	else for num_source= 0 , n-1 do S(num_source, *) =        real_part(  fft( SF(num_source, * ) , 1 ,/double) )
	;it should be possible to get rid of the loops.	


DONE:

return
end
  