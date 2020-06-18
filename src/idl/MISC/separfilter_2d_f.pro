;-------------------------------------------------------------------------------
;
; NAME:
;	separfilter_2d_f.pro
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
;	separfilter_2d_f,  X, Stats, Param, filter_type, S, MASK = mask
;
;
; INPUTS:
;	X : tx*ty*m matrix containing the observed data maps of size tx*ty in m channels
;
;	Stats : structure grouping spectral covariances computed from the data.
;			this is an output of the statistics computation step.
;			required fields used here are:
;				Stats.frings :  2*q  array specifying frequency bins
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
;   MASK : tx*ty boolean array with MASK(i,j) = 1 indicating valid samples in X(*, i,j)
;								  MASK(i,j) = 0 indicating invalid/masked samples in X(*, i,j)
;			
;			THE INPUT DATA X SHOULD BE MASKED IF NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO CALLING THE PROCEDURE.
;			ESTIMATED SOURCE SIGNALS ARE MULTIPLIED BY THE MASK BEFORE OUTPUT 
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	S :  tx*ty*n array containing the estimated signals from n sources over maps of size tx*ty
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		makefilter
;
; EXAMPLE:
;	separfilter_2d_f,  X, Stats, Param, filter_type, S, MASK = mask
;
; MODIFICATION HISTORY:
;   Yassir Moudden and Jean-Francois Cardoso, 2005
;-------------------------------------------------------------------------------   
pro separfilter_2d_f, X, Stats, Param, filter_type, S, MASK = mask

	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 5 $
	then begin 
        print, 'CALLING SEQUENCE: separfilter_2d_f, X, Stats, Param, filter_type, S, MASK = mask'
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
	tx = double((size(X))(1))		;data maps have size tx*ty
	ty = double((size(X))(2))

	;building the filter response for the q frequency bands
	makefilter, Param, filter_type, F
	
	
	;computing 2d FFT of the data, assuming equal sampling frequencies along x and y axes
	XF = dcomplexarr(tx, ty, m)
	for num_channel=0 , m-1 do XF(* , * , num_channel ) = (tx*ty)  * fft( X( * , * , num_channel ) , -1,/double )
	XF = transpose( reform(XF, tx*ty, m) )
	;is centering useful here? 
	; - mean( X(*,*, num_channel))
	;useless recomputation of XF since it was already computed once in data2stats_1d_f
	;we should try to keep it somewhere

	;computing a matrix containing the values of |f| in the 2d fft arrays
	fx = dindgen(tx , 1)/tx
	fx = fx - double(fx gt 0.5 ) 
	fy = dindgen(1 , ty)/ty 
	fy = fy - double(fy gt 0.5)
	fmap = sqrt(   (fx * fx) # (dblarr(1, tx )+1.)  +    (dblarr(ty, 1) +1. ) # (fy * fy)   )  
	
	
	; bandwise filtering of the data maps
	tol = 0.						;possible extra overlapping between bins
	SF = 0.*dcomplexarr(n, tx*ty)		;fourier transform of the reconstructed maps
	poids = 0.* dblarr(n, tx*ty)		;weighting matrix so that reconstructed fourier modes at bin edges are averaged contributions from both sides.

	for num_q = 0, q-1 do begin
		fstart = double( 	Stats.frings(0,num_q) )
		fstop  = double(    Stats.frings(1,num_q) )			
		f_okay   = WHERE(  (  fmap GE (fstart-tol)) AND ( fmap LE (fstop+tol)  )  )
		SF(*,f_okay) = SF(*, f_okay) + F(*,*,num_q) # XF(*,f_okay)				
		poids(*,f_okay) = poids(*,f_okay) +1.
	endfor 
	
	
	;STILL TO DO
	;what if some fourier modes do not fall in the specified modes?
	;do something like:
	;fx_not_in_the_bands = 
	;SF(*,fx_not_in_the_bands) = SF(*, fx_not_in_the_bands) + invert ( transpose(A)# A , /double) # transpose(A) # XF(*,fx_not_in_the_bands) 
	;poids(*,fx_not_in_the_bands) = poids(*,fx_not_in_the_bands) +1.
	
	
	poids = poids + (poids LT 1.)
	SF = SF/poids
  
	;inverse fourier transform to get the signals in the time representation
	S = 0. * dblarr(tx, ty, n)
	if mask_on $
	then for num_source= 0 , n-1 do S(* , * , num_source ) = mask * real_part(  fft( reform( SF(num_source, * ),tx, ty) , 1, /double ) )$   
	else for num_source= 0 , n-1 do S(* , * , num_source ) =        real_part(  fft( reform( SF(num_source, * ),tx, ty)  , 1, /double) )
	;it should be possible to get rid of the loops.	

DONE:

return
end

