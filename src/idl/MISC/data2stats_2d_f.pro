;-------------------------------------------------------------------------------
;+
; NAME:
;	data2stats_2d_f.pro
;	
; PURPOSE:
;	computes mean covariance matrices in the fourier representation over 
;   specified spectral rings 
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	data2stats_2d_f, X, freq_rings, Stats, MASK = mask
;
; INPUTS:
;	X : tx*ty*m matrix containing the observed data maps of size tx*ty in m channels
;   freq_rings : 2*q matrix containing the starting and ending reduced frequencies of each of 
;				q frequency rings. Each specified reduced frequency is between [0, 0.5[.
;
; OPTIONAL INPUTS:
;   MASK : tx*ty boolean array with MASK(i,j) = 1 indicating valid samples in X(*, i,j)
;								  MASK(i,j) = 0 indicating invalid/masked samples in X(*, i,j)
;			
;			THE INPUT DATA X SHOULD BE MASKED IF NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO COMPUTING THE STATISTICS.
;			IT IS THEN NOT USEFUL TO HAVE mask AS AN INPUT HERE. HOWEVER WE KEEP IT TO HAVE A UNIFORMITY
;			BETWEEN ALL data2stats PROCEDURES.
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
;				Stats.fmean : 1*q array containing the mean frequency vector norm over ring q
;				Stats.fwidth: 1*q array containing the width in frequency vector norm of ring q
;				Stats.frings:  2*q array equal to freq_rings
;	
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	data2stats_2d_f, X, freq_rings, Stats, MASK = mask
;
; MODIFICATION HISTORY:
;		Written: Yassir Moudden & Jean-Francois Cardoso 2005.
;
; COMMENTS:
;   - future versions wight want to include the possibility for different sampling frequencies
;   along the x and y axis of the 2D data maps X
;
;-------------------------------------------------------------------------------   

pro data2stats_2d_f, X, freq_rings, Stats, MASK = mask

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: data2stats_2d_f, X, freq_rings, Stats, MASK = mask'
        goto, DONE
	end
	
	if keyword_set(mask) eq 0 $					;rather useless since the fourier transform is not local!
	then mask_on = 0 $
	else mask_on = 1

	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	;recovering a few parameter values
   	q = (size(freq_rings))(2)		;number of spectral bands/of covariance matrices to be computed
	tx = double((size(X))(1))		;data maps have size tx*ty
	ty = double((size(X))(2))
	m = double((size(X))(3))		;number of channels
		
	;------------------------------------------------------------------------
	;create anonymous structure in which to group the summary covariance statistics
	Stats = {covmat:dblarr(m,m,q), weight:dblarr(1,q), fmean:dblarr(1,q), fwidth:dblarr(1,q), frings:dblarr(2, q)}
  
	;computing 2d FFT of the data, assuming equal sampling frequencies along x and y axes
	XF = dcomplexarr(tx, ty, m)
	for num_channel=0. , m-1. do XF(* , * , num_channel ) = (tx*ty)  * fft( X( * , * , num_channel ) - mean( X(*,*,num_channel)), -1, /double )
	XF = transpose( reform(XF, tx*ty, m) )
		

	;computing a matrix containing the values of |f| in the 2d fft arrays

	fx = dindgen(tx , 1)/tx
	fx = fx - double(fx gt 0.5 ) 
	fy = dindgen(1 , ty)/ty 
	fy = fy - double(fy gt 0.5)
	fmap = sqrt(   (fx * fx) # (dblarr(1, tx )+1.)  +    (dblarr(ty, 1) +1. ) # (fy * fy)   )  
	
	; computing the spectral covariances
	tol = 0.	;possible extra overlapping between bins

	for num_fring = 0.,q-1. do begin
		
		fstart = double( freq_rings(0,num_fring) )
		fstop  = double( freq_rings(1,num_fring) ) 
		f_okay   = WHERE(  (  fmap GE (fstart-tol)) AND ( fmap LE (fstop+tol)  )  )


		fring_weight = double(  (size(f_okay))(1)  )
		
		if fring_weight EQ 0.$
		then begin
			print, 'Error : Some specified frequency rings are empty.'
			goto, DONE
		endif
			
		XF_okay = XF(*,f_okay)
		for num_channel = 0., m-1. do XF_okay(num_channel, *) = XF_okay(num_channel, *) - mean( XF_okay(num_channel, *) )
	
		Stats.covmat(*,*,num_fring) = real_part( XF_okay # transpose(conj(XF_okay)) ) / fring_weight 
		Stats.weight(num_fring) = fring_weight 
		;for stationary maps, the statistical weight on each covariance is the number of fourier modes in each spectral ring
		
		Stats.fmean(num_fring) = (fstart+fstop)/2. 
		Stats.fwidth(num_fring) = (fstop-fstart)   

	endfor

	totalweight = total( Stats.weight, /double)
	Stats.weight = Stats.weight / totalweight
	;normalizing the weights
	Stats.frings = freq_rings
		

DONE:

return
end
