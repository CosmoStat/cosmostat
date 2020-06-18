;-------------------------------------------------------------------------------
;+
; NAME:
;	data2stats_1d_f.pro
;	
; PURPOSE:
;	computes mean covariance matrices in the fourier representation over 
;   specified spectral bands 
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	data2stats_1d_f, X, freq_bins, Stats, MASK = mask
;
; INPUTS:
;	X : m*T matrix containing the observed data in m channels during "time" T
;   freq_bins : 2*q matrix containing the starting and ending reduced frequencies of each of 
;				q frequency intervals. Each specified reduced frequency is between [0, 0.5[
;				and the qth average covariance matrix is computed over intervalle 
;				]-freq_bins(1, q), -freq_bins(0, q)]U[freq_bins(0, q), freq_bins(1, q)[ to have 
;				real statistics.
;
; OPTIONAL INPUTS:
;   MASK : 1*T boolean array with MASK(0, t) = 1 indicating valid samples in X(*, t)
;								  MASK(0, t) = 0 indicating invalid/masked samples in X(*, t)
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
;				Stats.covmat :  m*m*q 3D array of q m*m covaraince matrices
;								q is the number of bands or scales, etc.
;								m is the number of data channels
;
;				Stats.weight :  statistical significance weights of the covariance matrices
;
;			other fields for graphical purposes, source estimation, etc :
;		
;				Stats.fmean : 1*q array containing the central frequency of each freq_bin
;				Stats.fwidth: 1*q array containing the frequency width of each  freq_bin
;				Stats.fbins:  2*q array equal to freq_bins
;	
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	data2stats_1d_f, X, freq_bins, Stats, MASK = mask
;
; MODIFICATION HISTORY:
;		Written: Yassir Moudden 2005.
;-------------------------------------------------------------------------------   

pro data2stats_1d_f, X, freq_bins, Stats, MASK = mask

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: data2stats_1d_f, X, freq_bins, Stats, MASK = mask'
        goto, DONE
	end
	
	if keyword_set(mask) EQ 0 $				;the fourier transform is not local
	then mask_on = 0 $						;this is then quite useless
	else mask_on = 1						;however we keep for the homogeneity of the code

	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	;recovering a few parameter values
   	q = (size(freq_bins))(2)		;number of spectral bands/of covariance matrices to be computed
	T = double((size(X))(2))		;number of data samples
	m = double((size(X))(1))		;number of channels
  
	;------------------------------------------------------------------------
	;create anonymous structure in which to group the summary covariance statistics
	Stats = {covmat:dblarr(m,m,q), weight:dblarr(1,q), fmean:dblarr(1,q), fwidth:dblarr(1,q), fbins:dblarr(2, q)}
  
	;computing FFT of the data
	XF = dcomplexarr(m, T)
	for num_channel=0.,  m-1. do XF(num_channel, * ) = T * fft(  X(num_channel, * ) - mean( X(num_channel, * ) ), -1, /double )
	;it should be possible to escape the loop!
	;centering is because we are computing covariance matrices further down.
      
	;computing a matrix containing the values of |f| in the 1d fft arrays i.e. a frequency index
	fx = dindgen(T , 1)/T
	fx = abs( fx - double(fx gt 0.5 ) ) 

	; computing the spectral covariances
	tol = 0.	;possible extra overlapping between bins
	
	for num_fbin = 0.,q-1. do begin
   
		fstart = double( freq_bins(0,num_fbin) )
		fstop  = double( freq_bins(1,num_fbin) ) 
		fx_okay   = WHERE(  (  fx GE (fstart-tol)) AND ( fx LE (fstop+tol)  )  )

		fbin_weight = double(  (size(fx_okay))(1)  )
		
		if fbin_weight EQ 0.$
		then begin
			print, 'Error : Some specified frequency bins are empty.'
			goto, DONE
		endif
			
		XF_okay = XF(*,fx_okay)
		for num_channel = 0., m-1. do XF_okay(num_channel, *) = XF_okay(num_channel, *) - mean( XF_okay(num_channel, *) )
	
		Stats.covmat(*,*,num_fbin) = real_part( XF_okay # transpose(conj(XF_okay)) ) / fbin_weight 
		Stats.weight(num_fbin) = fbin_weight 
		;for stationary maps, the statistical weight on each covariance is the number of fourier modes in each spectral band
		
		Stats.fmean(num_fbin) = (fstart+fstop)/2. 
		Stats.fwidth(num_fbin) = (fstop-fstart)   

	endfor

	totalweight = total( Stats.weight, /double)
	Stats.weight = Stats.weight / totalweight
	;normalizing the weights
	Stats.fbins = freq_bins
	;------------------------------------------------------------------------


DONE:

return
end
