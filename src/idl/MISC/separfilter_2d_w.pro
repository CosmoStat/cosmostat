;-------------------------------------------------------------------------------
;
; NAME:
;	separfilter_2d_w.pro
;	
; PURPOSE:
;	Filters the data X using the optimized separating filter, in order to get estimates
;   of the component sources S as in X = AS + B
;
; EXPLANATION:
;   Further versions should allow to use differrent filters in different scales
;
; CALLING SEQUENCE:
;	separfilter_2d_w,  X, Stats, Param, filter_type, S, MASK = mask
;
;
; INPUTS:
;	X : tx*ty*m matrix containing the observed data maps of size tx*ty in m channels
;
;	Stats : structure grouping spectral covariances computed from the data.
;			this is an output of the statistics computation step.
;		NOT USED HERE BUT WE KEEP IT FORUNIFORMITY WITH OTHER separfilter PROCEDURES
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
;   MASK : tx*ty boolean array with MASK(i,j) = 1 indicating valid samples in X( i,j, *)
;								  MASK(i,j) = 0 indicating invalid/masked samples in X( i,j, *)
;
;
;			
;			THE INPUT DATA X SHOULD BE MASKED IF NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO CALLING THE PROCEDURE.
;			ESTIMATED SOURCE SIGNALS ARE MULTIPLIED BY THE MASK BEFORE OUTPUT 
;
;
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
;	separfilter_2d_w,  X, Stats, Param, filter_type, S, MASK = mask
;
; MODIFICATION HISTORY:
;   Yassir Moudden , 2005
;-------------------------------------------------------------------------------   
pro separfilter_2d_w, X, Stats, Param, filter_type, S, MASK = mask

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
	tx = double((size(X))(1))		;data maps have size tx*ty
	ty = double((size(X))(2))

	;building the filter response for the q scales
	makefilter, Param, filter_type, F
	
	;computing a trous wavelet transform of the data
	XW = dcomplexarr(tx, ty, m, q)
	for num_channel = 0 , m-1 do begin
		 atwt2d, X(*,*,num_channel), temp_out, Nscale = q
		 XW(*,*,num_channel, *) = temp_out
	end	 
	;THIS IS CLEARLY USELESS IF THE XW WERE SAVED TO DISK WHILE EXECUTING 
	;data2stats_2d_w. NEEDS TO BE MADE POSSIBLE IN FUTURE VERSIONS.



	;------------------------------------------------------------------------
	; scalewise filtering of the data maps
	;and inverse wavelet transform to get the spatial maps
	;this is a simple pixelwise sum of the different filtered maps
	SWq = dblarr(n, tx*ty)

	for num_scale = 0, q-1 do begin
		Xq = transpose( reform( XW(*, *, *, num_scale), tx*ty, m ) )
		SWq = SWq + F(*,*,num_scale) # Xq
	endfor

	
	S = dblarr(tx, ty, n)
	if mask_on $
	then for num_source= 0 , n-1 do S(* , * , num_source ) = mask * real_part( reform( SWq(num_source, * ),tx, ty) ) $   
	else for num_source= 0 , n-1 do S(* , * , num_source ) =        real_part( reform( SWq(num_source, * ),tx, ty) )
	;it should be possible to get rid of the loops.	

DONE:

return
end







	

	
