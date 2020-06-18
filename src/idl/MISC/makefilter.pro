;-------------------------------------------------------------------------------
;
; NAME:
;	makefilter.pro
;	
; PURPOSE:
;	Compute a separating filter in each of the q bands or scales, etc   
;
; EXPLANATION:
;   Future versions should include some sort of non linear filtering technique.
;
; CALLING SEQUENCE:
;	makefilter, Param, filter_type, F
;
;
; INPUTS:
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
;					NB "passband" is somewhat useless!!
;
; OPTIONAL INPUTS:
;  
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	F :  n*m*q matrix containing the n*m separating matrices in each of the q bands
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		
;
; EXAMPLE:
;	makefilter, Param, filter_type, F
;
; MODIFICATION HISTORY:
; Yassir Moudden and Jean-Francois Cardoso ,  2005
;-------------------------------------------------------------------------------   
pro	makefilter, Param, filter_type, F

	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: makefilter, Param, filter_type, F'
        goto, DONE
	end
	
	;checking that the filter type is supported meaning either "wiener" or "pinv" or "passband" or "pinvA"
	if ( (not strcmp(filter_type, 'wiener')) AND (not strcmp(filter_type , 'pinv')) AND $
		(not strcmp(filter_type , 'pinvA')) AND (not strcmp(filter_type, 'passband'))  ) $
	then begin
		print, 'Filter type not supported. Should be either wiener, pinv, pinvA or passband.'
		goto, DONE
	endif
	
	;checking that the parameter type is supported meaning either 1 or 2  
	if (Param.type NE 1) AND (Param.type NE 2)  then begin
		print, 'Incorrect parameter type. Type should be either 1 or 2.'
		goto, DONE
	endif
	;------------------------------------------------------------------------

	
	;------------------------------------------------------------------------
	;recovering a few model parameters
	A = Param.mixmat 
	Rp = Param.source 
	Rn = Param.noise  
	
	m = double( (size(A))(1) )
	n = double( (size(A))(2) )
	q = double( (size(Rp))(2) )
	
	F = dblarr(n,m,q)
	
	;in the case where parameter type is 1, extend the white noise to all bands  
	if (Param.type EQ 1) then Rn = Rn # (0.* dblarr(1,q) + 1. )

	;compute the filter matrix in the different spectral bands
	iRn = 0.* dblarr(m,m)  
	iRp = 0.* dblarr(n,n)   
	
	for num_q = 0, q-1 do begin
		iRn = diag_matrix( 1.0 / Rn(*,num_q))  ;validity of these matrix inverses?
		iRp = diag_matrix( 1.0 / Rp(*,num_q))  
	    
		case filter_type of
	
			'wiener' : F(*,*,num_q)   = invert ( transpose(A)#iRn# A + iRp, /double) # transpose(A) # iRn  
	
			'pinv' : F(*,*,num_q)   = invert ( transpose(A)#iRn# A , /double) # transpose(A) # iRn 
			
			'pinvA' : F(*,*,num_q)   = invert ( transpose(A)# A , /double) # transpose(A)
	
			'passband': F(*,*,num_q) = identity(m)  
	
		endcase
    
	endfor
  
  
  
DONE: 
return
  
end

