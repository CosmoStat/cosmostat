;-------------------------------------------------------------------------------
;+
; NAME:
;	data2stats_sph_f.pro
;	
; PURPOSE:
;	computes mean covariance matrices in the fourier (spectral harmonics) representation over 
;   specified spectral bands e.g. l in [l_min, l_max]
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	data2stats_sph_f, X, l_bins, Stats, MASK = mask
;
; INPUTS:
;	X : is an array of strings containing the file names for the m observed mixtures in Healpix NESTED format
;		the maps should be centered.
;
;   l_bins : 2*q matrix containing the starting and ending multipole number l of each of 
;				q spectral harmonics domains.
;
; OPTIONAL INPUTS:
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
;				Stats.lmean : 1*q array containing the mean frequency vector norm over ring q
;				Stats.lwidth: 1*q array containing the width in frequency vector norm of ring q
;				Stats.lbins:  2*q array equal to l_bins
;	
;
; DEPENDENCIES:
;		THe Healpix package for handling spherical maps and computing forward and reverse 
;		spherical harmonics transforms is used.
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	data2stats_sph_f, X, l_bins, Stats, MASK = mask
;
; MODIFICATION HISTORY:
;		Written: Yassir Moudden 2005.
;
; COMMENTS:
;   
;
;-------------------------------------------------------------------------------   

pro data2stats_sph_f, X, l_bins, Stats, MASK = mask

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: data2stats_sph_f, X, l_bins, Stats, MASK = mask'
        goto, DONE
	end
	
	if keyword_set(MASK) eq 0 $			;rather useless since the fourier transform is not local!
	then mask_on = 0 $
	else mask_on = 1
	
	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	;recovering a few parameter values
   	q = (size(l_bins))(2)			;number of spectral harmonics bands, number of covariance matrices to be computed
	m = double((size(X))(1))		;number of channels
		
	;------------------------------------------------------------------------
	;create anonymous structure in which to group the summary covariance statistics
	Stats = {covmat:dblarr(m,m,q), weight:dblarr(1,q), lmean:dblarr(1,q), lwidth:dblarr(1,q), lbins:dblarr(2, q)}

	


	;computing the spherical harmonics transform of each of the data maps
	for num_map = 0.,m-1. do begin
		name = X(num_map)
		;read_fits_map, name, tmp
		
		;the following generates a file used by anafast from the Healpix package
		command = 'command.dat'
		openw,com,command,/get_lun
		printf,com,'infile='+name
		printf,com,'outfile=!'+strcompress('power_'+name,  /remove_all)  
		printf,com,'nlmax='+string( ceil( l_bins[1,q-1] ) )
		printf,com,'simul_type='+string(1)
		printf,com,'outfile_alms=!'+strcompress( 'alm_'+name,  /remove_all)			; this file contains the spherical harmonics transform of the map
		free_lun,com																; in file name. we keep it for use in the filtering step after the 
																					; model parameter optimization procedure
																					; This way, the spherical harmonics transform is computed only once.
		spawn,'anafast '+command  
		;how about centering the maps before computing the alms?
		;centering should be done beforehand in the observed data maps
		
	endfor
     

	;computing the source covariance matrix in each of the q multipole bands
	for num_q = 0.,q-1. do begin
		lstart = floor( l_bins(0,num_q) ) 
		lstop  = ceil( l_bins(1,num_q) ) 
		
		lm2index2,lstart,0,indexmin				
		lm2index2,lstop+1.,0,indexmax			
		indexmax = indexmax - 1.			
		
		l_0 = lstart + findgen(lstop-lstart + 1)
		lm2index2, l_0, 0.*l_0, l_0_index			;recovering the index corresponding to the al0 for l in the specified range
		
		lsize_bis = 2.* ( double(indexmax) - double( indexmin) +1. ) 
		lsize = lsize_bis  -  (lstop-lstart + 1.)
		
		alm_okay = dcomplexarr(m , lsize_bis)
    
		for num_map=0.,m-1. do begin
			name = X(num_map)
			fits2alm,index,alm_tmp,strcompress( 'alm_'+ name,  /remove_all)
			alm_okay[num_map,*] = complex([alm_tmp(indexmin:indexmax, 0) , conj( alm_tmp(indexmin:indexmax, 0) ) ]  , [ alm_tmp(indexmin:indexmax, 1) , conj( alm_tmp(indexmin:indexmax, 1) ) ] )
			alm_okay[num_map,l_0_index] = 0.
			alm_okay[num_map,*] = alm_okay[num_map,*] - total( alm_okay[num_map,*], /double )/lsize								
			alm_okay[num_map,l_0_index] = 0.
		endfor
     
		Stats.covmat[*,*,num_q] = real_part( (alm_okay ) # transpose(  conj(alm_okay )) ) / lsize 
		
		;for stationary maps, the statistical weight on each covariance is the number of sph. harm.  modes in each spectral ring
		Stats.weight[num_q] = lsize		
				
		Stats.lmean[num_q] = (lstart+lstop)/2. 
		Stats.lwidth[num_q] = (lstop-lstart)   
	

	endfor

		totalweight = total( Stats.weight,/double)
		Stats.weight = Stats.weight / totalweight
		
		Stats.lbins =  l_bins


DONE:

return
end



;-------------------------------------------------------------------------------

pro lm2index2,l,m,index

index = ((long(l)*(long(l) + 1l)) / 2l) + m

end




