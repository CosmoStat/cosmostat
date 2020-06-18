;-------------------------------------------------------------------------------
;
; NAME:
;	separfilter_sph_f.pro
;	
; PURPOSE:
;	Filters the data X using the optimized separating filter, in order to get estimates
;   of the component sources S as in X = AS + B
;
; EXPLANATION:
;   Further versions should allow to use differrent filters in different bands, especially in the case
;   where the set of frequency bins is not a partition of the frequency intervalle 
;   (this is always the case for data mapped to the sphere: the high frequencies are always cut off)
;
; CALLING SEQUENCE:
;	separfilter_sph_f,  X, Stats, Param, filter_type, S, MASK = mask
;
;
; INPUTS:
;	X : is an array of strings containing the file names for the m observed mixtures in Healpix format
;		the maps should be centered
;
;	Stats : structure grouping spectral covariances computed from the data.
;			this is an output of the statistics computation step.
;			required fields used here are:
;				Stats.lbins :  2*q  array specifying spherical harmonic bins
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
;   MASK = mask: filename, a boolean Healpix map containing 0. to indicate invalid pixels and 1. to indicate valid pixels 
;			
;			THE INPUT DATA X SHOULD HAVE BEEN MASKED PRIOR TO CALLING THIS PROCEDURE IF
;			NECESSARY (MULTIPLIED BY THE MASK) PRIOR TO COMPUTING THE STATISTICS.
;			IT IS THEN NOT USEFUL TO HAVE mask AS AN INPUT HERE. HOWEVER WE KEEP IT TO HAVE A UNIFORMITY
;			BETWEEN ALL data2stats PROCEDURES.
;
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
;   This procedures assumes the alm's of the data maps computed with anafast while calling procedure data2stats_sph_f were saved to disk
;   with file names as in data2stats_sph_f.
;
; PROCEDURES USED:
;		makefilter
;
; EXAMPLE:
;	separfilter_sph_f,  X, Stats, Param, filter_type, S, MASK = mask
;
; MODIFICATION HISTORY:
;   Yassir Moudden , 2005
;-------------------------------------------------------------------------------   
pro separfilter_sph_f, X, Stats, Param, filter_type, S, MASK = mask

	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 5 $
	then begin 
        print, 'CALLING SEQUENCE: separfilter_sph_f, X, Stats, Param, filter_type, S, MASK = mask'
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

	;building the filter response for the q frequency bands
	makefilter, Param, filter_type, F

	nlmax= ceil( Stats.lbins(1,q-1) )		;same as in the call to Stats2data_sph_f 
	read_fits_map, X(0), tmp				
	nbside = npix2nside( (size(tmp))(1) )   ;same as in the call to Stats2data_sph_f

	if mask_on $
	then read_fits_map, mask, masque				;reading mask file for masking on exit
	
	;------------------------------------------------------------------------
	;bandwise filtering of the data maps
	;given the procedures available with healpix, we first loop on the reconstructed map number
	;and then on the frequency band.

	S = strarr(n)
	source_map = 'rec_sph_map_f_'
	
	
	name = X(0)
	fits2alm, index, almarr, strcompress( 'alm_' + name ,  /remove_all)  
	;these are the alms of the masked maps computed during the call to the procedure sphmaps2stats_f  
	

	nb_alms =  (size(almarr))(1)
	alm_S = dcomplexarr(  nb_alms   )
	alm_S_tmp = almarr
	;initializing	
	

	for num_source = 0, n-1 do begin			;loop over the reconstructed map numbers
	
		alm_S(*) = 0.

		S(num_source) = strcompress( source_map + string(floor(num_source)) + '.fits',  /remove_all)  
		poids  = 0.*dblarr(  nb_alms )
	
		for num_q = 0, q-1 do begin				;loop over the frequency bands
		
			lstart = floor(  Stats.lbins(0,num_q)  )
			lstop  = ceil ( Stats.lbins(1,num_q)   )		
		
			lm2index2,lstart,0,indexmin				
			lm2index2,lstop+1,0,indexmax			
			indexmax = indexmax - 1.				
			
			lsize = double( indexmax) - double( indexmin) +1.

			alm_okay = dcomplexarr(m,lsize)
   			
			for num_map=0,m-1 do begin
				name = X(num_map)
				fits2alm,index,alm_tmp,strcompress( 'alm_'+ name,  /remove_all)
				alm_okay[num_map,*] = complex(alm_tmp(indexmin:indexmax, 0), alm_tmp(indexmin:indexmax, 1))
			endfor
	
			alm_S(indexmin:indexmax) = alm_S(indexmin:indexmax)  + F(num_source,*,num_q) # alm_okay    ;is centering of alm_okay necessary here?
			poids(indexmin:indexmax) = poids(indexmin:indexmax) +1.
  
	 
		endfor
   
	poids = poids + (poids LT 1.)
	alm_S = alm_S / poids
	alm_S_tmp(*,0) = real_part(alm_S)
	alm_S_tmp(*,1) = imaginary(alm_S)
   
   
	alm2fits, index, alm_S_tmp, 'alm_S.fits'
   
   
	;synthesize map from the alms and save it
    command3 = 'command3.dat'
    openw, com3,command3,/get_lun
    printf,com3,"infile=''"
    printf,com3,'outfile=!'+S(num_source)
    printf,com3,'almsfile='+'alm_S.fits'		
	printf,com3,'nsmax='+string(floor( nbside) )			;same as initial data map
    printf,com3,'nlmax='+string(floor( nlmax ) )			;same as initial data map
    printf,com3,'fwhm_arcmin =0'
    printf,com3,'iseed=-1'
    free_lun,com3
	spawn,'synfast '+command3

	if mask_on then begin
		read_fits_map, S(num_source), temp_source			;multiplication on exit by the initial mask and reordering
		temp_source = reorder(temp_source,in='ring',out='nested')
		temp_source = real_part( temp_source) *masque
		write_fits_map, S(num_source), temp_source, /nested
	endif else begin
		read_fits_map, S(num_source), temp_source			;reorder the map on output
		temp_source = reorder(temp_source,in='ring',out='nested')
		write_fits_map, S(num_source), temp_source, /nested
	endelse	

endfor


DONE:

return
end


