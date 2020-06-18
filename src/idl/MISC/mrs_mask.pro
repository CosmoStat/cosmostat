;-------------------------------------------------------------------------------
;+
; NAME:
;		mrs_mask.pro
;	
; PURPOSE:
;		When gaps exist in a signal or a map, some wavelet coefficients located 
;		outside the initial mask are affected. The extent of the influence of the
;		mask depends on scale. The purpose of this function is to apply the specified wavelet
;		transform to the specified mask and to return a mask on each scale where 1s  
;		correspond to valid coefficients (i.e. coefficient which are contaminated by the mask but
;		below some threshold) and 0s correspond to contaminated coefficients. 
;
;		Implemented for three different topologies and "two" different 
;		transforms (ie undecimated a trous algorithm or orthogonal transform ): the undecimated
;		transform is the one used in mrs_smica whereas the orthogonal transform is used in mrs_jade.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;		mrs_mask, mask, topology, wt_type, nb_scales, mask_out,  nlmax = nlmax
;
; INPUTS:
;   mask :  either a length T array in the '1D' case, 
;			or a tx*ty array in the flat '2D' case,
;			or a spherical mask in HEALPIX NESTED FORMAT in the 'Sphere' case.
;
;   topology:   string specifying the topology of the mask
;				this should be either '1D' or '2D' or 'Sphere'
;				This is clearly redundant information, but makes things more convenient.
;				Any incoherence between the specified 'topology' and the structure of the input data leads to an error.
;
;   nb_scales : a scalar specifying a number of wavelet scales (including the smoothed array)
;
;   wt_type : string specifying the wavelet transform type
;			this should be either  'atrous' or 'ortho' respectively for the undecimated 
;			a trous wavelet transform (and extensions in different topologies) or the orthogonal wavelet transform (and extensions in the differenet topologies)
;
; OPTIONAL INPUTS:
;
;   nlmax : this is used in the undecimated spherical wavelet transform so when the specified topology is 'Sphere' and the 
;			specified transform is 'atrous'.
;			nlmax is not optional in the 'Sphere' AND 'atrous' case.
;			NB:The same value of nlmax should be used as in the spherical wavelet transform of the data maps.
; 
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	mask_out :  if wt_type = 'atrous', this is either 
;				an nb_scales*T array in the 1D case, 
;				or a tx*ty*nb_scales array in the flat 2D case,
;				or npix*nb_scales array of nb_scales spherical masks in healpix format where npix 
;				is the size of the initial mask in healpix format
;
;				if wt_type = 'ortho', this is either a length T array in the 1D case, 
;				or a tx*ty array in the flat 2D case,
;				or an array the same sizeas the initial mask in Healpix nested format 
;
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		the mrs package is used
;		the healpix package is used
;
; EXAMPLE:
;
;   mrs_mask, mask, 'sphere', 'atrous', 5, mask_out, nlmax = 512
;   
; MODIFICATION HISTORY:
;	Yassir Moudden , 2005
;-------------------------------------------------------------------------------   

pro mrs_mask, mask, topology, wt_type, nb_scales, mask_out,  nlmax = nlmax

COMMON MR1ENV


	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters and there coherence
	
	if N_PARAMS() LT 5 $
	then begin 
        print, 'CALLING SEQUENCE: mrs_mask, mask, topology, wt_type, nb_scales, mask_out,  nlmax = nlmax'
        goto, DONE
	end
		
	;making sure there are no useless dimensions in the specified mask
	mask = double(reform(mask))
	
	case topology of
		'1D'		: begin
						size_mask = size(mask) 
						if (size_mask(0) NE 1)$
						then begin
							print, 'error : dimensions of specified mask are incompatible with specified topology'
							goto, DONE
						end

						T = double(size_mask(1))
						
					  end	

		'2D'		: begin
						size_mask = size(mask) 
						if (size_mask(0) NE 2)$
						then begin
							print, 'error : dimensions of specified mask are incompatible with specified topology'
							goto, DONE
						end
		
						tx = double(size_mask(1))		;data maps have size tx*ty
						ty = double(size_mask(2))
					  end	
					
		'Sphere'	: begin
						size_mask = size(mask) 
						if (size_mask(0) NE 1)$
						then begin
							print, 'error : dimensions of specified mask are incompatible with specified topology'
							goto, DONE
						end
						
						npix = (size_mask)(1)
						nside = npix2nside(npix)

					  end
					  
		else		: begin
						print, 'error : specified topology is incorrect'
						goto, DONE
					  end	
					  
	endcase
	
	case wt_type of 
	
		'atrous'	: begin
						if ( topology eq 'Sphere' ) AND (keyword_set( nlmax) eq 0 ) $
						then begin
							print, 'error: nlmax is not an optional input in this context (ie Sphere and atrous)'
							goto, DONE
						endif	
					  end	
														
		'ortho'		: begin
					  end	
					  
		else		: begin
						print, 'error : specified transform is incorrect'
						goto, DONE
					  end

	endcase
	
	;------------------------------------------------------------------------


	;----------------------------------------------------------------------------------------
	;computing the wavelet transform of the mask and determining a mask on each scale
	;----------------------------------------------------------------------------------------

	case wt_type of 
	
	'atrous' : begin
					case topology of
					
					'1D'		: begin
										mask_out = dblarr(nb_scales, T)
										atwt1d, double( not( double( mask) )) , mask_out, Nscale = nb_scales
										;on each scale, the nul coefficients in  mask_out which are outside the
										;initial mask indicate a valid coefficient in the transformed data
										for num_scale = 0, nb_scales-2 do mask_out(num_scale, *) = (   (mask_out(num_scale, *) eq 0.) AND ( mask eq 1.)  )	
										mask_out(nb_scales-1, *) = mask_out(nb_scales-2, *)
										;at this point, 1s in mask_out indicate valid coefficients and 0s indicate invalid coefficients
								
								  end		;of the '1D' and 'atrous' case
									
					'2D'		: begin
										mask_out = dblarr( tx, ty, nb_scales)
										atwt2d, double(  not( double( mask)) ) , mask_out, Nscale = nb_scales
										; on each scale, the nul pixels in the mask_out which are outside the
										;initial mask indicate a valid coefficient in the transformed data
										for num_scale = 0, nb_scales-1 do mask_out(*, *, num_scale) = (   (mask_out(*, *, num_scale) eq 0) AND ( mask eq 1)  )	
										mask_out(nb_scales-1, *) = mask_out(nb_scales-2, *)
										;at this point, 1s in mask_out indicate valid coefficients and 0s indicate invalid coefficients
							
								  end		;of the '2D' and 'atrous' case
														
					'Sphere'	: begin
										mrs_wttrans, double( not( double( mask) )), wave_not_mask, NbrScale = nb_scales, lmax = nlmax
										mask_out = wave_not_mask.coef
										cumul_mask =  double( mask)
										for num_scale = 0, nb_scales - 2. do begin
											tmp_mask =  wave_not_mask.coef(*,num_scale)
											tmp_mask_max = max( tmp_mask^2 ) 
											tmp_mask = not( double(   abs(tmp_mask) GT  1e-2 *tmp_mask_max ) )
											cumul_mask = double( cumul_mask AND  tmp_mask)
											mask_out(*,num_scale) = cumul_mask
										endfor
										mask_out(*, nb_scales -1 ) = mask_out(*, nb_scales - 2)
					
								  end		;of the 'Sphere' and 'atrous' case

					endcase
				end
				
	'ortho' : begin
					case topology of
					
					'1D'		:  begin								;utiliser une ondelette a support compact dans la representation initiale
									mask_out = dblarr(T)
									temp_mask = double( not( double(mask) ))
									length_temp_mask = T
									
									if keyword_set(mr1ok) then begin					;if package mre is available
										
										wt_options = strcompress('-t15 -L ' + '-n'+ string(nb_scales)  )
										mr1d_trans, temp_mask , wave_temp_mask , Opt = wt_options

										for num_scale = 0, nb_scales-2 do begin
											start_index = wave_temp_mask.from(num_scale)
											end_index = wave_temp_mask.to(num_scale)

										
										if (end_index eq (2*start_index -1) ) $
											then begin
												tmp_mask = temp_mask( 2.*indgen( length_temp_mask/2.) ) OR temp_mask( 2.*indgen( length_temp_mask/2.) + 1. )
												mask_out(start_index:end_index) = tmp_mask OR (abs( wave_temp_mask.coef(start_index:end_index)) gt 1.e-4)
											
												length_temp_mask = length_temp_mask/2.
												temp_mask = mask_out(start_index:end_index) 
											
											endif else begin
												tmp_mask = temp_mask( 2.*indgen( floor(length_temp_mask/2.)) +1  )
												mask_out(start_index:end_index) = tmp_mask OR (abs( wave_temp_mask.coef(start_index:end_index)) gt 1.e-4)
										
												length_temp_mask = ceil(length_temp_mask/2.)
												temp_mask = temp_mask( 2.*indgen( length_temp_mask )  )	
												temp_mask(0:length_temp_mask-2) = temp_mask(0:length_temp_mask-2) OR  mask_out(start_index:end_index)
												temp_mask(1:length_temp_mask-1) = temp_mask(1:length_temp_mask-1) OR  mask_out(start_index:end_index)
									
											endelse 
	
										endfor
									
										start_index = wave_temp_mask.from(nb_scales - 1 )
										end_index = wave_temp_mask.to(nb_scales -1)

										mask_out( start_index:end_index ) = temp_mask
										mask_out =double( not( double(mask_out) )) 
										;at this point, 1s in mask_out indicate valid pixels and 0s indicate invalid pixels
									
									endif else begin			;use procedures in the mrs package
									
										wave_temp_mask = bwt01_lift(temp_mask, nb_scales-1)
										
										temp_n = T
										wave_temp_mask_index=lonarr(nb_scales+1)
										for num_band = 0, nb_scales-1 do begin 
											wave_temp_mask_index(num_band)=temp_n
											temp_mn= temp_n mod 2
											temp_n=temp_n/2+temp_mn    
										end   

										
										for num_scale = 0, nb_scales-2 do begin
																					
											start_index = wave_temp_mask_index(num_scale+1)
											end_index = wave_temp_mask_index(num_scale)-1
										
											if (end_index eq (2*start_index -1) ) $
											then begin
												tmp_mask = temp_mask( 2.*indgen( length_temp_mask/2.) ) OR temp_mask( 2.*indgen( length_temp_mask/2.) + 1. )
												mask_out(start_index:end_index) = tmp_mask OR (abs( wave_temp_mask(start_index:end_index)) gt 1.e-4)
											
												length_temp_mask = length_temp_mask/2.
												temp_mask = mask_out(start_index:end_index) 
											
											endif else begin
												tmp_mask = temp_mask( 2.*indgen( floor(length_temp_mask/2.)) +1  )
												mask_out(start_index:end_index) = tmp_mask OR (abs( wave_temp_mask(start_index:end_index)) gt 1.e-4)
										
												length_temp_mask = ceil(length_temp_mask/2.)
												temp_mask = temp_mask( 2.*indgen( length_temp_mask )  )	
												temp_mask(0:length_temp_mask-2) = temp_mask(0:length_temp_mask-2) OR  mask_out(start_index:end_index)
												temp_mask(1:length_temp_mask-1) = temp_mask(1:length_temp_mask-1) OR  mask_out(start_index:end_index)
									
											endelse 
	
										endfor
									
										start_index = wave_temp_mask_index(nb_scales  )
										end_index = wave_temp_mask_index(nb_scales - 1) -1
										
										mask_out( start_index:end_index ) = temp_mask
										mask_out =double( not( double(mask_out) )) 
										;at this point, 1s in mask_out indicate valid pixels and 0s indicate invalid pixels

									endelse
									
									
									
								   end			;of the '1D' and 'ortho' case 
								
									
					'2D'		:   begin								;utiliser une ondelette a support compact dans la representation initiale
									mask_out = dblarr(tx, ty)
									temp_mask = double( not( double( mask) )  )
									tx_temp_mask = tx
									ty_temp_mask = ty

									if keyword_set(mr1ok) then begin					;if package mre is available

										wt_options = '-t14 -L ' + strcompress('-n'+ string(nb_scales) )		;using an L2 normalized orthogonal wavelet transform 
										mr_transform, temp_mask, wave_temp_mask, Opt=wt_options	
									
									endif else begin									;use procedures in mrs
										wave_temp_mask = bwt01_lift(temp_mask, nb_scales-1)
									endelse		
												
									for num_scale = 0, nb_scales-2 do begin

										
										if floor(tx_temp_mask/2.) eq tx_temp_mask/2. then begin
										
												if floor(ty_temp_mask/2.) eq ty_temp_mask/2. then begin		;tx_temp_mask is pair and ty_temp_mask is pair
																
														low_x =	findgen(  tx_temp_mask/2.  )
														high_x = findgen(  tx_temp_mask/2.  ) + tx_temp_mask/2.
														low_y =	findgen( ty_temp_mask/2. )   
														high_y = findgen( ty_temp_mask/2.  ) + ty_temp_mask/2.
														
														tmp_mask = (temp_mask(2.*indgen( tx_temp_mask/2.), *))(*  ,   2.*indgen( ty_temp_mask/2.)  ) $
																OR (temp_mask( 2.*indgen( tx_temp_mask/2.) + 1.,*))(*,  2.*indgen( ty_temp_mask/2.) ) $
																		OR (temp_mask( 2.*indgen( tx_temp_mask/2.) ,*))(*,  2.*indgen( ty_temp_mask/2.) +1.) $
																			OR (temp_mask( 2.*indgen( tx_temp_mask/2.) + 1.,*))(*,  2.*indgen( ty_temp_mask/2.) +1. ) 
												
														temp_matrix = tmp_mask $
																	OR (abs( (wave_temp_mask(high_x ,*))(*, high_y) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(low_x , *))(*, high_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(high_x ,*))(*,low_y) ) gt 1.e-4)
														put_submatrix, mask_out, high_x, high_y, temp_matrix
														put_submatrix, mask_out, low_x, high_y, temp_matrix
														put_submatrix, mask_out, high_x, low_y, temp_matrix

														temp_mask = (mask_out(high_x ,*))(*,  high_y) 
														tx_temp_mask = tx_temp_mask/2.
														ty_temp_mask = ty_temp_mask/2.
												
												endif else begin			;tx_temp_mask is pair and ty_temp_mask is impair
												
														low_x =	findgen( tx_temp_mask/2.  )
														high_x = findgen( floor( tx_temp_mask/2.)  ) + tx_temp_mask/2. 
														low_y =	findgen( floor( ty_temp_mask/2.)  )   
														high_y = findgen( floor( ty_temp_mask/2.)  ) + ceil( ty_temp_mask/2.)

														tmp_mask = (temp_mask( 2.*indgen( tx_temp_mask/2.),*))(*,  2.*indgen(floor( ty_temp_mask/2.)) +1.) $
																	OR (temp_mask( 2.*indgen( tx_temp_mask/2.) + 1.,*))(*,  2.*indgen(floor( ty_temp_mask/2.)) +1. )

														temp_matrix = tmp_mask $
																	OR (abs( (wave_temp_mask(high_x ,*))(*, high_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(low_x ,*))(*, high_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(high_x,*))(*, low_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(high_x ,*))(*,low_y+1.) ) gt 1.e-4)
														
														put_submatrix, mask_out, high_x, high_y, temp_matrix
														put_submatrix, mask_out, low_x, high_y, temp_matrix
														put_submatrix, mask_out, high_x, low_y, temp_matrix
														temp_matrix = temp_matrix OR (mask_out(high_x ,*))(*, low_y+1)
														put_submatrix, mask_out, high_x, low_y+1., temp_matrix
														
														temp_mask = (temp_mask( 2.*indgen( tx_temp_mask/2.) ,*))(*,  2.*indgen(ceil( ty_temp_mask/2.)) -1.) $
																	OR (temp_mask( 2.*indgen( tx_temp_mask/2.) + 1.,*))(*,  2.*indgen(ceil( ty_temp_mask/2.)) -1. ) 
														
														tx_temp_mask = tx_temp_mask/2.
														ty_temp_mask = ceil( ty_temp_mask/2.)
											
														temp_mask(0: tx_temp_mask-1 , 0: ty_temp_mask -2 ) = temp_mask(0: tx_temp_mask-1 , 0: ty_temp_mask -2 )$
																											 OR ( mask_out( high_x,*))(*, high_y)
														temp_mask(0: tx_temp_mask-1 , 1: ty_temp_mask -1 ) = temp_mask(0: tx_temp_mask-1 , 1: ty_temp_mask -1 ) $
																											 OR  (mask_out(high_x ,*))(*, high_y)
														
												endelse
										
										endif else begin 
												
												if floor(ty_temp_mask/2.) eq ty_temp_mask/2. then begin		;tx_temp_mask is impair and ty_temp_mask is pair
	
														low_x =	findgen( floor( tx_temp_mask/2.)  )
														high_x = findgen( floor( tx_temp_mask/2.)  ) + ceil( tx_temp_mask/2.) 
														low_y =	findgen( ty_temp_mask/2.  )   
														high_y = findgen( ty_temp_mask/2.  ) + ty_temp_mask/2.
														
														tmp_mask = (temp_mask(2.*indgen(floor( tx_temp_mask/2.)) +1.,*))(*,  2.*indgen( ty_temp_mask/2.) ) $
																	OR (temp_mask(2.*indgen(floor( tx_temp_mask/2.)) +1.,*))(*,  2.*indgen( ty_temp_mask/2.) + 1. ) 

														temp_matrix = tmp_mask $
																	OR (abs( (wave_temp_mask(high_x ,*))(*, high_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(high_x ,*))(*, low_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(low_x,*))(*, high_y ) ) gt 1.e-4)$
																	OR (abs( (wave_temp_mask(low_x+1 ,*))(*,high_y) ) gt 1.e-4)
														
														put_submatrix, mask_out, high_x, high_y, temp_matrix
														put_submatrix, mask_out, high_x, low_y, temp_matrix
														put_submatrix, mask_out, low_x, high_y, temp_matrix
														temp_matrix = temp_matrix OR (mask_out(low_x+1. ,*))(*, high_y)
														put_submatrix, mask_out, low_x+1., high_y, temp_matrix
														
														
														temp_mask = (temp_mask(   2.*indgen(ceil( tx_temp_mask/2.)) -1.,*))(*, 2.*indgen( ty_temp_mask/2.)) $
																	OR (temp_mask(  2.*indgen(ceil( tx_temp_mask/2.)) -1. ,*))(*, 2.*indgen( ty_temp_mask/2.) + 1.) 
														
														ty_temp_mask = ty_temp_mask/2.
														tx_temp_mask = ceil( tx_temp_mask/2.)
														
														temp_mask( 0: tx_temp_mask -2, 0: ty_temp_mask-1  ) = temp_mask( 0: tx_temp_mask -2 , 0: ty_temp_mask-1 )$
																											 OR  (mask_out(  high_x,*))(*, high_y )
														temp_mask( 1: tx_temp_mask -1 , 0: ty_temp_mask-1 ) = temp_mask( 1: tx_temp_mask -1 , 0: ty_temp_mask-1 ) $
																											 OR ( mask_out(  high_x,*))(*, high_y )
														
														
													endif else begin											;tx_temp_mask is impair and ty_temp_mask is impair

														low_x =	findgen( floor( tx_temp_mask/2.)  )
														high_x = findgen( floor( tx_temp_mask/2.)  ) + ceil( tx_temp_mask/2.) 
														low_y =	findgen( floor( ty_temp_mask/2.)  )   
														high_y = findgen( floor( ty_temp_mask/2.)  ) + ceil( ty_temp_mask/2.)
														
														tmp_mask = ( temp_mask(  2.*indgen(floor( tx_temp_mask/2.)) +1.,*))(*,  2.*indgen(floor( ty_temp_mask/2.)) +1.) 

														temp_matrix = tmp_mask OR (abs( (wave_temp_mask( high_x ,*))(*, high_y ) )  gt 1.e-4 ) $
																					OR (abs( (wave_temp_mask(low_x ,*))(*, high_y ) ) gt 1.e-4) OR (abs( (wave_temp_mask(low_x+1. ,*))(*, high_y ) ) gt 1.e-4)$
																					OR (abs( (wave_temp_mask(high_x,*))(*, low_y ) ) gt 1.e-4) OR (abs( (wave_temp_mask(high_x ,*))(*, low_y + 1. ) ) gt 1.e-4)
														put_submatrix, mask_out, high_x, high_y, temp_matrix
														put_submatrix, mask_out, high_x, low_y, temp_matrix
														put_submatrix, mask_out, low_x, high_y, temp_matrix
														
														temp_matrix1 = temp_matrix OR (mask_out(low_x+1. ,*))(*, high_y)
														put_submatrix, mask_out, low_x+1., high_y, temp_matrix1
														temp_matrix = temp_matrix OR (mask_out(high_x ,*))(*, low_y+1)
														put_submatrix, mask_out, high_x, low_y+1., temp_matrix
														
														temp_mask = (temp_mask( 2.*indgen( ceil( tx_temp_mask/2.))-1.,*))(*,  2.*indgen(ceil( ty_temp_mask/2.)) -1.) 
																
														tx_temp_mask = ceil( tx_temp_mask/2.)
														ty_temp_mask = ceil( ty_temp_mask/2.)
																																								
														temp_mask(0: tx_temp_mask-2 , 0: ty_temp_mask -2 ) = temp_mask(0: tx_temp_mask-2 , 0: ty_temp_mask -2 )$
																											 OR  (mask_out( high_x ,*))(*, high_y)
														temp_mask(0: tx_temp_mask-2 , 1: ty_temp_mask -1 ) = temp_mask(0: tx_temp_mask-2 , 1: ty_temp_mask -1 ) $
																											 OR  (mask_out( high_x ,*))(*, high_y)
														temp_mask(1: tx_temp_mask-1 , 0: ty_temp_mask -2 ) = temp_mask(1: tx_temp_mask-1 , 0: ty_temp_mask -2 )$
																											 OR  (mask_out( high_x ,*))(*, high_y)
														temp_mask(1: tx_temp_mask-1 , 1: ty_temp_mask -1 ) = temp_mask(1: tx_temp_mask-1 , 1: ty_temp_mask -1 ) $
																											 OR  (mask_out( high_x ,*))(*, high_y)
														
													endelse
										endelse		
												 
									endfor
					
									mask_out( 0: tx_temp_mask-1., 0: ty_temp_mask-1. ) = temp_mask
									mask_out = double( not( double(mask_out) )) 
									;at this point, 1s in mask_out indicate valid pixels and 0s indicate invalid pixels
	
								   end			;of the '2D' and 'ortho' case 
								   
 					'Sphere'	: 	begin		;cette procedure sous estime le mask total sur la sphere
												;on ne tient pas compte de l'extension du filtre ondelette a chaque echelle de part 
												;et d'autre de la frontiere entre deux faces healpix.
									
									mask_out_cube = fltarr(nside, nside, 12)
									temp_mask = mask
									get_all_faces, temp_mask, temp_mask_cube
									
									for num_face = 0, 11 do begin
										mrs_mask, temp_mask_cube(*,*, num_face), '2D', 'ortho', nb_scales, mask_out_face
										mask_out_cube(*,*, num_face) = mask_out_face
										;at this point, 1s in mask_out_cube indicate valid pixels
										;and 0s indicate invalid pixels
									endfor
										 
									put_all_faces, mask_out_cube, mask_out
									
								   end			;of the 'sphere' and 'ortho' case 
								   
									
					endcase
				end
				
	endcase

DONE:



return
end 





;---------------------------------------------------------------
;a few utilities ....
;---------------------------------------------------------------

pro put_submatrix, big_matrix, ind_x, ind_y, sub_matrix

;purpose is to put the elements specified in sub matrix 
;at the coordinates specified by all couples of one entry in ind_x and one entry in ind_y

nb_x = (size(ind_x) )(1)
nb_y = (size(ind_y) )(1)

nb_x_bis = (size(sub_matrix) )(1)
nb_y_bis = (size(sub_matrix) )(2)
if (nb_x ne nb_x_bis) or (nb_y ne nb_y_bis) then begin
		print, 'ERROR : incompatible matrix dimensions'
		goto, DONE
endif

for num_x = 0, nb_x -1 do $
	for num_y = 0, nb_y -1 do $
			big_matrix( ind_x(num_x), ind_y(num_y) ) = sub_matrix(num_x, num_y)

DONE:

end





