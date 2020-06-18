;-------------------------------------------------------------------------------
;+
; NAME:
;		mrs_fastica.pro
;	
; PURPOSE:
;		Apply the ICA method FASTICA on data in different settings : 
;			- the multichannel data may consist of either 1D time series, 2D flat images or spherical maps  
;			- a mask can be specified to indicate missing or invalid pixels. 
;		
;		
;
; EXPLANATION:
;   The components to be separated are all assumed to be iid in the specified representation;
;
; CALLING SEQUENCE:
;	mrs_fastica, data, topology, nb_sources, sources, demixingmat, domain = domain, mask = mask, nb_scales=nb_scales
;
; INPUTS:
;   data :  either an m*T array in the 1D case, 
;			or a tx*ty*m array in the flat 2D case,
;			or an array of strings giving the filenames of the m spherical data maps in Healpix NESTED format.
;
;   topology:   string specifying the topology of the maps in the multichannel data to be processed 
;				this should be either '1D' or '2D' or 'Sphere'
;				This is clearly redundant information but makes things simpler.
;				Any incoherence between the specified 'topology' and the structure of the input data leads to an error.
;
;
;   nb_sources: number of independent sources to be recovered from the data.
;
;
; OPTIONAL INPUTS:
;	domain :	string specifying the representation in which the statistics should be computed
;				DEFAULT : domain = 'initial'
;
;	mask :  either a length T array in the 1D case, 
;			or a tx*ty array in the flat 2D case,
;			or a string giving the filename of a spherical map in Healpix format
;			The specified mask should be the same size as one of the data maps.
;			mask is an array of 0 and 1 where 0 indicates an invalid data sample, and 1 indicates a valid data sample.
;			IF A MASK IS SPECIFIED, THE DATA HAS TO BE MULTIPLIED BY THE MASK PRIOR TO CALLING THE MRS_FASTICA ROUTINE.
; 
;   nb_scales : number of scales in the wavelet transform
;				default is nb_scales = 4 although there is no verification that this is a valid number of scales
;   
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
;   sources :   either an nb_sources*T array in the 1D case, 
;				or a tx*ty*nb_sources array in the flat 2D case,
;				or an array of strings giving the filenames of the nb_sources spherical reconstructed source maps in Healpix NESTED format.
;   demixingmat : an nb_sources * m matrix used to estimate the source processes somewhat according to 'sources = demixingmat#data'.
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
;   mrs_fastica, data, topology, nb_sources, sources, demixingmat, domain = domain, mask = mask, nb_scales=nb_scales
;   
; MODIFICATION HISTORY:
;	Jerome Bobin and Yassir Moudden , 2005
;-------------------------------------------------------------------------------   
;	

pro mrs_fastica, data, topology, nb_sources, sources, demixingmat, domain = domain, mask = mask, nb_scales=nb_scales


COMMON MR1ENV


	;-----------------------------------------------------------------------
	;checking for proper number of i/o parameters and there coherence
	
	if N_PARAMS() LT 5 $
	then begin 
        print, 'CALLING SEQUENCE: mrs_fastica, data, topology, nb_sources, sources, demixingmat, domain = domain, mask = mask, nb_scales=nb_scales'
        goto, DONE
	end
	

	if keyword_set(domain) EQ 0 then domain = 'initial'		;setting default value of domain
	
	
	case topology of
		'1D'		: begin
						size_data = size(data) 
						if (size_data(0) NE 2)$
						then begin
							print, 'error : specified data array incompatible with specified topology'
							goto, DONE
						endif
						
						m = size_data(1)
						T = size_data(2)
						if nb_sources gt m $
						then begin
							print, 'error : specified number of sources is larger than number of channels'
							print, 'assuming the number of sources is equal to the number of channels'
							nb_sources = m	
						endif
						
					  end	

		'2D'		: begin
						size_data = size(data) 
						if (size_data(0) NE 3)$
						then begin
							print, 'error : specified data array incompatible with specified topology'
							goto, DONE
						end
						m = size_data(3)
						tx = double(size_data(1))		;data maps have size tx*ty
						ty = double(size_data(2))
						
						if nb_sources gt m $
						then begin
							print, 'error : specified number of sources is larger than number of channels'
							print, 'assuming the number of sources is equal to the number of channels'
							nb_sources = m	
						endif

					  end	
					
		'Sphere'	: begin
						size_data = size(data)
						size_size_data = size_data(0)+3
						if (size_data(0) NE 1) or (size_data(size_size_data -2 ) NE 7)  $		;is data an array of strings??
						then begin
							print, 'error : specified data array incompatible with specified topology'
							goto, DONE
						end
						m = size_data(1)
						
						if nb_sources gt m $
						then begin
							print, 'error : specified number of sources is larger than number of channels'
							print, 'assuming the number of sources is equal to the number of channels'
							nb_sources = m	
						endif
					  end
					  
		else		: begin
						print, 'error : specified topology is incorrect'
						goto, DONE
					  end	
					  
	endcase

	;------------------------------------------------------------------------



	;----------------------------------------------------------------------------------------
	;separation in the specified topology and domain
	;----------------------------------------------------------------------------------------
	
	case domain of 
	
	'initial' : begin
					case topology of
					
					'1D'		: begin
									if keyword_set(mask)$
									then data_ok = data(*, where( mask eq 1) )$
									else data_ok = data
	
									fastica, data_ok, nb_sources, A_est, demixingmat, sources
									sources = demixingmat # data 
								  end

					'2D'		: begin
									data_temp = transpose( reform(data, tx*ty, m ) )
									if keyword_set(mask)$
									then data_ok = data_temp(*, where( transpose(reform(mask, tx*ty,1))  eq 1) )$
									else data_ok = data_temp
	
									fastica, data_ok, nb_sources, A_est, demixingmat, sources
									sources = demixingmat # data_temp
									sources = reform(transpose(sources), tx, ty, nb_sources ) 
								  end
							
					'Sphere'	:  begin
									sources = strarr(nb_sources)
									source_map = 'sources_fastica_sph_'
									
									read_fits_map, data(0), map
									T = (size(map))(1)
									data_temp = dblarr(m, T)
									data_temp(0,*) = map
									for num_map = 1,m-1 do begin 
										read_fits_map, data(num_map), map
										data_temp(num_map,*) = map
									endfor
									 
									if keyword_set(mask) then begin
										read_fits_map, mask, mask_map
										data_ok = data_temp(*, where( mask_map eq 1) )
									endif else data_ok = data_temp
									
									fastica, data_ok, nb_sources, A_est, demixingmat, sources_temp
									sources_temp = demixingmat # data_temp 
									
									for num_source = 0, nb_sources -1 do begin
											sources(num_source) = strcompress(source_map + string(num_source) +'.fits', /remove_all )
											write_fits_map, sources(num_source),sources_temp(num_source, *), /nested
									endfor
								end		

					endcase
				end
				
	'wavelet' : begin					
					case topology of
					
					'1D'		:   begin
					
									if keyword_set(nb_scales) EQ 0 then nb_scales =4	;setting default number of scales
																						;no verification that this is a valid number of scales
									
									if keyword_set(mr1ok) then begin					;if package mre is available

										wt_options = strcompress('-t15 -L ' + '-n'+ string(nb_scales)  )		;using an L2 normalized orthogonal wavelet transform 
								
										data_wave = data
										for num_map = 0, m-1 do begin
											sig_in = data(num_map, *)
											mr1d_trans, sig_in ,  sig_out , Opt = wt_options		;using an orthogonal 1D wavelet transform
											data_wave(num_map, *) = sig_out.coef - mean(sig_out.coef)
										endfor
									
									endif else begin									;else, use procedures in the mrs package
									
										data_wave = data
										for num_map = 0, m-1 do begin
											sig_in = data(num_map, *)
											sig_out = bwt01_lift(sig_in, nb_scales-1)
											data_wave(num_map, *) = sig_out - mean(sig_out)
										endfor
										
									endelse
									
									
									if keyword_set(mask)$
									then begin 
										mrs_mask, mask, topology, 'ortho', nb_scales, mask_out
										data_ok = data_wave(*, where( mask_out eq 1) )
									endif else data_ok = data_wave
	
									fastica, data_ok, nb_sources, A_est, demixingmat, sources
																		
									;the sources may be reconstructed directly in the initial representation ...
									if 1 then begin 
									sources = demixingmat # data   
									endif
									
									; .. or in the wavelet representation as one may want to do some non linear thresholding before reconstruction
									;in which case, this should be done here before reconstruction using mr_recons or using the inverse of bwt01_lift
									;TO BE COMPLETED
									if 0 then begin									
									sources = demixingmat # data_wave
									for num_source = 0, nb_sources-1 do begin
										;begin thresholding to be inserted here ...
										;end thresholding
										if keyword_set(mr1ok) then sig_out.coef = sources(num_source,*)	$	;uses the structure sig_out that came as an output of the transform step
										else sig_out = sources(num_source,*)
										if keyword_set(mr1ok) then mr1d_recons, sig_out, rec_sig $
										else rec_sig = BWT01_INVERSE( sig_out, nb_scales-1)
										sources(num_source,*) = rec_sig
									endfor
									endif
									
									end		;of the 1D case

	
									


					'2D'		:   begin
									
									if keyword_set(nb_scales) EQ 0 then nb_scales =4	;setting default number of scales
																						;no verification that this is a valid number of scales


									if keyword_set(mr1ok) then begin					;if package mre is available

										MRFILE = 'xxfastica.mr'							;necessary if the procedure mr_Recons is to be used.
																					;otherwise, this is useless.
										data_temp = transpose( reform(data, tx*ty, m ) )
										data_wave = data_temp
										wt_options = strcompress('-t14 -L ' + '-n'+ string(nb_scales) )		;using an L2 normalized orthogonal wavelet transform 
										for num_map = 0, m-1 do begin
											map_in = data(*,*,num_map)
											map_out = map_in
											mr_transform, map_in, map_out, opt=wt_options, MR_File_Name=MRFILE
											data_wave(num_map, *) = map_out - mean(map_out)
										endfor
									
									endif else begin									;else, use procedures in the mrs package
									
										data_temp = transpose( reform(data, tx*ty, m ) )
										data_wave = data_temp
										
										for num_map = 0, m-1 do begin
											map_in = data(*,*,num_map)
											map_out = map_in
											map_out = bwt01_lift(map_in, nb_scales-1)
											data_wave(num_map, *) = map_out - mean(map_out)
										endfor
																											
									endelse
									
									
									if keyword_set(mask)$
									then begin 
										mrs_mask, mask, topology, 'ortho', nb_scales, mask_out
										data_ok = data_wave(*, where( transpose(reform(mask_out, tx*ty,1))  eq 1) )
									endif else data_ok = data_wave							
										
									fastica, data_ok, nb_sources, A_est, demixingmat, sources
																		
									;the sources may be reconstructed directly in the initial representation ...
									if 1 then begin 
									sources = demixingmat # data_temp   
									sources = reform(transpose(sources), tx, ty, nb_sources ) 
									endif
									
									; .. or in the wavelet representation as one may want to do some non linear thresholding before reconstruction
									;in which case, this should be done here before reconstruction using mr_reconsor using the inverse of bwt01_lift
									;TO BE COMPLETED
									if 0 then begin									
									map_out = readfits(MRFILE,Head_heur)			;simply to recover the header of the mrfile .mr, to be used in the inverse wavelet transform
									sources = demixingmat # data_wave
									sources = reform(transpose(sources), tx, ty, nb_sources ) 
									for num_source = 0, nb_sources-1 do begin
										map_in = sources(*,*,num_source)
										;begin thresholding to be inserted here ...
										;end thresholding
										map_out = map_in
										if keyword_set(mr1ok) then mr_recons, map_in, map_out, Header=Head_heur $
										else map_out = BWT01_INVERSE( map_in, nb_scales-1)
										sources(*,*,num_source) = map_out 
									end
									endif
									
									end		;of the 2D case
					

					'Sphere'	: 	begin
									sources = strarr(nb_sources)
									source_map = 'sources_fastica_sph_'
									if keyword_set(nb_scales) EQ 0 then nb_scales =4	;setting default number of scales
																						;no verification that this is a valid number of scales

									read_fits_map, data(0), map
									T = double( (size(map))(1) )
									nside = npix2nside(T)
									data_temp = dblarr(m, T)
									data_temp(0,*) = map
									for num_map = 1,m-1 do begin 
										read_fits_map, data(num_map), map
										data_temp(num_map,*) = map
									endfor
									
									data_wave = data_temp
									;wt_options = strcompress(' -n'+ string(nb_scales) )		;using an L2 normalized orthogonal wavelet transform 
									wt_options = strcompress(' ' )
									for num_map = 0, m-1 do begin
										map_in = reform( data_temp(num_map,*) )
										mrs_owttrans, map_in, map_out, NbrScale=nb_scales, opt=wt_options
										data_wave(num_map, *) = map_out.coef - mean(map_out.coef)		;map_out.coef is an nside*nside*12 array
									endfor
									
								
									if keyword_set(mask)$
									then begin 
										read_fits_map, mask, mask_map
										mrs_mask, mask_map, topology, 'ortho', nb_scales, mask_out
										data_ok = data_wave(*, where( mask_out eq 1) )
									endif else data_ok = data_wave							

									fastica, data_ok, nb_sources, A_est, demixingmat, sources_temp
																											
									;the sources may be reconstructed directly in the initial representation ...
									if 1 then begin 
									sources_temp = demixingmat # data_temp 
									for num_source = 0, nb_sources -1 do begin
										sources(num_source) = strcompress(source_map + string(num_source) +'.fits', /remove_all )
										write_fits_map, sources(num_source), sources_temp(num_source, *), /nested
									endfor
									endif

									; .. or in the wavelet representation as one may want to do some non linear thresholding before reconstruction
									;in which case, this should be done here before reconstruction using mr_recons or the inverse of bwt01_lift
									;TO BE COMPLETED
									if 0 then begin									
									sources_temp = demixingmat # data_wave
									for num_source = 0, nb_sources-1 do begin
										map_out.coef(*)  = sources_temp(num_source, *)			;map_out.coef is an nside*nside*12 idl array
										;begin thresholding to be inserted here ...
										;end thresholding
										mrs_owtrec, map_out, rec_map
										sources_temp(num_source, *) = rec_map 
										sources(num_source) = strcompress(source_map + string(num_source) +'.fits', /remove_all )
										write_fits_map, sources(num_source), sources_temp(num_source, *), /nested
									endfor
									endif
									
									end		;of the SPHERE case
																			
					endcase
				end
				
	endcase
	
	;------------------------------------------------------------

DONE : 

return
end 
