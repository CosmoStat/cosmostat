;+
; NAME:
;        WMAP9_DATA_routines
;
; PURPOSE:
;   Routines to launch SPSR on WMAP 9 year data
;
; CALLING:
;     PTFreeMap=WMAP9_PTSRCDATA_launch_SPSR(band,ctlg,beam_dir,map_dir,galmask_dir,nsigmafit=nsigmafit,Niter=Niter, timer=timer,verbose=verbose,INIT_SPSR= INIT_SPSR, RES_SPSR= RES_SPSR)

; INPUTS:
;     band -- WMAP frequency band [0-4]
;     ctlg -- Catalog of point sources, array containing colatitude and longitude
;     beam_dir  -- Where the beams are located
;     map_dir  -- Where the maps are located
;     galmask_dir  -- Where the galactic mask is located
;   
; OUTPUTS:
;     PTFreeMap -- Healpix map after point source removal
;
; KEYWORDS:
;     	nsigmafit		: Size of the fitting region in number of sigmas  [default: 3.5]
;     	Niter			: Number of iteration in the algorithm [default: 1]
;     	timer			: Print time spent in the various steps of the algorithm [default: 0]	
;     	verbose			: Verbosity [default: 0]
;	INIT_SPSR		: Structure containing all variables setting the initial state of the algorithm.  [default: all variables in the initial state set to 0]
;	RES_SPSR		: Structure containing all variables giving the final state of the algorithm. [default:no output]
;
; EXTERNAL CALLS:
;       create_all_axisymmetric_2dbeam_sphere, mrs_sparse_pointsource_removal 
;
; EXAMPLE:
;	Add 100 iterations of SPSR for WMAP channel K on a previous result PREV_SPSR to get NEXT_SPSR structure and PTSub map
;	PTSub= WMAP9_PTSRCDATA_launch_SPSR(0,ctlg,beam_dir,map_dir,galmask_dir,  Niter=100,INIT_SPSR= PREV_SPSR, RES_SPSR= NEXT_SPSR)

; HISTORY:
;	Written: Florent Sureau, September 2013
;--------------------------------------------------------------------------------------------------------

;-----------------------------------------------------------------------------------------------------
;Legendre transform of the Radial transfer function
function WMAP9_PTSRCDATA_create_axisymmetric_2dbeam_sphere,transfer_function,radial_distance
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	lmax=(size(transfer_function,/dim))[0]-1
	recfact=dblarr(lmax+1,3)
	den=dindgen(lmax+1)+1d
	recfact[*,0]=(2d * dindgen(lmax+1)+1d)/den
	recfact[*,1]=indgen(lmax+1)/den
	recfact[*,2]=(2d * dindgen(lmax+1)+1d)
	
	val_axi= radial_distance[*]*0.
	nd=(size(radial_distance[*],/dim))[0]
	RD=[double(radial_distance[*]),min(radial_distance[*]),max(radial_distance[*])] ;the last two are used for associated legendre polynomial integration (but not useful for spherical harmonics)


	wtrans= transfer_function * ((2d * dindgen(lmax+1)+1.)/(4d * !dpi)) ;;Legendre poly used in SPh Harm. are normalized according to Ll(1)=(2l+1)/(4d * !dpi) whereas here we assume Ll(1)=1
	;init_recursion
	Pkm1=dblarr(nd+2)+1d
	Pk= RD[*]
	val_axi= wtrans[0]*Pkm1+ wtrans[1]*Pk
	int_val=(RD[nd+1]-RD[nd])*wtrans[0] ;integration for l=0

	;recursion
	for k=2,lmax do begin
		Pkp1= RD[*] * Pk * recfact[k-1,0] - Pkm1  * recfact[k-1,1]
		val_axi= val_axi + wtrans[k]*Pkp1
		int_val = int_val +wtrans[k-1] * (( Pkp1[nd+1] - Pkp1[nd] ) - (Pkm1[nd+1] - Pkm1[nd]))/recfact[k-1,2] 
		Pkm1=Pk
		Pk=Pkp1
	endfor
	Pkp1= RD[*] * Pk * recfact[lmax,0] - Pkm1  * recfact[lmax,1]
	int_val = int_val +wtrans[lmax] * (( Pkp1[nd+1] - Pkp1[nd] ) - (Pkm1[nd+1] - Pkm1[nd]))/recfact[lmax,2] 
	val_axi=reform(val_axi[0:nd-1],size(radial_distance,/dim))
	return,val_axi
end

;-----------------------------------------------------------------------------------------------------
;Get Experimental Beams
function WMAP9_PTSRCDATA_get_exp_beam,da_index, beam_dir= beam_dir
	fname=["wmap_symm_beam_profile_K1_9yr_v5.txt","wmap_symm_beam_profile_Ka1_9yr_v5.txt","wmap_symm_beam_profile_Q1_9yr_v5.txt","wmap_symm_beam_profile_Q2_9yr_v5.txt", "wmap_symm_beam_profile_V1_9yr_v5.txt","wmap_symm_beam_profile_V2_9yr_v5.txt","wmap_symm_beam_profile_W1_9yr_v5.txt","wmap_symm_beam_profile_W2_9yr_v5.txt","wmap_symm_beam_profile_W3_9yr_v5.txt", "wmap_symm_beam_profile_W4_9yr_v5.txt"]
	ascii_tmpl={VERSION:1.0, DATASTART:10L, DELIMITER:32B, MISSINGVALUE:!Values.F_NAN, COMMENTSYMBOL:'',FIELDCOUNT:2L, FIELDTYPES:[4L,4L], FIELDNAMES:['field1','field2'], FIELDLOCATIONS:[5L,18L],FIELDGROUPS:[0L,1L]}
	bb=read_ascii(beam_dir+'/'+fname[da_index],templ=ascii_tmpl)
	beam_struc={x:bb.field1,y:bb.field2}
	return,beam_struc
end


;-----------------------------------------------------------------------------------------------------
;Get Transfer Function
function WMAP9_PTSRCDATA_get_transfer_function,da_index, beam_dir= beam_dir
	fname=["wmap_ampl_bl_K1_9yr_v5.txt","wmap_ampl_bl_Ka1_9yr_v5.txt","wmap_ampl_bl_Q1_9yr_v5.txt","wmap_ampl_bl_Q2_9yr_v5.txt","wmap_ampl_bl_V1_9yr_v5.txt","wmap_ampl_bl_V2_9yr_v5.txt", "wmap_ampl_bl_W1_9yr_v5.txt","wmap_ampl_bl_W2_9yr_v5.txt","wmap_ampl_bl_W3_9yr_v5.txt","wmap_ampl_bl_W4_9yr_v5.txt"]
	ascii_tmpl={VERSION:1.0, DATASTART:8L, DELIMITER:32B, MISSINGVALUE:!Values.F_NAN, COMMENTSYMBOL:'',FIELDCOUNT:2L, FIELDTYPES:[3L,4L], FIELDNAMES:['field1','field2'], FIELDLOCATIONS:[5L,11L],FIELDGROUPS:[0L,1L]}
	bb=read_ascii(beam_dir +'/'+fname[da_index],templ=ascii_tmpl)
	beam_struc={l:bb.field1,f:bb.field2}
	return,beam_struc
end


;-----------------------------------------------------------------------------------------------------
;Get Radial profile from transfer function
function WMAP9_PTSRCDATA_get_beam_from_transfer_function,da_index,normalize=normalize,sig_Eq=sig_Eq,normfact=normfact,sampling=sampling, beam_dir= beam_dir, map_dir= map_dir
	TF=WMAP9_PTSRCDATA_get_transfer_function(da_index, beam_dir= beam_dir)
	EXP_BEAM=WMAP9_PTSRCDATA_get_exp_beam(da_index, beam_dir= beam_dir)
	
	dd_k1= cos(EXP_BEAM.x*!dpi/180d)
	
	TBEAM=WMAP9_PTSRCDATA_create_axisymmetric_2dbeam_sphere(TF.f,dd_k1)
	if keyword_seT(normalize) then TBEAM=TBEAM / TBEAM[0] else TBEAM=TBEAM / total(TBEAM,/double) * total(EXP_BEAM.Y,/double)
	nbeam=(size(TBEAM,/dim))[0]
	
	if (N_ELEMENTS(sig_Eq) gt 0) then begin
		normfact=0d
		for k=0,nbeam-2 do normfact=normfact+abs((EXP_BEAM.x)[k]-(EXP_BEAM.x)[k+1])*TBEAM[k]*abs((EXP_BEAM.x)[k])
		normfact=normfact+abs((EXP_BEAM.x)[nbeam-1]-(EXP_BEAM.x)[nbeam-2])*TBEAM[nbeam-1]*abs((EXP_BEAM.x)[nbeam-1]) ;Use integrated PSF to compute equivalent sigma
		sigma_est=sqrt(normfact) 
		normfact=sqrt(2d * !dpi * normfact) ;Beam size
		nels=fix(sig_Eq * sigma_est/abs((EXP_BEAM.x)[1]-(EXP_BEAM.x)[0]))
	endif else nels=nbeam
	
	RadialProfile=TBEAM[0:nels-1]
	sampling=(EXP_BEAM.x)[0:nels-1]
	return,RadialProfile
end

;-----------------------------------------------------------------------------------------------------
function WMAP9_PTSRCDATA_get_band_beam_from_transfer_function,band,map_data,stdpix,beam_freq, nsigmafit = nsigmafit,normfact=normfact,sampling=sampling, beam_dir= beam_dir,map_dir=map_dir

	CASE band of 
		0:	begin
				beam_freq=WMAP9_PTSRCDATA_get_beam_from_transfer_function(0,/normalize,sig_eq=nsigmafit,normfact=normfact,sampling=sampling, beam_dir= beam_dir)
				nobs=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_K1_v5.fits'))[*,1]
				h=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_K1_v5.fits')
				map_data=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_K1_v5.fits'))[*,0] * 1e3 ;in muK
				stdpix=sqrt(fxpar(h,"SIGMA0_I")^2./nobs)*1e3 ;in muK
			end
		1:	begin
				beam_freq =WMAP9_PTSRCDATA_get_beam_from_transfer_function(1,/normalize,sig_eq=nsigmafit,normfact=normfact,sampling=sampling, beam_dir= beam_dir) 
				nobs=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_Ka1_v5.fits'))[*,1]  
				h=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_Ka1_v5.fits') 
				map_data=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_Ka1_v5.fits'))[*,0] * 1e3 ;in muK
				stdpix=sqrt(fxpar(h,"SIGMA0_I")^2./nobs)*1e3 
			end
		2: 	begin
				beam_da1=WMAP9_PTSRCDATA_get_beam_from_transfer_function(2,/normalize,sig_eq=nsigmafit,normfact=normfact1,sampling=xsampling1, beam_dir= beam_dir) 
				beam_da2=WMAP9_PTSRCDATA_get_beam_from_transfer_function(3,/normalize,sig_eq=nsigmafit,normfact=normfact2,sampling=xsampling2, beam_dir= beam_dir) 
				sampling=0.5*(xsampling1+ xsampling2) 
				beam_freq = 0.5*(beam_da1+ beam_da2) 
				map_data1=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_Q1_v5.fits'))[*,0] * 1e3  ;in muK
				nobs1=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_Q1_v5.fits'))[*,1] 
				h1=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_Q1_v5.fits') 
				map_data2=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_Q2_v5.fits'))[*,0] * 1e3  ;in muK
				nobs2=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_Q2_v5.fits'))[*,1] 
				h2=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_Q2_v5.fits')
				map_data=0.5*(map_data1+ map_data2) 
				stdpix=0.5*sqrt(fxpar(h1,"SIGMA0_I")^2./nobs1+fxpar(h2,"SIGMA0_I")^2./nobs2)*1e3 
			end
		3:	begin
				beam_da1=WMAP9_PTSRCDATA_get_beam_from_transfer_function(4,/normalize,sig_eq=nsigmafit,normfact=normfact1,sampling=xsampling1, beam_dir= beam_dir) 
				beam_da2=WMAP9_PTSRCDATA_get_beam_from_transfer_function(5,/normalize,sig_eq=nsigmafit,normfact=normfact2,sampling=xsampling2, beam_dir= beam_dir) 
				beam_freq = 0.5*(beam_da1+ beam_da2) 
				sampling=0.5*(xsampling1+ xsampling2) 
				map_data1=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_V1_v5.fits'))[*,0] * 1e3  ;in muK
				nobs1=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_V1_v5.fits'))[*,1]  
				h1=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_V1_v5.fits') 
				map_data2=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_V2_v5.fits'))[*,0] * 1e3  ;in muK
				nobs2=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_V2_v5.fits'))[*,1]  
				h2=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_V2_v5.fits') 
				map_data=0.5*(map_data1+ map_data2) 
				stdpix=0.5*sqrt(fxpar(h1,"SIGMA0_I")^2./nobs1+fxpar(h2,"SIGMA0_I")^2./nobs2)*1e3
	
			end
		4:	begin
				beam_da1=WMAP9_PTSRCDATA_get_beam_from_transfer_function(6,/normalize,sig_eq=nsigmafit,normfact=normfact1,sampling=xsampling1, beam_dir= beam_dir) 
				beam_da2=WMAP9_PTSRCDATA_get_beam_from_transfer_function(7,/normalize,sig_eq=nsigmafit,normfact=normfact2,sampling=xsampling2, beam_dir= beam_dir) 
				beam_da3=WMAP9_PTSRCDATA_get_beam_from_transfer_function(8,/normalize,sig_eq=nsigmafit,normfact=normfact1,sampling=xsampling3, beam_dir= beam_dir) 
				beam_da4=WMAP9_PTSRCDATA_get_beam_from_transfer_function(9,/normalize,sig_eq=nsigmafit,normfact=normfact2,sampling=xsampling4, beam_dir= beam_dir) 
				beam_freq = 0.25*(beam_da1+ beam_da2+beam_da3+ beam_da4) 
				sampling=0.25*(xsampling1+ xsampling2+ xsampling3+ xsampling4) 
				map_data1=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W1_v5.fits'))[*,0] * 1e3 ;in muK
				nobs1=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W1_v5.fits'))[*,1]  
				h1=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_W1_v5.fits') 
				map_data2=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W2_v5.fits'))[*,0] * 1e3  ;in muK
				nobs2=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W2_v5.fits'))[*,1]  
				h2=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_W2_v5.fits') 
				map_data3=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W3_v5.fits'))[*,0] * 1e3  ;in muK
				nobs3=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W3_v5.fits'))[*,1]  
				h3=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_W3_v5.fits') 
				map_data4=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W4_v5.fits'))[*,0] * 1e3 ;in muK
				nobs4=(rims(map_dir +'/wmap_deconv_imap_r9_9yr_W4_v5.fits'))[*,1]  
				h4=headfits(map_dir +'/wmap_deconv_imap_r9_9yr_W4_v5.fits') 
				map_data=0.25*(map_data1+ map_data2+map_data3+ map_data4) 
				stdpix=0.25*sqrt(fxpar(h1,"SIGMA0_I")^2./nobs1+fxpar(h2,"SIGMA0_I")^2./nobs2+fxpar(h3,"SIGMA0_I")^2./nobs3+fxpar(h4,"SIGMA0_I")^2./nobs4)*1e3
			end
		else: begin
			print,"ERROR: Band should be in [0-4]"
			error_code=1
			goto,END_WMAP9_BAND_BEAM
		end
	end
	error_code=0
	
	END_WMAP9_BAND_BEAM:
		return,error_code
end

;-----------------------------------------------------------------------------------------------------
function dist_sphere,colat1,long1,colat2,long2
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS
	return,acos(cos(colat1)*cos(colat2)+sin(colat1)*sin(colat2)*cos(long1-long2))
end

;-----------------------------------------------------------------------------------------------------
function WMAP9_PTSRCDATA_get_ctlg,ctlg_dir,excl_mask=excl_mask,incl_mask=incl_mask,verbose=verbose

	;We are only interestd in field glon and glat
	ascii_tmpl={VERSION:1.0, DATASTART:56L, DELIMITER:0B, MISSINGVALUE:!Values.F_NAN, COMMENTSYMBOL:'',FIELDCOUNT:3L, FIELDTYPES:[4L,4L,0L], FIELDNAMES:['glon','glat','field3'], FIELDLOCATIONS:[17L,24L,31L],FIELDGROUPS:[0L,1L,2L]}
	ctlg=read_ascii(ctlg_dir+"/wmap_ptsrc_catalog_9yr_v5.txt",templ= ascii_tmpl)
	npt1=(size(ctlg.glon,/dim))[0]
	ascii_tmpl2={VERSION:1.0, DATASTART:39L, DELIMITER:0B, MISSINGVALUE:!Values.F_NAN, COMMENTSYMBOL:'',FIELDCOUNT:3L, FIELDTYPES:[4L,4L,0L], FIELDNAMES:['glon','glat','field3'], FIELDLOCATIONS:[17L,24L,31L],FIELDGROUPS:[0L,1L,2L]}
	ctlg2=read_ascii(ctlg_dir+"/wmap_ptsrc_catalog_cmb_free_9yr_v5.txt",templ= ascii_tmpl2)
	npt2=(size(ctlg2.glon,/dim))[0]
	d2r=!dpi/180d
	
	;Merge catalogs
	;Attempt to get a common catalogue
	ctlg_all=dblarr(npt1+npt2,2)
	ctlg_all[0:npt1-1,0]=(90d - ctlg.glat)*d2r
	ctlg_all[0:npt1-1,1]= ctlg.glon*d2r

	nside=1024ul
	mask_map=dblarr(nside* nside*12ul,3)
	incl_mask =dblarr(nside* nside*12ul)
	mask_map[*,2]=-1

	;Check we do not have two point sources in the same pixel
	k=0ul
	while (k lt npt1) do begin 
		ang2pix_nest, nside, ctlg_all[k,0], ctlg_all[k,1], ipnest 
		if keyword_set(verbose) then if(mask_map[ipnest,2] gt 0) then print,"For point source ",k," Already point ",mask_map[ipnest,2] 
		mask_map[ipnest,0]= ctlg_all[k,0] 
		mask_map[ipnest,1]= ctlg_all[k,1] 
		mask_map[ipnest,2]=k  
		incl_mask[ipnest]=1
		k=k+1ul 
	endwhile
	if keyword_set(verbose) then print,"First Catalog Merged: NPT=",total(incl_mask,/double)

	;Check we do not have a source in the second catalog that is not 
	excl_mask =dblarr(nside* nside*12ul)
	radius= 0.22 * d2r ; WMAP W band FWHM as separation distance
	k=0ul
	off=npt1
	while (k lt npt2) do begin 
		if keyword_set(verbose) then print,"Process ",k," out of ",npt2
		colat_cur=(90d - (ctlg2.glat)[k])*d2r
		lon_cur=(ctlg2.glon)[k]*d2r
		ang2vec, colat_cur, lon_cur,vector
		ang2pix_nest, nside, colat_cur, lon_cur, ipnest
		query_disc , 1024,vector, radius, Listpix, Npix,/nested
		lst=where(mask_map[Listpix,2] gt 0,cnt)
		if(cnt gt 0) then begin ;We have a source close to the current one
			pix_comp=Listpix[lst]
			same_pt=0
			for k2=0ul, cnt-1 do begin
				dd=dist_sphere(mask_map[pix_comp[k2],0],mask_map[pix_comp[k2],1], colat_cur, lon_cur) 
				if(dd lt radius) then begin
					same_pt=1
					if keyword_set(verbose) then print,"Point ",k, " vs Point ", mask_map[pix_comp[k2],2], " (", colat_cur, lon_cur," vs ",mask_map[pix_comp[k2],0],mask_map[pix_comp[k2],1]," distance=",dd,")"
					break
				endif
			endfor
			if(same_pt eq 0) then begin ;Distance was not significant 
				ctlg_all[off,0]= colat_cur
				ctlg_all[off,1]= lon_cur
				incl_mask[ipnest]= 3
				off=off+1ul
			endif else begin ;Exclude point
				excl_mask[listPix]=1
				excl_mask[ipnest]=-10
			endelse
		endif else begin
			ctlg_all[off,0]=colat_cur
			ctlg_all[off,1]=lon_cur
			incl_mask[ipnest]= 2
			off=off+1ul
		endelse
		k=k+1ul
	endwhile

	ctlg_all= ctlg_all[0:off-1,*]
	return, ctlg_all	
end

function WMAP9_PTSRCDATA_launch_SPSR,band,ctlg,beam_dir,map_dir,galmask_dir,nsigmafit=nsigmafit,Niter=Niter, timer=timer,verbose=verbose,INIT_SPSR= INIT_SPSR, RES_SPSR= RES_SPSR

	if(N_ELEMENTS(nsigmafit) eq 0) then nsigmafit=3.5
	normfact=0
	err_code=WMAP9_PTSRCDATA_get_band_beam_from_transfer_function(band,map_data,stdpix,beam_freq, nsigmafit = nsigmafit,normfact=normfact,sampling=sampling, beam_dir= beam_dir, map_dir= map_dir)
	if(err_code eq 1) then goto, END_WMAP9_DATA

	;Create Beam Structure 
	npt=(size(ctlg,/dim))[0]
	d2r=!dpi/180d 
	theta_array=ctlg[*,0]
	phi_array=ctlg[*,1]
	nside=512ul
	sampling_rad = sampling*d2r
	beam_struc=create_all_axisymmetric_2dbeam_sphere(nside, sampling_rad, beam_freq,theta_array,phi_array,normfact=normfact)
	GalMask =(rims(galmask_dir +'/wmap_point_source_catalog_mask_r9_9yr_v5.fits'))[*,0]
	Tab_Lmax_algo=[600,800,900,1300,1536] ;maximal multipole (depending on the nside parameter, no more than 3nside

	nsidemap=sqrt((size(map_data,/dim))[0]/12d)
	struc_beam_ptsrc=beam_struc
	struc_beam_ptsrc.amp= 1d ;Beam not equal to 1 at 0, rather integrates to 1

	nptsrc=(size(STRUC_BEAM_PTSRC.amp,/dim))[0]

	if(Niter gt 0) then PTSub=mrs_sparse_pointsource_removal( map_data, GalMask, struc_beam_ptsrc, stdpix,  Niter=Niter, lmax= Tab_Lmax_algo[band],timer=timer,verbose=verbose,INIT_SPSR= INIT_SPSR, RES_SPSR= RES_SPSR)

	if(nsigmafit lt 5d) then begin
		;Subtract point sources fitted at 3.5 sigma in a 5 sigma region
		err_code=WMAP9_PTSRCDATA_get_band_beam_from_transfer_function(band,map_data,stdpix,beam_freq_5sig, nsigmafit = 5d,normfact=normfact,sampling=sampling_5sig, beam_dir= beam_dir, map_dir= map_dir)
		beam_struc_5sig=create_all_axisymmetric_2dbeam_sphere(nside, sampling_rad, beam_freq,theta_array,phi_array,normfact=normfact)
		for k=0, nptsrc-1 do beam_struc_5sig.amp[k]= RES_SPSR.xnp1[k]*max(*(beam_struc.PSF[k]))/max(*(beam_struc_5sig.PSF[k]))
		PT2Sub= PTSub*0.
		convol_map_local_2dbeam_sphere, PT2Sub, struc_beam_ptsrc
		PTSub=map_data-PT2Sub
	endif

END_WMAP9_DATA:

	return,PTSub

end
