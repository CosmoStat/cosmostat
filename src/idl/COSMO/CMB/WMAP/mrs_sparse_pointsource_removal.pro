;+
; NAME:
;        mrs_sparse_pointsource_removal
;
; PURPOSE:
;   Remove Estimated Point sources in a map with the SPSR algorithm  
;
; CALLING:
;     PTFreeMap =  mrs_remove_point_sources(Map, GalMask, BeamInfo,  Niter=Niter, lmax=lmax,timer=timer,verbose=verbose, INIT_SPSR = INIT_SPSR, RES_SPSR = RES_SPSR)
;
; INPUTS:
;     Map -- Healpix map to process (nested format)
;     GalMask -- Healpix map containing the mask indicating with 0 values where the extended compact sources need be estimated (nested format)
;     BeamInfo  -- Pointer array of structure containing info on the local beams (generated using create_all_axisymmetric_2dbeam_sphere)
;    
; OUTPUTS:
;     PTFreeMap -- Healpix map after point source removal
;
; KEYWORDS:
;     Niter			: Number of iteration in the algorithm [default: 1]
;     lmax			: Maximal multipole for the band limited diffuse component [default: min([2*nside,4200])]
;	  timer			: Print time spent in the various steps of the algorithm [default: 0]	
;	  verbose		: Verbosity [default: 0]
;	  INIT_SPSR		: Structure containing all variables setting the initial state of the algorithm.  [default: all variables in the initial state set to 0]
;	  RES_SPSR		: Structure containing all variables giving the final state of the algorithm. [default:no output]
;
; EXTERNAL CALLS:
;       query_disc, query_polygon (healpix software)
;       mrs_almtrans, mrs_almrec, mrs_wttrans, mrs_wtrec (iSAP software)
;		tag_exist (IDL ASTROLIB, also in iSAP)
;
; EXAMPLE:
;	Add 100 iterations of SPSR on a previous result PREV_SPSR to get NEXT_SPSR structure and PTSub map
;	PTSub=mrs_sparse_pointsource_removal( Map, GalMask, BeamInfo,std_map,  Niter=100,INIT_SPSR= PREV_SPSR, RES_SPSR= NEXT_SPSR)
;         
; HISTORY:
;	Written: Florent Sureau, September 2013
;--------------------------------------------------------------------------------------------------------


;-----------------------------------------------------------------------------------------------------
function create_ell_region_sphere,nside,sigma_x,sigma_y,theta_c,phi_c,size_radius=size_radius,orient=orient,polygon=polygon,nsigma=nsigma
	;Create structure containing list of pixels and distance to center in an elliptical region around point source in orthographic projection, with size given by nsigma and sigma_x / sigma_y
	;Theta_c and phi_c are in radians and corresponds to colatitude and longitude in galactic coordinates
	;sigma_x and sigma_y are the gaussian standard deviations, with orientation given by the parameter orient (in degree)
	;Size_radius is in degree and corresponds to the initial size of the search, either elliptical region or polygonal region (if keyword polygonn is set)
	;Then this region is refined using nsigma
	
	
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	if not keyword_set(nsigma) then nsigma=3.5
	if not keyword_set(size_radius) then size_rad=max([nsigma*sigma_x,nsigma*sigma_y]) $
	else size_rad=size_radius ;some distance containing most of power for Gaussian
	if not keyword_set(amp) then amp=0. ;will contain the amplitude of point source
	d2r=double(!dpi)/180d
	sigma_xrad=double(sigma_x[0])*d2r ;supposedly measured sigma in projection grid
	sigma_yrad=double(sigma_y[0])*d2r ;supposedly measured sigma in projection grid
	rFWHMx=sigma_xrad*sqrt(8d*alog(2d))
	if (N_ELEMENTS(orient) gt 0) then begin
		ort_ang=double(orient)*d2r
		cang=cos(ort_ang)
		sang=sin(ort_ang)
	endif
	ang2vec , theta_c, phi_c, VectorCent	
	;Simple polygon
	if keyword_set(polygon) then begin
		ctc=cos(double(theta_c))
		stc=sin(double(theta_c))
		cpc=cos(double(phi_c))
		spc=sin(double(phi_c))
		dtheta=[ctc*cpc,ctc*spc,-stc]
		dphi=[-spc,cpc,0.]
		if keyword_set(orient) then begin
			xfact=double(nsigma * sigma_xrad)
			yfact=double(nsigma * sigma_yrad)
			dxr=xfact*cang[0]- yfact*sang[0]
			dyr=xfact*sang[0]+ yfact*cang[0]
			DX1=double(VectorCent)+dxr[0]*dtheta+dyr[0]*dphi
			DX3=double(VectorCent)-dxr[0]*dtheta-dyr[0]*dphi
			dxr=-xfact*cang[0]- yfact*sang[0]
			dyr=-xfact*sang[0]+ yfact*cang[0]
			DX2=double(VectorCent)+dxr[0]*dtheta+dyr[0]*dphi
			DX4=double(VectorCent)-dxr[0]*dtheta-dyr[0]*dphi
		endif else begin
			dxr=double(dtheta[*]*nsigma * sigma_xrad)
			dyr=double(dphi[*]*nsigma * sigma_yrad)
			DX1=double(VectorCent+dxr+dyr)
			DX2=double(VectorCent+dxr-dyr)
			DX3=double(VectorCent-dxr-dyr)
			DX4=double(VectorCent-dxr+dyr)
		endelse		
		query_polygon , Nside, [[DX1[0],DX2[0],DX3[0],DX4[0]],[DX1[1],DX2[1],DX3[1],DX4[1]],[DX1[2],DX2[2],DX3[2],DX4[2]]], List, NpixList,/nested,/inclusive
	endif else query_disc , Nside, VectorCent, size_rad, List,NpixList,/nested,/deg,/inclusive ;or disk region
	
	;Refine search
	List=List[sort(List)]
	dd=dblarr(NpixList)
	dx=dblarr(NpixList)
	dy=dblarr(NpixList)
	sinC=sin(double(theta_c))
	cosC=cos(double(theta_c))
	for n=0,NpixList-1 do begin
		pix2ang_nest, nside, List[n], theta_p, phi_p
		sinP=sin(double(theta_p))
		dx[n]=sinP*cosC*cos(double(phi_p-phi_c))-sinC*cos(double(theta_p))
		dy[n]=sinP*sin(double(phi_p-phi_c))
		if keyword_set(orient) then begin
			temp=dx[n]*cang[0]+dy[n]*sang[0]
			dy[n]=-dx[n]*sang[0]+dy[n]*cang[0]
			dx[n]=temp
		endif
		dd[n]=(dx[n]/sigma_xrad)^2. + (dy[n]/sigma_yrad)^2.
	endfor
	ll=where(dd le double(nsigma*nsigma),NpixList)
	List=List[ll]
	dx=dx[ll]
	dy=dy[ll]
	dd=dd[ll]			
	distpix=sqrt(dx*dx+dy*dy) ;distance in projected plane
	struc_ell_region={List:List, nelem: NpixList, nside: nside,distpix:distpix,normdist: dd}
	return, struc_ell_region
end


;-----------------------------------------------------------------------------------------------------
function create_local_gaussian_2dbeam_sphere,nside,sigma_x,sigma_y,theta_c,phi_c,size_radius=size_radius,amp=amp,orient=orient,polygon=polygon,nsigma=nsigma,nonorm=nonorm
	;Create structure containing an elliptical gaussian beam, described by sigma_x and sigma_y (in degree), centered at theta_c and phi_c (in radians), with angle given by orient (degree)
	;Normalize to 1 by default the structure (for the algorithm)
	
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	if not keyword_set(amp) then amp=0. ;will contain the amplitude of point source
	d2r=double(!dpi)/180d
	area_pix=4d*double(!dpi)/(double(nside)*double(nside)*12d)
	sigma_xrad=double(sigma_x[0])*d2r ;supposedly measured sigma in projection grid
	sigma_yrad=double(sigma_y[0])*d2r ;supposedly measured sigma in projection grid
	normfact= 1d/(2d*double(!dpi)*sigma_xrad*sigma_yrad)*area_pix

	struc_ell_region=create_ell_region_sphere(nside,sigma_x,sigma_y,theta_c,phi_c,size_radius=size_radius,orient=orient,polygon=polygon,nsigma=nsigma)
	List=(struc_ell_region.List)
	NpixList= struc_ell_region.nelem

	PSF=dblarr(NpixList)
	Beam=dblarr(NpixList)
	for n=0,NpixList-1 do begin
		PSF[n]=normfact*exp(double(-0.5* (struc_ell_region.normdist)[n]))
		Beam[n]=PSF[n] ;symmetric
	endfor
	if not keyword_set(nonorm) then abstot=(total(abs(PSF[*]),/double)+1e-8) $
	else abstot=1d
	PSF=PSF/abstot
	Beam=Beam/abstot
	normfact= normfact*abstot
	struc_beam={List:List, PSF: PSF, Beam:Beam, nelem: NpixList, nside: nside, amp:amp,normfact:normfact, abstot: abstot}
	return,struc_beam
end

;-----------------------------------------------------------------------------------------------------
function create_all_gaussian_2dbeam_sphere,nside,FWHMx,FWHMy,theta_array,phi_array,size_radius=size_radius,orient=orient,polygon=polygon,nsigma=nsigma,nonorm=nonorm

	;Create point source structure with elliptical gaussian PSFs: 
	;nside: Healpix parameter nside (Npix=12*Nside^2)
	;FWHMx, FWHMy: full width at half maximum in the x/y direction (in arcminutes) local orthographic projection, with potential rotation angle given by orient 
	;theta (galactic colatitude: theta=0 galactic plane)/ phi (galactic longitude) in radians
	;
	;Optional keyword:
	;size_radius: size of radius for elliptical region (if set, nsigma can define it)
	;nsigma: number of sigma to define the region
	;orient: rotation angle in between the local coordinate system after orthographic projection and the major axis of the elliptical PSF (degree)
	;polygon: region defined using a polygon, not a circular region (useful if FHWMx and FWHMy differs a lot, i.e. very elongated point source region)
	;nonorm : if set, the PSF is not normalized to 1 (PSF normalization required by the algorithm).
	
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	if((size(theta_array))[0] gt 0) then nb=(size(theta_array,/dim))[0] $
	else nb=1
	List=ptrarr(nb,/allocate_heap)
	PSF=ptrarr(nb,/allocate_heap)
	Beam=ptrarr(nb,/allocate_heap)
	nelem=dblarr(nb)
	amp=dblarr(nb)
	normfact=dblarr(nb)
	abstot=dblarr(nb)
	rFWHMx=FWHMx/60d ;FWHMx -> deg
	rFWHMy=FWHMy/60d ;FWHMx -> deg	
	szrfx=(size(rFWHMx,/dim))[0]
	szrfy=(size(rFWHMy,/dim))[0]
	
	if((szrfx ne 0) && (szrfx ne nb)) then print,"Problem with FWHMx size: not 0 nor ",nb, " :",szrfx
	sigma_x=double(rFWHMx)/sqrt(8d*alog(2d)) 
	if((szrfy ne 0) && (szrfy ne nb)) then print,"Problem with FWHMy size: not 0 nor ",nb, " :",szrfy
	sigma_y=double(rFWHMy)/sqrt(8d*alog(2d))	
	
	if not keyword_set(orient) then orient=0 
		
	for k=0,nb-1 do begin
		;print,k+1," out of ",nb
		if((nb gt 1) && (szrfx gt 0)) then sgx=sigma_x[k] else sgx=sigma_x
		if((nb gt 1) && (szrfy gt 0)) then sgy=sigma_y[k] else sgy=sigma_y
		if((nb gt 1) && keyword_set(orient)) then ort=orient[k] else ort=orient
		struc_beam=create_local_gaussian_2dbeam_sphere(nside,sgx,sgy,theta_array[k],phi_array[k],size_radius=size_radius,orient=ort,polygon=polygon,nsigma=nsigma,nonorm=nonorm) 
		nelem[k]=struc_beam.nelem
		*(List[k])=struc_beam.List
		*(PSF[k])=struc_beam.PSF
		*(Beam[k])=struc_beam.Beam
		amp[k]=struc_beam.amp
		normfact[k]=struc_beam.normfact
		abstot[k]=struc_beam.abstot
	endfor
	
	struc_beam={List:List, PSF: PSF, Beam:Beam,nelem:nelem, nside: nside,amp: amp,normfact:normfact,nsigma:nsigma,abstot: abstot}
		
	return,struc_beam
end


;-----------------------------------------------------------------------------------------------------
;Create Axisymmetric beam
function create_local_axisymmetric_2dbeam_sphere,nside,RadialSampling,RadialProfile,theta_c,phi_c,amp=amp,nonorm=nonorm,normfact=normfact
	;Create structure containing an axisymmetric beam, described by RadialSampling (rad) and RadialProfile, centered at theta_c and phi_c (in radians)
	;normalize by default the structure
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	if not keyword_set(amp) then amp=0. ;will contain the amplitude of point source
	d2r=double(!dpi)/180d
	area_pix=4d*double(!dpi)/(double(nside)*double(nside)*12d)
	sprof=(size(RadialProfile,/dim))[0]
	dist_max=max(abs(RadialSampling[*]))*180d / !dpi;size in degrees

	if (N_ELEMENTS(normfact) lt 1) then begin
		normfact=0d
		for k=0,sprof-2 do normfact=normfact+abs(RadialSampling[k]-RadialSampling[k+1])*TBEAM[k]*abs(RadialSampling[k])
		normfact=normfact+abs(RadialSampling[sprof-1]-RadialSampling[sprof-2])*TBEAM[sprof-1]*abs(RadialSampling[sprof-1]) ;in degrees
		normfact=sqrt(2d * !dpi * normfact) ;Beam size
	endif

	struc_ell_region=create_ell_region_sphere(nside,1d,1d,theta_c,phi_c,nsigma=dist_max)
	List=(struc_ell_region.List)
	NpixList= struc_ell_region.nelem

	PSF=dblarr(NpixList)
	sindex=sort(struc_ell_region.distpix)
	new_dist=(struc_ell_region.distpix)[sindex]
	outpsf=spline(RadialSampling,RadialProfile,new_dist,/double)
	PSF[sindex]=outpsf
	lst_low=where(struc_ell_region.distpix lt RadialSampling[0],cnt_low)
	if(cnt_low gt 0) then PSF[lst_low]=RadialProfile[0]
	lst_hi=where(struc_ell_region.distpix gt RadialSampling[sprof-1],cnt_hi)
	if(cnt_hi gt 0) then PSF[lst_hi]=0d
	if not keyword_set(nonorm) then abstot=(total(abs(PSF[*]),/double)+1e-8) $
	else abstot=1d
	PSF=PSF/abstot
	Beam=PSF
	normfact= normfact*abstot
	struc_beam={List:List, PSF: PSF, Beam:Beam, nelem: NpixList, nside: nside, amp:amp,normfact:normfact}
	return,struc_beam
end



;-----------------------------------------------------------------------------------------------------
function create_all_axisymmetric_2dbeam_sphere,nside,RadialSampling,RadialProfile,theta_array,phi_array,nonorm=nonorm,normfact=normfact

	;Create point source structure with elliptical gaussian PSFs: 
	;nside: Healpix parameter nside (Npix=12*Nside^2)
	;theta (galactic colatitude: theta=0 galactic plane)/ phi (galactic longitude) in radians
	;
	;Optional keyword:
	;nsigma: number of sigma to define the region
	;nonorm : if set, the PSF is not normalized to 1 (normalization required by the algorithm).
	
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	if((size(theta_array))[0] gt 0) then nb=(size(theta_array,/dim))[0] $
	else nb=1
	List=ptrarr(nb,/allocate_heap)
	PSF=ptrarr(nb,/allocate_heap)
	Beam=ptrarr(nb,/allocate_heap)
	nelem=dblarr(nb)
	amp=dblarr(nb)
	normfact2=dblarr(nb)
			
	for k=0,nb-1 do begin
		;print,k+1," out of ",nb
		struc_beam=create_local_axisymmetric_2dbeam_sphere(nside,RadialSampling,RadialProfile,theta_array[k],phi_array[k],nonorm=nonorm,normfact=normfact) 		
		nelem[k]=struc_beam.nelem
		*(List[k])=struc_beam.List
		*(PSF[k])=struc_beam.PSF
		*(Beam[k])=struc_beam.Beam
		amp[k]=struc_beam.amp
		normfact2[k]=struc_beam.normfact
	endfor
	max_dist=max(abs(RadialSampling[*]))
	struc_beam={List:List, PSF: PSF, Beam:Beam,nelem:nelem, nside: nside,amp: amp,normfact:normfact2,nsigma:max_dist}
		
	return,struc_beam
end



;-----------------------------------------------------------------------------------------------------
pro convol_map_local_2dbeam_sphere,map_in,struc_beam
	;add to a map convolved point source, with beam structure struc_beam containing the amplitude of the sources and the psfs
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS
		
	nsidemap=struc_beam.nside
	if(nsidemap ne struc_beam.nside) then begin
		print,"Convol_LG2D: No match between map nside: ",nsidemap," and nside for beam: ",struc_beam.nside
		goto,END_CMG2S
	endif
	map_in=map_in*0d	
	nb=(size(struc_beam.nelem,/dim))[0]
	for kb=0,nb-1 do begin
		nelem=struc_beam.nelem[kb]
		for i=0,nelem-1 do begin
			map_in[(*(struc_beam.List[kb]))[i]]=map_in[(*(struc_beam.List[kb]))[i]]+(*(struc_beam.PSF[kb]))[i]*struc_beam.amp[kb]
		endfor
	endfor

	END_CMG2S:
	
end

;-----------------------------------------------------------------------------------------------------
pro backconvol_map_local_2dbeam_sphere,map_in,struc_beam
	;get the transpose of the convolution: get the amplitude by back-convolving the map with the beam

	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS
		
	nsidemap=struc_beam.nside
	if(nsidemap ne struc_beam.nside) then begin
		print,"Convol_LG2D: No match between map nside: ",nsidemap," and nside for beam: ",struc_beam.nside
		goto,END_BCMG2S
	endif
	
	nb=(size(struc_beam.nelem,/dim))[0]
	
	for kb=0,nb-1 do begin ;For point source kb
		nelem=struc_beam.nelem[kb]
		amp_value=0d
		for i=0,nelem-1 do begin ;Take everything in the mask
			amp_value=amp_value +map_in[(*(struc_beam.List[kb]))[i]] * (*(struc_beam.Beam[kb]))[i]
		endfor
		struc_beam.amp[kb]=amp_value
	endfor
	
	END_BCMG2S:
	
end

;-----------------------------------------------------------------------------------------------------
pro free_struc_beam,struc2free

	print,ntags(struc2free)
	ptr_free,struc2free.List
	ptr_free,struc2free.PSF
	ptr_free,struc2free.Beam

end


;-----------------------------------------------------------------------------------------------------
;SOFT THRESHOLDING in IMAGE SPACE
pro soft_threshold_im,im,thr,thr_im,fnorm=fnorm
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	ind=where(abs(im[*]) gt Thr,cnt)
	thr_im=im*0.
	if(cnt gt 0) then begin
		k=0l
		while(k lt cnt) do begin
			cur_ind=ind[k]
			thr_im[cur_ind]=im[cur_ind]-sign(im[cur_ind])*Thr
			k=k+1l
		endwhile
	endif
	if (N_ELEMENTS(fnorm) gt 0) then fnorm=total(abs(double(thr_im[*])),/double)

end


;-----------------------------------------------------------------------------------------------------
;SOFT THRESHOLDING in spherical harmonics
pro soft_threshold_ALM,im,Thr,thr_im,lmax=lmax,fnorm=fnorm
	
	;Note if we want to use soft thresholding on the real form of spherical harmonics:
	; For m=0, this is equivalent to thresholding the coefficient itself
	; Beware For m<>0, we need to multiply by the square root of 2 in front of imaginary and real parts.
	; In that case:  alm_temp[k,*]=(ALM.alm[k,*]*sqrt(2)-sign(ALM.alm[k,*])*Thr[k])/sqrt(2) ie we only need to divide the Thr[k] by sqrt(2)
	
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS
	npix=(size(im[*],/dimensions))[0]
	nside=sqrt(npix/12d)
	
	if not keyword_set(Lmax) then Lmax=min([4200,2*nside])
	nalms=(ulong(Lmax)+1ul)*(ulong(Lmax)+2ul)/2ul

	;THRESHOLDING OF ALMS DIRECTLY
	mrs_almtrans, im, ALM, lmax=lmax,/norm

	alm_temp=ALM.alm*0.
	Thr_real=dblarr(nalms)+abs(Thr)
	
	;Threshold the norm of complex spherical harmonics			
	k=0ul
	while(k lt nalms) do begin
		nrm2=sqrt(ALM.alm[k,0]^2d + ALM.alm[k,1]^2d)
		if(nrm2 gt Thr_real[k]) then begin
			alm_temp[k,0]=ALM.alm[k,0]*(1d - Thr_real[k]/nrm2)
			alm_temp[k,1]=ALM.alm[k,1]*(1d - Thr_real[k]/nrm2)
		endif
		k=k+1ul
	endwhile	
	
	if(N_ELEMENTS(fnorm) gt 0) then fnorm=total(abs(double(ALM.alm[*])),/double)
	ALM.alm=alm_temp

	mrs_almrec, ALM, thr_im
end

;-----------------------------------------------------------------------------------------------------
;Get inverse cumulative distribution function for chi2 r.v. with N degrees of freedom, using the incomplete gamma function (or regularized gamma functions, see abramowitz and stegun, 6.5.1->6.5.5) and linear interpolation
function tail_chi2,Nfree,alpha=alpha,prec=prec,verb=verb
;alpha is the desired lower probability (0.05 by default)
;asymptotic properties of Chi2-> N(k,2k)
	if not keyword_set(alpha) then alpha=0.05
	alpham=1.-alpha
	if not keyword_set(prec) then prec =(alpha)/1000.
	max_dicho=Nfree
	if((alpha lt 0)||(alpha gt 1)) then begin
		print,"tail_chi2: incorrect alpha: ", alpha
		propval=-1
		goto,END_tail_chi2
	endif
	scale= double(Nfree)/2.
	
	while(igamma(scale,max_dicho/2.) lt (alpham)) do  max_dicho= max_dicho*2.
	prop_val_cdf=igamma(scale,max_dicho/2.)
	prop_cdf=[0d,prop_val_cdf]
	prop_x= [0d,double(max_dicho)]
	prop_val= max_dicho
	while(abs(prop_val_cdf-alpham) gt prec) do begin
		prop_val =(INTERPOL(prop_x, prop_cdf, alpham))[0]
		prop_val_cdf=igamma(scale,prop_val/2.)	
		if keyword_set(verb) then begin
			print,"prop_val",prop_val
			print,"prop_val_cdf",prop_val_cdf
			print, prop_cdf, prop_x	
		endif
		if(prop_val_cdf gt alpham) then begin
			prop_cdf[1]= prop_val_cdf
			prop_x[1]=prop_val
		endif else begin
			prop_cdf[0]= prop_val_cdf
			prop_x[0]=prop_val
		endelse
	endwhile 

END_tail_chi2:
		return, prop_val
end

;-----------------------------------------------------------------------------------------------------
  
function PD_SPSR,map,Galmask,BeamInfo,std_map,niter=niter,thrcst=thrcst,lmax=lmax,timer=timer,verbose=verbose,save_struc=save_struc, init_ptsrc= init_ptsrc, init_diffuse= init_diffuse, NbrScale_MCA= NbrScale_MCA
;SPSR Algorithm 
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS
	
	;----------------------
	;initialization:
	
	;Map parameters
	if keyword_set(timer) then ref_time=systime(1)
	Npix=(size(map))[1]	
	nside=sqrt(double(Npix)/12d)	
	if not keyword_set(lmax) then lmax=min([2*nside,4200])
	if not keyword_set(thrcst) then thrcst =0.01
	
	;Noise parameters
	if (N_ELEMENTS(std_map) eq Npix) then invstdpix_norm= 1d / std_map else  invstdpix_norm=dblarr(Npix)+1d	;Sigma^{-1/2}	
	max_invstd_factor=max(invstdpix_norm)
	invstdpix_norm= invstdpix_norm /max_invstd_factor ;ensure that we have a norm <=1 for the diagonale std matrix
	alpha=0.05 ;parameter for the noise (1-alpha = percentile of noise)
	noise_thr=sqrt(tail_chi2(Npix,alpha=alpha,prec=1e-5)) /max_invstd_factor

	;Stats on map to set the hyperparameters
	mx_im=max(map[*])
	mrs_almtrans, map*invstdpix_norm, ALM, lmax=lmax, /norm 
	mx_trans=max(abs(ALM.ALM[*])) 

	if keyword_set(verbose) then print,"mx_trans=",mx_trans

	;Algorithm variables
	theta=1d
	;hyperparameters sigma and tau
	hpar=dblarr(2)+0.55
	if not keyword_set(niter) then niter =1 
		
	ptsrc_binmask= map*0.
	nptsrc=(size(BeamInfo.nelem,/dim))[0]
	for k=0,nptsrc-1 do for i=0, BeamInfo.nelem[k]-1 do ptsrc_binmask[(*(BeamInfo.List[k]))[i]]=1.

	;Initialize Extended sources parameters
	if not keyword_set(NbrScale_MCA) then NbrScale_MCA= 4
	lst_mca=where(Galmask eq 0,cnt_mca)
	if(cnt_mca lt 1) then begin
		print,"No Galmask"
		goto,END_PROX
	endif 
	Thr_mca=dblarr(Nbrscale_MCA)
	for k=0, Nbrscale_MCA -1 do Thr_mca[k]= thrcst*mx_trans

	;initialize algo parameters
	xn2=0
	xb2=0
	xnp2=0
	yn2=0
	yn1=0
	im_invnoise=map * hpar[0]* invstdpix_norm
	if keyword_set(save_struc) then begin
		if keyword_set(verbose) then print,"Initialize maps with save structure"
		xn0=save_struc.xnp0 ;diffuse background
		xn1=save_struc.xnp1 ;point source amplitude
		norm_l1_im=total(abs(save_struc.xnp1),/double)
		xb0=save_struc.xnp0+theta*(save_struc.xnp0-save_struc.xn0)
		xb1=save_struc.xnp1+theta*(save_struc.xnp1-save_struc.xn1)
		xn2=save_struc.xnp2 ;extended compact sources
		xb2=save_struc.xnp2+theta*(save_struc.xnp2-save_struc.xn2)
		yn=save_struc.yn ;dual variable
	endif else begin
		if (N_ELEMENTS(init_diffuse) ne Npix) then xn0 = map*0. else xn0 = init_diffuse ;diffuse background
		if (N_ELEMENTS(init_ptsrc) eq nptsrc) then xn1= init_ptsrc else xn1=dblarr(nptsrc)
		xb0=xn0
		xb1=xn1
		xb2=xn2
		yn=0.
	endelse		
	BeamInfo.amp=xn1
	BMnp1=map*0.
	convol_map_local_2dbeam_sphere,BMnp1, BeamInfo

	if keyword_set(verbose) then begin
		mrs_almtrans, xn0, ALM, lmax=lmax, /norm 
		norm_l1_trans=total(abs(sqrt(ALM.ALM[*,0]^2.+ALM.ALM[*,1]^2.)),/double)
		im_l2=norm(map[*]-xn0-BMnp1,l=2,/double)^2.
		print,"l1_norm start=",norm_l1_trans,"l2_norm start=", im_l2
	endif	

	;MAIN FOR LOOP FOR THE ALGORITHM
	for it=0, niter-1 do begin
		if keyword_set(timer) then begin
			start_it_time=systime(1)
			print,"Iteration "+strcompress(string(it),/remove_all)," start: ",start_it_time-ref_time
		endif

		if keyword_set(save_struc) then totit=it+save_struc.nit else totit=it
		if keyword_set(verbose) then print,"Nit=",it,"/",niter," (tot=",totit,") theta=",theta," hpar=",hpar
		
		;----------------------
		;start by computing y^{n+1}
		BeamInfo.amp=xb1
		rb_P=map*0.
		convol_map_local_2dbeam_sphere,rb_P, BeamInfo
			
		r=yn+hpar[0]*(xb0+xb2+rb_P)* invstdpix_norm  ;t^n+ K \bar{xi}^n

		if keyword_set(timer) then begin
			it_time=systime(1)
			print,"Residual: ",it_time-ref_time," (",it_time-start_it_time,")"
		endif
			
		;----------------------
		;Compute  prox_(sigma F*) (r)
		res=r-im_invnoise
		nrm_ball=norm(res,l=2,/double) 
		print,"nrmball=",nrm_ball/hpar[0], "vs ", noise_thr
		if(nrm_ball le (hpar[0]* noise_thr)) then yn=yn*0. else yn=res*(1.-hpar[0]* noise_thr/nrm_ball)								
					
		if keyword_set(verbose) then print,"MAX YN=",max(yn[*])
		rnp0=xn0-hpar[1]*yn*invstdpix_norm  ;gradient step for diffuse component
		backconvol_map_local_2dbeam_sphere,yn, BeamInfo 	
		rb_y= BeamInfo.amp ;gradient step for ptsrc flux 
		rnp2=xn2-hpar[1]*yn*invstdpix_norm* (1d - Galmask)  ;gradient step for Extended sources

		if keyword_set(verbose) then print,"max backconvol",max(BeamInfo.amp[*])
		rnp1=xn1-hpar[1]*rb_y*invstdpix_norm 

		if keyword_set(timer) then begin
			it_time=systime(1)
			print,"Prox Fidelity: ",it_time-ref_time," (",it_time-start_it_time,")"
		endif
			
		;----------------------
		;compute prox for diffuse component:
		;Soft thresholding in transform domain (ALM)
		Thr_trans=thrcst*mx_trans
		if keyword_set(verbose) then print,"ST TRANS with Thr=",Thr_trans," cst=",thrcst
		if keyword_set(verbose) then fnorm=0
		soft_threshold_ALM, rnp0,(Thr_trans*hpar[1]),xnp0,lmax=lmax,fnorm=fnorm
		if keyword_set(verbose) then print,"l1 norm",fnorm
		;----------------------
		;compute prox for ptsrc : Positivity constraint
		xnp1=rnp1
		lst_p=where(xnp1 lt 0d,cnt_p)
		if(cnt_p gt 0) then xnp1[lst_p]=0d
		if keyword_set(verbose) then print,"xnp1",max(abs(xnp1[*]-xn1))
						
		;----------------------
		;compute prox for extended compact sources : Soft thresholding in wavelet domain
		mrs_wttrans,rnp2,rnp2_wtc,lmax=lmax,NbrScale=Nbrscale_MCA
		for k=0, Nbrscale_MCA-2 do begin
			soft_threshold_im, rnp2_wtc.coef[*,k],(Thr_mca[k]*hpar[1]),temp,fnorm=fnorm
			rnp2_wtc.coef[*,k]=temp
		end
		rnp2_wtc.coef[*, Nbrscale_MCA-1]=0
		mrs_wtrec, rnp2_wtc,temp
		xnp2=temp* (1d - Galmask)
		
		if keyword_set(timer) then begin
			it_time=systime(1)
			print,"Prox constraint: ",it_time-ref_time," (",it_time-start_it_time,")"
		endif
			
		xb1=xnp1+theta*(xnp1-xn1)
		xb0=xnp0+theta*(xnp0-xn0)
		xb2=xnp2+theta*(xnp2-xn2)
	
		if keyword_set(verbose) then	print,"Max Flux:",max(abs(xnp1[*]))
			
		BeamInfo.amp =xnp1
		BMnp1=map*0.
		convol_map_local_2dbeam_sphere,BMnp1, BeamInfo ;Estimate of ptsrc map  

		if keyword_set(timer) then begin
			it_time=systime(1)
			print,"Convol end: ",it_time-ref_time," (",it_time-start_it_time,")"
		endif
					
		if(it lt niter-1) then begin ;Update primal variables	
			xn0=xnp0
			xn1=xnp1
			xn2=xnp2
		endif
	endfor
		
	if keyword_set(save_struc) then totit=it+save_struc.nit else totit=it

	RES_SPSR = {xn0: xn0, xn1: xn1,xn2:xn2, xnp0: xnp0, xnp1: xnp1, xnp2:xnp2,yn: yn,BMnp1: BMnp1, nit: totit}

	END_PROX:
	 
	return, RES_SPSR
end


;-----------------------------------------------------------------------------------------------------

function mrs_sparse_pointsource_removal, Map, GalMask, BeamInfo,std_map,  Niter=Niter, lmax=lmax,timer=timer,verbose=verbose,INIT_SPSR= INIT_SPSR, RES_SPSR= RES_SPSR

;Check compatibility in between all inputs
MapNpix=(size(Map,/dim))[0]
MaskNpix=(size(GalMask,/dim))[0]
StdNpix=(size(std_map,/dim))[0]
PSFreeMap=0

if(MapNpix ne MaskNpix) then begin
	print,"Galmask not compatible with Map (",MaskNpix," vs ",MapNpix,")"
	goto,END_SPSR
endif
if(StdNpix ne MaskNpix) then begin
	print,"Std Map not compatible with Map (",StdNpix," vs ",MaskNpix,")"
	goto,END_SPSR
endif
NsMap= npix2nside(MapNpix)

if(NsMap ne BeamInfo.nside) then begin
	print,"BeamInfo not compatible with Map (nside: ",BeamInfo.nside," vs ",NsMap," )"
	goto,END_SPSR
endif

if keyword_set(Niter) then if(Niter le 0) then begin
	print,"Niter not correct -> set to 0"
	Niter=0
endif
if keyword_set(lmax) then if(lmax le 0) then begin
	print,"lmax not correct -> set to 0"
	lmax=0
endif
if keyword_set(INIT_SPSR) then begin
	if not tag_exist(INIT_SPSR,'nit') OR not tag_exist(INIT_SPSR,'xn0') OR not tag_exist(INIT_SPSR,'xnp0') OR not tag_exist(INIT_SPSR,'xn1') OR not tag_exist(INIT_SPSR,'xnp1') OR not tag_exist(INIT_SPSR,'xn2') OR not tag_exist(INIT_SPSR,'xnp2') OR not tag_exist(INIT_SPSR,'yn') then begin
		print,"INIT_SPSR structure not compatible with SPSR"
		goto,END_SPSR
	endif
endif

;Launch main program
RES_SPSR = PD_SPSR(Map,  GalMask, BeamInfo,std_map,Niter=Niter,thrcst=0.01,lmax=lmax,timer=timer,verbose=verbose,save_struc= INIT_SPSR)

;Point source Free Map
PSFreeMap=Map-RES_SPSR.bmnp1

END_SPSR:
	return, PSFreeMap

end


;-----------------------------------------------------------------------------------------------------



