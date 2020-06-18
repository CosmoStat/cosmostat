;+
; NAME:  
;       mrs_gls_inpainting
;
; PURPOSE:
;        To compute the GLS solution for real spherical harmonic multipoles using multistep gradient descent (Nesterov scheme).
 ;        This routine is optimized for low l (l < 50) inpainting. It has been test with map of size nside=32.
;        All tests have been done using CMB maps in muK units.
;function mrs_gls_inpainting, Dat, Mask, Niter=Niter, lmax=lmax, cl=cl,  AlmRes= AlmRes,  verb=verb

; CALLING:
;     InpaintMap= mrs_gls_inpainting(Imag, Mask, Niter=Niter, lmax=lmax,  AlmRes = AlmRes, Cl=Cl) 
;
; INPUTS:
;     Imag -- IDL 1D array: Input Healpix image to be inpainted 
;     
; OUTPUTS:
;     InpaintMap -- IDL 1D array: Output inpainted map   
;          
; INPUT KEYWORDS:
;      Niter: int: number of iterations used in the reconstruction
;      Lmax      : Number of spherical harmonics computed in the decomposition (default is 50)
;      Cl            : Theoretical Cl  (used by Energy constraint)
;
; OUTPUT KEYWORDS:
;     AlmRes: Alm coefficients of the inpainted map.
;    
; EXAMPLE:
;     InpMap =  mrs_gls_inpainting(Map, Mask, Cl=Cl)
;     mrs_tv, InpMap
;
; HISTORY:
;       Written : Florent Sureau,   2012.
;-
;-----------------------------------------------------------------


function degrade_ideal,obs,final_nside,idealbeam=idealbeam
	;To degrade a map or several maps "obs" to nside "final_nside"
	;Optional return of the ideal beam (going from 1 to 0 between 2*nside<= l <= 3 Nside)
 
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	nside_obs=sqrt((size(obs,/dim))[0] / 12d)
	
	if(nside_obs lt final_nside) then begin
		print,"Degrade only"
		degraded_obs=0
		goto, END_DEGRADE
	endif

	if((size(obs))[0] eq 1) then nmaps=1 $
	else nmaps=(size(obs,/dim))[1] 
	
	beamdata=dblarr(3*final_nside+1)
	idealbeam=getidealbeam(beamdata, lmin=2*final_nside, lmax=3*final_nside, /tozero)
	
	degraded_obs=dblarr(ulong(final_nside)*ulong(final_nside)*12ul,nmaps)

	for km=0,nmaps-1 do begin
	 	mrs_almtrans,obs[*,km],obs_alms,/tab,lmax=3*final_nside 
		for kl=0,3*final_nside do obs_alms.alm[kl,*,*]=obs_alms.alm[kl,*,*] *idealbeam[kl]
		obs_alms.nside=final_nside
		mrs_almrec,obs_alms,deg_obs
		degraded_obs[*,km]=temporary(deg_obs)
	endfor
END_DEGRADE:
	return,degraded_obs
end


;------------------------------------------------------------
; Algos for GLS estimates
function GLS_svd,obs,almmatrix,lmin,lmax,sigma_noise=sigma_noise,pspec_alm=pspec_alm,pseudo_inv=pseudo_inv,singular=singular
	;To compute the GLS solution for real spherical harmonic multipoles using svds
	;
	;-Inputs
	;"obs": observed signal
	;"almmatrix": real spherical harmonics matrix restricted to observed pixels
	;"lmin": min multipole to estimate
	;"lmax": max multipole to estimate
	;"sigma_noise" = spatial standard deviation of the (detector noise)
	;"pspec_alm" =theoretical power spectrum of signal 
	;
	;-Options: 
	;"singular": flag to specify if the covariance matrix is singular
	
	;-Optional Return:
	;"pseudo-inverse": pseudo-inverse ;can be directly applied to any signal provided (with same mask, same noise, same theoretical power spectrum) to obtain the solution
	;
	;-Note: Some tricks are involved, based on the structure of the problem (covar=spatial part + alm part), mainly using Woodbury identity to perform inversion of the covariance matrix. For lmax=32 or 48, the time is reasonable
	COMPILE_OPT DEFINT32, STRICTARR,STRICTARRSUBS

	x=-1
	;check enough multipoles in power spectrun
	lmax_pspec =(size(pspec_alm,/dim))[0] 
	if(lmax gt lmax_pspec) then begin
		print,"Error : trying to inpaint till multipole ",strcompress(string(lmax),/remove_all)," but only powerspectrum till ",strcompress(string(lmax_pspec),/remove_all)," provided"
		goto,END_GLS_SVD
	endif	

	nobs=(size(obs,/dim))[0]
	if((size(almmatrix,/dim))[1] ne nobs) then begin	
		print,"Error : real sph. harmonics matrix not computed at every point "
		goto,END_GLS_SVD
	endif
	if not keyword_set(nit) then nit=1ul else nit=ulong(nit)

	;min and max indices for the estimation part
	index_min=(lmin)^2.
	index_max=(lmax+1)^2.-1

	;Initialization
	cov_space= sigma_noise* sigma_noise
	cov_l= pspec_alm
	cov_l[lmin:lmax]=0 ;only this part will be estimated

	nalms2est=index_max-index_min+1

	if keyword_set(singular) then cov_l[lmin:lmax]=1.

	;almmatrix which contains only the Ylms corresponding to the alms to be estimated
	sub_matrix=dblarr(nalms2est,nobs)
	for k=0,nobs-1 do sub_matrix[0:nalms2est-1, k]= almmatrix[index_min:index_max, k]

	;Include only relevant Ylms and weight them 
	index_excluded=where(cov_l ne 0,cnt_incl)
	if(cnt_incl gt 0) then begin
		total_excl=0
		for kl=0, cnt_incl-1 do total_excl=total_excl+2*index_excluded[kl]+1
		alm_matrix_nonest=dblarr(total_excl,nobs)
		offset=0ul	
		;Weight the Ylm by sqrt(Cov_l)
		for kl=0, cnt_incl-1 do begin
			index_min2=(index_excluded[kl])^2
			index_max2=(index_excluded[kl]+1)^2-1
			for kp=0,nobs-1 do alm_matrix_nonest[offset:offset+(index_max2-index_min2),kp]=almmatrix[index_min2:index_max2,kp]
			offset=offset+(index_max2-index_min2+1)
		endfor
	endif 

	;use woodbury matrix identity for inversion->do not have to invert a [nobs,nobs] complex matrix, but only a [nalms,nalms] matrix with nalms number of multipoles non estimated and a easy to invert diagonal [nobs,nobs] matrix 
	;S contains spherical harmonics not estimated
	;(C_sp + S##C_alm##S^*)^{-1}=C_sp^{-1} -C_sp^{-1}##S##(C_alm^{-1} +S^*## C_sp^{-1} ## S)^{-1} ## S^* ## C_sp^{-1}
	if(cnt_incl gt 0) then begin
		weight_sp_sub_matrix= alm_matrix_nonest
		for k=0,nobs-1 do weight_sp_sub_matrix[*, k]= alm_matrix_nonest[*, k]/sigma_noise[k] ; C_sp^{-1/2} ## S
		inv_interm_bef=reform(transpose(weight_sp_sub_matrix)##weight_sp_sub_matrix) ;S^*## C_sp^{-1} ## S
		weight_sp_sub_matrix=0
		offset=0
		inv_interm= temporary(inv_interm_bef)
		for kl=0, cnt_incl-1 do begin
			index_min2=(index_excluded[kl])^2
			index_max2=(index_excluded[kl]+1)^2-1
			for kp= offset,offset+(index_max2-index_min2) do  inv_interm[kp,kp]= inv_interm[kp,kp]+1d / cov_l[index_excluded[kl]];C_alm^{-1} +S^*## C_sp ## S
			offset=offset+(index_max2-index_min2+1)
		endfor
		svdc, inv_interm,d_int,u_int,v_int
		inv_interm=0.
		print,"Cond number for covar inversion:",max(d_int[*])/min(d_int[*])
		matrix_inv_1= temporary(v_int)##diag_matrix(1d / d_int)##transpose(temporary(u_int))	
		;compute produc_matrix_inv1
		matrix_inv= alm_matrix_nonest##temporary(matrix_inv_1)##transpose(alm_matrix_nonest)
		;Get the inverse of (At \Simga ^-1 A) with A matrix containing the amls to estimate
		matrix_1=dblarr(nalms2est,nobs)
		for k=0,nobs-1 do matrix_1[*,k]=sub_matrix[*,k]/cov_space[k] ;matrix_1=C_sp^{-1} A
		Prematrix=transpose(sub_matrix)##matrix_1 - transpose(matrix_1)##matrix_inv##matrix_1
	endif else begin
		matrix_1=dblarr(nalms2est,nobs)
		for k=0,nobs-1 do matrix_1[*,k]=sub_matrix[*,k]/cov_space[k] ;matrix_1=C_sp^{-1} A
		Prematrix=transpose(temporary(sub_matrix))##matrix_1
	endelse
	alm_matrix_nonest=0
	sub_matrix=0		
	;Get the inverse of (At \Simga ^-1 A) with A matrix containing the alms to estimate
	svdc,temporary(Prematrix),d_pre,u_pre,v_pre
	print,"Cond number for Pseudo-inverse inversion:",max(d_pre[*])/min(d_pre[*])
	inv_prematrix=temporary(v_pre)##diag_matrix(1d / d_pre)##transpose(temporary(u_pre))
	pseudo_inv=temporary(inv_prematrix)##temporary(transpose(matrix_1))
	pseudo_inv=float(pseudo_inv) ;necessary because of space lims
	if(cnt_incl gt 0) then begin
		 matrix_2=fltarr(nobs,nobs)
		 for k=0,nobs-1  do matrix_2[k,*]= -matrix_inv[k,*]/cov_space[k] 
	 	 for k=0,nobs-1  do matrix_2[k,k]= 1d + matrix_2[k,k]
	 	 matrix_inv=0
		 ;A^t \Sigma^{-1}=transpose(matrix_1)##matrix_2=A^t C_sp^{-1} ## (Id - matrix_inv##C_sp^{-1})
		 pseudo_inv=pseudo_inv##temporary(matrix_2) 
	endif 
	
	x=reform(pseudo_inv##obs)	
END_GLS_SVD:
	sol=dblarr(index_max+1)
	sol[index_min:index_max]= x
	return,sol
end


function CG_GLS_inpainting,obs,almmatrix,lmin,lmax,sigma_noise= sigma_noise,pspec_alm=pspec_alm,nit=nit,norm_matrix=norm_matrix
	;-To compute the GLS solution for real spherical harmonic multipoles using conjugate gradients
	;-Inputs
	;"obs": observed signal
	;"almmatrix": real spherical harmonics matrix restricted to observed pixels
	;"lmin": min multipole to estimate
	;"lmax": max multipole to estimate
	;"sigma_noise" = spatial standard deviation of the (detector noise)
	;"pspec_alm" =theoretical power spectrum of signal 
	;"nit": number of iterations
	;
	;-options: 
	;"norm_matrix": should be set if a mask is used (take into account the loss of power for multipoles to be estimated)
	;-optional return: 

	x=-1
	;check enough multipoles in power spectrun
	lmax_pspec =(size(pspec_alm,/dim))[0] 
	if(lmax gt lmax_pspec) then begin
		print,"Error : trying to inpaint till multipole ",strcompress(string(lmax),/remove_all)," but only powerspectrum till ",strcompress(string(lmax_pspec),/remove_all)," provided"
		goto,END_CG
	endif	

	nobs=(size(obs,/dim))[0]
	if((size(almmatrix,/dim))[1] ne nobs) then begin	
		print,"Error : real sph. harmonics matrix not computed at every point "
		goto,END_CG
	endif
	if not keyword_set(nit) then nit=200ul else nit=ulong(nit)

	;min and max indices for the estimation part
	index_min=(lmin)^2.
	index_max=(lmax+1)^2.-1

	;Initialization
	cov_space= sigma_noise* sigma_noise
	cov_l= pspec_alm
	cov_l[lmin:lmax]=0 ;this part will be estimated

	;almmatrix which contains only the Ylms corresponding to the alms to be estimated
	sub_matrix=dblarr(index_max-index_min+1,nobs)
	for k=0,nobs-1 do sub_matrix[0:index_max-index_min, k]= almmatrix[index_min:index_max, k]

	if keyword_set(norm_matrix) then begin ;orthonormalize the gram matrix
		gram=transpose(sub_matrix)##(sub_matrix)
		svdc,gram,d_svd,u_svd,v_svd
		;print,"eigenvalues of gram estimation matrix=",d_svd
		weight_matrix=u_svd##diag_matrix(1d / sqrt(d_svd))##transpose(u_svd)
		sub_matrix= sub_matrix##weight_matrix ;we are now looking for the solution in the whole space spanned by the reweighted alms
	endif


	;Include only relevant Ylms and weight them 
	index_included=where(cov_l ne 0,cnt_incl)
	if(cnt_incl gt 0) then begin
		total_alms=0
		for kl=0, cnt_incl-1 do total_alms=total_alms+2*index_included[kl]+1
		alm_matrix_cov=dblarr(total_alms,nobs)
		alm_matrix_nonest=dblarr(total_alms,nobs)
		offset=0ul	
		;Weight the Ylm by sqrt(Cov_l)
		for kl=0, cnt_incl-1 do begin
			index_min2=(index_included[kl])^2
			index_max2=(index_included[kl]+1)^2-1
			for kp=0,nobs-1 do alm_matrix_nonest[offset:offset+(index_max2-index_min2),kp]=almmatrix[index_min2:index_max2,kp]
			for kp=0,nobs-1 do alm_matrix_cov[offset:offset+(index_max2-index_min2),kp]=almmatrix[index_min2:index_max2,kp]*sqrt(cov_l[index_included[kl]])
			offset=offset+(index_max2-index_min2+1)
		endfor
	endif else for kp=0,nobs-1 do alm_matrix_cov= dblarr(2,nobs)
	augm_vector=dblarr(2*nobs)
	augm_vector[0:nobs-1]=obs ;observation in first nobs, lagrange multipliers in the second part
	init=dblarr(nobs)

	;Initialize 
	conju=augm_vector	
	conjv=augm_vector
	conju_norm=norm(conju,l=2)
	nconju= conju/conju_norm

	;At u1
	;first project the lagrange multipliers into the orthogonal of the range of S
	proj=nconju[nobs:2*nobs-1]-reform(sub_matrix##reform(transpose(sub_matrix)##nconju[nobs:2*nobs-1]))
	conjv[0:nobs-1]=nconju[0:nobs-1]+proj
	;Then project the observations into the orthogonal of the range of S and multiply with covariance
	proj=nconju[0:nobs-1]-reform(sub_matrix##reform(transpose(sub_matrix)##nconju[0:nobs-1]))
	conjv[nobs:2*nobs-1]=reform(alm_matrix_cov##reform(transpose(alm_matrix_cov)##proj))+proj*cov_space
	conjv_norm=norm(conjv,l=2)
	nconjv=conjv/conjv_norm
	conjw=nconjv
	x=0
	phibar= conju_norm
	rhobar=conjv_norm
	it=1ul
	notconvergence =1
	while(notconvergence) do begin
		;A v_l-alpha_l u_l
		proj=nconjv[nobs:2*nobs-1]-reform(sub_matrix##reform(transpose(sub_matrix)##nconjv[nobs:2*nobs-1]))
		conju[0:nobs-1]=reform(alm_matrix_cov##reform(transpose(alm_matrix_cov)##proj))+proj*cov_space+nconjv[0:nobs-1]-conjv_norm*nconju[0:nobs-1]
		conju[nobs:2*nobs-1]=nconjv[0:nobs-1]-reform(sub_matrix##reform(transpose(sub_matrix)##nconjv[0:nobs-1]))-conjv_norm*nconju[nobs:2*nobs-1];
		conju_norm=norm(conju,l=2)
		nconju= conju/conju_norm
		;At u_l+1 - beta_l+1 v_l
		proj=nconju[nobs:2*nobs-1]-reform(sub_matrix##reform(transpose(sub_matrix)##nconju[nobs:2*nobs-1]))
		conjv[0:nobs-1]=nconju[0:nobs-1]+proj-conju_norm*nconjv[0:nobs-1]
		proj=nconju[0:nobs-1]-reform(sub_matrix##reform(transpose(sub_matrix)##nconju[0:nobs-1]))
		conjv[nobs:2*nobs-1]=reform(alm_matrix_cov##reform(transpose(alm_matrix_cov)##proj))+proj*cov_space-conju_norm*nconjv[nobs:2*nobs-1]
		conjv_norm=norm(conjv,l=2)
		nconjv= conjv/conjv_norm
		rho=sqrt(rhobar*rhobar+conju_norm*conju_norm)
		c=rhobar/rho
		s=conju_norm/rho
		theta=s*conjv_norm
		rhobar=-c*conjv_norm
		phi=c*phibar
		phibar=s*phibar
		xp=x+(phi/rho)*conjw
		conjwp=nconjv - (theta/rho) * conjw
		print,"Conj grad it:",strcompress(string(it),/remove_all),", norm of resitual=",strcompress(string(phibar),/remove_all)
		it=it+1ul	
		if(it gt nit) then notconvergence=0	
		x=xp
		conjw=conjwp
	endwhile	
	
	END_CG:
	if keyword_set(norm_matrix) then begin
		solw=reform(transpose(sub_matrix)##x[0:nobs-1]) ;solution in weighted frame
		x_w=reform(weight_matrix##solw)
	endif else x_w=reform(transpose(sub_matrix)##x[0:nobs-1])
		sol=dblarr(index_max+1)
		sol[index_min:index_max]= x_w
		return,sol

end


function conjugate_gradient_matrix_lsqr_cov_alm_and_space,vector,alm_matrix,cov_l=cov_l,cov_space=cov_space,nit=nit
	;subroutine to invert a covariance matrix with both elements in space and spherical harmonics
	;called by the gradient descent algos
	nobs=(size(vector,/dim))[0]
	init=dblarr(nobs)

	;Initialize 
	conju=vector	
	conju_norm=norm(conju,l=2)
	nconju= conju/conju_norm
	
	;Include only relevant Ylms
	index_included=where(cov_l gt 0,cnt_incl)
	if(cnt_incl gt 0) then begin
		total_alms=0
		for kl=0, cnt_incl-1 do total_alms=total_alms+2*index_included[kl]+1
		alm_matrix_incl=dblarr(total_alms,nobs)
		offset=0ul	
		;Weight the Ylm by sqrt(Cov_l)
		for kl=0, cnt_incl-1 do begin
			index_min=(index_included[kl])^2
			index_max=(index_included[kl]+1)^2-1
			for kp=0,nobs-1 do alm_matrix_incl[offset:offset+index_max-index_min,kp]= alm_matrix[index_min:index_max,kp]*sqrt(cov_l[index_included[kl]])
			offset=offset+(index_max-index_min+1)
		endfor
	endif else for kp=0,nobs-1 do alm_matrix_incl= dblarr(2,nobs)

	cov_alm=reform(alm_matrix_incl##transpose(alm_matrix_incl))
	
	conjv= cov_alm##nconju+nconju*cov_space	
	conjv_norm=norm(conjv,l=2)
	nconjv=conjv/conjv_norm
	conjw=nconjv
	x=0
	phibar= conju_norm
	rhobar=conjv_norm
	it=1ul
	notconvergence =1
	bk_norm2=0
	dk=conjv*0.
	dk_norm2=0.
	theta=0.
	while(notconvergence) do begin
		conju=cov_alm##nconjv+ nconjv*cov_space-conjv_norm*nconju	
		conju_norm=norm(conju,l=2)
		nconju= conju/conju_norm
		rho=sqrt(rhobar*rhobar+conju_norm*conju_norm)
		;we have bk_norm2<=||A||_F^2<=dk_norm2
		bk_norm2=bk_norm2+conjv_norm*conjv_norm+conju_norm*conju_norm
		dk=1d / rho * (nconjv - theta*dk)
		dk_norm2=dk_norm2+norm(dk[*],l=2)^2.
		conjv=cov_alm##nconju+nconju*cov_space-conju_norm*nconjv
		conjv_norm=norm(conjv,l=2)
		nconjv= conjv/conjv_norm
		c=rhobar/rho
		s=conju_norm/rho		
		theta=s*conjv_norm
		rhobar=-c*conjv_norm
		phi=c*phibar
		phibar=s*phibar
		xp=x+(phi/rho)*conjw
		conjwp=nconjv - (theta/rho) * conjw
		it=it+1ul	
		if(it gt nit) then notconvergence=0
		Atrk=phibar*conjv_norm*abs(c)
		x=xp
		conjw=conjwp
		test=Atrk / (phibar*bk_norm2)
	endwhile
	print,"Conj grad it:",strcompress(string(it),/remove_all),", norm of resitual=",strcompress(string(phibar),/remove_all)," , test=",test

	return,x
end

function first_order_gradient_descent_GLS_inpainting,obs,almmatrix,lmin,lmax,sigma_noise= sigma_noise,pspec_alm=pspec_alm,nit=nit,conjnit=conjnit,CG=CG

	;-To compute the GLS solution for real spherical harmonic multipoles using gradient descent
	;-Inputs
	;"obs": observed signal
	;"almmatrix": real spherical harmonics matrix restricted to observed pixels
	;"lmin": min multipole to estimate
	;"lmax": max multipole to estimate
	;"sigma_noise" = spatial standard deviation of the (detector noise)
	;"pspec_alm" =theoretical power spectrum of signal 
	;"nit": number of iterations
	;
	;-options: 
	;"CG": should be set if a conjugate gradient step is used for covariance matrix inversion (default: clever SVD inversion that should be used for small scales problems)
	;"nitconj": if the "CG" flag is used, number of iterations of CG for inversion (should be high enough)
	;-optional return: 

	lmax_pspec =(size(pspec_alm,/dim))[0] 
	if(lmax gt lmax_pspec) then begin
		print,"Error : trying to inpaint till multipole ",strcompress(string(lmax),/remove_all)," but only powerspectrum till ",strcompress(string(lmax_pspec),/remove_all)," provided"
		goto,END_GRADIENT
	endif	

	index_min=(lmin)^2.
	index_max=(lmax+1)^2.-1
	ak=dblarr(index_max-index_min+1)

	nobs=(size(obs,/dim))[0]
	if((size(almmatrix,/dim))[1] ne nobs) then begin	
		print,"Error : real sph. harmonics matrix not computed at every point "
		goto,END_GRADIENT
	endif


	;Initialization
	cov_space= sigma_noise* sigma_noise
	cov_l= pspec_alm
	cov_l[lmin:lmax]=0

	DEF_ALM_FAST=0
	DEF_ALM_NITER=0
	notconvergence=1
	it=1ul
	if not keyword_set(nit) then nit=1ul else nit=ulong(nit)
	if not keyword_set(nitconj) then nitconj =40ul else nit=ulong(nit)
	;almmatrix which contains only the Ylms corresponding to the alms to be estimated
	sub_matrix=dblarr(index_max-index_min+1,nobs)
	for k=0,nobs-1 do sub_matrix[0:index_max-index_min, k]= almmatrix[index_min:index_max, k]
	cov_sub_matrix=transpose(sub_matrix)##sub_matrix
	svdc,temporary(cov_sub_matrix),d_sub,u_sub,v_sub	
	print,"Cond number for submatrix inversion:",max(d_sub[*])/min(d_sub[*])
	
	;Include only relevant Ylms
	index_included=where(cov_l ne 0,cnt_incl)
	if(cnt_incl gt 0) then begin
		total_alms=0
		for kl=0, cnt_incl-1 do total_alms=total_alms+2*index_included[kl]+1
		alm_matrix_incl=dblarr(total_alms,nobs)
		alm_matrix_weight=dblarr(total_alms,nobs)
		offset=0ul	
		;Weight the Ylm by sqrt(Cov_l)
		for kl=0, cnt_incl-1 do begin
			index_min2=(index_included[kl])^2
			index_max2=(index_included[kl]+1)^2-1
			for kp=0,nobs-1 do alm_matrix_incl[offset:offset+index_max2-index_min2,kp]= almmatrix[index_min2:index_max2,kp]
			for kp=0,nobs-1 do alm_matrix_weight[offset:offset+index_max2-index_min2,kp]= almmatrix[index_min2:index_max2,kp]*sqrt(cov_l[index_included[kl]])
			offset=offset+(index_max2-index_min2+1)
		endfor
	endif else for kp=0,nobs-1 do alm_matrix_incl= dblarr(2,nobs)
	cov_alm=alm_matrix_weight##transpose(alm_matrix_weight)
	alm_matrix_weight=0
	
	;we need to know the norm of S^* Sigma^-1 S
	;easy to show that >=lambda_min(AAt+YclYt)>=min(lambda_min(AAt)+lambda_min(BBt))
	min_eig_space=min(cov_space[*])
	min_eig_alm=0 ;;cov_alm is singular
	max_inv=max(d_sub[*]) / (min_eig_space + min_eig_alm)	
	max_mu=2d/max_inv ;maximal gradient step
	print,max_mu

	if not keyword_set(cg) then begin
		;use woodbury matrix identity for inversion->do not have to invert a [nobs,nobs] matrix, but only a [nalms,nalms] matrix with nalms number of multipoles non estimated
		;S contains spherical harmonics not estimated
		;(C_sp + S##C_alm##S^*)^{-1}=C_sp^{-1} -C_sp^{-1}##S##(C_alm^{-1} +S^*## C_sp^{-1} ## S)^{-1} ## S^* ## C_sp^{-1}
		if(cnt_incl gt 0) then begin
			weight_sp_sub_matrix= alm_matrix_incl
			for k=0,nobs-1 do weight_sp_sub_matrix[*, k]= alm_matrix_incl[*, k]/sigma_noise[k] ; C_sp^{-1/2} ## S
			inv_interm=reform(transpose(weight_sp_sub_matrix)##weight_sp_sub_matrix) ;S^*## C_sp^{-1} ## S
			weight_sp_sub_matrix=0
			offset=0
			for kl=0, cnt_incl-1 do begin
				index_min2=(index_included[kl])^2
				index_max2=(index_included[kl]+1)^2-1
				for kp= offset,offset+(index_max2-index_min2) do  inv_interm[kp,kp]= inv_interm[kp,kp]+1d / cov_l[index_included[kl]];C_alm^{-1} +S^*## C_sp ## S
				offset=offset+(index_max2-index_min2+1)
			endfor
			svdc, temporary(inv_interm),d_int,u_int,v_int
			print,"Cond number for covar inversion:",max(d_int[*])/min(d_int[*])
			matrix_inv_1= temporary(v_int)##diag_matrix(1d / d_int)##transpose(temporary(u_int))	
			inv_covar= alm_matrix_incl##temporary(matrix_inv_1)##transpose(alm_matrix_incl)
			ALM_MATRIX_INCL=0
			for k1=0,nobs-1 do for k2=0,nobs-1 do inv_covar[k1,k2]=-inv_covar[k1,k2]/cov_space[k1]/cov_space[k2]
			for k1=0,nobs-1 do 	inv_covar[k1,k1]=inv_covar[k1,k1]+1d/cov_space[k1]
			;Now check if necessary
	;		check=inv_covar##diag_matrix(cov_space)+ inv_covar##reform(alm_matrix_incl_weighted##reform(transpose(alm_matrix_incl_weighted)))
		endif else inv_cov=diag_matrix(1d/cov_space)
	endif

	ak=(reform(transpose(almmatrix)##obs))[index_min:index_max]*0.
	while(notconvergence) do begin
		res= obs-reform(sub_matrix##ak)
		print,"it",it,"residual=",total(res*res,/double)
		if keyword_set(CG) then uk=conjugate_gradient_matrix_lsqr_cov_alm_and_space(res, almmatrix,cov_l=cov_l,cov_space=cov_space,nit=conjnit) $
		else uk=inv_covar##res
		res_inv=temporary(res)-(uk*cov_space+cov_alm##uk) 
		print,"it",it,"residual inversion=",total(res_inv*res_inv,/double)
		res_inv=0
		ak=ak+ max_mu*0.9 *reform(transpose(sub_matrix)##uk)
		it=it+1ul	
		if(it gt nit) then notconvergence=0	
	endwhile

	END_GRADIENT:
		sol=dblarr(index_max+1)
		sol[index_min:index_max]= ak
		return,sol
end

function mstep_first_order_gradient_descent_GLS_inpainting,obs,almmatrix,lmin,lmax,sigma_noise= sigma_noise,pspec_alm=pspec_alm,nit=nit,conjnit=conjnit,cg=cg,singular=singular

	;-To compute the GLS solution for real spherical harmonic multipoles using multistep gradient descent (Nesterov scheme)
	;-Inputs
	;"obs": observed signal
	;"almmatrix": real spherical harmonics matrix restricted to observed pixels
	;"lmin": min multipole to estimate
	;"lmax": max multipole to estimate
	;"sigma_noise" = spatial standard deviation of the (detector noise)
	;"pspec_alm" =theoretical power spectrum of signal 
	;"nit": number of iterations
	;
	;-options: 
	;"CG": should be set if a conjugate gradient step is used for covariance matrix inversion (default: clever SVD inversion that should be used for small scales problems)
	;"nitconj": if the "CG" flag is used, number of iterations of CG for inversion (should be high enough)
	;-optional return: 


	lmax_pspec =(size(pspec_alm,/dim))[0] 
	if(lmax gt lmax_pspec) then begin
		print,"Error : trying to inpaint till multipole ",strcompress(string(lmax),/remove_all)," but only powerspectrum till ",strcompress(string(lmax_pspec),/remove_all)," provided"
		goto,END_GRADIENT
	endif	

	index_min=(lmin)^2.
	index_max=(lmax+1)^2.-1
	ak=dblarr(index_max-index_min+1)

	nobs=(size(obs,/dim))[0]
	if((size(almmatrix,/dim))[1] ne nobs) then begin	
		print,"Error : real sph. harmonics matrix not computed at every point "
		goto,END_GRADIENT
	endif

	;Initialization
	cov_space= sigma_noise* sigma_noise
	cov_l= pspec_alm

	alm_excl=where(pspec_alm ne 0,cnt_excl)

	DEF_ALM_FAST=0
	DEF_ALM_NITER=0
	notconvergence=1
	it=1ul
	if not keyword_set(nit) then nit=1ul else nit=ulong(nit)
	if not keyword_set(nitconj) then nitconj =40ul else nit=ulong(nit)
	;almmatrix which contains only the Ylms corresponding to the alms to be estimated
	sub_matrix=dblarr(index_max-index_min+1,nobs)
	for k=0,nobs-1 do sub_matrix[0:index_max-index_min, k]= almmatrix[index_min:index_max, k]	
	cov_sub_matrix=transpose(sub_matrix)##sub_matrix
	svdc,temporary(cov_sub_matrix),d_sub,u_sub,v_sub	
	print,"Cond number for submatrix inversion:",max(d_sub[*])/min(d_sub[*])
		
	if keyword_set(singular) then cov_l[lmin:lmax]=1d $
	else cov_l[lmin:lmax]=0

	;Include only relevant Ylms
	index_included=where(cov_l ne 0,cnt_incl)
	if(cnt_incl gt 0) then begin
		total_alms=0
		for kl=0, cnt_incl-1 do total_alms=total_alms+2*index_included[kl]+1
		alm_matrix_incl=dblarr(total_alms,nobs)
		alm_matrix_incl_weighted=dblarr(total_alms,nobs)
		offset=0ul	
		for kl=0, cnt_incl-1 do begin
			index_min2=(index_included[kl])^2
			index_max2=(index_included[kl]+1)^2-1
			for kp=0,nobs-1 do alm_matrix_incl[offset:offset+index_max2-index_min2,kp]= almmatrix[index_min2:index_max2,kp]
		;Weight the Ylm by sqrt(Cov_l)
			for kp=0,nobs-1 do alm_matrix_incl_weighted[offset:offset+index_max2-index_min2,kp]= almmatrix[index_min2:index_max2,kp]*sqrt(cov_l[index_included[kl]])
			offset=offset+(index_max2-index_min2+1)
		endfor
	endif else begin
		for kp=0,nobs-1 do alm_matrix_incl= dblarr(2,nobs)
		for kp=0,nobs-1 do alm_matrix_incl_weighted = dblarr(2,nobs)
	endelse
	cov_alm= alm_matrix_incl_weighted##transpose(alm_matrix_incl_weighted)
	ALM_MATRIX_INCL_WEIGHTED=0
	;we need to know the norm of S^* Sigma^-1 S
	;easy to show that >=lambda_min(AAt+YclYt)>=min(lambda_min(AAt)+lambda_min(BBt))
	min_eig_space=min(cov_space[*])
	min_eig_alm=0 ;cov_alm is singular
	max_inv=max(d_sub[*]) / (min_eig_space+min_eig_alm)	
	max_mu=1d/max_inv ;maximal gradient step
	print,max_mu
	print,index_min,index_max

	;multi-step parameters
	alphak=double(indgen(nit+1)+1) / 2d
	tauk=2d / double(indgen(nit+1)+3)

	ak=(reform(transpose(almmatrix)##obs))[index_min:index_max]
	bk=ak*0.
	gk=ak*0.

	if not keyword_set(cg) then begin
		;use woodbury matrix identity for inversion->do not have to invert a [nobs,nobs] matrix, but only a [nalms,nalms] matrix with nalms number of multipoles non estimated
		;S contains spherical harmonics not estimated
		;(C_sp + S##C_alm##S^*)^{-1}=C_sp^{-1} -C_sp^{-1}##S##(C_alm^{-1} +S^*## C_sp^{-1} ## S)^{-1} ## S^* ## C_sp^{-1}
		if(cnt_incl gt 0) then begin
			weight_sp_sub_matrix= alm_matrix_incl
			for k=0,nobs-1 do weight_sp_sub_matrix[*, k]= alm_matrix_incl[*, k]/sigma_noise[k] ; C_sp^{-1/2} ## S
			inv_interm=reform(transpose(weight_sp_sub_matrix)##weight_sp_sub_matrix) ;S^*## C_sp^{-1} ## S
			weight_sp_sub_matrix=0

			offset=0
			for kl=0, cnt_incl-1 do begin
				index_min2=(index_included[kl])^2
				index_max2=(index_included[kl]+1)^2-1
				for kp= offset,offset+(index_max2-index_min2) do  inv_interm[kp,kp]= inv_interm[kp,kp]+1d / cov_l[index_included[kl]];C_alm^{-1} +S^*## C_sp ## S
				offset=offset+(index_max2-index_min2+1)
			endfor
			svdc, temporary(inv_interm),d_int,u_int,v_int
			print,"Cond number for covar inversion:",max(d_int[*])/min(d_int[*])
			matrix_inv_1= temporary(v_int)##diag_matrix(1d / d_int)##transpose(temporary(u_int))	
			inv_covar= alm_matrix_incl##temporary(matrix_inv_1)##transpose(alm_matrix_incl)
			ALM_MATRIX_INCL=0
			for k1=0,nobs-1 do for k2=0,nobs-1 do inv_covar[k1,k2]=-inv_covar[k1,k2]/cov_space[k1]/cov_space[k2]
			for k1=0,nobs-1 do 	inv_covar[k1,k1]=inv_covar[k1,k1]+1d/cov_space[k1]
			;Now check if necessary
	;		check=inv_covar##diag_matrix(cov_space)+ inv_covar##reform(alm_matrix_incl_weighted##reform(transpose(alm_matrix_incl_weighted)))
		endif else inv_cov=diag_matrix(1d/cov_space)
	endif
	mu_grad=max_mu*0.9
	print,"mu_grad=",mu_grad	

	while(notconvergence) do begin
		res= obs-reform(sub_matrix##bk)
		if keyword_set(CG) then uk=conjugate_gradient_matrix_lsqr_cov_alm_and_space(res, almmatrix,cov_l=cov_l,cov_space=cov_space,nit=conjnit) $
		else uk=inv_covar##res
		res_inv=temporary(res)-(uk*cov_space+cov_alm##uk) 
		print,"it",it,"residual inversion=",total(res_inv*res_inv,/double)		
		res_inv=0
		ak=bk+  mu_grad*reform(transpose(sub_matrix)##uk) ;yk=TQ(xk)=xk-Grad f(x_k)
		gk=gk+  mu_grad*alphak[it-1]*reform(transpose(sub_matrix)##uk) ;weighted gradient accumulation: zk
		bk=(1d - tauk[it-1])*ak[*] +tauk[it-1]*gk[*] ;xk+1=tau_k z_k + (1- tau_k) yk
		res= obs-reform(sub_matrix##ak)
		print,"residual=",total(res*res,/double)
		it=it+1ul	
		if(it gt nit) then notconvergence=0	
	endwhile

	END_GRADIENT:
		sol=dblarr(index_max+1)
		sol[index_min:index_max]= ak
		return,sol
end


;------------------------------------------------------------
;Real Spherical harmonics routines
pro precompute_spherical_harmonics_factors,lmax,Fm,logval
	;Routine called by "compute_real_spherical_harmonics"
	;Start by computing the Legendre polynomials using recurrence
	;note Rm(l,theta)=(-1)^m *sqrt((2l+1)/(4*!dpi) * fact(l-m)/fact(l+m)) * Plm(cos(theta))
	;Then Rm(l+1,theta)= Rm(l,theta)*Fm(l)*cos(theta)-Rm(l-1)* Fm(l)/Fm(l-1)
	Fm=dblarr(lmax+1,lmax+1,2)
	for km=0, lmax do begin	
		F_old=1.
		for kl=km,lmax do begin
			Fm[kl,km,0]=sqrt(2d *kl +3d)/sqrt(kl+km+1d) * sqrt(2d *kl + 1d)/sqrt(kl - km + 1d) ;FM that is used for kl+1
			Fm[kl,km,1]=Fm[kl,km,0]/F_old
			F_old=Fm[kl,km,0]
		endfor
	endfor
	logval=dblarr(lmax+1) ;this corresponds to log(Rl(l,theta))-l*log(sin(theta))
	new_fact=0d
	fourpi=(4d * !dpi)
	logval[0]=-0.5*alog(fourpi)
	for kl=1,lmax do begin 
		new_fact=new_fact+alog(1d - 1d/(2d * kl))
		logval[kl]=0.5*(alog((2d * kl +1d)/fourpi) + new_fact)
	endfor
end


pro real2complex_sph_list,real_list,complex_tab,tab=tab
	;This function transforms real spherical harmonics into complex spherical harmonics, either the usual tab if tab is set or the usual list
	lmax=sqrt((size(real_list,/dim))[0])-1
	if((size(real_list))[0] eq 1) then nb_pix=0 $
	else nb_pix=(size(real_list,/dim))[1]
	if keyword_set(tab) then begin
		if(nb_pix gt 0) then complex_tab=dblarr(lmax+1,lmax+1,2,nb_pix) $
		else complex_tab=dblarr(lmax+1,lmax+1,2)
		for kl=0,lmax do begin
			offset_l=ulong(kl)^2.
			complex_tab[kl,0,0,*]= real_list[offset_l,*]
			for km=1,kl do begin
				offset_m=offset_l+2*km-1
				complex_tab[kl,km,0,*]= real_list[offset_m,*]/sqrt(2d)
				complex_tab[kl,km,1,*]= -real_list[offset_m +1,*]/sqrt(2d) ;comes from the definition of real/complex spherical harmonics
			endfor
		endfor
	endif else begin
		if(nb_pix gt 0) then complex_tab=dblarr((lmax+1ul)*(lmax+2ul)/2ul,2,nb_pix) $
		else complex_tab=dblarr((lmax+1ul)*(lmax+2ul)/2ul,2)
		for kl=0,lmax do begin
			offset_l=ulong(kl)^2.
			new_offset=kl*(kl+1ul)/2ul
			complex_tab[new_offset,0,*]= real_list[offset_l,*]
			for km=1,kl do begin
				offset_m=offset_l+2*km-1
				complex_tab[new_offset+km,0,*]= real_list[offset_m,*]/sqrt(2d)
			 	complex_tab[new_offset+km,1,*]= -real_list[offset_m +1,*]/sqrt(2d) ;comes from the definition of real/complex spherical harmonics
			endfor
		endfor
	endelse
end

pro real2pspec_sph_list,real_list,pspec
	;This function computes power spectrum from real spherical harmonics 
	lmax=sqrt((size(real_list,/dim))[0])-1
	if((size(real_list))[0] eq 1) then nb_pix=0 $
	else nb_pix=(size(real_list,/dim))[1]

	if(nb_pix gt 0) then complex_tab=dblarr(lmax+1,lmax+1,2,nb_pix) $
	else complex_tab=dblarr(lmax+1,lmax+1,2)
	pspec=dblarr(lmax+1)
	for kl=0,lmax do begin
		index_min=ulong(kl)^2.
		index_max=ulong(kl+1)^2.-1.
		pspec[kl]=total(real_list[index_min: index_max]*real_list[index_min: index_max])/(2d * kl +1.)
	endfor
	
end

function compute_real_spherical_harmonics,nside,lmax,coords,recursive_phi=recursive_phi,normalize=normalize, weightring= weightring
	
	;- Function used to compute the real spherical harmonic matrix
	;-Inputs
	;"nside": nside of the maps (give the number of lines)
	;"lmax": max multipole to compute (give the number of columns = (lmax+1)^2.)
	;
	;-options: 
	;"recursive_phi": compute the multiple sine and cosine recursively (no real gain in time)
	;"normalize": the real spherical harmonics are normalized (in the sense that multipole 0 = sqrt(4pi/Nelem) for all elements)
	;"weightring": the real spherical harmonics are weighted by an function of colatitude for better accuracy (coming from HEALPIX)
	;-optional return: 

	;NOTE: Considering that IDL cannot go further in precision than double
	;and that we use recursions to compute the weighted associated Legendre polynomials
	;this function should only be used for small multipoles, i.e. do not go higher than for instance 64
	;this is adapted from Healpix way of computing the ALMs (see ylmgen.h, by Reinecke)
	;Essentially, recursive computation of the associated legendre polynomials require to pre-compute some factors, which is done in precompute_spherical_harmonics_factors

	;coords=[npix,2]
	precompute_spherical_harmonics_factors,lmax,Fm,logval
	nb_pix=(size(coords,/dim))[0]
	alm_matrix=dblarr((lmax+1)^2.,nb_pix)
	fourpi=(4d * !dpi)
	npixtot=ulong(nside)^2ul * 12ul

	;sort theta:
	srt_theta=SORT(coords[*,0]) ;sorted indices
	theta_sorted=coords[srt_theta,0]	
	list_theta_uq=UNIQ(theta_sorted) ;list of indices in coords containing uniq indices 
	nb_theta_uq=(size(list_theta_uq,/dim))[0] ; nb of unique thetas
	val_theta_uq= theta_sorted[list_theta_uq] ;corresponding value of theta
	theta_sorted=0
	index_theta=ulonarr(nb_pix) ;index in the unique list for each original theta
	index_theta[srt_theta[0:list_theta_uq[0]]]= 0
	for k=1, nb_theta_uq-1 do index_theta[srt_theta[list_theta_uq[k-1]+1:list_theta_uq[k]]]=k	
	precision = 15 ;15 digits precision
	small_nb=10^(-precision)
	big_nb=10^(precision)

	;sort phi
	srt_phi=SORT(coords[*,1]) ;sorted indices
	phi_sorted=coords[srt_phi,1]	
	list_phi_uq=UNIQ(phi_sorted) ;list of indices in coords containing uniq indices 
	nb_phi_uq=(size(list_phi_uq,/dim))[0] ; nb of unique thetas
	val_phi_uq= phi_sorted[list_phi_uq] ;corresponding value of theta
	phi_sorted=0
	index_phi=ulonarr(nb_pix) ;index in the unique list for each original theta
	index_phi[srt_phi[0: list_phi_uq[0]]]= 0
	for k=1, nb_phi_uq-1 do index_phi[srt_phi[list_phi_uq[k-1]+1: list_phi_uq[k]]]=k	

	;Recursive computation of Rm(l,theta)
	large_exp=40. ;precision of 2^-40
	fsmall=2d^(-large_exp) ;we will factor 
	fbig=2d^(large_exp)
	minscale=-15.
	cf=dblarr(15)
	for k=0,14 do cf[k]=2d^((k+minscale)*large_exp) ;scales are between 2^(-600) and 2^(400)

	cos_theta=cos(double(val_theta_uq))
	sin_theta=sin(double(val_theta_uq))
	alm_factor_theta=dblarr((lmax+1)*(lmax+2)/2,nb_theta_uq)
	for kt=0, nb_theta_uq -1 do begin
		lsin=alog(sin_theta[kt])
		;Go through all ms	
		for km=0,lmax do begin
			lfactor=logval[km]+double(km)*lsin
			;follows what has been done for healpix
			scale=floor(lfactor/large_exp/alog(2d))-minscale
			if(scale lt 0) then corfact= 0 else corfact=cf[scale] ;factorize cf[scale]
			Rmm=exp(lfactor-(scale+minscale)*large_exp*alog(2.)) ;factor the necessary corfact in Rmm
			if (km mod 2 eq 1) then Rmm=-Rmm			
			prev2=0d
			prev1=Rmm
			curl=km
			cond_true=1
			;discard all values < small_nb
			while(cond_true) do begin
				if(abs(Rmm*corfact)> small_nb) then break ;We start to reach the necessary precision
				curl=curl+1
				if(curl gt lmax) then break
				Rmm=prev1*Fm[curl-1,km,0]*cos_theta[kt]-prev2*Fm[curl-1,km,1] ;recursion
				prev2=prev1
				prev1=Rmm
				if(abs(Rmm*corfact)> small_nb) then begin ;We start to reach the necessary precision
					cond_true=0
					break
				endif
				while(abs(prev1) >fbig) do begin ;we need to increase the scale
					prev2=prev2*fsmall
					prev1=prev1*fsmall
					scale=scale+1
					if(scale lt 0) then corfact= 0 else corfact=cf[scale] ;factorize cf[scale]
				endwhile
			endwhile
			if(curl le lmax) then begin
				;print,"for sin_theta =",strcompress(string(sin_theta[kt]),/remove_all),", start at l,m=",curl,km
				prev1=prev1*corfact
				prev2=prev2*corfact
				index_l=ulong(curl) * ulong(curl +1) / 2ul + km
				alm_factor_theta[index_l,kt]=prev1
				for kl2=curl+1,lmax do begin
					index_l=ulong(kl2) * ulong(kl2+1) / 2ul + km
					Rmm=prev1*Fm[kl2-1,km,0]*cos_theta[kt]-prev2*Fm[kl2-1,km,1] ;recursion
					alm_factor_theta[index_l,kt]=Rmm
					prev2=prev1
					prev1=Rmm				
				endfor
			endif
		endfor				
	endfor
	if keyword_set(weightring) then begin
		spawn,"echo $HEALPIX",HEALPIX_DIR
		ring_weights=1d + (mrdfits(HEALPiX_DIR+"/data/weight_ring_n"+string(nside,format='(I05)')+".fits",1)).TEMPERATURE_WEIGHTS
		nrings=4*Nside-1
		theta_ring=dblarr(nrings)
		for k=0,Nside-2 do theta_ring[k]=acos(1-(k+1d)^2d/(3d * nside*nside))
		for k=Nside-1,2*Nside-1 do theta_ring[k]=acos(4d/3d - 2d*(k+1d)/(3d * nside))
		for k=2*Nside,nrings-1 do theta_ring[k]=!dpi -theta_ring[nrings-1-k] 
		
		for kt=0,nb_theta_uq -1 do begin
			ring_nb=where(abs(val_theta_uq[kt]-theta_ring) lt 1e-4,cnt)
			if(cnt lt 0) then begin
				print,"Cannot find which ring the point belongs to "
			endif else begin
				if(ring_nb gt 2*Nside-1)then ring_nb=nrings-1-ring_nb
;				print,kt,val_theta_uq[kt],ring_nb,ring_weights[ring_nb]
				for ka=0,(lmax+1)*(lmax+2)/2-1 do alm_factor_theta[ka,kt]=alm_factor_theta[ka,kt]* ring_weights[ring_nb]
			endelse
		endfor
	endif

	;Compute the exp(imphi)	
	alm_factor_phi=dblarr((lmax+1), nb_phi_uq,2) ;[m,phi,0:1] real
	if keyword_set(recursive_phi) then begin	;recursive computation of exp(imphi)
		alm_factor_phi[0,*,0]=1 ;case m=0
		cos_phi=cos(double(val_phi_uq))
		sin_phi=sin(double(val_phi_uq))

		alm_factor_phi[1,*,0]=cos_phi ;case m=1, real
		alm_factor_phi[1,*,1]=sin_phi ;case m=1, imag

		for kp=0, nb_phi_uq -1 do begin
			prev1_c=cos_phi[kp]
			prev2_c=1.
			prev1_s=sin_phi[kp]
			prev2_s=0.
			for km=2,lmax do begin
				alm_factor_phi[km,kp,0]=2d * cos_phi[kp] * prev1_c - prev2_c
				alm_factor_phi[km,kp,1]=2d * cos_phi[kp] * prev1_s - prev2_s
				prev2_c=prev1_c
				prev1_c=alm_factor_phi[km,kp,0]
				prev2_s=prev1_s
				prev1_s=alm_factor_phi[km,kp,1]
			endfor
		endfor	
	endif else begin
		for kp=0, nb_phi_uq -1 do for km=0,lmax do begin
			alm_factor_phi[km,kp,0]=cos(km * val_phi_uq[kp])
			alm_factor_phi[km,kp,1]=sin(km * val_phi_uq[kp])
		endfor
	endelse
	alm_factor_phi[1:lmax,*,*]=alm_factor_phi[1:lmax,*,*]*sqrt(2d) ;real spherical harmonics
	if keyword_set(normalize) then alm_factor_phi= alm_factor_phi*sqrt(fourpi/npixtot)

	;Now we can build the spherical harmonic matrix
	for kp=0, nb_pix -1 do begin
		for kl=0,lmax do begin
			offset_l=ulong(kl)^2.
			offset_old=kl*(kl+1)/2
			alm_matrix[offset_l,kp]= alm_factor_theta[offset_old,index_theta[kp]]*alm_factor_phi[0,index_phi[kp],0];cosine
			for km=1,kl do begin
				offset_m=offset_l+2*km-1
				alm_matrix[offset_m,kp]= alm_factor_theta[offset_old +km,index_theta[kp]]*alm_factor_phi[km,index_phi[kp],0];cosine
				alm_matrix[offset_m+1,kp]= alm_factor_theta[offset_old +km,index_theta[kp]]*alm_factor_phi[km,index_phi[kp],1];sine
			endfor
		endfor
	endfor
	return,alm_matrix
end

function compute_masked_spherical_harmonics,mask,lmax
	;Compute spherical harmonics up to l=lmax sampled outside of the max
	npix=(size(mask,/dim))[0]
	nside=sqrt(npix / 12d)
	lst_unmask=where(mask ne 0, cnt_unmasked)
	if(cnt_unmasked gt 0) then begin
		coords=dblarr(cnt_unmasked ,2)
		el=0ul
		while (el lt cnt_unmasked) do begin	
			pix2ang_nest,nside, lst_unmask[el],theta,phi 
			coords[el,*]=[theta,phi] 
			el=el+1ul
		endwhile
		alm_matrix_mask=compute_real_spherical_harmonics(nside,lmax,coords,/normalize)
	endif else begin
		print,"ERROR: No point outside of the mask"
		alm_matrix_mask=-1
	endelse
	
	return,	alm_matrix_mask

end

;===================================================================================================


function mrs_gls_inpainting, Dat, Mask, Niter=Niter, lmax=lmax, cl=cl,  AlmRes= AlmRes,  verb=verb

if not keyword_set(CL) then begin
   print, "Error: Cl keyword is required."
   InpMapGrad = -1
   goto, DONE
end

if not keyword_set(lmax) then lmax=10
if not keyword_set(Niter) then Niter =100
 
analysis_lmax = lmax
nside = gnside(Dat)
lmaxMat = 3 * analysis_lmax
almmatrix_cut=compute_masked_spherical_harmonics(Mask , lmaxMat)
ClD3 = cl[0:lmaxMat]
    
ind = where(Mask EQ 1)
PixOk = Dat[ind]
SigmaNoiseMap = Dat * 0 + 0.001
 
fake_m=dblarr(ulong(nside)*ulong(nside)*12ul)+1.
almmatrix=compute_masked_spherical_harmonics(fake_m , analysis_lmax)
 
ak_mstep_cut = mstep_first_order_gradient_descent_GLS_inpainting(PixOk,  almmatrix_cut, 2,  analysis_lmax, sigma_noise=SigmaNoiseMap[ind], pspec_alm=ClD3,nit=500)
InpMapGrad = reform(almmatrix##ak_mstep_cut)

real2complex_sph_list, ak_mstep_cut,  agrad, /tab
mrs_almtrans, Dat, a, /tab, lmax=analysis_lmax,/norm ;BEWARE we use NORMALIZED spherical harmonics
a.alm = agrad
a.alm = a.alm / a.NormVal
a.norm = 0
AlmRes = a 

DONE:

 return, InpMapGrad
end

