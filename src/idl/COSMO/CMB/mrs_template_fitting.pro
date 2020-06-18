;+
; NAME:
;        mrs_template_fitting
;
; PURPOSE:
;	It aims at solving a template fitting problem on the sphere
;
; CALLING:
;
;	mrs_template_fitting,map,template,BS=BS,NbrScale=NbrScale,alpha = alpha
;    
; INPUT:
;     	map -- IDL array of healpix map: Input image 
;	  	template -- IDL array of healpix map : Template image
;
; OUTPUT:
;     	Image obtained by fitting (i.e. regression coefficients x template)
;
; INPUT KEYWORDS:
;		NbrScale -- scalar: number of wavelet scales used for analysis - default is 4
;		BS -- scalar: minimal patch size used to compute the noise variance - default is 16
;		denoise -- scalar : if set, it does not perform fitting when noise dominates at level denoise x sigma_noise (e.g. put denoise = 3 or more)
;		mask -- IDL array : binary mask to neglect some pixel when computing the regression coefficient at large scale
;       scale_coarse -- scalar : scale at which the large scale regression coefficient should be estimated
;		positive -- binary : if set the regression coefficients are forced to be positive
;
; OUTPUT KEYWORDS:
;		alpha -- IDL array of healpix map: regression coefficients per scale
;
; EXAMPLE:
;       fImage = mrs_template_fitting(Imag,template,BS=8,NbrScale=5,alpha = alpha)
;         
; HISTORY:
;	Written: Jérôme Bobin, 2012
;	September, 2012 File creation
;
;---------------------------------------------------------------------------------------------------------------------------------------------------------




function mrs_template_fitting,map,template,BS=BS,NbrScale=NbrScale,alpha = alpha,denoise=denoise,mask=mask,scale_coarse=scale_coarse,positive=positive, WT_map=WT_map, WT_template=WT_template, TabCorrel=TabCorrel

if not keyword_set(NbrScale) then NbrScale = 4
if not keyword_set(BS) then BS = 16
if keyword_set(denoise) then begin
	print,'Fitting with thresholding'
	thrd = denoise
endif else begin
	thrd = 0.	
endelse

if keyword_set(WT_map) then cmap = WT_map.coef $
else begin
   mrs_wttrans,map, WT_map,Nbr = NbrScale
   cmap = WT_map.coef
end

if keyword_set(WT_template) then ctemp = WT_template.coef $
else begin
   mrs_wttrans, template, WT_template,Nbr = NbrScale
   ctemp = WT_template.coef
end
out = WT_template
nside= gnside(map)

BSl = BS

alpha = 0.*ctemp
TabCorrel = alpha

for ll=0,NbrScale-2 do begin  ;---- All scales

	noise_level = mad(cmap[*,ll])

	for f=0,11 do begin   ;---- Face
	
		fmap = get_one_face(cmap[*,ll],f)
		ftemp = get_one_face(ctemp[*,ll],f)
		falpha = 0.*ftemp
		fcorrel = 0.*ftemp
		nb_bin = nside/BSl
		
		for bx = 0,nb_bin-1. do begin ;--- scan bx
		for by = 0,nb_bin-1. do begin
		
			pmap = fmap[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]
			ptemp = ftemp[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]
			
			pmap = pmap - mean(pmap)
			ptemp = ptemp - mean(ptemp)
			
			ind = where(pmap gt thrd*noise_level,count_map)
			ind = where(ptemp gt thrd*noise_level,count_temp)
			
			if (count_map gt 0.) AND (count_temp gt 0.) then begin
			     falpha[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = total(pmap*ptemp)/total(ptemp*ptemp)
			     fcorrel[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = correlate(pmap, ptemp)
			     end
		endfor
		endfor ;--- patch loops

		
		temp = alpha[*,ll]
		put_one_face,temp,falpha,f
		alpha[*,ll] = temp
		
		cor = TabCorrel[*,ll] 
		put_one_face,cor, fcorrel, f
	   TabCorrel[*,ll] = cor
	endfor   ;---- Face
	
	BSl = 2.*BSl

endfor   ;---- All scales

if not keyword_set(scale_coarse) then begin
	fmap = cmap[*,NbrScale-1] - mean(cmap[*,NbrScale-1])
	ftemp = ctemp[*,NbrScale-1] - mean(ctemp[*,NbrScale-1])
endif else begin
	fmap = cmap[*,NbrScale-1-scale_coarse] - mean(cmap[*,NbrScale-1-scale_coarse])
	ftemp = ctemp[*,NbrScale-1-scale_coarse] - mean(ctemp[*,NbrScale-1-scale_coarse])
endelse

if not keyword_set(mask) then begin
	alpha[*,NbrScale-1] = total(fmap*ftemp)/total(ftemp*ftemp)
endif else begin
	ind = where(mask eq 1.)
	fmap = fmap[ind]-mean(fmap[ind])
	ftemp = ftemp[ind] - mean(ftemp[ind])
	
endelse

alpha[*,NbrScale-1] = total(fmap*ftemp)/total(ftemp*ftemp)
TabCorrel[*,NbrScale-1] = correlate(fmap, ftemp)

if keyword_set(positive) then alpha = alpha*(alpha gt 0.)

for scale = 0,NbrScale - 2 do begin 
	
	BSl = BS*(2.^scale)
	enside = nside/BSl
			
	f = reform(alpha[*,scale]) 
	filter = getidealbeam(0.*dblarr(2.*nside),lmin=0.,lmax=min([8.*enside,min([2.*nside-100.,16.*enside])]),/tozero)
			
	; fout = convol_maps(f,filter)
	mrs_convol, f, filter, fout
	
	alpha[*,scale] = fout
			
endfor

out.coef = alpha*ctemp

mrs_wtrec,out,tmap

return, tmap

end
