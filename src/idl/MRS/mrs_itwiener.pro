;+
; NAME:
;        mrs_ItWiener
;
; PURPOSE:
;	Iterative Wiener or Cole filtering of an image on the sphere (Healpix pixel NESTED representation)
;
; CALLING:
;
;	mrs_ItWiener,Imag,Powspec,fImag,Noise=Noise,fNoise = fNoise,niter=niter,cole=cole,StartWiener=StartWiener,Filter=Filter,VarNoise=VarNoise,NbrScale=NbrScale,BS=BS,Frg=Frg,RMSmoothing = RMSmoothing,tol = tol,verbose=verbose
;    
; INPUT:
;     	Imag -- IDL array of healpix map: Input image be filtered 
;	  	Powspec -- theoretical power spectrum of the image
;
; OUTPUT:
;     	fImag -- filtered image
;
; INPUT KEYWORDS:
;		Imag -- IDL array of healpix map: Input image be filtered
;	    Noise -- IDL -- array of healpix map: Noise realization from which the noise statistics are measured
;		niter -- scalar: number of iterations - default is 100
;		cole -- if set, applies a Cole filter instead of the Wiener filter - default is no
;		StartWiener -- if set, the algorithm is initialized with the global Wiener filter (i.e. in the spherical harmonics)
;		NbrScale -- scalar: number of wavelet scales used for analysis - default is 4
;		BS -- scalar: minimal patch size used to compute the noise variance - default is 16
;		Frg -- if set, also estimates the contribution coming from extra foreground residuals - default is no
;		RMSmoothing --- if set, applies a smoothing of the estimated noise variance maps - default is no
;		tol -- convergence precision - default is 1e-6
;		verbose -- scalar, verbose mode - default is no 
;
; INPUT/OUTPUT:
;  		filter - Wiener filter
;       VarNoise -- IDL array of healpix map per wavelet scale - Noise variance per wavelet scale - if not set, it is estimated from the noise realization Noise
;		lmax : int = maximum l value in the Spherical Harmonic Space (Healpix)
;
; OUTPUT KEYWORDS:
;		fNoise -- IDL array of healpix map: filtered noise realization
;
; EXAMPLE:
;       Filter an image with 5 scales. The result is stored in fImag
;               mrs_ItWiener,Imag,Powspec,fImag,Noise=Noise,NbrScale=5
;         
; HISTORY:
;	Written: Jérôme Bobin, 2012
;	February, 2012 File creation
;
;---------------------------------------------------------------------------------------------------------------------------------------------------------

;###################################################################################

function convol_maps,map,beam,lmax=lmax,enside=enside

lmax = n_elements(beam)-1
mrs_almtrans,map,alm,/tab,lmax=lmax
for l=1,Alm.lmax do  alm.alm[l,*,*] = alm.alm[l,*,*] *  beam[l]
if keyword_set(enside) then alm.nside=enside
mrs_almrec, alm, rec
return,rec

end


function WT_Apply,r,alpha,wfilter,NbrScale=NbrScale,BS=BS,lmax=lmax

nside=double(gnside(r))

;--- transform the maps

mrs_wttrans,r,out,NbrScale=NbrScale,lmax=lmax
mapc = out.coef

BSl = BS

;--- Wiener for all scales

for s=0,NbrScale-2. do begin ;--- We do not process the coarse scale
	
	c = mapc[*,s]
	cw = wfilter[*,s]
		
	nb_bin = double(nside)/double(BSl)

	for facenum=0,11 do begin ;--scan all faces
	
		f = get_one_face(c,facenum)
		fw = get_one_face(cw,facenum)

		for bx = 0,nb_bin-1. do begin ;--- scan bx
			for by = 0,nb_bin-1. do begin
			
				p = f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]
				w = mean(fw[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])
				
				p = w*p ;--- we should record w

				f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = p

			endfor
		endfor ;--- patch loops
	
		put_one_face,c,f,facenum

	endfor ;--- face loop
	
	mapc[*,s] = alpha[s]*c

	BSl = 2.*BSl
	
endfor ;--- Scale loop

out.coef = mapc
mrs_wtrec,out,wr

return,wr

end

;###################################################################################


function mrs_bandpass_sigma, Cl, Filter,  lmax=lmax, npix=npix, NormVal= NormVal

if not keyword_set(NormVal) then  NormVal = 1.
vs = size(filter)
lm = vs[1]
vs = size(Cl)
lm1 = vs[1]
if not keyword_set(lmax) then lmax = min([lm, lm1]) -1
 
 CoefN = NormVal* NormVal
 IndL = lindgen(LMAX+1)
IndL = 2. * IndL + 1.
Sig = total( IndL * Cl[0:Lmax]  * Filter[0:Lmax]^2 / CoefN)
Sig = Sig / (4. * !DPI)
return, sqrt(Sig)
end

;###################################################################################

function GetNoiseVariance,cn,BS,nside


NbrScale = (size(cn))[2]
noise = cn
BSl = BS

for s=0,NbrScale-1 do begin

	temp = cn[*,s]

	for facenum=0,11 do begin ;--scan all faces
	
		fn = get_one_face(temp,facenum)
		gn = fn
		
		nb_bin = double(nside)/double(BSl)
				
		for bx = 0,nb_bin-1. do begin ;--- scan bx
			for by = 0,nb_bin-1. do begin
			
				pn = fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]	
				gn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = sigma(pn)^2. ;--- Attention aux normalisations !!
				
			endfor
				
		endfor ;--- patch loops
	
		put_one_face,temp,gn,facenum

	endfor ;--- face loop

	noise[*,s] = temp
	
	BSl = 2.*BSl

endfor

return, noise

end

;###################################################################################

function GetFrgVariance,c,cn,BS,nside,Cl,T

NbrScale = (size(cn))[2]
noise = cn

Tall = 0.*T

BSl = BS

for s=0,NbrScale-2 do begin

	tempn = cn[*,s]	
	temp = c[*,s]
	
	CMBVar = mrs_bandpass_sigma(Cl, T[*,s])^2.

	for facenum=0,11 do begin ;--scan all faces
	
		f = get_one_face(temp,facenum)
		fn = get_one_face(tempn,facenum)
		gn = fn
		
		nb_bin = double(nside)/double(BSl)
				
		for bx = 0,nb_bin-1. do begin ;--- scan bx
			for by = 0,nb_bin-1. do begin
			
				p = f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]
				pn = fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]	
				gn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = max([sigma(pn)^2.,sigma(p)^2. - CMBVar]) ;--- Attention aux normalisations !!
				
			endfor
				
		endfor ;--- patch loops
	
		put_one_face,temp,gn,facenum

	endfor ;--- face loop

	noise[*,s] = temp
	
	Tall = Tall + T[*,s]
	
	BSl = 2.*BSl

endfor

s = NbrScale-1

CMBVar = mrs_bandpass_sigma(Cl, 1-Tall)^2.

tempn = cn[*,s]	
	temp = c[*,s]
	
	for facenum=0,11 do begin ;--scan all faces
	
		f = get_one_face(temp,facenum)
		fn = get_one_face(tempn,facenum)
		gn = fn
		
		nb_bin = double(nside)/double(BSl)
				
		for bx = 0,nb_bin-1. do begin ;--- scan bx
			for by = 0,nb_bin-1. do begin
			
				p = f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]
				pn = fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]	
				gn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = max([sigma(pn)^2.,sigma(p)^2. - CMBVar]) ;--- Attention aux normalisations !!
				
			endfor
				
		endfor ;--- patch loops
	
		put_one_face,temp,gn,facenum

	endfor ;--- face loop

	noise[*,s] = temp

return, noise

end

;###################################################################################

pro mrs_ItWiener,Imag,Powspec,fImag,Noise=Noise,fNoise=fNoise,niter=niter,cole=cole,StartWiener=StartWiener,Filter=Filter,VarNoise=VarNoise,NbrScale=NbrScale,BS=BS,Frg=Frg,RMSmoothing = RMSmoothing,tol = tol,verbose=verbose

;--- Initialization

if not keyword_set(tol) then tol = 1e-6
if not keyword_set(niter) then niter = 100.
if not keyword_set(NbrScale) then NbrScale = 4


nside=double(gnside(Imag))
BSl = BS
lmax = n_elements(Powspec)-1.

;--- Get the noise covariance matrix

x0 = Imag

if keyword_set(Noise) then begin

	if keyword_set(verbose) then print,'1- Get the noise covariance matrix'

	mrs_wttrans,Noise,out,NbrScale=NbrScale,lmax=floor(lmax)
	noisec = out.coef
	
	xn0 = Noise
	
	if keyword_set(StartWiener) then begin
	
		if keyword_set(verbose) then print,'Starting point is the Wiener solution'
		
		pn = mrs_powspec(Noise) ;--- It should be filtered
		w = Powspec/(Powspec + pn)*getidealbeam(0.*dblarr(lmax+1),lmin=2800,lmax=3200)
	
		x0 = convol_maps(Imag,w) 
		xn0 = convol_maps(Noise,w)
	
	endif
	
	T = out.TabPsi
	Tall = reform(0.*T[*,0])
	
	if keyword_set(Frg) then begin
	
		if keyword_set(verbose) then print,'Include foreground in the noise covariance matrix'
			
		mrs_wttrans,Imag,out,NbrScale=NbrScale,lmax=floor(lmax)
		c = out.coef

		VarNoise = GetFrgVariance(c,noisec,BSl,nside,Powspec,T)	
		c = 0
		
	endif else begin
	
		VarNoise = GetNoiseVariance(noisec,BSl,nside)
	
	endelse
	
	noisec = 0

endif else begin

	mrs_wttrans,Imag,out,NbrScale=NbrScale,lmax=floor(lmax)
	T = out.TabPsi
	Tall = reform(0.*T[*,0])

	if not keyword_set(VarNoise) then begin
		
		print,'Be careful, you should provide at least a noise realization or the noise variance per wavelet scale . . . '
		fImag = 0
		goto,DONE
		
	endif

endelse

NormVal = double(n_elements(Imag)/(4.*!DPI))

BSl = BS
SigmaI = 1./VarNoise

;--- Initialization

fx = x0
if keyword_set(Noise) then fxn = xn0

filter = 0.*dblarr(n_elements(Powspec))
alpha = 0.*dblarr(NbrScale)

BSl = BS

for s=0,NbrScale-2 do begin

	if keyword_set(RMSmoothing) then begin
	
		if keyword_set(verbose) then print,'Variance map smoothing'
		
		enside = double(nside)/BSl
		temp = mrs_resize(SigmaI[*,s],nside=enside)   ;--- "zero padding" in sph. harmonics
		mrs_almtrans,temp,alm,/tab
		alm.nside = nside
		mrs_almrec,alm,temp
		SigmaI[*,s] = reform(temp)
		
	endif

	BSl = 2.*BS

	alpha[s] = 1./max(SigmaI[*,s])
	if keyword_set(Cole) then filter = filter + sqrt(Powspec/(Powspec + alpha[s]/NormVal))*T[*,s] else filter = filter + Powspec/(Powspec + alpha[s]/NormVal)*T[*,s]
	Tall = Tall + T[*,s]

endfor

s = NbrScale-1
alpha[s] = 1./max(SigmaI[*,s])

if keyword_set(Cole) then filter = filter + sqrt(Powspec/(Powspec + alpha[s]/NormVal))*(1. - Tall) else filter = filter + Powspec/(Powspec + alpha[s]/NormVal)*(1. - Tall)

;--- Regularizing the filter - we do not filter the first 20 l et filter out large values of l

;if not keyword_set(filter) then begin
filter[0:20] = 1. 
filter = filter*getidealbeam(0.*filter,lmin=floor(0.8*lmax),lmax=floor(0.95*lmax),/tozero)
;endif

out = 0

fy = fx
if keyword_set(Noise) then fyn = fxn

;--- Main loop

lambda = 1.
dlambda = lambda/niter
tk = 1.

srold = 1.
sx = sigma(Imag)
if keyword_set(Noise) then sxn = sigma(Noise)

for ll=0,niter-1 do begin

	;--- Compute the gradient of f2 - STEP 1
	
	r = (Imag - fy)
	if keyword_set(Noise) then rn = (Noise - fyn)
	
	sr = sigma(r)/sx;
	if keyword_set(Noise) then srn = sigma(rn)/sxn;
	
	com = 'It. # '+strn(ll+1)+' - Relative residual norm : '+strn(sr)+ ' / Variation : '+strn(abs(sr-srold))
	if keyword_set(verbose) then print,com

	;--- Compute the gradient of f2 - STEP 2
	
	fy = fy + WT_Apply(r,alpha,SigmaI,NbrScale=NbrScale,BS=BS)
	if keyword_set(Noise) then fyn = fyn + WT_Apply(rn,alpha,SigmaI,NbrScale=NbrScale,BS=BS)
	
	;--- Applying the prox of f1 (filtering)

	fxp = convol_maps(fy,filter)
	if keyword_set(Noise) then fxnp = convol_maps(fyn,filter)
	
	;--- Multistep
	
	tkp = 0.5*(1 + sqrt(1+4*tk^2));
    beta = (tk - 1)/tkp;
    tk = tkp;
    
   	fy = fxp + beta*(fxp - fx);
   	if keyword_set(Noise) then fyn = fxnp + beta*(fxnp - fxn);
    
    fx = fxp
    if keyword_set(Noise) then fxn = fxnp

	;--- Stopping criterion

	if ll gt 10. then if abs(sr - srold) lt tol then goto, DONE0
	
	srold = sr

endfor

DONE0:

fImag = fx
if keyword_set(Noise) then fNoise = fxn

DONE:

end
