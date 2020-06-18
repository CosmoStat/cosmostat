;+
; NAME:
;        mrs_pwtrec
;
; PURPOSE:
;	Compute the inverse pyramidal wavelet transform on the sphere.
;
; CALLING:
;
;     mrs_pwtrec, Trans, Rec, filter=filter        
;    
; INPUT:
;     Trans -- IDL structures with the following fields:  
;                    NbrScale : int = number of scales 
;                       nside : int = Healpix nside parameter, only present with healpix input image
;                        npix : long = Number of pixels
;                      Scale1 : finest scale (highest frequencies). A IDL array of healpix map or IDL structure of a Glesp map
;                      Scale2 : Second scale    
;                      Scalej : j th scale
;                        ...
;                      ScaleJ : with J = NbrScale, coarsest resolution
;                       lmax  : int = nlmax parameter at the first scale
;		   Tab_lmax[NbrScale] : int array Tab_lmax[j] = lmax at scale j+1, j=0...NbrScale-1
;		  Tab_nside[NbrScale] : int array Tab_nside[j] = nside parameter of the scale j+1, j=0..NbrScale-1 (Healpix input map)
;														 nx number of rings, Glesp parameter of the scale j+1, j=0..NbrScale-1 (Glesp input map)
;					 UseGLESP : int = 1 if the input image was in Glesp format or keyword Healpix_with_Glesp used, otherwise 0
;				    MeyerWave : int = 1 if the keyword MeyerWave used, otherwise 0
;					  DifInSH : int = 1 if the keyword DifInSH used, otherwise 0
;					  	   nx : int = number of rings at the first scale, Glesp parameter, only present with glesp input image
;						   np : int = max number of pixel on a ring at the first scale, Glesp parameter, only present with glesp input image
;
; OUTPUT:
;     Imag -- IDL array of healpix map or IDL structure of a Glesp map: reconstructed image from the wavelet coefficients   
;
; KEYWORDS:
;      filter : Use filters for the reconstructions. If this keyword is not set, the reconstructed image
;               is obtained by a simple addition of all wavelet scales. When computing direct transform
;				(mrs_pwttrans function), if the keywords DifInSH or MeyerWave were setted, filter is 
;				automatically used.
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;   	synfast (healpix software)
;   	alm_product2 (idl)
;   	compute_g (idl)
;   	compute_h (idl)
;   	compute_gtilde (idl)
;   	compute_htilde (idl)
;
; EXAMPLE:
;       Compute the inverse pyramidal wavelet transform:
;        The result is stored in Imag 
;               mrs_pwtrec, Trans, Imag 
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2005
;	February, 2005 File creation
;------------------------------------------------------------------------


pro mrs_glesp_pwtrec, Trans, Imag, filter=filter 

 
if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_pwtrec, Trans, Imag, filter=filter'
        goto, DONE
        end

NbrScale = Trans.NbrScale
npix =  Trans.npix
nx = Trans.nx
np = Trans.np
nlmax = Trans.lmax

if Trans.MeyerWave EQ 1 then filter=1
if Trans.DifInSH   EQ 1 then filter=1
 
;nside = Trans.nside
NbrScale = Trans.NbrScale

lmax = Trans.lmax

ech = 2^(NbrScale-2)
 
if keyword_set(filter) then begin

 HHScale = mrs_wtget(Trans, 0)
 ;  LScale = mrs_wtget(Trans, NbrScale-1) 
   ;nside = Trans.Tab_nside[NbrScale-1]
   ;alm_trans, LScale, ALM_LScale, nlmax=lmax, nx=nx,np=np, index=index
   ;LScale.lmax = (Lscale.nx-1)/2
   mrs_almtrans,HHScale,ALM_HHScale,lmax=(HHscale.nx-1)/2;nlmax
   
   for j =(NbrScale-1),1,-1 do begin
       LScale = mrs_wtget(Trans, j)
       mrs_almtrans,LScale,alm_l
       ALM_HHScale.alm = ALM_HHScale.alm + alm_l.alm
    endfor
    mrs_almrec, ALM_HHScale,Imag, nx=nx,np=np;,nlmax=nlmax



endif else begin 
    Imag = mrs_wtget(Trans, NbrScale-1)  
    for j =(NbrScale-2),0,-1 do begin
        nlmax = Trans.Tab_lmax[j]
        nx = Trans.Tab_nside[j]
        ; print,j,nlmax,nx
	;alm_trans, Imag, ALM_LScale, nlmax=nlmax/2, index=index2
 	mrs_almtrans, Imag, ALM_LScale, lmax=(imag.nx-1)/2;nlmax/2-1;, index=index2
 	
	;alm_itrans, ALM_LScale, LReScale, nx=nx,np=np*2 nlmax=nlmax/2, index=index2
        WScale = mrs_wtget(Trans, j)
	mrs_almrec,ALM_LScale,LReScale,nx=wscale.nx,np=wscale.np
	
	;hs,wscale
	;hs,lrescale
 	Imag = LreScale
	Imag.t_sky = LReScale.t_sky + WScale.t_sky
	;imag.nx = LReScale.nx
    endfor
end

 
DONE:

END


;============================================================================================

pro mrs_pwtrec, Trans, Imag, filter=filter 

COMMON MR1ENV
 
if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_pwtrec, Trans, Imag, filter=filter'
        goto, DONE
        end

if Trans.UseGLESP EQ 1 then mrs_glesp_pwtrec, Trans, Imag, filter=filter $
else BEGIN  ; HEALPIX

NbrScale = Trans.NbrScale
npix =  Trans.npix
nside = Trans.nside
nlmax = Trans.lmax

if Trans.MeyerWave EQ 1 then filter=1
if Trans.DifInSH   EQ 1 then filter=1
 
nside = Trans.nside
NbrScale = Trans.NbrScale

lmax = Trans.lmax

ech = 2^(NbrScale-2)
 
if keyword_set(filter) then begin
	LScale = mrs_wtget(Trans, NbrScale-1) 
	nside = Trans.Tab_nside[NbrScale-1]
	if keyword_set(HealpixCXX) then begin
		alm_cxxtrans, LScale, ALM_LScale, nlmax=lmax, nside=nside, index=index
	end else begin
		alm_trans, LScale, ALM_LScale, nlmax=lmax, nside=nside, index=index
	end
	
	for j =(NbrScale-2),0,-1 do begin
		nlmax = Trans.Tab_lmax[j]
		if Trans.MeyerWave EQ 0 then begin
			compute_gtilde, lmax, ech, gtilde
			compute_htilde, lmax, ech, htilde
		end else hgmey, lmax, j, htilde, gtilde, dif=dif
		
		alm_product2, ALM_LScale, htilde, ConvH
		WScale = mrs_wtget(Trans, j)
      
		index = 0
		if keyword_set(HealpixCXX) then begin
			alm_cxxtrans, WScale, ALM_WScale, nside=nside, nlmax=lmax, index=index
		end else begin
			alm_trans, WScale, ALM_WScale, nside=nside, nlmax=lmax, index=index
		end
		alm_product2, ALM_WScale, gtilde, ConvG
		ALM_LScale = ConvH + ConvG
		ech = ech/2
    endfor
    nside = Trans.nside
    if keyword_set(HealpixCXX) then begin
    	alm_cxxitrans, ALM_LScale, Imag, nside=nside, nlmax=lmax, index=index
    end else begin
    	alm_itrans, ALM_LScale, Imag, nside=nside, nlmax=lmax, index=index
    end
end else begin
    Imag = mrs_wtget(Trans, NbrScale-1)  
    for j=NbrScale-2,0,-1 do begin;	j=(NbrScale-2),0,-1
        nlmax = Trans.Tab_lmax[j]
        nside = Trans.Tab_nside[j]
        if keyword_set(HealpixCXX) then begin
			alm_cxxtrans, Imag, ALM_LScale, nlmax=nlmax/2, index=index2
 			alm_cxxitrans, ALM_LScale, LReScale, nside=nside, nlmax=nlmax/2, index=index2
        endif else begin
			alm_trans, Imag, ALM_LScale, nlmax=nlmax/2, index=index2
 			alm_itrans, ALM_LScale, LReScale, nside=nside, nlmax=nlmax/2, index=index2
		endelse
	
		WScale = mrs_wtget(Trans, j)
 		Imag = LReScale + WScale
    endfor
end

END ; HEALPIX

DONE:

END

