;+
; NAME:
;        GMCA2D
;
; PURPOSE:
;	 Blind source separation using the fast-GMCA algorithm
;
; CALLING:
;     GMCA2D,  Data, NbrSources, Sources, MixMat, NbrIter=NbrIter,  COL_CMB=col_cmb,  NbrScale=NbrScale
;
; INPUTS:
;    Data -- IDL array[*,*,0:Nx-1]  = multichannel data, Nx being the number of channels
;    NbrSources  -- int: Number of sources 
; 
; OUTPUTS:
;   Sources-- IDL array of healpix map [*,*,0:Ns-1]  = estimated sources
;   MixMat = Mixing matrix
;
; KEYWORDS:
;        NBRSCALE  -- LONG: Number of scales of the wavelet transform. Default is 4.
;        NbrIter -- LONG: Number of iterations. Default is 40.
;        COL_CMB -- IDL array [0:Nx-1]: if set, the first column of the matrix is fixed, and is equal to COL_CMB
; 
; EXTERNAL CALLS:
;         bwt01_lift (written by N. Aghanim and O. Forni)
;
; HISTORY:
;	Written: Jerome Bobin, 2008
;-
;-----------------------------------------------------------------

pro  gmca2d,  Data,  ns,  s,  AA,   NbrIter= NbrIter,  COL_CMB=col_cmb, NbrScale=NbrScale, FirstThreshold=FirstThreshold, LastThreshold=LastThreshold, Norm= Norm
COMMON MR1ENV

if N_PARAMS() LT 4  then begin 
        print, 'CALLING SEQUENCE: GMCA2D,  Data,  NbrSources, Sources, MixMat, NbrIter=NbrIter,  COL_CMB=col_cmb,  NbrScale=NbrScale'
        goto, DONE
        end

X = Data

if not keyword_set(NbrIter) then NbrIter =40
nit_max=NbrIter
if not keyword_set(FirstThreshold) then  FirstThreshold=15.
if not defined(LastThreshold) then  LastThreshold =2.
ts_max = FirstThreshold
ts_min = LastThreshold

if not keyword_set(NbrScale) then NbrScale=4
Opt='-t14 -L -n' + strc(NbrScale)

nx = size(x)
n = nx(1) ; side
nc = nx(3) ; nb. de comp.

if   keyword_set(Norm) then  for i=0,nc-1 do X[*,*,i] = (X[*,*,i]-mean(X[*,*,i])) / sigma( X[*,*,i])

;#### Donnees

AA = randomn(seed,nc,ns) 

if keyword_set(COL_CMB) then begin
            AA[*,0] = col_cmb
            AA[*,1:ns-1] = randomn(1,nc,ns-1)
endif

; Normalize the matrix
for pp=0,ns-1 do AA[*,pp] = AA[*,pp]/norm(AA(*,pp))	
if keyword_set(mr1ok) then begin				;if package mre is available
         print, 	'using mre'
         mr_transform, x[*,*,0], coefs1, opt=Opt
end else  coefs1 = bwt01_lift( x[*,*,0],  NbrScale-1)
 
nc1 = double(n)*double(n)		    
xcoeffs = dblarr(nc,nc1)
xcoeffs[*,*]=0.0
coeffs = dblarr(ns,nc1)
coeffs[*,*]=0.0
		    
xcoeffs[0,*] = reform(coefs1,1,nc1)
		    
for tt=1,nc-1 do begin
    if keyword_set(mr1ok) then  mr_transform, x[*,*,tt],coefs1,opt=Opt   $
   	else  coefs1 = bwt01_lift(x[*,*,tt],  NbrScale-1)
    xcoeffs[tt,*] = reform(coefs1,1,nc1)	          
endfor

;#### Seuil en k*sigma_ns avec k decroissant lineairement
ts = ts_max
dts = (ts_max - ts_min)/(nit_max - 1.0)
nit = 0.0

;#### Corps du programme fast_gmca
while (nit lt nit_max) do begin
                        nit = nit + 1.0
                        scoeffs = invert(transpose(AA)#AA,/DOUBLE)#transpose(AA)#xcoeffs
                        for pp =0,ns-1 do begin                                        
                                        thrd = ts * mad(scoeffs(pp,*))
                                        scoeffs[pp,*] = scoeffs[pp,*]*(abs(scoeffs[pp,*]) gt thrd)
                        endfor
                        n_v = dblarr(1,ns)
                        for ll=0,ns-1 do n_v(ll) = norm(scoeffs(ll,*))
                        if keyword_set(COL_CMB) then n_v(0) = 0.0
                        indix = where(n_v gt 1e-7)
                        if (max(n_v) gt 1e-7) then begin
                                        tA = xcoeffs#transpose(scoeffs(indix,*))#invert(scoeffs(indix,*)#transpose(scoeffs(indix,*)),/DOUBLE)
                                        AA(*,indix) = tA
                        endif
                        for pp=0,ns-1 do begin
                                   na = norm(AA(*,pp))
                                   if (na gt 1e-7) then AA(*,pp) = AA(*,pp)/na
                        endfor
                        ts = ts - dts
endwhile

;---- Final Inversion
x_r = dblarr(nc,n*n)
for pp=0,nc-1 do x_r[pp,*] = reform(X[*,*,pp],1,n*n)

invA = invert(transpose(AA)#AA,/DOUBLE)#transpose(AA)
s_r = invA#x_r
s = dblarr(n,n,ns)
s[*,*,*]=0.0
for ll=0,ns-1 do s[*,*,ll] = reform(s_r[ll,*],n,n)

if keyword_set(COL_CMB) then  s[*,*,0] = s[*,*,0]/sigma(COL_CMB)*sigma(AA[*,0])

DONE:

end


