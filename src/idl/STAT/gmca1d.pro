;+
; NAME:
;        GMCA1D
;
; PURPOSE:
;	 Blind source separation using the fast-GMCA algorithm
;
; CALLING:
;     GMCA1D,  Data, NbrSources, Sources, MixMat, NbrIter=NbrIter,  COL_CMB=col_cmb,  NbrScale=NbrScale
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
; 
; EXTERNAL CALLS:
;         bwt01_lift (written by N. Aghanim and O. Forni)
;
; HISTORY:
;	Written: Jerome Bobin, 2008
;-
;-----------------------------------------------------------------

function  gmca1d,  Data,  ns,  AA=AA,   NbrIter= NbrIter,   NbrScale=NbrScale, FirstThreshold=FirstThreshold, LastThreshold=LastThreshold, Norm=Norm
COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
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
Opt='-t16 -L -n' + strc(NbrScale)

nx = size(x)
n = nx(1) ; side
nc = nx(2) ; nb. de comp.

if   keyword_set(Norm) then  for i=0,nc-1 do X[*,i] = (X[*,i]-mean(X[*,i]))  ; / sigma( X[*,i])

;#### Donnees

AA = randomn(seed,nc,ns) 

; Normalize the matrix
for pp=0,ns-1 do AA[*,pp] = AA[*,pp]/ norm(AA(*,pp))	

mr1d_trans, x[*,0], w, opt=Opt

coefs1 = w.COEF

nc1 = double(n)    
xcoeffs = dblarr(nc,nc1)
xcoeffs[*]=0.0
coeffs = dblarr(ns,nc1)
coeffs[*]=0.0
		    
xcoeffs[0,*] = reform(coefs1,1,nc1)
		    
for tt=1,nc-1 do begin
     mr1d_trans, x[*,tt], w, opt=Opt    
     coefs1 = w.COEF
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
x_r = dblarr(nc,n)
for pp=0,nc-1 do x_r[pp,*] = reform(X[*,pp],1,n)

invA = invert(transpose(AA)#AA,/DOUBLE)#transpose(AA)
s_r = invA#x_r
s = dblarr(n,ns)
s[*,*]=0.0
for ll=0,ns-1 do s[*,ll] = reform(s_r[ll,*],n)

DONE:
return, s
end


