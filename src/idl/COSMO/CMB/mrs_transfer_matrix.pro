; NAME:
;        MRS_TRANSFER_MATRIX
;
; PURPOSE:
;	Computes the radiation transfer matrix using CLASS. The CMB
;	power spectrum computed with this matrix is in units of
;	l*(l+1)/2pi uK^2
;
; CALLING:
;       mrs_transfer_matrix, T, k, l, parameters=parameters, TCMB=TCMB, beam_FWHM=fwhm, healpix_NSIDE=nside, matmask=matmask, kmin=kmin, kmax=kmax, nklog=nklog, lmax=lmax
;    
; OUTPUTS:
;       T: IDL 2D array = Radiation transfer matrix
;       k: IDL 1D array = Wave number sampling of the primordial power
;           spectrum in Mpc^-1
;       l: IDL 1D array = Multipole sampling of the CMB power spectrum
;
; INPUT KEYWORDS:
;       parameters:    string = Path to the CLASS input file (.ini) for custom parameters, default is none, Planck PR1 best fit parameters are used
;       TCMB:          float  = Temperature of the CMB for normalisation of the power spectrum, default is 2.7255e6 uK
;       beam_FWHM:     float  = FWHM of the instrumental or effective beam to apply to the CMB, default is none
;       healpix_NSIDE: float  = Healpix resolution to apply the corresponding window size to the power spectrum, default is none
;       matmask:       string = Path to the Master matrix fits file computed with mrs_matmask, default none
;       kmin:          float  = Minimum scale for the primordial power spectrum sampling, default 1e-4 Mpc^-1
;       kmax:          float  = Maximum scale for the primordial power spectrum sampling, default 0.3 Mpc^-1
;       nklog:         int    = Number of points in the primordial power spectrum sampling, default is 512
;       lmax:          int    = Maximum multipole l computed by the transfer matrix 
;
; EXTERNAL CALLS:
;       cmbclass
;
; EXAMPLE:
;       mrs_transfer_matrix,t,k,l,beam_FWHM=5, healpix_NSIDE=2048, matmask='masterMatrix.fits'
;       pk = 2.215084e-9*(k/0.05)^(0.96235340-1)
;       cl = T ## pk 
;         
; HISTORY:
;	Written:  Francois Lanusse, May 2013
;-

pro mrs_transfer_matrix, T, kl, l, parameters=parameters, TCMB=TCMB, beam_FWHM=fwhm, healpix_NSIDE=nside, matmask=matmask, kmin=kmin, kmax=kmax, nklog=nklog, lmax=lmax

if not keyword_set(TCMB) then TCMB=2.7255e6

if keyword_set(parameters) then spawn, 'cmbclass -I'+parameters+' -T xx' else spawn, 'cmbclass -p -T xx'

if not keyword_set(lmax) then lmax=3994

; Extracting lensing contribution
readcol,'cl.dat',l0,cl0
readcol,'cl_lensed.dat',l1,cl1
lensing_contrib = cl1/cl0

k = mrdfits('xx_transmat_ks.fits')
l = mrdfits('xx_transmat_ls.fits')
trans = mrdfits('xx_transmat.fits')

delete, 'xx_transmat_ks.fits'
delete, 'xx_transmat_ls.fits'
delete, 'xx_transmat.fits'

T =  4.0*!pi*transpose(double(trans))^2

nk = size(k,/n_elem)
nl= size(l,/n_elem)

nl = min([nl, size(l0,/n_elem),size(l1,/n_elem),lmax]) 

T = T[*,0:nl-1]
lensing_contrib = lensing_contrib[0:nl-1]

dk = dblarr(nk)
for i = 1,nk-2 do dk(i) = k(i+1)-k(i-1)
dk = dk*0.5
dk[0]    = k[1]-k[0]
dk[nk-1] = k[nk-1] - k[nk-2]

for ll=0,nl-1 do T[*,ll]  *= dk/k

; Apply lensing
for kk=0,nk -1 do T[kk,*] *= lensing_contrib

; Applying beam and healpix window
beam = l*0+1.0

if keyword_set(fwhm) then begin
   b = getbeam(fwhm=fwhm)
   beam *= b[2:*]
endif

if keyword_set(nside) then begin
   hpw = healpixwindow(nside)
   beam *= hpw[2:*,0]
endif

for ll=0,nl-1 do T[*,ll] *= beam[ll]^2

; Applying master Matrix 
if keyword_set(matmask) then begin
   mat = rim(matmask)

   ;; Remove monopole and dipole contributions
   s = size(mat) 
   lmax = min([nl, s[1]-2])
   mat = double(mat[2:lmax+1,2:lmax+1])

   ;; Apply the matmask to the transfer matrix
   T = transpose(mat) ## T[*,0:lmax-1]
   l = l[0:lmax-1]
endif


; Apply inerpolation matrix to fit the k grid
if not keyword_set(kmin) then kmin = 1e-4
if not keyword_set(kmax) then kmax = 3e-1
if not keyword_set(nklog) then nklog=512

Nklog = 512
kl =  kmin*exp( (1.0/(float(Nklog-1)))*alog(kmax/kmin)*findgen(Nklog) )
interpMat = dblarr(Nklog,nk)
for i=0,nk-1 do begin
   if k[i] le kmin then interpMat[0,i] = 1 else if k[i] ge kmax then interpMat[nklog-1,i] = 1 else begin
      indmin = max(where(kl lt k[i]))>0
      dmin = k[i] - kl[indmin]
      dmax = kl[indmin +1] - k[i]
      interpMat[indmin,i] = dmax/(dmin+dmax)
      interpMat[indmin+1,i] = dmin/(dmin+dmax)
   endelse
endfor

l = l[0:nl-1]

for kk=0,nk -1 do T[kk,*] *= TCMB^2*l*(l+1.0)/(2.0*!pi)

T = T ## interpmat

end
