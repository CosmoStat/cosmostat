; NAME:
;        MRS_PRISM
;
; PURPOSE:
;	Reconstructs the primordial power spectrum using the PRISM
;	algorithm.
;
; CALLING:
;       pk = mrs_prism(Cl, noise=Nl, TransferMat=mat,k=k,l=l, Opt=opt)
;
; INPUTS:
;       Cl: IDL 1D array = measured CMB auto-powerspectrum, by default
;       in units of l*(l+1)/2pi uK^2 starting at l=2
;
;       Nl: IDL 1D array = Estimated noise power spectrum [optional]
;
;       TransferMat: IDL 2D array = Radiative transfer matrix, if not
;           provided a default matrix for the Planck best fit
;           parameters will be used 
; OUTPUTS:
;       pk: IDL 1D array = Estimated primordial power spectrum
;
; INPUT KEYWORDS:
;       noise: IDL 1D array = Estimated noise power spectrum, default none
;
;       TransferMat: IDL 2D array = Radiative transfer matrix, if not
;           provided a default matrix for the Planck best fit
;           parameters will be used 
;           
;       k: IDL 1D array = Wavenumber sampling of the primordial power
;       spectrum if a custom transfer matrix is provided
;
;       l: IDL 1D array = Multipole sampling of the input Cl if a
;       custom transfer matrix is provided
;
;       opt: string = Options of the prism executable:
;         [-T type_of_filters]
;               1: Biorthogonal 7/9 filters 
;               2: Daubechies filter 4 
;               3: Biorthogonal 2/6 Haar filters 
;               4: Biorthogonal 2/10 Haar filters 
;               5: Odegard 9/7 filters 
;               6: 5/3 filter 
;               7: Battle-Lemarie filters (2 vanishing moments) 
;               8: Battle-Lemarie filters (4 vanishing moments) 
;               9: Battle-Lemarie filters (6 vanishing moments) 
;               10: User's filters 
;               11: Haar filter 
;               12: 3/5 filter 
;               13: 4/4 Linar spline filters 
;               14: Undefined sub-band filters 
;              default is Battle-Lemarie filters (4 vanishing moments)
;
;          [-i number_of_iterations]
;              Maximum number of iterations
;              default is 200.
;
;          [-I ExtIter]
;               Number of reweighting iterations 
;               default is 3
;
;          [-G RegulParam]
;               Regularization parameter 
;               default is 0.000000
;
;          [-n number_of_scales]
;              number of scales used in the multiresolution transform.
;               default is 9
;
;          [-s nsigma]
;              Thresholding at nsigma * SigmaNoise
;              default is  5.           
;
;          [-v]                                                                                                                
;              Verbose. Default is no. 
;
; EXTERNAL CALLS:
;       prism
;
; EXAMPLE:
;       ; Compute a transfer matrix with mrs_transfer_matrix
;       mrs_transfer_matrix,t,k,l,beam_FWHM=5, healpix_NSIDE=2048, matmask='masterMatrix.fits'
;       ; Reconstruct a power spectrum with
;       pk = mrs_prism(Cl, noise=Nl, TransferMat=t,k=k,l=l, opt='-v')
;         
; HISTORY:
;	Written:  Francois Lanusse, June 2014
;-

function mrs_prism,Cl, noise=Nl, TransferMat=mat,k=k,l=l, Opt=opt

if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE:  Pk = mrs_prism(Cl, noise=Nl, TransferMat=mat,k=k,l=l, Opt=opt)'
        pk=-1.
        goto, DONE
        end


if not keyword_set(Nl)  then Nl = cl*0
if not keyword_set(mat) then mrs_transfer_matrix,mat,k,l
if not keyword_set(Opt) then Opt = ' '

; Outputing files 
writefits,'xx_cl.fits', Cl
writefits,'xx_transmat.fits',transpose(mat*2.4e-9) ; Just a normalisation of the matrix to avoid numrical problems
writefits,'xx_ks.fits',k
writefits,'xx_ls.fits',l
writefits,'xx_noise.fits',Nl

; Launching prism
com = 'prism ' + Opt + ' xx_transmat.fits xx_ls.fits xx_ks.fits xx_noise.fits xx_cl.fits xx_pk.fits'

spawn, com

pk = mrdfits('xx_pk.fits')*2.4e-9

delete, 'xx_cl.fits'
delete, 'xx_transmat.fits'
delete, 'xx_ks.fits'
delete, 'xx_ls.fits'
delete, 'xx_noise.fits'
delete, 'xx_pk.fits'

DONE:

return, pk

end
