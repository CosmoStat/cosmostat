;+
; NAME:
;        mrs_wiener.
;
; PURPOSE:
;    Perform wiener filtering in the spherical harmonic space.
;       and Alm_output = Alm_input * WienerFilter
;          where  WienerFilter =  P^* S / [ | P P^* |^2 S + N ]
;                 P = instrumental beam (i.e. PSF, by default P = 1)
;                 S = A priori Signal Power Spectrum (default, power spectrum of the data)
;                 N = Noise power spectrum
;   If the keyword ALM is set (it should be a structure (see mrs_almtrans) containg the ALM coeff of the data, 
;   then the first parameter is not used and the alm are not calculated in this routine 
;   The Wiener is optimal in the sens of the Least mean square error of the reconstructed map. 
;   the power spectrum of the reconstructed map is biased (i.e. PowSpectrum(WienerMap) = PowSpectrum(RealMap) / (1 + N / (P^2S) ))
;     
;   If the Cole keyword is set, then the Wiener filtering is replaced by the Cole filtering and the power spectrum of 
;   the reconstructed map is unbiased, but the estimation is not optimal anymore for the east mean square error criterion.
;   The Cole filter is:
;             ColeFilter =  [ S / [ | P P^* |^2 S + N ] ]^0.5
;
;   If the keyword DataPrior is set, then the vector given by SignalPrior corresponds to an a priori on the power spectrum of the signal multiplied by P^2
;
; CALLING:
;     mrs_wiener, Imag, NoiseSpectrum, Recons, SignalPrior=SignalPrior, alm=alm, Spec1D=Spec1D, WienerFilter=WienerFilter, lmax=lmax, Psf=Psf, Cole=Cole, bias=bias
;
; INPUTS:
;     Imag -- IDL array of healpix map or GLESP structure: Input image to be denoised 
;     NoiseSpectrum, -- float or IDL 1D array: variance or power spectrum of the noise
;
; OUTPUTS:
;     Recons -- IDL array of healpix map or GLESP structure: Output denoised  image 
;
; KEYWORDS:
;     SignalPrior -- input IDL 1D array: 1D power spectrum of the expected power spectrum
;                                  By default, it is estimated from the data 
;     DataPrior -- -- scalar: if set, then the vector given by SignalPrior corresponds to an a priori on the power spectrum of the signal multiplied by Psf^2
;     Psf  -- input IDL 1D array: Instrumental beam (i.e. PSF[l] = Spherical Harmonics a_{l,0} = ... = a_{l,l+1} of the 
;                                 intrumental beam (i.e. Point Spread Funciton).
;  
;     alm -- input/output ALM structure: IF this keyword is set, the input image is not used, and the ALM given by this
;                                 keyword are used. The denoised ALM are stored in the structure
;     Spec1D -- output IDL 1D array: Estimated power spectrum from the denoised ALM
;     WienerFilter -- output IDL 1D array: Wiener or Cole filter window
;     lmax -- Maximum l used in the calculation of the ALM. This kewyord is not used if the keyword ALM is set.
;     Cole -- scalar: if set, then the window filter is not the Wiener filter, but the Cole filter which
;     bias -- output scalar or IDL 1D array: Estimated power spectrum bias
;
; EXTERNAL CALLS:
;       mrs_alm2spec, mrs_almrec,mrs_almtrans
;
; EXAMPLE:
;       
;         
; HISTORY:
;	Written:  Jean-Luc Starck, 2011
;	December, 2011 File creation
;--------------------------------------------------------------------------------------------------------

;=====================================================================

pro mrs_wiener, Imag, NoiseSpectrum, Recons, SignalPrior=SignalPrior, DataPrior=DataPrior, alm=alm, Spec1D=Spec1D, WienerFilter=WienerFilter,lmax=lmax, Psf=Psf, Cole=Cole, bias=bias, NoiseRea=NoiseRea,  WienerNoiseRea=WienerNoiseRea

   if not keyword_set(Psf) then Psf=1.
   if not keyword_set(alm) then mrs_almtrans, imag, alm, /tab, /complex, /norm,lmax=lmax $
   else begin
     if alm.COMPLEX_ALM NE 0 or alm.norm NE 1 or alm.tab NE 1 then begin
       print, 'Error: the input ALM must be computed with options "/tab,/norm" and not "/psp" or "/complex"'
       goto, DONE
    end
   end
   
   PSNoiseRea = 0
   if keyword_set(NoiseRea)  and not keyword_set(NoiseSpectrum)  then begin
      PSNoiseRea = 1
      mrs_almtrans, NoiseRea, AlmNoise, /tab, /complex, /norm,lmax=lmax 
      NoiseSpectrum = mrs_alm2spec(AlmNoise)
   end
   
   if not keyword_set(SignalPrior) then begin
      PSig = mrs_alm2spec(ALM)
      SignalPrior = PSig - NoiseSpectrum
      ind = where(SignalPrior lt 0, c)
      if c GT 0 then SignalPrior[ind] = 0
      DataPrior=1
   end
   H2 = Psf^2.
   if not keyword_set(Cole) then begin
     if keyword_set(DataPrior) then WienerFilter =  Psf*SignalPrior / (H2*(SignalPrior  +  NoiseSpectrum)) $
     else WienerFilter =  Psf*SignalPrior / (H2*SignalPrior  +  NoiseSpectrum) 
   end else begin
     if keyword_set(DataPrior) then WienerFilter = sqrt( SignalPrior / (H2*(SignalPrior  +  NoiseSpectrum) )) $
     else WienerFilter = sqrt( SignalPrior / (H2*SignalPrior  +  NoiseSpectrum) )
   end
   WienerFilter[0] = 1
   for i=0,alm.lmax do  alm.alm[i,*,*] = alm.alm[i,*,*] * WienerFilter[i]

   Spec1D = mrs_alm2spec(ALM)
   mrs_almrec,alm, Recons
   
   if keyword_set(NoiseRea) then begin
      if PSNoiseRea  EQ 0 then mrs_almtrans, NoiseRea, AlmNoise, /tab, /complex, /norm,lmax=lmax 
      for i=0,alm.lmax do  AlmNoise.alm[i,*,*] = AlmNoise.alm[i,*,*] * WienerFilter[i]
      mrs_almrec,  AlmNoise, WienerNoiseRea
   end
   
   bias = SignalPrior
   bias[*]=1.
   if not keyword_set(Cole) then begin
     if keyword_set(DataPrior) then bias = SignalPrior / (SignalPrior + NoiseSpectrum) $
     else bias = H2*SignalPrior / (H2*SignalPrior + NoiseSpectrum)
   end
DONE:

END

;=====================================================================
