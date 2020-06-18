FUNCTION CROSS_SPECTRUM, x, y, DX=dx, DY=dy, DT=dt, N_SEG=n_seg, $ 
                         POIS=pois, COH=coh, CO_ERR=co_err, TLAG=tlag, LAG_ERR=lag_err, $
                         PHASE=phase, PH_ERR=ph_err, F_C=f_c, DF=df, BINF=binf, BIN_N=bin_n, $
                         P_1=binP_1, P_2=binP_2, CTR=ctr, MINBIN=minbin,  $
                         F_L = f_l, F_U=f_u, S_1=binS_1, S_2=binS_2, RMS=rms, $
                         P1ERR=P_1err, P2ERR=P_2err, COVAR=covar, COV_ERR=cov_err, $
                         ENDMATCH=endmatch, CHATTER=chatter, MEAN_1=mean_1, MEAN_2=mean_2, $
                         NN=nn, ALL_COV=covar_all, ALL_ERR_COV=cov_all_err

; ----------------------------------------------------------
;+
; NAME:
;       CROSS_SPECTRUM
;
; PURPOSE:
;       Calculate cross-spectral products of two time series
;
; AUTHOR:
;       Simon Vaughan
;
; CALLING SEQUENCE:
;       C = CROSS_SPECTRUM(x,y)
;
; INPUTS:
;       x       - (float array) first time series, evenly sampled
;       y       - (float array) second time series, same sampling
;
; OPTIONAL INPUTS:
;       dx      - (float array) 'errors' on series 1 data
;       dy      - (float array) 'errors' on series 2 data
;       dt      - (float) sampling period (same for x and y)
;       n_seg   - (integer) number of points in each segment 
;       pois    - (logical) if set then data are Poisson
;       binf    - (float) frequency binning factor 
;       ctr     - (logical) indicates data are in ct/bin (cf. ct/s)
;       minbin  - (logical) minimum points per bin (for log-freq rebinning)
;       rms     - (logical) whether to use 1/mean
;                           normalisation on FTs
;       endmatch - (logical) should we end-match each segment? (default=FALSE)
;       CHATTER  - (integer) more or less feedback to screen?
;
; OUTPUTS:
;       C       - (complex array) Complex cross-spectrum of two series
;       
; OPTIONAL OUTPUTS:
;       coh     - (float array) coherence of two series
;       co_err  - (float array) error on coherence
;       phase   - (float array) phase lag (-pi, +pi)
;       ph_err  - (float array) error on phase lag
;       tlag    - (float array) time lag
;       lag_err - (float array) error on time lag
;       covar   - (float array) covariance spectrum 
;       cov_err - (float array) error on covariance 
;       cov_all - (float) covariance computed over all frequencies
;       f_c     - (float array) frequency (centre of bin)
;       f_l     - (float array) lower bound of frequency bin
;       f_u     - (float array) upper bound of frequency bin
;       df      - (float) resolution of periodograms
;       P_1     - (float array) auto-periodogram of series 1 (x)
;       P_2     - (float array) auto-periodogram of series 2 (y)
;       bin_n   - (int array) number of periodogram points in each bin
;       nn      - (float array) bias correction for Poisson noise
;
; DETAILS:
;       Computes the cross spectral properties of two
;       evenly sampled time series x and y. The cross-spectrum
;       of two time series is the complex-valued product of the
;       Fourier transforms of the two series. If X = FT(x) and
;       Y = FT(y) then the cross-spectrum of x and y is:
;
;         C_xy(f) = X'(f) * Y(f)
;
;       where the ' represents the complex conjugate.
;       The cross-spectrum of a series with itself is the
;       auto-spectrum (usually called the power spectrum). 
;       In reality we only have discrete, finite
;       time series x_i and y_i, with DFTs X_j and Y_j, where j
;       represents the Fourier frequencies. In order
;       to estimate the cross-spectrum we calculate the
;       cross-periodogram from the DFTs of the two series
;
;       C_xy(j) = X(j)' * Y(j).
;
;       The cross-spectrum (or cross-periodogram) are
;       complex valued. The magnitude of this is used to
;       calculate the 'coherence spectrum' and the phase is used
;       to calculate the 'phase spectrum'.
;       The coherence estimated from the the normalised square
;       of the amplitude of the complex cross-periodogram:
;
;         coh(j) = |<C(j)>|^2 / <|X(j)|^2> <|Y(j)|^2>
;
;       where <.> denotes the ensemble average (or expectation).
;       The denominator is the product of the two (auto-) periodgrams
;       for x and y.
;       The coherence is a measure of the linear correlation
;       between the two series and ranges from 0 to 1, where
;       0 means no correlation and 1 means perfect correlation.
;
;       The phase lag is the phase of the complex cross spectrum.
;
;           phi(j) = arg( <C(j)> )
;
;       and is usually defined on the interval (-pi,pi).
;       The phase lag can be simply converted into a time lag
;
;           tlag(j) = phi(j) / 2*pi*f(j)
;
;       and records the average delay betwen the two series
;       as a function of Fourier frequency. This can be positive
;       or negative. A positive lag indicates x leading y.
;       Note that the phase is limited to (-pi,pi) which can
;       cause a wrap-around effect when the true lag is greater
;       than pi in magnitude, ie. |time lag| > 1/2f. Therefore
;       a lead of 3/4f (3pi/2) could be measured as a lag of
;       1/4f (pi/2).
;
;       The cross-perdiogram is relatively straightforward to
;       calculate. The phase and coherence require
;       averaging over realisations in order to approach the
;       expectation <.>. This routine handles this in two
;       ways.
;
;       The time series of length N is broken into equal length
;       segments of length N_SEG, where N_SEG is specified on input.
;       Each of these segments is considered an independent
;       realisation of the underlying process. The cross-periodgram
;       for each segement are averaged together to provide phase and
;       coherence estimates at each Fourier frequency. The Fourier
;       frequencies for a time series of length N_SEG with sampling
;       interval dT are therefore f(j) = j/(N*dT) with
;       j=1,2,...,N_SEG/2. 
;
;       In addition to this the cross-spectrum may
;       be averaged over consecutive Fourier frequencies. This
;       obviously lowers the resolution of the resultant spectrum but
;       improves the quality of the estimate (decreases sampling
;       fluctuations). The degree of frequency averaging is specified
;       with the BINF input keyword. Using BINF > 1 the routine
;       performs even frequency binning, averging together BINF
;       Fourier frequencies. With BINF < -1 the routine uses
;       logarithmic frequency binning, with each bin spanning a factor
;       of |BINF|. E.g. with BINF = -1.2 each bin spans f_j -> 1.2f_j
;
;       If the POIS keyword is set the data are assumed
;       to be Poisson data in counts/sec units with dT in
;       units of seconds. The coherence and its error are
;       then calculated using equation 8 of Vaughan & Nowak
;       which applies a Poisson noise correction.
;       
;       If the 'errors' (dx, dy) are supplied for the time series 
;       then the noise terms in the coherence correction are derived
;       from these instead of the usual sqrt(counts) rule.
;
;       By default:
;       FTs are not normalised before averaging cross spectra
;       over segments - if there are large variances in rms and mean 
;       level it may be more approprate to renormalise according to 
;       e.g. mean level before averaging over segments.
;       If RMS keyword is TRUE then:
;       Each FFT is normalised by the mean rate of the data segment. 
;       This will effectively divide out any rms-flux relation. The FFTs
;       are normalised by SQRT(2.0*dt/N_seg/<x>^2) so that when squared
;       they give a periodogram with [rms/mean]^2 Hz^-1 fractional units. 
;
;       [18/07/2013] An additional output is the "co-rms". This was 
;       introduced by Wilkinson & Uttley (2009). The idea is that if
;       we have two simultaneous time series, x(t) a "subject" and 
;       y(t) a "reference", we can extract the power density of x that
;       is covariant with y. This can be seen in the Fourier domain
;       (Uttley et al. 2011; Cassatella et al. 2012). It is essentially
;       
;          covar = coh_xy * <|X|^2 - |N_x|^2> * (f_max - f_min)
;          
;       i.e. the PSD of x, integrated over some frequency band to give the
;       total power, multiplied by its coherence with y. If the two series 
;       are perfectly coherent the "covar" is just the PSD of x. If x is
;       partially coherent with y it reveals only the covariance power.
;       We use the noise-corrected coherence and noise-subjected PSDs.
;       (See references for more details.)
;       
;       We then take the sqrt of this to give the "coherent rms" or co-rms.
;       This is given for each frequency as COVAR (vector) with error COV_ERR. 
;       We also return the value estimated by averaging over all available
;       frequencies in the COV_ALL (and COV_ERR_ALL) outputs.
;       By computing COVAR for different "subject" time series x_i using the
;       same "reference" time series y one can estimate, at a given frequency,
;       the co-rms of each x_i. If the "reference" time series has high S/N
;       this will give smaller uncertainties than the rms, to which it is
;       closely related.
;
;       For more details see:
;         Bendat J. S., Piersol A. G., 1986, 
;            Random Data: Analysis & Measurement Procedures
;         Vaughan B., Nowak M., 1997; ApJ, 474, L43   
;         Wilkinson & Uttley 2009, MNRAS, 397, 666
;         Uttley et al. 2011, MNRAS, 414, L60 
;         Cassatella et al. 2012, MNRAS, 427, 2985 
;
;       Example usage
;       IDL> coh = CROSS_SPECTRUM(x, y, DT=2.0, F_C=f, /POIS)
;
; HISTORY:
;       01/02/07  - v1.0 - first working version
;                           major re-write and bug fix
;       30/04/09  - v2.1 - added COUNT keyword and functionality
;       25/02/10  - v2.2 - several minor bugs fixed. Added MINBIN
;                           input
;       21/03/11  - v3.0 - rewrote binning routine to use REGROUP.
;       19/12/11  - v3.1 - added S_1, S_2 keywords as output;
;                          added errors (dx, dy) as input for noise 
;                          correction terms
;       21/02/12  - v3.2 - Change COUNT keyword to CTR to avoid confusion 
;                          with 'count' as used in WHERE().
;                          Output PHASE and PH_ERR for phase data, as well
;                          as TLAG and LAG_ERR for time lag data.
;                          Output COVAR (and COV_ERR) for 'covariance'
;                          spectrum as defined in e.g. Uttley et al. (2011)
;       29/02/12  - v3.3 - replace TOTAL(...,1)/n statements with the more 
;                          transparent MEAN(...,DIMENSION=1). 
;                          Added functionality for the RMS keyword.
;       26/11/12  - v3.4 - bug fix. [Thanks to DB.]
;       15/01/13  - v3.5 - bug fix. (Skip ERR_1 arrays when POIS not set.)
;       25/04/13  - v3.6 - added ENDMATCH input
;       08/07/13  - v3.7 - added CHATTER keyword and MEAN_1/MEAN_2 output
;                           added NN output
;       19/07/13  - v3.8 - added COV_ALL output. Made correction to COV_ERR
;                           output
;
; USES:
;       REGROUP
;
; NOTES: 
;       No correction for deadtime in Fourier transforms.
;;
;-
; ----------------------------------------------------------

; options for compilation (recommended by RSI)

  COMPILE_OPT idl2

; watch out for errors

  ON_ERROR, 3

; ----------------------------------------------------------
; Check the arguments 

  n = N_ELEMENTS(x)

; is the sampling period dt supplied? (default = 1.0)

  if (N_ELEMENTS(dt) eq 0) then dt=1.0

; is the number of segments supplied? (default = 20)

  if (N_ELEMENTS(n_seg) eq 0) then n_seg = n / 20

; check that N_SEG is integer type

  type = SIZE(n_seg, /TYPE)
  if (type ne 2 and type ne 3) then begin
      PRINT,'** N_SEG not an integer in CROSS_SPECTRUM'
      RETURN, !NULL
  endif

; is it large enough?

  if (n_seg lt 16) then begin
      PRINT,'** N_SEG too small in CROSS_SPECTRUM.'
      RETURN, !NULL
  endif

; set BINF to default (=1) if not set

  if N_ELEMENTS(binf) eq 0 then binf=1

; is MINBIN set?

  if (not KEYWORD_SET(minbin)) then minbin = 1

; define the level out output to screen

  IF NOT KEYWORD_SET(chatter) THEN chatter = 0

  IF (chatter ge 0) THEN BEGIN
    PRINT,"-- CROSS_SPECTRUM v3.5"
  ENDIF
  
; ----------------------------------------------------------
; Check the input data

; is time series x data array well-defined?

  if (n lt n_seg) then begin
      PRINT,'** Not enough data in CROSS_SPECTRUM.'
      RETURN, !NULL
  endif
  
  if (TOTAL(FINITE(x)) lt n) THEN BEGIN
     PRINT, '** Non-finite values in X'
     RETURN, !NULL
  ENDIF

; check y has same number of elements

  if (N_ELEMENTS(y) ne n) then begin
      PRINT,'** X and Y of different lengths in CROSS_SPECTRUM.'
      RETURN, !NULL
  endif

  if (TOTAL(FINITE(y)) lt n) THEN BEGIN
     PRINT, '** Non-finite values in Y'
     RETURN, !NULL
  ENDIF

; check the size of the error arrays

  if KEYWORD_SET(pois) then begin
      n_err = N_ELEMENTS(dx)
      if (n_err eq 0) then err_1 = SQRT(ABS(x/dt))
      if (n_err eq 1) then err_1 = MAKE_ARRAY(n, VALUE=dx)
      if (n_err eq n) then err_1 = dx
      if (n_err gt 1 and n_err lt n) then begin
          PRINT, '** DX not 0, 1 or N in length.'
          RETURN, !NULL
      endif 

      n_err = N_ELEMENTS(dy)
      if (n_err eq 0) then err_2 = SQRT(ABS(y/dt))
      if (n_err eq 1) then err_2 = MAKE_ARRAY(n, VALUE=dy)
      if (n_err eq n) then err_2 = dy
      if (n_err gt 1 and n_err lt n) then begin
          PRINT, '** DY not 0, 1 or N in length.'
          RETURN, !NULL
      endif 
  endif

  if KEYWORD_SET(ctr) then begin
    err_1 = SQRT(ABS(X))
    err_2 = SQRT(ABS(Y))
  endif 

; ----------------------------------------------------------
; Break the time series into M segments of length N_seg each

  m = n/n_seg
  n_cut = m*n_seg
  x_1 = REFORM(x[0:n_cut-1], n_seg, m)
  x_2 = REFORM(y[0:n_cut-1], n_seg, m)
  if KEYWORD_SET(pois) then begin
    err_1 = REFORM(err_1[0:n_cut-1], n_seg, m)
    err_2 = REFORM(err_2[0:n_cut-1], n_seg, m)
  ENDIF

  IF (chatter ge 1) THEN BEGIN
    PRINT, '-- Segment length:', n_seg
    PRINT, '-- No. segments: ', m
    PRINT, '-- Raw data points ', n, '; Used data points: ', n_cut
  ENDIF

; no. +ve frequencies

  nf = n_seg/2                      

; frequency resolution

  df = 1.0/(dt*n_seg)

; frequency array: f_j = j / (N*dt) where j = 1, 2, ..., N_seg/2 

  f = (FINDGEN(nf)+1) * df 
  IF (chatter ge 1) THEN BEGIN
    PRINT,'-- Number of frequencies:', nf
    PRINT,'-- Resolution ', df, ' and range: ', MIN(f), MAX(f)
  ENDIF

; ----------------------------------------------------------
; Calculate the auto- and cross-periodograms

; calculate mean (DC component) for each segment

  mean_1 = MEAN(x_1, DIMENSION=1, /DOUBLE) 
  mean_2 = MEAN(x_2, DIMENSION=1, /DOUBLE)

; subtract mean from each segment

  x_1 = TEMPORARY(x_1) - TRANSPOSE(REBIN(mean_1, m, n_seg, /SAMPLE))
  x_2 = TEMPORARY(x_2) - TRANSPOSE(REBIN(mean_2, m, n_seg, /SAMPLE))

; end-matching ?

  IF (KEYWORD_SET(endmatch)) THEN begin
    IF (endmatch ne 0) THEN BEGIN
      IF (chatter ge 1) THEN PRINT,'-- Applying End-Matching.'
      FOR i = 0, m-1 DO BEGIN
        x_1[*, i] = ENDMATCH(REFORM(x_1[*, i]))
        x_2[*, i] = ENDMATCH(REFORM(x_2[*, i]))
      ENDFOR
    ENDIF
  ENDIF

; Calculate the Fourier transform of mean-subtracted data
; (NB: 1/N factor when for forward DFT)

  dft_1 = FFT(x_1, 1, DIMENSION=1, /DOUBLE)
  dft_2 = FFT(x_2, 1, DIMENSION=1, /DOUBLE)

; extract positive (non-zero) frequencies only

  dft_1 = TEMPORARY(dft_1[1:nf, *])
  dft_2 = TEMPORARY(dft_2[1:nf, *])
  
; normalise the DFTs

  norm = SQRT(2.0D * dt / FLOAT(n_seg))
  norm1 = norm
  norm2 = norm

; Fractional rms normalisation, when /RMS is used

  if KEYWORD_SET(rms) then begin

    IF (MIN(mean_1) eq 0) OR (MIN(mean_2) eq 0) THEN BEGIN
      PRINT, '** Warning: mean[i] = 0 in CROSS_SPECTRUM'
    ENDIF
    means_1 = TRANSPOSE(REBIN(mean_1, m, nf, /sample))
    means_2 = TRANSPOSE(REBIN(mean_2, m, nf, /sample))
    norm1 = norm1 / means_1
    norm2 = norm2 / means_2
  endif
  
  
  dft_1 = TEMPORARY(dft_1) * norm1 
  dft_2 = TEMPORARY(dft_2) * norm2
  
; Quality check of the FFT data
  
  if (TOTAL(FINITE(dft_1)) lt N_ELEMENTS(dft_1)) THEN BEGIN
     PRINT, '** Non-finite values in DFT_1'
   ENDIF
  if (TOTAL(FINITE(dft_2)) lt N_ELEMENTS(dft_2)) THEN BEGIN
     PRINT, '** Non-finite values in DFT_2'
  ENDIF
  
; square the DFT: power = |DFT|^2

  P_1 = ABS(dft_1)^2
  P_2 = ABS(dft_2)^2

; Calculate the complex cross periodogram C = X' * Y

  C = CONJ(dft_1) * dft_2

; ----------------------------------------------------------
; Average over the ensamble of time series segments and
; adjacent frequencies

; average the second-order quantities: C, P_1, P_2
; over the ensemble of M segments 

  if (m gt 1) then begin
      binC = MEAN(C, DIMENSION=2, /DOUBLE) 
      binP_1 = MEAN(P_1, DIMENSION=2, /DOUBLE)
      binP_2 = MEAN(P_2, DIMENSION=2, /DOUBLE)
  endif else begin
      binC = C
      binP_1 = P_1
      binP_2 = P_2
  endelse

; bin over BINF adjacent frequency bins 
; remember that C is complex, so bin real and imaginary parts

  if (binf gt 1) then begin
      binwidth = binf * df
  endif else begin
      binwidth = binf
  endelse

  if (binf gt 1 or binf lt 0) then begin
      Re = REAL_PART(binC)
      Im = IMAGINARY(binC)
      bin_Re = REGROUP(Re, f, DX=df, BINWIDTH=binwidth, MINBIN=minbin, BIN_X=bin_f, $
                       BIN_N=bin_n, BIN_L=f_l, BIN_U=f_u)
      bin_Im = REGROUP(Im, f, DX=df, BINWIDTH=binwidth, MINBIN=minbin)
      binC = COMPLEX(bin_Re, bin_Im)
      binP_1 = REGROUP(binP_1, f, DX=df, BINWIDTH=binwidth, MINBIN=minbin)
      binP_2 = REGROUP(binP_2, f, DX=df, BINWIDTH=binwidth, MINBIN=minbin)
  endif else begin
      bin_n = MAKE_ARRAY(N_ELEMENTS(binC), /INTEGER, VALUE=1)
      bin_f = f
      f_u = f + df*0.5
      f_l = f - df*0.5
  endelse

; Number of periodogram points per binned data point
; from M segments and BIN_N(i) adjacent frequencies

  bin_n = bin_n * m

; calculate coherence |<C>|^2 / <|X|^2><|Y|^2>

  coh = ABS(binC)^2 / (binP_1 * binP_2)

; calculate the phase lag (phase of complex cross spectrum)

  phase = ATAN(binC, /PHASE)

; ----------------------------------------------------------
; Calculate errors from standard sampling formulae
; See Vaughan & Nowak (1997) and table 9.6 of Bendat & Pierson (1986)

; calculate error on coherence

  co_err = SQRT(2.0 / bin_n) * coh * (1.0-coh) / SQRT(coh)  

; calculate the phase lag error (see Bendat & Piersol 1986)

  ph_err = SQRT(1.0 / bin_n) * SQRT( ABS((1.0-coh)/(2.0*coh)) )

; convert from phase to time lag

  tlag = phase / (2.0*!pi*bin_f)
  lag_err = ph_err / (2.0*!pi*bin_f)

; if no POIS correction then set binS_1 = binP_1 etc. [26/11/2012 - SV]

  binS_1 = binP_1
  binS_2 = binP_2

; ----------------------------------------------------------
; if POIS keyword set then apply the Poisson filtering to coherence
; see Vaughan & Nowak (1997)

  N_1 = 0.0
  N_2 = 0.0
  nn = 0.0
  if KEYWORD_SET(pois) then begin

; Calculate expected Poisson rates: N_1, N_2
; where N_1 = <x> * N_seg/dT and so on.
;
; NB: This is the case because we did not normalise the periodograms.
; If we used the "rms" norm, 2*dT/(<x>^2*N_seg), where x is in ct/s,
; the Poisson noise level would be N = 2*dT*<err^2>/<x>^2, or in the 
; Poisson case where x is in units of ct/s with no background, N = 2/<x>.
; If we do not re-normalise (i.e. norm of 1) the
; Poisson noise level is N = <err^2> * N_seg.

      meanerrsq_1 = TOTAL(err_1^2, 1, /DOUBLE)
      meanerrsq_2 = TOTAL(err_2^2, 1, /DOUBLE)
      IF (SIZE(norm1, /N_DIMENSION) gt 0) THEN BEGIN
        norm1 = MEAN(norm1, DIMENSION=1, /DOUBLE)
        norm2 = MEAN(norm2, DIMENSION=1, /DOUBLE)
      ENDIF
      N_1 = meanerrsq_1 * norm1^2
      N_2 = meanerrsq_2 * norm2^2
      N_1 = MEAN(N_1, /DOUBLE)
      N_2 = MEAN(N_2, /DOUBLE) 

; Calculate ensemble averaged noise-subtracted spectra: <S_1(j)> 
; S_1(j) = |X(j)|^2 - |N(j)|^2 where N is DFT(noise)
; so <S_1(j)>  = <|X(j)|^2> - <|N(j)|^2>

      if (m gt 1) then begin
          binS_1 = MEAN(P_1, DIMENSION=2, /DOUBLE) - N_1
          binS_2 = MEAN(P_2, DIMENSION=2, /DOUBLE) - N_2
      endif else begin
          binS_1 = P_1 - N_1
          binS_2 = P_2 - N_2
      endelse

; bin over BINF adjact frequency bins 

      if (binf ne 1) then begin
          binS_1 = REGROUP(binS_1, f, DX=df, BINWIDTH=binwidth, MINBIN=minbin)
          binS_2 = REGROUP(binS_2, f, DX=df, BINWIDTH=binwidth, MINBIN=minbin)
      endif

; apply Poisson filtering to coherence

      nn = (binS_1*N_2 + binS_2*N_1 + N_1*N_2) / bin_n
      coh = ( ABS(binC)^2 - nn ) / (binS_1 * binS_2)

; and calculate error on Poisson-corrected coherence
; (see equation 8 of Vaughan & Nowak 1997)

      err = 2.0 * nn^2 * bin_n / (ABS(binC)^2 - nn)^2 + $
            (N_1/binS_1)^2 + (N_2/binS_2)^2 + bin_n*co_err^2/coh^2
      co_err = ABS(coh)/SQRT(bin_n) * SQRT(err)

  endif

; ----------------------------------------------------------
; compute the 'covariance' spectrum
; Assume X (and therefore S_1) is a 'comparison' series and
; Y (and S_2) is a 'reference' band. We calculate the covariance
; of X with Y. This resembles the part of the power spectrum of 
; the comparison band that is coherent with the reference band.
 
  covar2 = coh * binS_1 * (f_u - f_l)
  covar = SQRT(covar2)
  
;  cov_err = (binS_1 * N_2 * (f_u - f_l)) + (binS_2 * N_1 * (f_u - f_l)) + (N_1 * N_2 * (f_u - f_l))
;  cov_err = SQRT(cov_err / binS_2 / (2*bin_n))

; Error formula (Phil Uttley - priv. comm. July 2013)

  cov_err = SQRT(nn * (f_u - f_l) / binS_2 / 2.0)

; ----------------------------------------------------------
; compute the covariance over all frequencies (a useful check)

  C_all = MEAN(C)
  P_1_all = MEAN(P_1)
  P_2_all = MEAN(P_2)
  N_1_all = MEAN(N_1)
  N_2_all = MEAN(N_2)
  nn_all = 0.0
  S_1_all = P_1_all 
  S_2_all = P_2_all 
  if KEYWORD_SET(pois) THEN BEGIN
    S_1_all = P_1_all - N_1_all
    S_2_all = P_2_all - N_2_all
    nn_all = (S_1_all*N_2_all + S_2_all*N_1_all + N_1_all*N_2_all) / (nf * m)
  ENDIF

  delta_f = (MAX(f_u) - MIN(f_l))
  covar_all = ( ABS(C_all)^2 - nn_all ) * delta_f / (S_2_all)
  covar_all = SQRT(covar_all)
  
  cov_all_err = (S_1_all * N_2_all * delta_f) + $
                (S_2_all * N_1_all * delta_f) + $
                (N_1_all * N_2_all * delta_f)
  cov_all_err = SQRT(cov_all_err / S_2_all / (2 * nf * m))

  cov_all_err = SQRT(nn_all * delta_f / S_2_all / 2.0)


; ----------------------------------------------------------
; calculate errors on output PSDs

  npb = SQRT(bin_n > 1)
  P_1err = binP_1 / npb
  P_2err = binP_2 / npb
  
; ----------------------------------------------------------
; return finished coherence and freqiencies

  f_c = bin_f
  RETURN, binC

END


