; NAME:
;        mrsp_threshold
;
; PURPOSE:
;		Threshold the decomposition coefficents of a polarized healpix map. Work on sevral decompositions see mrsp_trans.pro
;
; CALLING:
;		mrsp_threshold, Trans, NSigma=Nsigma, KillLastScale=KillLastScale
;
; INPUT / OUTPUT:
;
;     Trans -- IDL structures with the following fields (see mrsp_trans.pro):  
;			NBRSCALE  -- LONG: Number of scales of the wavelet transform
;			Nside -- Nside value
;			Npix -- Number of pixels 
;			DEC1 -- IDL structure: First component transformation (depends on the chosen transform)
;			DEC2 -- IDL structure: Second component transformation (depends on the chosen transform)
;			DEC3 -- IDL structure: Third component transformation (depends on the chosen transform)
;			ebdec -- int: set to 1 if an EB decomposiiton has been applied
;			lmax -- int: maximum used spherical harmonic (for isoptropic wavelet transform only)
;			TransChoice -- int: Transform number, 1 to 6.
;			TransTypeName -- strarray: array of transform names	
;
; INPUT KEYWORDS:
;		Nsigma: float = Level of thresholding (default is 3)
;		KillLastScale: if set, the last scale is set to zero
;
; EXTERNAL CALLS:
;         mrsp_wtget
;         mrsp_wtput
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of a vector field I with five scales
;       The result is stored in WT, then wavelet coefficients are threshold at 2 sigma;
;               mrsp_trans, Imag, WT, NbrScale=5, /UWT
;				mrsp_threshold, WT, NSigma=2
;         
; HISTORY:
;	Written:  Jean-Luc Starck, May 2008
;---------------------------------------------------------------------------------------------------------------------

pro mrsp_threshold, Trans, NSigma=Nsigma, Mad=Mad, KillLastScale=KillLastScale

if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: mrsp_threshold, Trans, Nsigma=Nsigma, KillLastScale=KillLastScale'
        goto, DONE
end
	 
if  not keyword_set(NSigma) then NSigma = 3

TabNbrBandPerScale = Trans.TabNbrBandPerScale
for c=0,2 do begin
for j1=0,Trans.NbrScale-2 do begin
for j2=0,TabNbrBandPerScale[j1]-1 do begin
  Band = mrsp_wtget(Trans, c, j1, BandNumber =j2)
  
  MadSigma = mad(Band) 
  ThresholdLevel = NSigma *  MadSigma
  ind = where( ABS(BAND) LT ThresholdLevel, N)
  if N LT 0 then N = 0
  P = float(N) / float(N_ELEMENTS(Band)) * 100.
  ; Band = mrs_absthreshold(Band, ThresholdLevel, soft=soft)
  if N GT 0 then Band[ind] = 0
  
  print, "BAND ", c+1, j1+1, j2+1, ": sigma = ", sigma(Band), ", MadSig = ", MadSigma, ",  P = ", P
  ; Band[*] = 0
  mrsp_wtput, Trans, Band, c, j1, BandNumber =j2
end
end
end


j = Trans.NbrScale-1
if keyword_set(KillLastScale) then begin
  for c=0,2 do begin
     Band = mrsp_wtget(Trans, c, j)
     Band[*] = 0
     mrsp_wtput, Trans, Band, c, j
  end
end

DONE:

END




