;+
; NAME:
;        mrs_master_powspec
;
; PURPOSE:
;   Computes the power spectrum of a map with missing data, 
;   using the HEALPix representation (nested data
;   representation by default). The MASTER method is used.
;  
;
; CALLING:
;     P = mrs_master_powspec( Imag, Mask,  PSCMB_Mask=PSCMB_Mask)
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image to be transformed 
;     Mask -- IDL array of healpix map: Input Mask of missing data (Mask[k] =0, if pixel k is missing) 
;    
; OUTPUTS:
;     P -- 1D IDL fltarr: Power Spectrum. P[k] = Mean(  POWSPECTRUM[*,l]  )
;
; OPTIONAL INPUT KEYWORDS:
;     PSCMB_Mask -- IDL array:   power  spectrum of the mask
;
;
; EXTERNAL CALLS:
;       mrs_almtrans.pro
;
; EXAMPLE:
;       Compute the power spectrum of an image. 
;               P = mrs_master_powspec(Imag, Mask) 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	December, 2005 File creation
;--------------------------------------------------------------------------------------------------------


;====================================================================

function mrs_master_powspec, Imag, Mask,  PSCMB_Mask=PSCMB_Mask, BIN=BIN, lmax=lmax

COMMON C_PLANCK

if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: p = mrs_master_powspec( Imag, Mask,  PSCMB_Mask=PSCMB_Mask)'
	Ret=-1
        goto, DONE
        end

  MASTER_BIN=0
  if keyword_set(BIN) then MASTER_BIN = 1

  if not keyword_set(PSCMB_Mask) then PSCMB_M = mrs_powspec(Imag*mask) $
  else  PSCMB_M=PSCMB_Mask

  ; MASTER CALL TO COMPUTE THE POWER SPECTRUM
   PS_Master=0
   NewMasterPSR = 0
   mll = 0
   invMll = 0
    ; mll = output  coupling matrice
    ; mbb = output  Binned  coupling matrice
    ; matp = output matrice to transform the Cl to the binned Cl
    ; matq = output matrice to the binned Cl to the CL (transpose of matp)
    ; lmax = input lmax
    if not keyword_set(lmax) then lmax = P_LMAX  
    ; nbins = number bins we want
    nbins = lmax
    ;vs = size(NewPS)
    ;nbins = vs[1]
    deltl = lmax / nbins
    m1024 = mask   ; mrs_resize(mask, nside=1024)
    ; MASTER_BIN = 0  ; if MASTER_BIN == 1 then applied a binning
    if MASTER_BIN EQ 0 then begin
      compute_mllp, mll, mbb, matp, matq, lmax, nbins, mask=m1024, filebinmax=filebinmax, filebinmin=filebinmin, deltl=deltl
      invMll = invert(mll, /double)
      PS_Master = invMll[0: lmax,0: lmax] # PSCMB_M
    end else begin
     DIR = getenv('ISAP')
     filebinmin = DIR+'/data/fmin'
     filebinmax = DIR+'/data/fmax'
     nbin=140
     compute_mllp, mll, mbb, matp, matq, lmax, nbins, mask=m1024, filebinmax=filebinmax, filebinmin=filebinmin
     invMll = invert(mll, /double)
     PS_Master = invMll[0: lmax,0: lmax] # PSCMB_M    ; Powspectrum non binne
     PS_MasterB = mbb[*,0: lmax] # PSCMB_M             ; Binned Powspectrum 
     PS_MasterRebin =  matq # mbb[*,0: lmax] # PSCMB_M  ; Rebin to 3001
     end
 
DONE:
 
  RETURN, PS_Master
end


