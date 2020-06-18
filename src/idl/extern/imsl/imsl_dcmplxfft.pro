; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_dcmplxfft.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
; NAME:
;       DCMPLXFFT
; PURPOSE:
;       compute a complex FFT using double precision.
;
; CATEGORY:
;       Signal Processing
; CALLING SEQUENCE:
;       DCMPLXFFT, r_in, i_in, r_out, i_out
; INPUTS:
;       r_in    Real part of the complex array to be transformed.
;       i_in    Imaginary part of the complex array to be transformed.
; KEYWORD PARAMETERS:
;       Backward If present and nonzero, then an inverse FFT is
;                computed.
; OUTPUTS:
;       r_out   Real part of the transformed array
;       i_out   Imaginary part of the transformed array
; SIDE EFFECTS:
;       None
; RESTRICTIONS:
;       None yet discovered
; PROCEDURE:
;       Use the real FFT provided by FFTCOMP() to compute a complex
;       FFT in double precision.
;
PRO imsl_Dcmplxfft, r_in, $
               i_in, $
               r_out, $
               i_out, $
               Backward=backward

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Convert input to double precision.
   r_in_cvt = double(r_in)
   i_in_cvt = double(i_in)
   n = IMSL_N_ELEMENTS(r_in_cvt)
   n2 = n/IMSL_2
   ;
   ;
   ; Transform input.
   r_in_t = imsl_fftcomp(r_in_cvt)
   i_in_t = imsl_fftcomp(i_in_cvt)
   ;
   ; Get space for returned arrays.
   r_out = dblarr(n, /NOZERO)
   i_out = dblarr(n, /NOZERO)
   ;
   ; If the length of the transform is less than 4 then compute
   ; it in a straight forward way.
   IF (n EQ 2) THEN BEGIN
      r_out = r_in_t
      i_out = i_in_t
      RETURN
   END
   IF (n EQ 3) THEN BEGIN
      r_out(0) = r_in_t(0)
      r_out(1) = r_in_t(1) - i_in_t(2)
      r_out(2) = r_in_t(1) + i_in_t(2)
      i_out(0) = i_in_t(0)
      i_out(1) = i_in_t(1) - r_in_t(2)
      i_out(2) = i_in_t(1) + r_in_t(2)
      RETURN
   END
   ;
   ; Length of sequence to be transformed is at least 4.
   ;
   ; Compute indexes needed to get the complex transform
   ; from the two real transforms.
   IF (n MOD 2) THEN idx = lindgen(n2) ELSE idx = lindgen(n2-1)
   ;
   n2p1 = n2+1
   idxp1 = idx+1
   idx1 = 2*idx+1
   idx2 = idx1+1
   ; The following two asignments are only usefull in the even case.
   ; It is here to avoid the 'IF' statement.  We do it anyway, and
   ; if n is not even, it gets overwritten by the correct value.
   r_out(n2) = r_in_t(n-1)
   i_out(n2) = i_in_t(n-1)
   r_out(0) = r_in_t(0)
   i_out(0) = i_in_t(0)
   IF KEYWORD_SET(backward) THEN BEGIN
      r_out(idxp1) = r_in_t(idx1) + i_in_t(idx2)
      r_out(n2p1:*) = reverse(r_in_t(idx1) - i_in_t(idx2))
      i_out(idxp1) =  i_in_t(idx1) - r_in_t(idx2)
      i_out(n2p1:*) = reverse(r_in_t(idx2) + i_in_t(idx1))
   END ELSE BEGIN
      r_out(idxp1) = r_in_t(idx1) - i_in_t(idx2)
      r_out(n2p1:*) = reverse(r_in_t(idx1) + i_in_t(idx2))
      i_out(idxp1) = r_in_t(idx2) + i_in_t(idx1)
      i_out(n2p1:*) = reverse(i_in_t(idx1) - r_in_t(idx2))
   END


END


