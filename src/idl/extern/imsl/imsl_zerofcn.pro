; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_zerofcn.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_zerofcn, f, $                   ;INPUT Scalar STRING
                double=double, $            ;INPUT Scalar ON/OFF flag
                n_roots=n_roots, $         ;INPUT Scalar LONG
                itmax=itmax, $             ;INPUT Scalar LONG
                eps=eps, $                 ;INPUT Scalar floating point
                err_abs=err_abs, $         ;INPUT Scalar floating point
                err_rel=err_rel, $         ;INPUT Scalar floating point
                eta=eta, $                 ;INPUT Scalar floating point
                xguess=xguess, $           ;INPUT 1-D array: floating point
                info=info                  ;OUTPUT 1-D array: LONG

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - N_ROOTS must be positive.
   ; - XGUESS must be a 1-D array with N_ROOTS elements (if N_ROOTS gt 1)
   ;   or either a scalar or 1-D array of length 1 if N_ROOTS == 1.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   IF (KEYWORD_SET(n_roots)) THEN BEGIN
      n_roots_cvt = IMSL_LONG(n_roots(0))
      IF (n_roots_cvt LT 1) THEN message, 'N_ROOTS must be positive.'
   END ELSE n_roots_cvt = IMSL_1
   IF (KEYWORD_SET(xguess)) THEN BEGIN
      size_xguess = IMSL_SIZE(xguess)
      IF (n_roots_cvt GT 1) THEN BEGIN
         IF ((size_xguess(0) NE 1) OR (N_ELEMENTS(xguess) NE n_roots_cvt)) $
           THEN message, 'XGUESS is not the correct size.'
      END ELSE BEGIN
         IF (N_ELEMENTS(xguess) NE n_roots_cvt) $
           THEN message, 'XGUESS is not the correct size.'
      END
   END
   ;
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(itmax)) THEN itmax_cvt = IMSL_LONG(itmax(0)) ELSE itmax_cvt = IMSL_LONG(100)
   ; Output LONG keyword(s)
   ;
   IF (ARG_PRESENT(info)) THEN info_spc = IMSL_LONARR(n_roots_cvt)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n_roots_cvt)
      IF (KEYWORD_SET(eps)) THEN eps_cvt = double(eps(0))
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = double(err_abs(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = double(err_rel(0))
      IF (KEYWORD_SET(eta)) THEN eta_cvt = double(eta(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess) $
        ELSE xguess_cvt = dblarr(n_roots_cvt)
      IF (n_roots_cvt EQ 1) THEN xguess_cvt = xguess_cvt(0)
   END ELSE BEGIN
      result = fltarr(n_roots_cvt)
      IF (KEYWORD_SET(eps)) THEN eps_cvt = float(eps(0))
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = float(err_abs(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = float(err_rel(0))
      IF (KEYWORD_SET(eta)) THEN eta_cvt = float(eta(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess) $
        ELSE xguess_cvt = fltarr(n_roots_cvt)
      IF (n_roots_cvt EQ 1) THEN xguess_cvt = xguess_cvt(0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_196, type, err_status, $
                              f, $
                              n_roots_cvt, $
                              itmax_cvt, $
                              eps_cvt, $
                              err_abs_cvt, $
                              err_rel_cvt, $
                              eta_cvt, $
                              xguess_cvt, $
                              info_spc, $
                              result

   IF (ARG_PRESENT(info)) THEN info = info_spc
 RETURN, result
END

