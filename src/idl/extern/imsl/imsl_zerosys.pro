; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_zerosys.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Zerosys, f, $          ;INPUT Scalar STRING
                n, $                       ;INPUT Scalar LONG
                double=double, $           ;INPUT Scalar ON/OFF flag
                jacobian=jacobian, $       ;INPUT Scalar STRING
                itmax=itmax, $             ;INPUT Scalar LONG
                err_rel=err_rel, $         ;INPUT Scalar floating point
                xguess=xguess, $           ;INPUT 1-D array: floating point
                fnorm=fnorm                ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - N is converted to TYP_MEMINT, and must be positive.
   ; - XGUESS must be a 1-D array with N_ROOTS elements (if N gt 1)
   ;   or either a scalar or 1-D array of length 1 if N == 1.
   ; - If JACOBIAN present, it must be a scalar string.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   n_cvt = IMSL_LONG(n(0))
   IF (n_cvt LT 1) THEN message, 'N must be positive.'
   IF (KEYWORD_SET(xguess)) THEN BEGIN
      size_xguess = IMSL_SIZE(xguess)
      IF (n_cvt GT 1) THEN BEGIN
         IF ((size_xguess(0) NE 1) OR (N_ELEMENTS(xguess) NE n_cvt)) $
           THEN message, 'XGUESS is not the correct size.'
      END ELSE BEGIN
         IF (N_ELEMENTS(xguess) NE n_cvt) $
           THEN message, 'XGUESS is not the correct size.'
      END
   END
   IF (KEYWORD_SET(jacobian)) THEN BEGIN
      size_jac = IMSL_SIZE(jacobian)
      IF ((N_ELEMENTS(jacobian) NE 1) OR $
          (size_jac(N_ELEMENTS(size_jac)-2) NE 7)) THEN $
        message, 'JACOBIAN must be a scalar string.'
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
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n_cvt)
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = double(err_rel(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess) $
        ELSE xguess_cvt = dblarr(n_cvt)
      IF (n_cvt EQ 1) THEN xguess_cvt = xguess_cvt(0)
      IF (ARG_PRESENT(fnorm)) THEN fnorm_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(n_cvt)
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = float(err_rel(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess) $
        ELSE xguess_cvt = fltarr(n_cvt)
      IF (n_cvt EQ 1) THEN xguess_cvt = xguess_cvt(0)
      IF (ARG_PRESENT(fnorm)) THEN fnorm_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_198, type, err_status, $
                              f, $
                              n_cvt, $
                              jacobian, $
                              itmax_cvt, $
                              err_rel_cvt, $
                              xguess_cvt, $
                              fnorm_spc, $
                              result

   IF (ARG_PRESENT(fnorm)) THEN fnorm = fnorm_spc
 RETURN, result
END

