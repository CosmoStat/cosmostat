; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_intfcnhyper.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_intfcnhyper, f, $                      ;INPUT Scalar STRING
                   a, $                           ;INPUT 1-D arry: floating point
                   b, $                           ;INPUT 1-D arry: floating point
                   double=double, $               ;INPUT Scalar ON/OFF flag
                   err_abs=err_abs, $             ;INPUT Scalar floating point
                   err_rel=err_rel, $             ;INPUT Scalar floating point
                   err_est=err_est, $             ;OUTPUT Scalar floating point
                   max_evals=max_evals            ;INPUT Scalar LONG
@imsl_init.pro   
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - A must be a 1-D array (with ndim elements)
   ; - B must be a 1-D array (with ndim elements)
   ;                       
   nargs = n_params()
   IF (nargs NE 3) THEN message, 'Incorrect number of arguments.'
   
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 1) THEN message, 'A must be a 1-D array.'
   ndim = IMSL_N_ELEMENTS(a)
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
   IF (N_ELEMENTS(b) NE ndim) THEN message, 'B must be the same size as A.'
   ;
   ;Select the precision to use.
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE ELSE type = TYP_FLOAT
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Long args/keywords
   IF (KEYWORD_SET(max_evals)) THEN max_evals_cvt = IMSL_LONG(max_evals(0))
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      tmp = imsl_machine(/double)
      err_abs_cvt = sqrt(tmp.(3))
      err_rel_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = double(err_abs(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = double(err_rel(0))
      result = double(0.0)
      a_cvt = double(a)
      b_cvt = double(b)
      IF (ARG_PRESENT(err_est)) THEN err_est_spc = double(0.0)
   END ELSE BEGIN
      tmp = imsl_machine(/float)
      err_abs_cvt = sqrt(tmp.(3))
      err_rel_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = float(err_abs(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = float(err_rel(0))
      result = float(0.0)
      a_cvt = float(a)
      b_cvt = float(b)
      IF (ARG_PRESENT(err_est)) THEN err_est_spc = float(0.0)
   END
   
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_151, type, err_status, $
                              f, $
                              a_cvt, $
                              b_cvt, $
                              ndim, $
                              max_evals_cvt, $
                              err_abs_cvt, $
                              err_est_spc, $
                              err_rel_cvt, $
                              result
   IF (ARG_PRESENT(err_est)) THEN err_est = err_est_spc
 RETURN, result
END

