; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_fminv.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_fminv, f, $                    ;INPUT Scalar STRING
                n, $                      ;INPUT Scalar LONG
                grad=grad, $              ;INPUT Scalar STRING
                double = double, $        ;INPUT Scalar ON/OFF flag
                max_evals=max_evals, $    ;INPUT Scalar LONG 
                itmax=itmax, $            ;INPUT Scalar LONG 
                max_grad=max_grad, $      ;INPUT Scalar LONG 
                n_digit=n_digit, $        ;INPUT Scalar LONG 
                ihess=ihess, $            ;INPUT Scalar LONG
                fscale=fscale, $          ;INPUT Scalar floating point
                max_step=max_step, $      ;INPUT Scalar floating point
                tol_grad=tol_grad, $      ;INPUT Scalar floating point
                tol_rfcn=tol_rfcn, $      ;INPUT Scalar floating point
                tol_step=tol_step, $      ;INPUT Scalar floating point
                xguess=xguess, $          ;INPUT Scalar floating point
                xscale=xscale, $          ;INPUT 1-D array floating point
                fvalue=fvalue             ;OUTPUT Scalar floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - If GRAD is supplied,  must be a scalar string.
   ; - TOL_GRAD is only allowed if GRAD is also supplied.
   ; - MAX_GRAD is only allowed if GRAD is also supplied.
   ; - if XGUESS is supplied, it must be a 1-D array of length n.
   ; - if XSCALE is supplied, it must be a 1-D array of length n.
   ;                       
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   n_cvt = IMSL_LONG(n(0))
   IF (KEYWORD_SET(grad)) THEN BEGIN
      size_grad = IMSL_SIZE(grad)
      IF ((N_ELEMENTS(grad) NE 1) OR (size_grad(N_ELEMENTS(size_grad)-2) NE 7)) THEN $
        message, 'GRAD must be a scalar string.'
   END ELSE BEGIN
      IF (KEYWORD_SET(tol_grad)) THEN message, $
        'TOL_GRAD is only valid when GRAD is also supplied.'
      IF (KEYWORD_SET(max_grad)) THEN message, $
        'MAX_GRAD is only valid when GRAD is also supplied.'
   END
   IF (KEYWORD_SET(xguess)) THEN BEGIN
      size_xguess = IMSL_SIZE(xguess)
      IF (size_xguess(0) NE 1) THEN message, 'XGUESS must be a 1-D array.'
      IF (size_xguess(1) NE n_cvt) THEN message, 'XGUESS is not the correct size.'
   END
   IF (KEYWORD_SET(xscale)) THEN BEGIN
      size_xscale = IMSL_SIZE(xscale)
      IF (size_xscale(0) NE 1) THEN message, 'XSCALE must be a 1-D array.'
      IF (size_xscale(1) NE n_cvt) THEN message, 'XSCALE is not the correct size.'
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
   IF (KEYWORD_SET(max_evals)) THEN max_evals_cvt = IMSL_LONG(max_evals(0)) $
     ELSE max_evals_cvt = IMSL_LONG(400)
   IF (KEYWORD_SET(max_grad)) THEN max_grad_cvt = IMSL_LONG(max_grad(0)) $
     ELSE max_grad_cvt = IMSL_LONG(400)
   IF (KEYWORD_SET(n_digit)) THEN n_digit_cvt = IMSL_LONG(n_digit(0))
   IF (KEYWORD_SET(ihess)) THEN ihess_cvt = IMSL_LONG(ihess(0))
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n_cvt)
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = double(tol_grad(0))
      IF (KEYWORD_SET(fscale)) THEN fscale_cvt = double(fscale(0))
      IF (KEYWORD_SET(max_step)) THEN max_step_cvt = double(max_step(0))
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = double(tol_grad(0))
      IF (KEYWORD_SET(tol_rfcn)) THEN tol_rfcn_cvt = double(tol_rfcn(0))
      IF (KEYWORD_SET(tol_step)) THEN tol_step_cvt = double(tol_step(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess)
      IF (KEYWORD_SET(xscale)) THEN xscale_cvt = double(xscale)
      IF (ARG_PRESENT(fvalue)) THEN fvalue_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(n_cvt)
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = float(tol_grad(0))
      IF (KEYWORD_SET(fscale)) THEN fscale_cvt = float(fscale(0))
      IF (KEYWORD_SET(max_step)) THEN max_step_cvt = float(max_step(0))
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = float(tol_grad(0))
      IF (KEYWORD_SET(tol_rfcn)) THEN tol_rfcn_cvt = float(tol_rfcn(0))
      IF (KEYWORD_SET(tol_step)) THEN tol_step_cvt = float(tol_step(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess)
      IF (KEYWORD_SET(xscale)) THEN xscale_cvt = float(xscale)
      IF (ARG_PRESENT(fvalue)) THEN fvalue_spc = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_143, type, err_status, $
                         f, $
                         n_cvt, $
                         grad, $
                         max_evals_cvt, $
                         itmax_cvt, $
                         max_grad_cvt, $
                         n_digit_cvt, $
                         fscale_cvt, $
                         ihess_cvt, $
                         max_step_cvt, $
                         tol_grad_cvt, $
                         tol_rfcn_cvt, $
                         tol_step_cvt, $
                         xguess_cvt, $
                         xscale_cvt, $
                         fvalue_spc, $
                         result
   
   IF (ARG_PRESENT(fvalue)) THEN fvalue = fvalue_spc
 RETURN, result
END

