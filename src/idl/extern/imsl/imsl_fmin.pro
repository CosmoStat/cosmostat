; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_fmin.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_fmin, f, $                     ;INPUT Scalar STRING
                a, $                      ;INPUT Scalar floating point
                b, $                      ;INPUT Scalar floating point
                g, $                      ;INPUT Scalar STRING
                double = double, $        ;INPUT Scalar ON/OFF flag
                err_abs=err_abs, $        ;INPUT Scalar floating point
                err_rel=err_rel, $        ;INPUT Scalar floating point
                fvalue=fvalue, $          ;OUTPUT Scalar floating point
                gvalue=gvalue, $          ;OUTPUT Scalar floating point
                max_evals=max_evals, $    ;INPUT Scalar LONG 
                max_fcn=max_fcn, $        ;INPUT Scalar LONG 
                step=step, $              ;INPUT Scalar floating point
                tol_grad=tol_grad, $      ;INPUT Scalar floating point
                xguess=xguess             ;INPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - If (nargs EQ 4) g must be a scalar string.
   ; - Keyword STEP is only valid if (argc LT 4)
   ; - Keyword ERR_ABS is only valid if (argc LT 4)
   ; - Keyword FVALUE is only valid if (argc EQ 4)
   ; - Keyword GVALUE is only valid if (argc EQ 4)
   ; - Keyword ERR_REL is only valid if (argc EQ 4)
   ; - Keyword TOL_GRAD is only valid if (argc EQ 4)
   ;                       
   nargs = n_params()
   IF ((nargs NE 3) AND (nargs NE 4)) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   IF (nargs EQ 4) THEN BEGIN
      size_g = IMSL_SIZE(g)
      IF ((N_ELEMENTS(g) NE 1) OR (size_g(N_ELEMENTS(size_g)-2) NE 7)) THEN $
        message, 'GRAD must be a scalar string.'
      IF (KEYWORD_SET(step)) THEN message, $
        'STEP is not valid when GRAD is also supplied.'
      IF (KEYWORD_SET(err_abs)) THEN message, $
        'ERR_ABS is not valid when GRAD is also supplied.'
   END ELSE BEGIN
      IF (KEYWORD_SET(fvalue)) THEN message, $
        'FVALUE is only valid when GRAD is also supplied.'
      IF (KEYWORD_SET(gvalue)) THEN message, $
        'GVALUE is only valid when GRAD is also supplied.'
      IF (KEYWORD_SET(err_rel)) THEN message, $
        'ERR_REL is only valid when GRAD is also supplied.'
      IF (KEYWORD_SET(tol_grad)) THEN message, $
        'TOL_GRAD is only valid when GRAD is also supplied.'
   END
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   size_a = IMSL_SIZE(a)
   size_b = IMSL_SIZE(b)
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(itmax)) THEN itmax_cvt = IMSL_LONG(itmax(0)) ELSE itmax_cvt = IMSL_LONG(100)
   IF (KEYWORD_SET(max_evals)) THEN max_evals_cvt = IMSL_LONG(max_evals(0)) $
     ELSE max_evals_cvt = IMSL_LONG(1000)
   IF (KEYWORD_SET(max_fcn)) THEN max_evals_cvt = IMSL_LONG(max_fcn(0)) $
     ELSE max_evals_cvt = IMSL_LONG(1000)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      a_cvt = double(a(0))
      b_cvt = double(b(0))
      fvalue_spc = double(0.0)
      gvalue_spc = double(0.0)
      err_abs_cvt = double(.0001)
      tmp = imsl_machine(/double)
      err_rel_cvt = sqrt(tmp.(3))
      tol_grad_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = double(err_abs(0))
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = double(tol_grad(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = double(err_rel(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = double(xguess(0)) $
        ELSE xguess_cvt = (b_cvt + a_cvt)/2.
      IF (KEYWORD_SET(step)) THEN step_cvt = double(step(0)) $
        ELSE step_cvt = double(1.0)
   END ELSE BEGIN
      result = float(0.0)
      a_cvt = float(a(0))
      b_cvt = float(b(0))
      fvalue_spc = float(0.0)
      gvalue_spc = float(0.0)
      err_abs_cvt = float(.0001)
      tmp = imsl_machine(/float)
      err_rel_cvt = sqrt(tmp.(3))
      tol_grad_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(err_abs)) THEN err_abs_cvt = float(err_abs(0))
      IF (KEYWORD_SET(tol_grad)) THEN tol_grad_cvt = float(tol_grad(0))
      IF (KEYWORD_SET(err_rel)) THEN err_rel_cvt = float(err_rel(0))
      IF (KEYWORD_SET(xguess)) THEN xguess_cvt = float(xguess(0)) $
        ELSE xguess_cvt = (b_cvt + a_cvt)/2.
      IF (KEYWORD_SET(step)) THEN step_cvt = float(step(0)) $
        ELSE step_cvt = float(1.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_142, type, err_status, $
                              f, $
                              a_cvt, $
                              b_cvt, $
                              g, $
                              err_abs_cvt, $
                              err_rel_cvt, $
                              fvalue_spc, $
                              gvalue_spc, $
                              max_evals_cvt, $
                              step_cvt, $
                              tol_grad_cvt, $
                              xguess_cvt, $
                              result
   
   IF (ARG_PRESENT(fvalue)) THEN fvalue = fvalue_spc
   IF (ARG_PRESENT(gvalue)) THEN gvalue = gvalue_spc
 RETURN, result
END

