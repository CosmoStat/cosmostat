; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_ode.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_ode, t, $                       ;INPUT 1-D array: floating point 
                y, $                       ;INPUT 1-D array: floating point 
                f, $                       ;INPUT Scalar STRING
                floor=floor, $             ;INPUT Scalar floating point
                hinit=hinit, $             ;INPUT Scalar floating point
                hmax=hmax, $               ;INPUT Scalar floating point
                hmin=hmin, $               ;INPUT Scalar floating point
                jacobian=jacobian, $       ;INPUT Scalar STRING
                max_evals=max_evals, $     ;INPUT Scalar LONG
                max_ord=max_ord, $         ;INPUT Scalar LONG
                max_steps=max_steps, $     ;INPUT Scalar LONG
                method=method, $           ;INPUT Scalar LONG
                miter=miter, $             ;INPUT Scalar LONG
                norm=norm, $               ;INPUT Scalar LONG
                n_evals=n_evals, $         ;OUTPUT Scalar LONG 
                n_jevals=n_jevals, $       ;OUTPUT Scalar LONG 
                n_steps=n_steps, $         ;OUTPUT Scalar LONG 
                r_k_v=r_k_v, $             ;INPUT Scalar ON/OFF flag
                scale=scale, $             ;INPUT Scalar floating point
                tolerance=tolerance, $     ;INPUT Scalar floating point
                double=double              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - t is a 1-D array. (n_t_vals = length.)
   ; - y is a 1-D array. (neq = length.)
   ; - f must be a scalar string.
   ; - If R_K_V is present, then the following keywords can't be present.
   ;      JACOBIAN, N_JEVALS, METHOD, MAX_ORD, MITER.
   ; - If JACOBIAN present, it must be a scalar string.
   ; - If JACOBIAN is not present, then N_JEVALS can't be present.
   ;                       
   nargs = n_params()
   IF (nargs NE 3) THEN message, 'Incorrect number of arguments.'
   size_t = IMSL_SIZE(t)
   size_y = IMSL_SIZE(y)
   IF (size_t(0) NE 1) THEN message, 'T must be a 1-D array.'
   n_t_vals = IMSL_LONG(size_t(1))
   IF (size_y(0) NE 1) THEN message, 'Y must be a 1-D array.'
   neq = IMSL_LONG(size_y(1))
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   IF (KEYWORD_SET(r_k_v)) THEN BEGIN
      err_str1 = 'The keywords R_K_V and '
      err_str2 = ' are mutually exclusive.'
      IF KEYWORD_SET(jacobian) THEN message, err_str1+'JACOBIAN'+err_str2
      IF KEYWORD_SET(n_jevals) THEN message, err_str1+'N_JEVALS'+err_str2
      IF KEYWORD_SET(method) THEN message, err_str1+'METHOD'+err_str2
      IF KEYWORD_SET(max_ord) THEN message, err_str1+'MAX_ORD'+err_str2
      IF KEYWORD_SET(miter) THEN message, err_str1+'MITER'+err_str2
   END
   IF (KEYWORD_SET(jacobian)) THEN BEGIN
   size_jac = IMSL_SIZE(jacobian)
   IF ((N_ELEMENTS(jacobian) NE 1) OR (size_jac(N_ELEMENTS(size_jac)-2) NE 7)) THEN $
     message, 'JACOBIAN must be a scalar string.'
   END ELSE BEGIN
      IF (ARG_PRESENT(n_jevals)) THEN $
        message, 'The keyword N_JEVALS is valid only if JACOBIAN is also present.'
   END
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_t(N_ELEMENTS(size_t)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(max_evals)) THEN max_evals_cvt = IMSL_LONG(max_evals(0))
   IF (KEYWORD_SET(max_ord)) THEN max_ord_cvt = IMSL_LONG(max_ord(0))
   IF (KEYWORD_SET(max_steps)) THEN max_steps_cvt = IMSL_LONG(max_steps(0)) $
     ELSE max_steps_cvt = IMSL_LONG(500)
   IF (KEYWORD_SET(method)) THEN method_cvt = IMSL_LONG(method(0))
   IF (KEYWORD_SET(miter)) THEN miter_cvt = IMSL_LONG(miter(0))
   IF (KEYWORD_SET(norm)) THEN norm_cvt = (IMSL_LONG(norm))(0)
   IF (KEYWORD_SET(r_k_v)) THEN r_k_v_cvt = IMSL_1
   ; Output LONG keyword(s)
   ;
   ; In order to make things a little more simple, we
   ; will alway ask for n_step, and n_evals.
   ; We will only return them to the user if they are requested.
   n_evals_spc = IMSL_LONARR(n_t_vals)
   n_steps_spc = IMSL_LONARR(n_t_vals)
   IF (ARG_PRESENT(n_jevals)) THEN n_jevals_spc = IMSL_LONARR(n_t_vals)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(neq, n_t_vals)
      t_cvt = double(t)
      y_cvt = double(y)
      IF (KEYWORD_SET(floor)) THEN floor_cvt = double(floor(0))
      IF (KEYWORD_SET(hinit)) THEN hinit_cvt = double(hinit(0))
      IF (KEYWORD_SET(hmax)) THEN hmax_cvt = double(hmax(0))
      IF (KEYWORD_SET(hmin)) THEN hmin_cvt = double(hmin(0)) ELSE hmin_cvt = double(0.0)
      IF (KEYWORD_SET(scale)) THEN scale_cvt = double(scale(0))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0)) $
        ELSE tolerance_cvt = double(0.001)
   END ELSE BEGIN
      result = fltarr(neq, n_t_vals)
      t_cvt = float(t)
      y_cvt = float(y)
      IF (KEYWORD_SET(floor)) THEN floor_cvt = float(floor(0))
      IF (KEYWORD_SET(hinit)) THEN hinit_cvt = float(hinit(0))
      IF (KEYWORD_SET(hmax)) THEN hmax_cvt = float(hmax(0))
      IF (KEYWORD_SET(hmin)) THEN hmin_cvt = float(hmin(0)) ELSE hmin_cvt = float(0.0)
      IF (KEYWORD_SET(scale)) THEN scale_cvt = float(scale(0))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0)) $
        ELSE tolerance_cvt = float(0.001)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_171, type, err_status, t_cvt, $
                              y_cvt, $
                              f, $
                              n_t_vals, $
                              neq, $
                              floor_cvt, $
                              hinit_cvt, $
                              hmax_cvt, $
                              hmin_cvt, $
                              jacobian, $
                              max_evals_cvt, $
                              max_ord_cvt, $
                              max_steps_cvt, $
                              method_cvt, $
                              miter_cvt, $
                              norm_cvt, $
                              n_evals_spc, $
                              n_jevals_spc, $
                              n_steps_spc, $
                              r_k_v_cvt, $
                              scale_cvt, $
                              tolerance_cvt, $
                              result
   
   IF (ARG_PRESENT(n_evals)) THEN n_evals = n_evals_spc
   IF (ARG_PRESENT(n_jevals)) THEN n_jevals = n_jevals_spc
   IF (ARG_PRESENT(n_steps)) THEN n_steps = n_steps_spc
 RETURN, result
END

