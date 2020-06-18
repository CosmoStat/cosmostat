; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_nonlinopt.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION IMSL_Nonlinopt,   f,  $                   ;INPUT Scalar STRING
                      n_parameters,  $        ;INPUT Scalar LONG
                      x,  $                   ;INPUT 2-D array floating point
                      y,  $                   ;INPUT 1-D array floating point
                      a_matrix=a_matrix,  $   ;INPUT 2-D array floating point
                      b=b,  $                 ;INPUT 1-D array floating point
                      xlb=xlb,  $             ;INPUT 1-D array floating point
                      xub=xub,  $             ;INPUT 1-D array floating point
                      meq=meq,  $             ;INPUT Scalar LONG 
                      weights=weights, $             ;INPUT 1-D array: floating point 
                      frequencies=frequencies, $     ;INPUT 1-D array: floating point 
                      theta_guess=theta_guess,  $    ;INPUT 1-D array floating point
                      max_sse_evals=max_sse_evals,  $  ;INPUT Scalar LONG 
                      jacobian=jacobian,  $          ;INPUT Scalar STRING
                      acc=acc,  $                    ;INPUT Scalar floating point
                      sse=sse,  $                    ;OUTPUT Scalar floating point
                      stop_info=stop_info,  $        ;OUTPUT Scalar LONG
                      num_active=num_active,  $      ;OUTPUT Scalar LONG
                      active_const=active_const,  $  ;OUTPUT 1-D array LONG
                      lagrange_mult=lagrange_mult, $ ;OUTPUT 1-D array floating point
                      predicted=predicted, $         ;OUTPUT 1-D array floating point
                      residual=residual, $           ;OUTPUT 1-D array floating point
                      double  =  double              ;INPUT Scalar ON/OFF flag
   
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - F must be a scalar string.
   ; - X must be a 2-D array. set nobs and n_indep based on A.
   ; - Y must be a 1-D array of length nobs
   ; - If supplied, THETA_GUESS must be a 1D array of size N_PARAMETERS.
   ; - If JACOBIAN is supplied, it must be a scalar string.
   ; - If supplied, A_MATRIX and B must be used together.
   ; - If supplied, A_MATRIX must be a 2D array of size N_CON by N_PARAMETERS. (set NCON here)
   ; - If supplied, B must be a 1D array of size N_CON
   ; - If supplied, XLB must be a 1D array of size N_PARAMETERS.
   ; - If supplied, XUB must be a 1D array of size N_PARAMETERS.
   ; - If supplied, FREQUENCIES must be a 1D array of size NOBS.
   ; - If supplied, WEIGHTS must be a 1D array of size NOBS.
   ;                       
   nargs = n_params()
   IF (nargs NE 4) THEN message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'

   n_parameters_cvt = IMSL_LONG(n_parameters(0))
   size_x = IMSL_SIZE(x)
   IF (size_x(0) EQ 1) THEN BEGIN
      nobs  =  size_x(1)
      n_indep  =  IMSL_1
   END ELSE IF (size_x(0) EQ 2) THEN BEGIN
      nobs  =  size_x(1)
      n_indep  =  size_x(2)
   END ELSE MESSAGE,  'X is not the correct size.'
   
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN message, 'Y must be a 1-D array.'
   IF (N_ELEMENTS(y) NE nobs) THEN message, 'Y is not the correct size'

   IF (KEYWORD_SET(theta_guess)) THEN BEGIN
      size_tg = IMSL_SIZE(theta_guess)
      IF (size_tg(0) NE 1) THEN message, 'THETA_GUESS must be a 1-D array.'
      IF (size_tg(1) NE n_parameters_cvt) THEN message, 'THETA_GUESS is not the correct size'
   END

   IF (KEYWORD_SET(jacobian)) THEN BEGIN
      size_jacobian = IMSL_SIZE(jacobian)
      IF ((N_ELEMENTS(jacobian) NE 1) OR (size_jacobian(N_ELEMENTS(size_jacobian)-2) NE 7)) THEN $
        message, 'JACOBIAN must be a scalar string.'
   END
   n_con  =  IMSL_0
   IF ((KEYWORD_SET(a_matrix)+KEYWORD_SET(B)) EQ 1) THEN $
     MESSAGE,  "Keywords A_MATRIX and B must be used together."
   IF (KEYWORD_SET(a_matrix)) THEN BEGIN
      size_a_matrix = IMSL_SIZE(a_matrix)
      IF (size_a_matrix(0) EQ 1) THEN BEGIN
         n_con  =  size_a_matrix(1)
         IF (n_parameters_cvt NE 1) THEN MESSAGE, "Keyword A_MATRIX is no the correct size."
      END ELSE IF (size_a_matrix(0) EQ 2) THEN BEGIN
         n_con  =  size_a_matrix(1)
         IF (n_parameters_cvt NE size_a_matrix(2)) THEN MESSAGE,  $
           "Keyword A_MATRIX is not the correct size."
      END ELSE MESSAGE, "Keyword A_MATRIX is not the correct size."
      size_b = IMSL_SIZE(b)
      IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
      IF (N_ELEMENTS(b) NE n_con) THEN message, 'B is not the correct size'

   END
   IF (KEYWORD_SET(xlb)) THEN BEGIN 
      size_xlb = IMSL_SIZE(xlb)
      IF (size_xlb(0) NE 1) THEN message, 'XLB must be a 1-D array.'
      IF (N_ELEMENTS(xlb) NE n_parameters_cvt) THEN MESSAGE,  'XLB is not the correct size'
   END
   IF (KEYWORD_SET(xub)) THEN BEGIN 
      size_xub = IMSL_SIZE(xub)
      IF (size_xub(0) NE 1) THEN message, 'XUB must be a 1-D array.'
      IF (N_ELEMENTS(xub) NE n_parameters_cvt) THEN MESSAGE,  'XUB is not the correct size'
   END
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF (size_weights(0) NE 1) THEN message, 'WEIGHTS must be a 1-D array.'
      IF (N_ELEMENTS(weights) NE n_parameters_cvt) THEN MESSAGE,  'WEIGHTS is not the correct size'
   END
   IF (KEYWORD_SET(frequencies)) THEN BEGIN 
      size_frequencies = IMSL_SIZE(frequencies)
      IF (size_frequencies(0) NE 1) THEN message, 'FREQUENCIES must be a 1-D array.'
      IF (N_ELEMENTS(frequencies) NE n_parameters_cvt) THEN MESSAGE,  'FREQUENCIES is not the correct size'
   END

   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input/Output LONG keyword(s)
   ;
   num_active_spc = IMSL_0
   IF (KEYWORD_SET(meq)) THEN meq_cvt = IMSL_LONG(meq(0)) ELSE meq_cvt = IMSL_0
   IF (KEYWORD_SET(max_sse_evals)) THEN max_sse_evals_cvt = IMSL_LONG(max_sse_evals(0))
   IF (n_con GT 0) THEN active_const_spc = IMSL_LONARR(n_con)
   IF (n_con GT 0) THEN num_active_spc = IMSL_0
   IF (ARG_PRESENT(stop_info)) THEN stop_info_spc = IMSL_0
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result  =  DBLARR(n_parameters_cvt)
      ; Input
      IF (n_indep GT 1) THEN x_cvt = TRANSPOSE(DOUBLE(x)) ELSE x_cvt = DOUBLE(x)
      y_cvt = DOUBLE(y)
      IF (KEYWORD_SET(theta_guess)) THEN theta_guess_cvt = double(theta_guess)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = double(frequencies)
      IF (KEYWORD_SET(xlb)) THEN xlb_cvt = double(xlb)
      IF (KEYWORD_SET(xub)) THEN xub_cvt = double(xub)
      IF (KEYWORD_SET(b)) THEN b_cvt = double(b)
      IF (KEYWORD_SET(a_matrix)) THEN BEGIN
         IF (n_parameters_cvt GT 1) THEN a_matrix_cvt  =  TRANSPOSE(DOUBLE(a_matrix)) $
           ELSE a   =   DOUBLE(a_matrix)
      END      
      IF (KEYWORD_SET(acc)) THEN acc_cvt  =  DOUBLE(acc(0))
      ; Output
      IF (n_con GT 0) THEN lagrange_m_spc = DBLARR(n_con)
      IF (ARG_PRESENT(predicted)) THEN predicted_spc = DBLARR(nobs)
      IF (ARG_PRESENT(residual)) THEN residual_spc = DBLARR(nobs)
      IF (ARG_PRESENT(sse)) THEN sse_spc = double(0.0)
   END ELSE BEGIN
      result  =  FLTARR(n_parameters_cvt)
      ; Input
      IF (n_indep GT 1) THEN x_cvt = TRANSPOSE(FLOAT(x)) ELSE x_cvt = FLOAT(x)
      y_cvt = FLOAT(y)
      IF (KEYWORD_SET(theta_guess)) THEN theta_guess_cvt = float(theta_guess)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      IF (KEYWORD_SET(frequencies)) THEN frequencies_cvt = float(frequencies)
      IF (KEYWORD_SET(xlb)) THEN xlb_cvt = float(xlb)
      IF (KEYWORD_SET(xub)) THEN xub_cvt = float(xub)
      IF (KEYWORD_SET(b)) THEN b_cvt = float(b)
      IF (KEYWORD_SET(a_matrix)) THEN BEGIN
         IF (n_parameters_cvt GT 1) THEN a_matrix_cvt  =  TRANSPOSE(FLOAT(a_matrix)) $
           ELSE a = FLOAT(a_matrix)
      END      
      IF (KEYWORD_SET(acc)) THEN acc_cvt  =  FLOAT(acc(0))
      ; Output
      IF (n_con GT 0)THEN  lagrange_m_spc = FLTARR(n_con)
      IF (ARG_PRESENT(predicted)) THEN predicted_spc = FLTARR(nobs)
      IF (ARG_PRESENT(residual)) THEN residual_spc = FLTARR(nobs)
      IF (ARG_PRESENT(sse)) THEN sse_spc = float(0.0)
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_283, type, err_status, $
                         f, $
                         jacobian, $
                         nobs, $
                         n_con, $
                         n_indep, $
                         n_parameters_cvt, $
                         meq_cvt, $
                         max_sse_evals_cvt, $
                         x_cvt, $
                         y_cvt, $
                         theta_guess_cvt, $
                         weights_cvt, $
                         frequencies_cvt, $
                         xlb_cvt, $
                         xub_cvt, $
                         a_matrix_cvt, $
                         b_cvt, $
                         acc_cvt, $
                         active_const_spc, $
                         num_active_spc, $
                         lagrange_mult_spc, $
                         stop_info_spc, $
                         predicted_spc, $
                         residual_spc, $
                         sse_spc, $
                         result

   IF (ARG_PRESENT(predicted)) THEN predicted = predicted_spc
   IF (ARG_PRESENT(residual)) THEN residual = residual_spc
   IF (ARG_PRESENT(sse)) THEN sse = sse_spc
   num_active = num_active_spc
   IF (num_active GT 0) THEN BEGIN 
     IF (ARG_PRESENT(active_const)) THEN $
       active_const = active_const_spc(0:num_active-1)
     IF (ARG_PRESENT(lagrange_mult)) THEN $
       lagrange_mult = lagrange_m_spc(0:num_active-1)
  END
 RETURN, result
END
