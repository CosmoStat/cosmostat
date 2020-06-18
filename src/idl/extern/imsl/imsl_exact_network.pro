; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_exact_network.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_exact_network, table, $     ;INPUT 2-D array: floating point 
           wk_params=wk_params, $         ;INPUT 1-D array: LONG
           approx_params=approx_params, $ ;INPUT 1-D array: floating point
           no_approx=no_approx, $         ;INPUT Scalar ON/OFF flag
           n_attempts=n_attempts,  $      ;OUTPUT Scalar LONG
           p_value=p_value,  $            ;OUTPUT Scalar floating point
           prob_table=prob_table,  $      ;OUTPUT Scalar floating point
           double=double                  ;INPUT Scalar ON/OFF flag
           
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - TABLE must be either a 1-D or 2-D array. (m x n)
   ;  - if WK_PARAMS is present, it must be a 1-D array of length 3.
   ;  - if APPROX_PARAMS is present, it must be a 1-D array of length 3.
   ;          
   nargs = n_params()
   IF (nargs NE 1)  THEN message, 'Incorrect number of arguments.'
     
   size_table = IMSL_LONG(size(table))
   IF ((size_table(0) NE 1) AND (size_table(0) NE 2)) THEN $
     message, 'TABLE must be a 1-D or 2-D array.'
   m = IMSL_LONG(size_table(1))
   IF (size_table(0) EQ 2) THEN n  =  IMSL_LONG(size_table(2)) ELSE n  =  IMSL_1
   IF (KEYWORD_SET(wk_params)) THEN BEGIN
      size_wk_params = IMSL_SIZE(wk_params)
      IF (size_wk_params(0) NE 1) THEN $
         message, "The input array, WK_PARAMS, must be 1-D"
      IF (size_wk_params(1) NE 3) THEN $ 
         message, "The length of the WK_PARAMS array  be 3."
   END
   If (KEYWORD_SET(approx_params)) THEN BEGIN
      size_approx_params = IMSL_SIZE(approx_params)
      IF (size_approx_params(0) NE 1) THEN $
         message, "The input array, APPROX_PARAMS, must be 1-D"
      IF (size_approx_params(1) NE 3) THEN $ 
         message, "The length of the APPROX_PARAMS array  be 3."
   END   
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_table(N_ELEMENTS(size_table)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   n_attempts_spc    =    IMSL_0
   workspace = IMSL_LONARR(4)
   IF (KEYWORD_SET(wk_params)) THEN workspace(0:2)  =  wk_params(0:2)  ELSE workspace(0:2)  =  [100, 3000, 10]
   IF (KEYWORD_SET(no_approx)) THEN no_approx_cvt = 1 ELSE no_approx_cvt = 0
   ; Floating point arguments and keywords
   ;
   IF (NOT KEYWORD_SET(approx_params)) THEN approx_params  =  [5., 80, 1]
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = DOUBLE(0.)
      p_value_spc = double(0.)
      prob_table_spc = double(0.)
   END ELSE BEGIN
      result = FLOAT(0.)
      p_value_spc = float(0.)
      prob_table_spc = float(0.)
   END
   table_cvt = imsl_cvt_arr(table, type)
   approx_params_cvt = imsl_cvt_arr(approx_params,     type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_254,    type,    err_status,    table_cvt,    m,    n,  $
     p_value_spc, prob_table_spc,  n_attempts_spc, $
     approx_params_cvt, workspace, no_approx_cvt, $
     result
   ;
   p_value = p_value_spc
   prob_table = prob_table_spc
   n_attempts = workspace(3)
   ;
   ; Return
   RETURN, result
END
