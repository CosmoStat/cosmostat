; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_norm1samp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_norm1samp, x, $                               ;INPUT 1-D array: floating point
                    Chi_sq_null_hyp=chi_sq_null_hyp, $ ;INPUT Scalar floating point
                    Conf_mean=conf_mean, $             ;INPUT Scalar floating point
                    T_null_hyp=t_null_hyp, $           ;INPUT Scalar floating point
                    Conf_var=conf_var, $               ;INPUT Scalar floating point
                    Double=double, $                   ;INPUT Scalar ON/OFF flag
                    Chi_sq_test=chi_sq_test, $         ;OUTPUT 1-D array: floating point
                    Ci_mean=ci_mean, $                 ;OUTPUT 1-D array: floating point
                    Ci_var=ci_var, $                   ;OUTPUT 1-D array: floating point
                    T_test=t_test, $                   ;OUTPUT 1-D array: floating point
                    Stdev=stdev                        ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN $
         message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN
      message, "X must be a 1-D array."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   nobs = size_x(N_ELEMENTS(size_x)-1)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      ; 
      ; Input
      x_cvt = double(x)
      IF (KEYWORD_SET(Conf_mean) EQ TRUE) THEN $
        Conf_mean_cvt = (double(Conf_mean))(0)
      IF (KEYWORD_SET(Conf_var) EQ TRUE) THEN $
        Conf_var_cvt = (double(Conf_var))(0)
      IF (KEYWORD_SET(T_null_hyp) EQ TRUE) THEN $
        T_null_hyp_cvt = (double(T_null_hyp))(0)
      IF (KEYWORD_SET(Chi_sq_null_hyp) EQ TRUE) THEN $
        Chi_sq_null_cvt = (double(Chi_sq_null_hyp))(0)
      ; 
      ; Output
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test_spc = dblarr(3)
      IF (ARG_PRESENT(Ci_mean) EQ TRUE) THEN $
        Ci_mean_spc = dblarr(2)
      IF (ARG_PRESENT(Ci_var) EQ TRUE) THEN $
        Ci_var_spc = dblarr(2)
      IF (ARG_PRESENT(Stdev) EQ TRUE) THEN $
        stdev_spc = double(0.0)
      IF (ARG_PRESENT(T_test) EQ TRUE) THEN $
        T_test_spc = dblarr(3)
   END ELSE BEGIN
      result = float(0.0)
      ; 
      ; Input
      x_cvt = float(x)
      IF (KEYWORD_SET(Conf_mean) EQ TRUE) THEN $
        Conf_mean_cvt = (float(Conf_mean))(0)
      IF (KEYWORD_SET(Conf_var) EQ TRUE) THEN $
        Conf_var_cvt = (float(Conf_var))(0)
      IF (KEYWORD_SET(T_null_hyp) EQ TRUE) THEN $
        T_null_hyp_cvt = (float(T_null_hyp))(0)
      IF (KEYWORD_SET(Chi_sq_null_hyp) EQ TRUE) THEN $
        Chi_sq_null_cvt = (float(Chi_sq_null_hyp))(0)
      ; 
      ; Output
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test_spc = fltarr(3)
      IF (ARG_PRESENT(Ci_mean) EQ TRUE) THEN $
        Ci_mean_spc = fltarr(2)
      IF (ARG_PRESENT(Ci_var) EQ TRUE) THEN $
        Ci_var_spc = fltarr(2)
      IF (ARG_PRESENT(Stdev) EQ TRUE) THEN $
        stdev_spc = float(0.0)
      IF (ARG_PRESENT(T_test) EQ TRUE) THEN $
        T_test_spc = fltarr(3)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_167, type, err_status, x_cvt, nobs, $
                          Chi_sq_null_cvt, $
                          Chi_sq_test_spc, $
                          Ci_mean_spc, $
                          Ci_var_mean_spc, $
                          Conf_mean_cvt, $
                          Conf_var_cvt, $
                          Stdev_spc, $
                          T_null_hyp_cvt, $
                          T_test_spc, $
                          result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
     Chi_sq_test = Chi_sq_test_spc
   IF (ARG_PRESENT(Ci_mean) EQ TRUE) THEN $
     Ci_mean = Ci_mean_spc
   IF (ARG_PRESENT(Ci_var) EQ TRUE) THEN $
     Ci_var = Ci_var_spc
   IF (ARG_PRESENT(Stdev) EQ TRUE) THEN $
     stdev = stdev_spc
   IF (ARG_PRESENT(T_test) EQ TRUE) THEN $
     T_test = T_test_spc
   ;
   ; Return.
   ;
   RETURN, result
END
