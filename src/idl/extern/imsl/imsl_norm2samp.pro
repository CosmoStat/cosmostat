; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_norm2samp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_norm2samp, x1, $                              ;INPUT 1-D array: floating point
                    x2, $                              ;INPUT 1-D array: floating point
                    Conf_mean=conf_mean, $             ;INPUT Scalar floating point
                    T_test_null_hyp=t_test_null_hyp, $ ;INPUT Scalar floating point
                    Conf_var=conf_var, $               ;INPUT Scalar floating point
                    Chi_sq_null_hyp=chi_sq_null_hyp, $ ;INPUT Scalar floating point
                    Double=double, $                   ;INPUT Scalar ON/OFF flag
                    Mean_x1 = mean_x1, $               ;OUTPUT 1-D array: floating point
                    Mean_x2=mean_x2, $                 ;OUTPUT 1-D array: floating point
                    Ci_diff_eq_var=ci_diff_eq_var, $   ;OUTPUT 1-D array: floating point
                    Ci_diff_ne_var=ci_diff_ne_var, $   ;OUTPUT 1-D array: floating point
                    T_test_eq_var=t_test_eq_var, $     ;OUTPUT 1-D array: floating point
                    T_test_ne_var=t_test_ne_var, $     ;OUTPUT 1-D array: floating point
                    Ci_comm_var=ci_comm_var, $         ;OUTPUT 1-D array: floating point
                    Chi_sq_test=chi_sq_test, $         ;OUTPUT 1-D array: floating point
                    Ci_ratio_var=ci_ratio_var, $       ;OUTPUT 1-D array: floating point  
                    F_test=f_test, $                   ;OUTPUT 1-D array: floating point
                    Pooled_var=pooled_var, $           ;OUTPUT Scalar floating point
                    Stdev_x1=stdev_x1, $               ;OUTPUT Scalar floating point
                    Stdev_x2=stdev_x2                  ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are simple arrays.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_x1 = IMSL_SIZE(x1)
   IF (size_x1(0) NE 1) THEN BEGIN
      message, "X1 must be a 1-D array."
   END
   size_x2 = IMSL_SIZE(x2)
   IF (size_x2(0) NE 1) THEN BEGIN
      message, "X2 must be a 1-D array."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x1(N_ELEMENTS(size_x1)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_x2(N_ELEMENTS(size_x2)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   nobs1 = size_x1(N_ELEMENTS(size_x1)-1)
   nobs2 = size_x2(N_ELEMENTS(size_x1)-1)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      ; 
      ; Input 
      x1_cvt = double(x1)
      x2_cvt = double(x2)
      IF (KEYWORD_SET(Conf_mean) EQ TRUE) THEN $
        Conf_mean_cvt = (double(Conf_mean))(0)
      IF (KEYWORD_SET(T_test_null_hyp) EQ TRUE) THEN $
        T_test_null_cvt = (double(T_test_null_hyp))(0)
      IF (KEYWORD_SET(Conf_var) EQ TRUE) THEN $
        Conf_var_cvt = (double(Conf_var))(0)
      IF (KEYWORD_SET(Chi_sq_null_hyp) EQ TRUE) THEN $
        Chi_sq_null_cvt = (double(Chi_sq_null_hyp))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(Mean_x1) EQ TRUE) THEN $
        Mean_x1_spc = double(0.0)
      IF (ARG_PRESENT(Mean_x2) EQ TRUE) THEN $
        Mean_x2_spc = double(0.0)
      IF (ARG_PRESENT(Ci_diff_eq_var) EQ TRUE) THEN $
        Ci_diff_eq_spc = dblarr(2)
      IF (ARG_PRESENT(Ci_diff_ne_var) EQ TRUE) THEN $
        Ci_diff_ne_spc = dblarr(2)
      IF (ARG_PRESENT(T_test_eq_var) EQ TRUE) THEN $
        T_test_eq_spc = dblarr(3)
      IF (ARG_PRESENT(T_test_ne_var) EQ TRUE) THEN $
        T_test_ne_spc = dblarr(3)
      IF (ARG_PRESENT(Pooled_var) EQ TRUE) THEN $
        Pooled_var_spc = double(0.0)
      IF (ARG_PRESENT(Ci_comm_var) EQ TRUE) THEN $
        Ci_comm_var_spc = dblarr(2)
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test_spc = dblarr(3)
      IF (ARG_PRESENT(Stdev_x1) EQ TRUE) THEN $
        stdev_x1_spc = double(0.0)
      IF (ARG_PRESENT(Stdev_x2) EQ TRUE) THEN $
        stdev_x2_spc = double(0.0)
      IF (ARG_PRESENT(Ci_ratio_var) EQ TRUE) THEN $
        Ci_ratio_spc = dblarr(2)
      IF (ARG_PRESENT(F_test) EQ TRUE) THEN $
        F_test_spc = dblarr(4)
   END ELSE BEGIN
      result = float(0.0)
      ; 
      ; Input 
      x1_cvt = float(x1)
      x2_cvt = float(x2)
      IF (KEYWORD_SET(Conf_mean) EQ TRUE) THEN $
        Conf_mean_cvt = (float(Conf_mean))(0)
      IF (KEYWORD_SET(T_test_null_hyp) EQ TRUE) THEN $
        T_test_null_cvt = (float(T_test_null_hyp))(0)
      IF (KEYWORD_SET(Conf_var) EQ TRUE) THEN $
        Conf_var_cvt = (float(Conf_var))(0)
      IF (KEYWORD_SET(Chi_sq_null_hyp) EQ TRUE) THEN $
        Chi_sq_null_cvt = (float(Chi_sq_null_hyp))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(Mean_x1) EQ TRUE) THEN $
        Mean_x1_spc = float(0.0)
      IF (ARG_PRESENT(Mean_x2) EQ TRUE) THEN $
        Mean_x2_spc = float(0.0)
      IF (ARG_PRESENT(Ci_diff_eq_var) EQ TRUE) THEN $
        Ci_diff_eq_spc = fltarr(2)
      IF (ARG_PRESENT(Ci_diff_ne_var) EQ TRUE) THEN $
        Ci_diff_ne_spc = fltarr(2)
      IF (ARG_PRESENT(T_test_eq_var) EQ TRUE) THEN $
        T_test_eq_spc = fltarr(3)
      IF (ARG_PRESENT(T_test_ne_var) EQ TRUE) THEN $
        T_test_ne_spc = fltarr(3)
      IF (ARG_PRESENT(Pooled_var) EQ TRUE) THEN $
        Pooled_var_spc = float(0.0)
      IF (ARG_PRESENT(Ci_comm_var) EQ TRUE) THEN $
        Ci_comm_var_spc = fltarr(2)
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test_spc = fltarr(3)
      IF (ARG_PRESENT(Stdev_x1) EQ TRUE) THEN $
        stdev_x1_spc = float(0.0)
      IF (ARG_PRESENT(Stdev_x2) EQ TRUE) THEN $
        stdev_x2_spc = float(0.0)
      IF (ARG_PRESENT(Ci_ratio_var) EQ TRUE) THEN $
        Ci_ratio_spc = fltarr(2)
      IF (ARG_PRESENT(F_test) EQ TRUE) THEN $
        F_test_spc = fltarr(4)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_168, type, err_status, x1_cvt, x2_cvt, nobs1, nobs2, $
                          Conf_mean_cvt, $
                          T_test_null_cvt, $
                          Conf_var_cvt, $
                          Chi_sq_null_cvt, $
                          Mean_x1_spc , $
                          Mean_x2_spc, $
                          Ci_diff_eq_spc, $
                          Ci_diff_ne_spc, $
                          T_test_eq_spc, $
                          T_test_ne_spc, $
                          Pooled_var_spc, $
                          Ci_comm_var_spc, $
                          Chi_sq_test_spc, $
                          Stdev_x1_spc, $
                          Stdev_x2_spc, $
                          Ci_ratio_spc, $
                          F_test_spc, $
                          result
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test=Chi_sq_test_spc
      IF (ARG_PRESENT(Ci_comm_var) EQ TRUE) THEN $
        Ci_comm_var=Ci_comm_var_spc
      IF (ARG_PRESENT(Ci_diff_eq_var) EQ TRUE) THEN $
        Ci_diff_eq_var=Ci_diff_eq_spc
      IF (ARG_PRESENT(Ci_diff_ne_var) EQ TRUE) THEN $
        Ci_diff_ne_var=Ci_diff_ne_spc
      IF (ARG_PRESENT(Ci_ratio_var) EQ TRUE) THEN $
        Ci_ratio_var=Ci_ratio_spc
      IF (ARG_PRESENT(F_test) EQ TRUE) THEN $
        F_test=F_test_spc
      IF (ARG_PRESENT(Mean_x1) EQ TRUE) THEN $
        Mean_x1=Mean_x1_spc
      IF (ARG_PRESENT(Mean_x2) EQ TRUE) THEN $
        Mean_x2=Mean_x2_spc
      IF (ARG_PRESENT(Pooled_var) EQ TRUE) THEN $
        Pooled_var=Pooled_var_spc
      IF (ARG_PRESENT(Stdev_x1) EQ TRUE) THEN $
        stdev_x1=stdev_x1_spc
      IF (ARG_PRESENT(Stdev_x2) EQ TRUE) THEN $
        stdev_x2=stdev_x2_spc
      IF (ARG_PRESENT(T_test_eq_var) EQ TRUE) THEN $
        T_test_eq_var=T_test_eq_spc
      IF (ARG_PRESENT(T_test_ne_var) EQ TRUE) THEN $
        T_test_ne_var=T_test_ne_spc
   ;
   ; Return.
   ;
   RETURN, result
END
