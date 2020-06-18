; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_multipredict.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_multipredict, info_v, $              ;INPUT structure
                       x, $                        ;INPUT 1-D or 2-D array: floating point
                       Double=double, $                ;INPUT Scalar ON/OFF flag
                       Weights=weights, $          ;INPUT 1-D array: floating point
                       Confidence=confidence, $    ;INPUT Scalar: floating point
                       Y=y, $                      ;INPUT 1-D array: floating point
                       Ci_scheffe=Ci_scheffe, $    ;OUTPUT 2-D array: floating point
                       Ci_ptw_pop_mean=Ci_ptw_pop_mean, $ ;OUTPUT 2-D array: floating point
                       Ci_ptw_new_samp=Ci_ptw_new_samp, $ ;OUTPUT 2-D array: floating point
                       leverage=leverage, $        ;OUTPUT 1-D array: floating point
                       Residual=residual, $        ;OUTPUT 1-D array: floating point
                       Std_Residual=std_residual, $;OUTPUT 1-D array: floating point
                       Del_Residual=del_residual, $;OUTPUT 1-D array: floating point
                       Cooks_d=cooks_d, $          ;OUTPUT 1-D array: floating point
                       Dffits=dffits               ;OUTPUT 1-D array: floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are 1-D or 2-D arrays.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF ((size_x(0) LT 0) OR (size_x(0) GT 2)) THEN BEGIN 
      message, "X must be a scalar, 1-D, or 2-D"
   END
   IF (size_x(0) EQ 0) THEN BEGIN
      n_predict = IMSL_1
      n_indep_vars = IMSL_1
   END ELSE IF (size_x(0) EQ 1) THEN BEGIN
      n_predict = size_x(1)
      n_indep_vars = IMSL_1
   END ELSE BEGIN     
      n_predict = size_x(1)
      n_indep_vars = size_x(2)
   END
   IF (KEYWORD_SET(weights) EQ TRUE) THEN BEGIN
     size_weights = IMSL_SIZE(weights)
     IF (size_weights(0) NE 1) THEN BEGIN
        message, "The weights array must be 1-D"
     END
     IF (size_weights(1) NE n_predict) THEN BEGIN
        message, "WEIGHTS is NOT the correct length."
     END
   END
  size_info_v = IMSL_SIZE(info_v)
  IF ((size_info_v(0) NE 1) OR (size_info_v(N_ELEMENTS(size_info_v)-2) NE 1)  ) THEN $
      message, "PREDICT_INFO must be a 1-D array of type BYTE"
  IF (size_info_v(1) LT 23) THEN $
    message, "PREDICT_INFO is not valid, it must be exactly as returned by MULTIREGRESS"
  ; 
  ; Make a number of checks concerning the keyword Y.
  IF (ARG_PRESENT(residual) AND (NOT KEYWORD_SET(Y))) THEN $
    message, "You must specify the keyword Y if you specify the keyword RESIDUAL"
  IF (ARG_PRESENT(std_residual) AND (NOT KEYWORD_SET(Y))) THEN $
    message, "You must specify the keyword Y if you specify the keyword STD_RESIDUAL"
  IF (ARG_PRESENT(del_residual) AND (NOT KEYWORD_SET(Y))) THEN $
    message, "You must specify the keyword Y if you specify the keyword DEL_RESIDUAL"
  IF (ARG_PRESENT(cooks_d) AND (NOT KEYWORD_SET(Y))) THEN $
    message, "You must specify the keyword Y if you specify the keyword COOKS_D"
  IF (ARG_PRESENT(dffits) AND (NOT KEYWORD_SET(Y))) THEN $
    message, "You must specify the keyword Y if you specify the keyword DFFITS"
  IF (KEYWORD_SET(y) EQ TRUE) THEN BEGIN
     IF (N_ELEMENTS(y) NE n_predict) THEN $
       message, "The number OF observed responses, as specified by the keyword Y, does NOT agree with the argument X"
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
   ; Result array
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      l_result = dblarr(n_predict)
      ; 
      ; Input 
      x_cvt = double(transpose(x))
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = double(weights)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (double(confidence))(0)
      IF (KEYWORD_SET(y) EQ TRUE) THEN $
        y_cvt = double(y)
      ; 
      ; Output 
      IF (ARG_PRESENT(ci_scheffe) EQ TRUE) THEN $
        ci_scheffe_spc = dblarr(n_predict, 2)
      IF (ARG_PRESENT(Ci_ptw_pop_mean) EQ TRUE) THEN $
        Ci_pop_mean_spc = dblarr(n_predict, 2)
      IF (ARG_PRESENT(Ci_ptw_new_samp) EQ TRUE) THEN $
        Ci_new_samp_spc = dblarr(n_predict, 2)
      IF (ARG_PRESENT(leverage) EQ TRUE) THEN $
        leverage_spc = dblarr(n_predict)
      IF (ARG_PRESENT(residual) EQ TRUE) THEN $
        residual_spc = dblarr(n_predict)
      IF (ARG_PRESENT(std_residual) EQ TRUE) THEN $
        std_resid_spc = dblarr(n_predict)
      IF (ARG_PRESENT(del_residual) EQ TRUE) THEN $
        del_resid_spc = dblarr(n_predict)
      IF (ARG_PRESENT(cooks_d) EQ TRUE) THEN $
        cooks_d_spc = dblarr(n_predict)
      IF (ARG_PRESENT(dffits) EQ TRUE) THEN $
        dffits_spc = dblarr(n_predict)
   END ELSE BEGIN
      l_result = fltarr(n_predict)
      ; 
      ; Input 
      x_cvt = float(transpose(x))
      IF (KEYWORD_SET(weights) EQ TRUE) THEN $
        weights_cvt = float(weights)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (float(confidence))(0)
      IF (KEYWORD_SET(y) EQ TRUE) THEN $
        y_cvt = float(y)
      ; 
      ; Output 
      IF (ARG_PRESENT(ci_scheffe) EQ TRUE) THEN $
        ci_scheffe_spc = fltarr(n_predict, 2)
      IF (ARG_PRESENT(Ci_ptw_pop_mean) EQ TRUE) THEN $
        Ci_pop_mean_spc = fltarr(n_predict, 2)
      IF (ARG_PRESENT(Ci_ptw_new_samp) EQ TRUE) THEN $
        Ci_new_samp_spc = fltarr(n_predict, 2)
      IF (ARG_PRESENT(leverage) EQ TRUE) THEN $
        leverage_spc = fltarr(n_predict)
      IF (ARG_PRESENT(residual) EQ TRUE) THEN $
        residual_spc = fltarr(n_predict)
      IF (ARG_PRESENT(std_residual) EQ TRUE) THEN $
        std_resid_spc = fltarr(n_predict)
      IF (ARG_PRESENT(del_residual) EQ TRUE) THEN $
        del_resid_spc = fltarr(n_predict)
      IF (ARG_PRESENT(cooks_d) EQ TRUE) THEN $
        cooks_d_spc = fltarr(n_predict)
      IF (ARG_PRESENT(dffits) EQ TRUE) THEN $
        dffits_spc = fltarr(n_predict)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_161, type, err_status, info_v, IMSL_N_ELEMENTS(info_v),  $
                    x_cvt, n_indep_vars,  n_predict, $
                    weights_cvt, $
                    confidence_cvt, $
                    y_cvt, $
                    ci_scheffe_spc, $
                    Ci_pop_mean_spc, $
                    Ci_new_samp_spc, $
                    leverage_spc, $
                    residual_spc, $
                    std_resid_spc, $
                    del_resid_spc, $
                    cooks_d_spc,$
                    dffits_spc, $
                    l_result
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(ci_scheffe) EQ TRUE) THEN $
        ci_scheffe = transpose(ci_scheffe_spc)
      IF (ARG_PRESENT(Ci_ptw_pop_mean) EQ TRUE) THEN $
        ci_ptw_pop_mean = transpose(ci_pop_mean_spc)
      IF (ARG_PRESENT(ci_ptw_new_samp) EQ TRUE) THEN $
        ci_ptw_new_samp = transpose(Ci_new_samp_spc)
      IF (ARG_PRESENT(cooks_d) EQ TRUE) THEN $
        cooks_d = cooks_d_spc
      IF (ARG_PRESENT(del_residual) EQ TRUE) THEN $
        del_residual = del_resid_spc
      IF (ARG_PRESENT(dffits) EQ TRUE) THEN $
        dffits = dffits_spc
      IF (ARG_PRESENT(leverage) EQ TRUE) THEN $
        leverage = leverage_spc
      IF (ARG_PRESENT(residual) EQ TRUE) THEN $
        residual = residual_spc
      IF (ARG_PRESENT(std_residual) EQ TRUE) THEN $
        std_residual = std_resid_spc
   ;
   ; Return.
   ;
   RETURN, l_result
END
   

