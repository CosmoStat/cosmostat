; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_regressors.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_regressors, x, $                              ;INPUT 2-D array: floating point
                     n_class, $                        ;INPUT Scalar LONG
                     n_continuous, $                   ;INPUT Scalar LONG
                     class_columns=class_columns, $    ;INPUT 1-D array: LONG
                     double=double, $                  ;INPUT Scalar ON/OFF flag
                     dummy_method=dummy_method, $      ;INPUT Scalar LONG
                     indices_effects=indices_effects, $;INPUT 1-D array: LONG
                     order=order, $                    ;INPUT Scalar LONG
                     var_effects=var_effects           ;INPUT 1-D array: LONG

@imsl_init.pro   
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; The following checks are performed.
   ;          x : must be an array.
   ;             -If it is a 1-D array, then 
   ;                nobs = first dimension
   ;              and (n_class+c_continous) must be 1.
   ;             -If it is a 2-D array, then 
   ;                nobs = first dimension
   ;              and second dimension must 
   ;              be equal to the sum of the second and 
   ;              third arguments.
   ; 
   ;   If CLASS_COLUMNS is present, should be a TYP_MEMINT array of 
   ;           length n_class (the third argument).
   ;              
   ; 
   ;   VAR_EFFECTS and INDICES_EFFECTS must be used together.
   ;             -VAR_EFFECTS must specify a 1-D array, and assign
   ;              n_effects = the # of elements of this array.
   ;             -INDICES_EFFECTS should be a 1-D array whose length
   ;              is equal to the sum of the elements of the 
   ;              VAR_EFFECTS array.
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN $
         message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   n_class_cvt = (IMSL_LONG(n_class))(0)
   n_continuous_cvt = (IMSL_LONG(n_continuous))(0)
   IF ((size_x(0) NE 1) AND (size_x(0) NE 2)) THEN BEGIN
      message, "X must be a 1-D or 2-D array."
   END
   IF (size_x(0) EQ 1) then BEGIN
      nobs = size_x(1)
      IF ((n_class_cvt + n_continuous_cvt) NE 1) THEN $
        message, "When X is a 1-D array, the sum of N_CLASS and N_CONTINUOUS must be equal to one"
   END ELSE BEGIN  ;X is 2-D
      nobs = size_x(1)
      IF ((n_class_cvt + n_continuous_cvt) NE size_x(2)) THEN $
        message, "The second dimension of X must equal the sum of N_CLASS and N_CONTINUOUS"
   END
   IF KEYWORD_SET(class_columns) THEN BEGIN
      IF (N_ELEMENTS(class_columns) NE n_class_cvt) THEN $
        message, "CLASS_COLUMNS must have N_CLASS elements."
   END
   IF ((KEYWORD_SET(var_effects)) OR  (KEYWORD_SET(indices_effects))) THEN BEGIN 
      i_tmp = KEYWORD_SET(var_effects) + KEYWORD_SET(indices_effects)
      IF (i_tmp EQ 1) THEN message, "VAR_EFFECTS and INDICES_EFFECTS must be used together."
      size_tmp = IMSL_SIZE(var_effects)
      IF (size_tmp(0) NE 1) THEN message, "var_EFFECTS must be a 1-D array."
      n_effects = size_tmp(1)
      var_effects_cvt = IMSL_LONG(var_effects)
      size_tmp = IMSL_SIZE(indices_effects)
      IF (size_tmp(0) NE 1) THEN message, "INDICES_EFFECTS must be a 1-D array."
      IF (size_tmp(1) NE total(var_effects_cvt)) THEN $
        message, "The number of elements in INDICES_EFFECTS must equal the sum " + $
        "of the elements OF VAR_EFFECTS"
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
   ; Input LONG keyword(s)
   ;
   IF (KEYWORD_SET(class_columns) EQ TRUE) THEN $
     class_cols_cvt = IMSL_LONG(class_columns)
   IF (KEYWORD_SET(dummy_method) EQ TRUE) THEN $
     dummy_cvt = IMSL_LONG(dummy_method)
   IF (KEYWORD_SET(indices_effects) EQ TRUE) THEN $
     indices_eff_cvt = IMSL_LONG(indices_effects)
   IF (KEYWORD_SET(order) EQ TRUE) THEN $
     order_cvt = IMSL_LONG(order)
   IF (KEYWORD_SET(var_effects) EQ TRUE) THEN $
     var_eff_cvt = IMSL_LONG(var_effects)
   ;
   ; Result: Scalar LONG
   result = IMSL_0
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input 
      x_cvt = double(transpose(x))
   END ELSE BEGIN
      ; 
      ; Input 
      x_cvt = float(transpose(x))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_184, type, err_status, x_cvt, $ 
                     n_class_cvt, $
                     n_continuous_cvt, $
                     nobs, $
                     n_effects, $
                     class_cols_cvt, $
                     dummy_cvt, $
                     indices_eff_cvt, $
                     order_cvt, $
                     var_eff_cvt, $
                     result
   ;
   ; Return.
   ;
   RETURN, transpose(result)
END
