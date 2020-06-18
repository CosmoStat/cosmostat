; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_random_arma.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_random_arma, n, $               ;INPUT Scalar LONG
                     nparams, $               ;INPUT 1-D array: LONG
                     arg3, $                  ;INPUT 1-D array: floating point
                     arg4, $                  ;INPUT 1-D array: floating point
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     accept_reject= accept_reject, $ ;INPUT Scalar ON/OFF flag
                     const=const, $           ;INPUT Scalar floating point
                     var_noise=var_noise, $   ;INPUT Scalar floating point
                     input_noise=input_noise, $   ;INPUT 1-D array: floating point
                     ar_lags=ar_lags, $       ;INPUT 1-D array: LONG
                     ma_lags=ma_lags, $       ;INPUT 1-D array: LONG
                     w_init=w_init, $         ;INPUT 1-D array: floating point
                     output_noise=output_noise    ;OUTPUT 1-D array: floating point
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - N greater than 0.
   ;   - NPARAMS a 2-element array. (p=nparams(0), q = nparams(1))
   ;   - If (NPARAMS(*) must be greater than or equal to 0.
   ;   - If (p != 0) and (q == 0), arg3 is AR, arg4 is not allowed.
   ;   - If (p == 0) and (q != 0), arg3 is MA, arg4 is not allowed.
   ;   - If (p != 0) and (q != 0), arg3 is AR, arg4 is MA.
   ;   - If (p == 0) and (q == 0), arg3 and arg4 are not allowed.
   ;   - If present, AR is a 1D array of length p.
   ;   - If present, MA is a 1D array of length q.
   ;   - If present, AR_LAGS must be a 1D array of length p.
   ;   - If present, MA_LAGS must be a 1D array of length q.
   ;   - If present, W_INIT must be a 1D array of length max(ar_lags).
   ;   - VAR_NOISE and INPUT_NOISE cannot be used together.
   ;   - OUTPUT_NOISE and INPUT_NOISE cannot be used together.
   ;   - if INPUT_NOISE is supplied, it must be a 1D array of length n+max(ma_lags).
   ;
   nargs = n_params()
   IF (nargs LT 2) THEN MESSAGE,  "Incorrect number of arguments."
   IF (n LE 0) THEN MESSAGE,  "N_OBSERVATIONS must be greater than 0."
   n_cvt = IMSL_LONG(n(0))
   size_nparams = IMSL_SIZE(nparams)
   IF (size_nparams(0) NE 1) THEN BEGIN
      MESSAGE,  "NPARAMS must be a 1-D array."
   END
   IF (N_ELEMENTS(nparams) NE 2) THEN MESSAGE,   "NPARAMS is not the correct length."
   p_cvt = IMSL_LONG(nparams(0))
   q_cvt = IMSL_LONG(nparams(1))
   IF (MAX(nparams) LT 0) THEN MESSAGE, "All elements of NPARAM must be nonnegative."
   ; Case with 2 positional arguments
   IF (nargs EQ 2) THEN BEGIN
      IF ((p_cvt NE 0) OR (q_cvt NE 0)) THEN MESSAGE,  "Incorrect number of arguments."
      ar  =  0
      ma  = 0
      size_ar = IMSL_SIZE(ar)
      size_ma = IMSL_SIZE(ma)
   END
   ; Cases with 3 positional arguments
   IF (nargs EQ 3) THEN BEGIN
      IF ((p_cvt EQ 0) AND  (q_cvt EQ 0)) THEN $
        MESSAGE,  "Incorrect number of arguments."
      IF ((p_cvt NE 0) AND  (q_cvt NE 0)) THEN $
        MESSAGE,  "Incorrect number of arguments."
      IF ((p_cvt EQ 0) AND  (q_cvt NE 0)) THEN BEGIN
         ma  =  arg3
         ar  =  0
         size_ma  = IMSL_SIZE(ma)
         size_ar  = IMSL_SIZE(ar)
         IF (size_ma(0) NE 1) THEN BEGIN
            message, "MA must be a 1-D array."
         END
         IF (N_ELEMENTS(ma) NE q_cvt) THEN MESSAGE, "..MA is not the correct length."
      END
      IF ((p_cvt NE 0) AND  (q_cvt EQ 0)) THEN BEGIN
         ar  =  arg3
         ma  =  0
         size_ar  = IMSL_SIZE(ar)
         size_ma  = IMSL_SIZE(ma)
         IF (size_ar(0) NE 1) THEN BEGIN
            message, "AR must be a 1-D array."
         END
         IF (N_ELEMENTS(ar) NE p_cvt) THEN MESSAGE, "AR is not the correct length."
      END
   END
   ; Cases with 4 positional arguments
   IF (nargs EQ 4) THEN BEGIN
      IF ((p_cvt  <  q_cvt) EQ 0) THEN $
        MESSAGE,  "Incorrect number of arguments."
      IF ((p_cvt NE 0) AND  (q_cvt NE 0)) THEN BEGIN
         ar  =  arg3
         ma  =  arg4
         size_ar  = IMSL_SIZE(ar)
         IF (size_ar(0) NE 1) THEN BEGIN
            message, "AR must be a 1-D array."
         END
         IF (N_ELEMENTS(ar) NE p_cvt) THEN MESSAGE, "AR is not the correct length."
         size_ma  = IMSL_SIZE(ma)
         IF (size_ma(0) NE 1) THEN BEGIN
            message, "MA must be a 1-D array."
         END
         IF (N_ELEMENTS(ma) NE q_cvt) THEN MESSAGE, "MA is not the correct length."
      END
   END
   
   IF (KEYWORD_SET(ar_lags)) THEN BEGIN
      size_ar_lags = IMSL_SIZE(ar_lags)
      IF ((p_cvt GT 1) AND (size_ar_lags(0) NE 1)) THEN BEGIN
         MESSAGE,  "AR_LAGS must be a 1-D array."
      END
      IF (N_ELEMENTS(ar_lags) NE p_cvt) THEN MESSAGE,  "AR_LAGS is not the correct length."   
      max_ar_lags = IMSL_LONG(MAX(ar_lags))
   END ELSE max_ar_lags = p_cvt  
   IF (KEYWORD_SET(ma_lags)) THEN BEGIN
      size_ma_lags = IMSL_SIZE(ma_lags)
      IF ((q_cvt GT 1) AND (size_ma_lags(0) NE 1)) THEN BEGIN
         MESSAGE,  "MA_LAGS must be a 1-D array."
      END
      IF (N_ELEMENTS(ma_lags) NE q_cvt) THEN MESSAGE,   "MA_LAGS is not the correct length."
      max_ma_lags = IMSL_LONG(MAX(ma_lags))
   END ELSE max_ma_lags = q_cvt  
   IF (KEYWORD_SET(w_init)) THEN BEGIN
      size_w_init  = IMSL_SIZE(w_init)
      IF ((max_ar_lags GT 1) AND (size_w_init(0) NE 1)) THEN BEGIN
         message, "W_INIT must be a 1-D array."
      END
      IF (N_ELEMENTS(w_init) NE max_ar_lags) THEN MESSAGE, "W_INIT is not the correct length."
   END
   IF ((KEYWORD_SET(input_noise) + ARG_PRESENT(output_noise)) EQ 2) THEN $
     MESSAGE, "Keywords INPUT_NOISE and OUTPUT_NOISE cannot be used together."
   IF (KEYWORD_SET(input_noise)) THEN BEGIN
      size_in = IMSL_SIZE(input_noise)
      IF (size_in(0) NE 1) THEN BEGIN
         message, "INPUT_NOISE must be a 1-D array."
      END
      IF (N_ELEMENTS(input_noise) NE n_cvt) THEN $
        MESSAGE,  "INPUT_NOISE is not the correct length."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_ma(N_ELEMENTS(size_ma)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_ar(N_ELEMENTS(size_ar)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(ar_lags) EQ TRUE) THEN ar_lags_cvt = IMSL_LONG(ar_lags)
   IF (KEYWORD_SET(ma_lags) EQ TRUE) THEN ma_lags_cvt = IMSL_LONG(ma_lags)
   IF (KEYWORD_SET(accept_reject) EQ TRUE) THEN acc_rej_cvt = IMSL_1
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      ar_cvt = DOUBLE(ar)
      ma_cvt  =  DOUBLE(ma)
      IF (KEYWORD_SET(const) EQ TRUE) THEN const_cvt = DOUBLE(const(0)) ELSE const_cvt = DOUBLE(0.0)
      IF (KEYWORD_SET(var_noise) EQ TRUE) THEN var_noise_cvt = DOUBLE(var_noise(0))
      IF (KEYWORD_SET(input_noise) EQ TRUE) THEN input_noise_cvt = DOUBLE(input_noise)
      IF (KEYWORD_SET(w_init) EQ TRUE) THEN w_init_cvt = DOUBLE(w_init)
      ; Output
      IF (ARG_PRESENT(output_noise)) THEN on_spc = DBLARR(n_cvt+max_ma_lags)      
      result = DBLARR(n_cvt)
   END ELSE BEGIN
      ; Input
      ar_cvt = FLOAT(ar)
      ma_cvt  =  FLOAT(ma)
      IF (KEYWORD_SET(const) EQ TRUE) THEN const_cvt = FLOAT(const(0)) ELSE const_cvt = FLOAT(0.0)
      IF (KEYWORD_SET(var_noise) EQ TRUE) THEN var_noise_cvt = FLOAT(var_noise(0))
      IF (KEYWORD_SET(input_noise) EQ TRUE) THEN input_noise_cvt = FLOAT(input_noise)
      IF (KEYWORD_SET(w_init) EQ TRUE) THEN w_init_cvt = FLOAT(w_init)
      ; Output
      IF (ARG_PRESENT(output_noise)) THEN on_spc = FLTARR(n_cvt+max_ma_lags)      
      result = FLTARR(n_cvt)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_268,  type,  err_status,  $
                                      n_cvt,   $
                                      p_cvt,   $
                                      ar_cvt,  $
                                      q_cvt,   $
                                      ma_cvt,  $
                                      const_cvt,  $
                                      var_noise_cvt,     $
                                      input_noise_cvt,   $
                                      w_init_cvt,        $
                                      ar_lags_cvt,  $
                                      ma_lags_cvt,  $
                                      acc_rej_cvt,  $
                                      on_spc,     $
                                      result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(output_noise) EQ TRUE) THEN output_noise = on_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
