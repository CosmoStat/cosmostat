; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_garch.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_garch, p, $                  ;INPUT Scalar LONG
                q, $                       ;INPUT Scalar LONG
                y, $                       ;INPUT 1-D array: floating point
                xguess, $                  ;INPUT 1-D array: floating point
                double=double, $           ;INPUT Scalar ON/OFF flag
                max_sigma=max_sigma, $     ;INPUT 1-D Scalar floating point
                var=var, $                 ;OUTPUT 2-D array: floating point
                log_likelihood=log_likelihood, $ ;OUTPUT Scalar floating point
                aic=aic                    ;OUTPUT Scalar floating point

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   Y must be a 1D array   
   ;   XGUESS must be a 1D array of length [p+q+1]
   ;   P+Q+1 must be positive.
   nargs = n_params()
   IF (nargs NE 4) THEN message, "Incorrect number of arguments."
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   m = IMSL_N_ELEMENTS(y)
   size_xguess = IMSL_SIZE(xguess)
   IF (size_xguess(0) NE 1) THEN BEGIN
      message, "XGUESS must be a 1-D array."
   END
   IF (N_ELEMENTS(xguess) NE (p+q+1)) THEN MESSAGE, "XGUESS is not the correct size."
   IF ((P+Q+1) LT 1) THEN MESSAGE, "P+Q+1 must be positive."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_xguess(N_ELEMENTS(size_xguess)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s) and arguments
   p_cvt = IMSL_LONG(p)
   q_cvt = IMSL_LONG(q)
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      y_cvt = double(y)
      xguess_cvt = double(xguess)
      IF (KEYWORD_SET(max_sigma) EQ TRUE) THEN max_sigma_cvt = DOUBLE(max_sigma) ELSE max_sigma_cvt = DOUBLE(10.0)
      var_spc = DBLARR(p+1+q, p+q+1)
      aic_spc = DOUBLE(0.0)
      llh_spc = DOUBLE(0.0)
      result = dblarr(p+q+1)
   END ELSE BEGIN
      ; Input
      y_cvt = float(y)
      xguess_cvt = float(xguess)
      IF (KEYWORD_SET(max_sigma) EQ TRUE) THEN max_sigma_cvt = float(max_sigma) ELSE max_sigma_cvt = float(10.0)
      var_spc = fltarr(p+1+q, p+q+1)
      aic_spc = float(0.0)
      llh_spc = float(0.0)
      result = fltarr(p+q+1)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_266,type,err_status, p_cvt, q_cvt, m, y_cvt, xguess_cvt, $
                           max_sigma_cvt, $
                           var_spc, $
                           llh_spc, $
                           aic_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(var) EQ TRUE) THEN var = TRANSPOSE(var_spc)
   IF (ARG_PRESENT(aic) EQ TRUE) THEN aic = aic_spc
   IF (ARG_PRESENT(log_likelihood) EQ TRUE) THEN log_likelihood = llh_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
