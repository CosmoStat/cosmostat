; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_signtest.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_signtest, x, $                     ;INPUT 1-D array: floating point
                   double=double, $         ;INPUT Scalar ON/OFF flag
                   percentage=percentage, $ ;INPUT Scalar floating point
                   percentile=percentile, $ ;INPUT Scalar floating point
                   n_pos_dev=n_pos_dev, $   ;OUTPUT Scalar LONG
                   n_zero_dev=n_zero_dev ;OUTPUT Scalar LONG
   
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 1) THEN BEGIN
      size_x = IMSL_SIZE(x)
      IF (size_x(0) NE 1) THEN BEGIN
         message, "X must be a 1-D array."
      END
   END ELSE message, "Incorrect number of arguments."
      
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
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(n_pos_dev) EQ TRUE) THEN n_pos_dev_spc = IMSL_0
   IF (ARG_PRESENT(n_zero_dev) EQ TRUE) THEN n_zero_dev_spc = IMSL_0
   ; 
   ; Input LONG argument(s)
   nobs = size_x(N_ELEMENTS(size_x)-1)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result
      ;
      result = double(0.0)
      ; 
      ; Input 
      x_cvt = double(x)
      IF (KEYWORD_SET(percentile) EQ TRUE) THEN $
        percentile_cvt = (double(percentile))(0)
      IF (KEYWORD_SET(percentage) EQ TRUE) THEN $
        percentage_cvt = (double(percentage))(0)
   END ELSE BEGIN
      ; Result
      ;
      result = float(0.0)
      ; 
      ; Input 
      x_cvt = float(x)
      IF (KEYWORD_SET(percentile) EQ TRUE) THEN $
        percentile_cvt = (float(percentile))(0)
      IF (KEYWORD_SET(percentage) EQ TRUE) THEN $
        percentage_cvt = (float(percentage))(0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_186, type, err_status, x_cvt, nobs, percentage_cvt, percentile_cvt, $
                         n_pos_dev_spc, n_zero_dev_spc, result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(n_pos_dev) EQ TRUE) THEN n_pos_dev=n_pos_dev_spc
   IF (ARG_PRESENT(n_zero_dev) EQ TRUE) THEN n_zero_dev=n_zero_dev_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
