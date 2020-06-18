; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_kolmogorov2.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_kolmogorov2, x, $                  ;INPUT 1-D array: floating point
                      y, $                       ;INPUT 1-D array: floating point
                     double=double, $            ;INPUT Scalar ON/OFF flag
                     differences=differences, $  ;OUTPUT 1-D array: floating point
                     nmissingx=nmissingx, $      ;OUTPUT Scalar LONG
                     nmissingy=nmissingy         ;OUTPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 1D array   
   ;   Y must be a 1D array   
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN
      message, "X must be a 1-D array."
   END
   nobsx = IMSL_N_ELEMENTS(x)
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   nobsy = IMSL_N_ELEMENTS(y)
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(nmissingx) EQ TRUE) THEN nmissingx_spc = IMSL_LONG(0)
   IF (ARG_PRESENT(nmissingy) EQ TRUE) THEN nmissingy_spc = IMSL_LONG(0)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt = double(x)
      y_cvt = double(y)
      ; Output
      IF (ARG_PRESENT(differences) EQ TRUE) THEN diff_spc = DBLARR(3)
      result = dblarr(3)
   END ELSE BEGIN
      ; Input
      x_cvt = float(x)
      y_cvt = float(y)
      ; Output
      IF (ARG_PRESENT(differences) EQ TRUE) THEN diff_spc = fltarr(3)
      result = fltarr(3)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_265,type,err_status, x_cvt,  y_cvt, nobsx, nobsy, $
                           diff_spc, $
                           nmissingx_spc, $
                           nmissingy_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(nmissingx) EQ TRUE) THEN nmissingx = nmissingx_spc
   IF (ARG_PRESENT(nmissingy) EQ TRUE) THEN nmissingy = nmissingy_spc
   IF (ARG_PRESENT(differences) EQ TRUE) THEN differences = diff_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
