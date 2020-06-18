; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_ktrends.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_ktrends, n, $                   ;INPUT 1-D array: LONG
                   y, $                       ;INPUT 1-D array: floating point
                     double=DOUBLE            ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   N must be a 1D array
   ;   Y must be a 1D array
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_n = IMSL_SIZE(n)
   IF (size_n(0) NE 1) THEN BEGIN
      message, "N must be a 1-D array."
   END
   ngroups = IMSL_N_ELEMENTS(n)
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   IF (N_ELEMENTS(y) NE TOTAL(n)) THEN MESSAGE, "Y is not the correct length."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   ; Output LONG keyword(s)/Args
   n_cvt = IMSL_LONG(n)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      y_cvt = double(y)
      result = dblarr(17)
   END ELSE BEGIN
      ; Input
      y_cvt = float(y)
      result = fltarr(17)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_259,type,err_status, n_cvt, y_cvt,  ngroups, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
