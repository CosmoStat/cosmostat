; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_nctrends.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_nctrends, x, $                  ;INPUT 1-D array: floating point
                     double=double, $         ;INPUT Scalar ON/OFF flag
                     fuzz=fuzz, $             ;INPUT 1-D Scalar floating point
                     nstat=nstat, $           ;OUTPUT 1-D array: LONG
                     nmissing=nmissing        ;OUTPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 1D array   
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN
      message, "X must be a 1-D array."
   END
   nobs = IMSL_N_ELEMENTS(x)
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
   ; Output LONG keyword(s)
   IF (ARG_PRESENT(nstat) EQ TRUE) THEN nstat_spc = IMSL_LONARR(6)
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing_spc = IMSL_LONG(0)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt = double(x)
      IF (KEYWORD_SET(fuzz) EQ TRUE) THEN fuzz_cvt = DOUBLE(fuzz) ELSE fuzz_cvt = DOUBLE(0.0)
      result = dblarr(3)
   END ELSE BEGIN
      ; Input
      x_cvt = float(x)
      IF (KEYWORD_SET(fuzz) EQ TRUE) THEN fuzz_cvt = FLOAT(fuzz) ELSE fuzz_cvt = FLOAT(0.0)
      result = fltarr(3)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_258,type,err_status, x_cvt,  nobs, $
                           fuzz_cvt, $
                           nstat_spc, $
                           nmissing_spc, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(nmissing) EQ TRUE) THEN nmissing = nmissing_spc
   IF (ARG_PRESENT(nstat) EQ TRUE) THEN nstat = nstat_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
