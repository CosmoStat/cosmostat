; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_partial_ac.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_partial_ac, cf, $                  ;INPUT 1-D array: floating point
                     double=DOUBLE               ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   CF must be a 1D array, length = (lagmax+1)
   ;
   nargs   =   N_PARAMS()
   IF (nargs NE 1) THEN message, "Incorrect number of arguments."
   size_cf = IMSL_SIZE(cf)
   IF (size_cf(0) NE 1) THEN BEGIN
      message, "CF must be a 1-D array."
   END
   lagmax = IMSL_N_ELEMENTS(cf)-1
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_cf(N_ELEMENTS(size_cf)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      cf_cvt = double(cf)
      result = dblarr(lagmax)
   END ELSE BEGIN
      ; Input
      cf_cvt = float(cf)
      result = fltarr(lagmax)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_262,type,err_status, cf_cvt,  lagmax, $
                           result
   ;
   ; Now copy over all output keywords results.
   ;
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
