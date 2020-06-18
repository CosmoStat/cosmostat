; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_lack_of_fit.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_lack_of_fit, nobs, $                   ;INPUT 1-D Scalar LONG
                     cf, $                           ;INPUT 1-D array: floating point
                     npfree, $                       ;INPUT 1-D Scalar LONG
                     lagmin=lagmin, $                ;INPUT Scalar LONG
                     double=DOUBLE                   ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   CF must be a 1D array
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN message, "Incorrect number of arguments."
   size_cf = IMSL_SIZE(cf)
   IF (size_cf(0) NE 1) THEN BEGIN
      message, "CF must be a 1-D array."
   END
   lagmax  =  IMSL_N_ELEMENTS(cf)-1
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
   ; Input LONG keyword(s) and arguments
   nobs_cvt = IMSL_LONG(nobs(0))
   npfree_cvt = IMSL_LONG(npfree(0))
   IF (KEYWORD_SET(lagmin) EQ TRUE) THEN lagmin_cvt = IMSL_LONG(lagmin) ELSE lagmin_cvt = IMSL_1
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      cf_cvt = double(cf)
      result = dblarr(2)
   END ELSE BEGIN
      ; Input
      cf_cvt = float(cf)
      result = fltarr(2)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_264,type,err_status, cf_cvt, nobs_cvt, lagmax, npfree_cvt, $
                           lagmin_cvt, $
                           result                       
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
