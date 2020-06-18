; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_faure_init.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_faure_init, ndim, $                  ;INPUT Scalar LONG
                  base=base, $                     ;INPUT Scalar LONG
                  skip=skip                        ;INPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input arguments are 1-D or 2-D arrays.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN $
     MESSAGE, "Incorrect number of arguments."
   ndim_cvt = IMSL_LONG(ndim(0))
   IF (ndim_cvt LT 1) THEN MESSAGE, "NDIM must be positive."
   IF (KEYWORD_SET(base) EQ TRUE) THEN base_cvt = IMSL_LONG(base(0)) 
   IF (KEYWORD_SET(skip) EQ TRUE) THEN skip_cvt = IMSL_LONG(skip(0))

   ; Get value for maxDigits
   err_flag = IMSL_0
   maxDigits = IMSL_0
   MATHSTAT_313, ndim_cvt, base_cvt, maxDigits, err_flag
   IF (err_flag NE 0) THEN MESSAGE, "BASE must be a prime at least as large as the dimension."
   ;
   ; Get space for result.
   nskip_spc = IMSL_0
   y_spc = IMSL_LONARR(ndim_cvt*maxDigits)
   maxdigits_spc = IMSL_0
   base_spc = IMSL_0
   digitsN_spc = IMSL_LONARR(maxDigits)
   c_spc = IMSL_LONARR(maxDigits*maxDigits)
   power_spc = IMSL_LONARR(ndim_cvt*2*maxDigits)
   scale_spc = DOUBLE(0.0)
   ;
   ; Call the system function.
   ;
   err_status = 0L
    MATHSTAT_314,  err_status, $
                ndim_cvt, $
                base_cvt, $
                skip_cvt, $
                nskip_spc, $
                y_spc, $
                maxdigits_spc, $
                base_spc, $
                digitsN_spc, $
                c_spc, $
                power_spc, $
                scale_spc
   ;
   ; Now copy over all output keywords results.
   ;
    result = {NSKIP:nskip_spc, $
               Y:y_spc, DIM:ndim_cvt, $
               MAXDIGITS:maxDigits_spc, $
               BASE:base_spc, $
               DIGITSN:digitsN_spc, $
               C:c_spc, $
               POWER:power_spc, $
               SCALE:scale_spc}
   ;
   ; Return.
   ;
   RETURN, result
END
   

