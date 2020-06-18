; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_rand_from_data.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_rand_from_data, n, $            ;INPUT Scalar LONG
                     x, $                     ;INPUT 2-D array: floating point
                     nn, $                    ;INPUT Scalar LONG
                     double=double            ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - N greater than 0.
   ;   - x is size (nsamp, ndim) if 2D, or size (nsamp) if 1D.
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN MESSAGE,  "Incorrect number of arguments."
   IF (n LE 0) THEN MESSAGE,  "N must be greater than 0."
   n_cvt = IMSL_LONG(n(0))

   size_x = IMSL_SIZE(x)
   IF ((size_x(0) NE 1) AND (size_x(0) NE 2)) THEN BEGIN
      message, "X must be a 1-D or 2-D array."
   END
   IF (size_x(0) EQ 1) then BEGIN
      nsamp = size_x(1)
      ndim = 1
   END ELSE BEGIN  ;X is 2-D
      nsamp = size_x(1)
      ndim = size_x(2)
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
   ; Input LONG
   nsamp_cvt = IMSL_LONG(nsamp)
   ndim_cvt = IMSL_LONG(ndim)
   nn_cvt = IMSL_LONG(nn(0))
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
     x_cvt = DOUBLE(TRANSPOSE(x))
     result = DBLARR(ndim_cvt, n_cvt)
   END ELSE BEGIN
     x_cvt = FLOAT(TRANSPOSE(x))
     result = FLTARR(ndim_cvt, n_cvt)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_304,  type,  err_status,  $
                                      n_cvt,   $
                                      nn_cvt,   $
                                      nsamp_cvt,  $
                                      ndim_cvt,   $
                                      x_cvt,  $
                                      result
   ;
   ; Now copy over all output keywords results.
   ;
   ; Return.
   ;
   if (n_cvt GT 1) THEN result = transpose(result)
   RETURN, result
END

                   
                   
                   

  
      

  
