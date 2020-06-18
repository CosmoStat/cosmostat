; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_difference.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_difference, z, $                       ;INPUT 1-D array: floating point
                     periods, $                      ;INPUT 1-D input array: LONG
                     double=double, $                ;INPUT Scalar ON/OFF flag
                     Exclude_first=exclude_first, $  ;INPUT Scalar ON/OFF flag
                     First_to_nan=first_to_nan, $    ;INPUT Scalar ON/OFF flag
                     Orders=orders, $                ;INPUT 1-D input array: LONG
                     Num_lost=num_lost               ;OUTPUT Scalar LONG output

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   check first arg is a simple array
   ;   check second arg is a simple array
   ;   if ORDERS is present, check that it specified a 1-D
   ;    array the same length as PERIODS.
   ;
   nargs = n_params()
   IF (nargs EQ 2) THEN BEGIN
      size_z = IMSL_SIZE(z)
      IF (size_z(0) NE 1) THEN BEGIN
         message, "Z must be a 1-D array."
      END
      size_periods = IMSL_SIZE(periods)
      IF (size_periods(0) NE 1) THEN BEGIN
         message, "PERIODS must be a 1-D array."
      END
      n_differences = size_periods(N_ELEMENTS(size_periods)-1)
   END ELSE message, "Incorrect number of arguments."

   IF (KEYWORD_SET(orders) EQ TRUE) THEN BEGIN
      size_orders = IMSL_SIZE(orders)
      IF (size_orders(0) NE 1) THEN  $
         message, "ORDERS must be a 1-D array."
      IF (size_orders(N_ELEMENTS(size_orders)-1) NE n_differences) THEN $
         message, "ORDERS must specify an array the that is the same length as PERIODS."
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_z(N_ELEMENTS(size_z)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_periods(N_ELEMENTS(size_periods)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(Orders) EQ TRUE) THEN orders_cvt = IMSL_LONG(orders)
   IF (KEYWORD_SET(Exclude_first) EQ TRUE) THEN Excl_first_cvt = IMSL_1
   IF (KEYWORD_SET(First_to_nan) EQ TRUE) THEN  First_nan_cvt = IMSL_1
   ;
   ; Output LONG keyword(s)
   ; We are always going to ask for num_lost, but only return if it was specified.
   num_lost_spc = IMSL_LONG(0)
   ;
   ; Input LONG argument(s)
   nobs = size_z(N_ELEMENTS(size_z)-1)
   periods_cvt = IMSL_LONG(periods)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      IF (N_ELEMENTS(z) NE 0) THEN z_cvt = double(z)
      ; Output
      result = dblarr(nobs)
   END ELSE BEGIN
      ; Input
      IF (N_ELEMENTS(z) NE 0) THEN z_cvt = float(z)
      ; Output
      result = fltarr(nobs)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_132,type,err_status, z_cvt, periods_cvt, nobs, n_differences, $
                           Excl_first_cvt, $
                           First_nan_cvt, $
                           Num_lost_spc, $
                           Orders_cvt, $
                           result
                           
   ;
   ; IF EXCLUDE_FIRST was set, then there is the possibility that we
   ; want to truncate the returned array by NUM_LOST elements.
   ;
   IF (KEYWORD_SET(Exclude_first) EQ TRUE) THEN result = result(0:((nobs-num_lost_spc-1) > 0))
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(Num_lost) EQ TRUE) THEN num_lost = num_lost_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
