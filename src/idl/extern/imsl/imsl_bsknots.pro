; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_bsknots.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_bsknots, xdata, $            ;INPUT 1-D array: floating point
               order=order, $           ;INPUT Scalar LONG
               optimum=optimum, $       ;INPUT Scalar ON/OFF flag
               itmax=itmax, $           ;INPUT Scalar LONG
               double=double            ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking.  
   ;  xdata must be a 1-D array.
   ;  If ORDER is present, it must be greater than 1.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN $
         message, "Incorrect number of arguments."
   ;
   size_xdata = IMSL_SIZE(xdata)
   IF (size_xdata(0) NE 1) THEN BEGIN 
      message, "XDATA must be a 1-D array."
   END
   nxdata = IMSL_LONG(size_xdata(1))
   IF (KEYWORD_SET(order)) THEN $
     IF (IMSL_LONG(order(0)) Le IMSL_1) THEN message, "ORDER must be greater than 1."
   ;
   ;ERROR CHECKING COMPLETE.
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_xdata(N_ELEMENTS(size_xdata)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(order) EQ TRUE) THEN $
     order_cvt = (IMSL_LONG(order))(0) ELSE order_cvt = IMSL_4
   IF (KEYWORD_SET(optimum) EQ TRUE) THEN $
     optimum_cvt = IMSL_1
   IF (KEYWORD_SET(itmax) EQ TRUE) THEN $
     itmax_cvt = (IMSL_LONG(itmax))(0)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      result = dblarr(nxdata+order_cvt)
      xdata_cvt = double(xdata)
   END ELSE BEGIN
      ; Result vector.
      result = fltarr(nxdata+order_cvt)
      xdata_cvt = float(xdata)
   END
      ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_114,   type, err_status, $
               xdata_cvt, $
               nxdata, $
               order_cvt, $
               optimum_cvt, $
               itmax_cvt, $
               result
   ;
   ; Return.
   ;
   RETURN, result
END
   

