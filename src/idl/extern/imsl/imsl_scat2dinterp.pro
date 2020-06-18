; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_scat2dinterp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_scat2dinterp, xydata, $           ;INPUT 2-D array: floating point
               fdata, $                 ;INPUT -D array: floating point
               xout, $                   ;INPUT -D array: floating point
               yout, $                   ;INPUT -D array: floating point
               double=double            ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking.  
   ;
   ; xydata
   ;   - Must be a 2-D array with first dimension 2.
   ;     Set (n = second dimension)
   ; fdata
   ;   - Must be a 1-D array of length (n).
   ; xout
   ;   - Must be a 1-D array. Set nxout = (its length)
   ; yout
   ;   - Must be a 1-D array. Set nyout = (its length)
   ;
   nargs = n_params()
   IF (nargs NE 4) THEN $
         message, "Incorrect number of arguments."
   ;
   size_xydata = IMSL_SIZE(xydata)
   size_fdata = IMSL_SIZE(fdata)
   size_xout = IMSL_SIZE(xout)
   size_yout = IMSL_SIZE(yout)

   IF (size_xydata(0) NE 2) THEN BEGIN 
      message, "XYDATA must be a 2-D array."
   END
   m = IMSL_LONG(size_xydata(1))
   n = IMSL_LONG(size_xydata(2))
   IF ((size_fdata(0) NE 1) OR (N_ELEMENTS(fdata) NE n)) THEN BEGIN 
      message, 'FDATA must be a 1-D arrayof length '+strtrim(n,1)+'.'
   END
   IF (size_xout(0) NE 1) THEN BEGIN 
      message, "XOUT must be a 1-D array."
   END
   nxout = IMSL_LONG(size_xout(1))
   IF (size_yout(0) NE 1) THEN BEGIN 
      message, "YOUT must be a 1-D array."
   END
   nyout = IMSL_LONG(size_yout(1))
   ;
   ;ERROR CHECKING COMPLETE.
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_xydata(N_ELEMENTS(size_xydata)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_fdata(N_ELEMENTS(size_fdata)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_xout(N_ELEMENTS(size_xout)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_yout(N_ELEMENTS(size_yout)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
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
      result = dblarr(nyout, nxout)
      xydata_cvt = double(xydata)
      fdata_cvt = double(fdata)
      xout_cvt = double(xout)
      yout_cvt = double(yout)
   END ELSE BEGIN
      ; Result vector.
      result = fltarr(nyout, nxout)
      xydata_cvt = float(xydata)
      fdata_cvt = float(fdata)
      xout_cvt = float(xout)
      yout_cvt = float(yout)
   END
      ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_185,   type, err_status, $
               xydata_cvt, $
               fdata_cvt, $
               xout_cvt, $
               yout_cvt, $
               n, $
               nxout, $
               nyout, $
               result
   ;
   ; Return.
   ;
   RETURN, transpose(result)
END
   

