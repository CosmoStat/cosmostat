; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_fcnlsq.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_fcnlsq, F, $                         ;INPUT Scalar STRING
                 nbasis, $                      ;INPUT Scalar LONG
                 xdata, $                       ;INPUT 1-D array: floating point
                 fdata, $                       ;INPUT 1-D array: floating point
                 weights=weights, $             ;INPUT 1-D array: floating point
                 double=double, $               ;INPUT Scalar ON/OFF flag
                 intercept=intercept, $         ;OUTPUT Scalar floating point
                 sse=sse                        ;OUTPUT Scalar floating point
                 

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; - F must be a scalar string.
   ; - nbasis is converted to TYP_MEMINT, and must be positive.
   ; - xdata must be a 1-D array (with ndata elements)
   ; - fdata must be a 1-D array (with ndata elements)
   ; - WEIGHTS must be a 1-D array with ndata elements.
   ;
   nargs = n_params()
   IF (nargs NE 4) THEN $
     message, 'Incorrect number of arguments.'
   size_f = IMSL_SIZE(f)
   IF ((N_ELEMENTS(f) NE 1) OR (size_f(N_ELEMENTS(size_f)-2) NE 7)) THEN $
     message, 'F must be a scalar string.'
   nbasis_cvt = IMSL_LONG(nbasis(0))
   IF (nbasis_cvt LT 1) THEN message, 'NBASIS must be positive.'
   size_xdata = IMSL_SIZE(xdata)
   IF (size_xdata(0) NE 1) THEN message, 'XDATA must be a 1-D array.' $
     ELSE ndata = IMSL_N_ELEMENTS(xdata)
   size_fdata = IMSL_SIZE(fdata)
   IF ((size_fdata(0) NE 1) OR (N_ELEMENTS(fdata) NE ndata)) THEN $
        message, 'FDATA must have the same dimensions as XDATA.'
   IF (KEYWORD_SET(weights)) THEN BEGIN 
      size_weights = IMSL_SIZE(weights)
      IF ((size_weights(0) NE 1) OR (N_ELEMENTS(weights) NE ndata)) THEN $
        message, 'WEIGHTS must have the same dimensions as XDATA.'
   END
   ;
   ;Select the precision to use.
   type = TYP_FLOAT
   IF (size_xdata(N_ELEMENTS(size_xdata)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_fdata(N_ELEMENTS(size_fdata)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE ELSE type = TYP_FLOAT
   ;
   ; Setup the parameters for the call to the system function.
   ;----------------------------------------------------------
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(nbasis_cvt)
      ;
      ; Input arguments and keyword(s)
      ;
      xdata_cvt = double(xdata)
      fdata_cvt = double(fdata)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = double(weights)
      ; Output keyword(s)
      ;
      IF (ARG_PRESENT(intercept)) THEN intercept_spc = double(0.0)
      IF (ARG_PRESENT(sse)) THEN sse_spc = double(0.0)
   END ELSE BEGIN
      result = fltarr(nbasis_cvt)
      ;
      ; Input arguments and keyword(s)
      ;
      xdata_cvt = float(xdata)
      fdata_cvt = float(fdata)
      IF (KEYWORD_SET(weights)) THEN weights_cvt = float(weights)
      ; Output keyword(s)
      ;
      IF (ARG_PRESENT(intercept)) THEN intercept_spc = float(0.0)
      IF (ARG_PRESENT(sse)) THEN sse_spc = float(0.0)
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_139, type, err_status, $
                 F, $
                 nbasis_cvt, $
                 xdata_cvt, $
                 fdata_cvt, $
                 ndata, $
                 weights_cvt, $
                 intercept_spc, $
                 sse_spc, $
                 result

   IF ARG_PRESENT(intercept) THEN intercept = intercept_spc
   IF ARG_PRESENT(sse)       THEN       sse = sse_spc

   ; return
   RETURN, result
END
