; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_spvalue.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_spvalue, data1, $         ;INPUT Scalar or 1-D array: floating point 
                    data2, $         ;INPUT Spline structure,
                                     ;      Scalar or 1-D array: floating point 
                    data3, $         ;INPUT Spline structure,
                    xderiv=xderiv, $ ;INPUT Scalar LONG
                    yderiv=yderiv    ;INPUT Scalar LONG

@imsl_init.pro   
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; Error checking.
   ;  (nargs EQ 2)
   ;    - We are evaluating a 1-D spline.
   ;    - DATA1 must be either a scalar or a 1-D array of length (nx).
   ;    - DATA2 must be a spline stucture.
   ;
   ;  (nargs EQ 3)
   ;    - We are evaluating a 2-D spline.
   ;    - DATA1 must be either a scalar or a 1-D array of length (nx).
   ;    - DATA2 must be either a scalar or a 1-D array of length (ny).
   ;    - DATA3 must be a spline stucture.
   ;
   nargs = n_params()
   IF ((nargs NE 2) AND (nargs NE 3)) THEN message, 'Incorrect number of arguments.'
   size_data1 = IMSL_SIZE(data1)
   IF (size_data1(0) GT 1) THEN $
     message, 'XDATA must be a either a scalar or a 1-D array.' $
     ELSE nx = IMSL_N_ELEMENTS(data1)
   
   IF (nargs EQ 2) THEN BEGIN
      data2_cvt = imsl_chk_spline_struct(data2, type, struct_is_ppoly)
      domain_dim = IMSL_1
   END ELSE BEGIN
      ; (nargs eq 3)
      size_data2 = IMSL_SIZE(data2)
      IF (size_data2(0) GT 1) THEN $
        message, 'YDATA must be a either a scalar or a 1-D array.' $
        ELSE ny = IMSL_N_ELEMENTS(data2)
      data3_cvt = imsl_chk_spline_struct(data3, type, struct_is_ppoly)
      domain_dim = IMSL_2
   END

   ; Error checking complete.
   ;
   ; Decide on what precision to use.
   ; The precision is dependent upon the precision of the data in the
   ; spline structure.  The procedure chk_spline_struct() returns the precision.
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(xderiv)) THEN xderiv_cvt = IMSL_LONG(xderiv(0)) ELSE xderiv_cvt = IMSL_0
   IF (KEYWORD_SET(yderiv)) THEN yderiv_cvt = IMSL_LONG(yderiv(0)) ELSE yderiv_cvt = IMSL_0
   ;
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      IF (nargs EQ 2) THEN BEGIN
         data1_cvt = double(data1)
         result = dblarr(nx)
      END ELSE BEGIN
         data1_cvt = double(data1)
         data2_cvt = double(data2)
         result = dblarr(ny, nx)
      END
   END ELSE BEGIN
      IF (nargs EQ 2) THEN BEGIN
         data1_cvt = float(data1)
         result = fltarr(nx)
      END ELSE BEGIN
         data1_cvt = float(data1)
         data2_cvt = float(data2)
         result = fltarr(ny, nx)
      END
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_190, type, err_status, data1_cvt, data2_cvt, data3_cvt, nx, ny, $
                              xderiv_cvt, yderiv_cvt, $
                              domain_dim, struct_is_ppoly, $
                              result
                           
   ; return
   IF (nargs EQ 3) THEN RETURN, transpose(result) ELSE RETURN, result
END

                   
                   
                   

  
      

  
