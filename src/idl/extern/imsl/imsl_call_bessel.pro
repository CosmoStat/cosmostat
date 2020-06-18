; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_call_bessel.pro#1 $   
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_Call_bessel, order, $             ;INPUT scalar floating point 
                z, $ ;INPUT scalar or 1-D array: floating point or complex
                fcn_name=fcn_name, $       ;INPUT Scalar string.
                double=double, $           ;INPUT Scalar ON/OFF flag
                sequence=sequence          ;INPUT Scalar LONG

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - ORDER
   ;         This argument defines the scalar order, thus it must be
   ;         a scalar.
   ;    - Z
   ;         There are two cases allowed here:
   ;         1. Z a scalar.
   ;         2. Z a 1-D array.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_order = IMSL_SIZE(order)
   IF (N_ELEMENTS(order) NE 1) THEN message, 'ORDER must be a scalar.'
   size_z = IMSL_SIZE(z)
   IF ((size_z(0) NE 0) AND (size_z(0) NE 1)) $
     THEN message, 'Z must be either a scalar or a 1-D array.'
   ;
   ; Decide on what precision to use.
   ; Select the data type to use.
   ; since the return type of this function is always TYP_COMPLEX,
   ; we are trying to decide whether to call the single precision
   ; complex C/Math routine, or the double precision complex C/Math routine
   ;
   type = TYP_COMPLEX
   IF (size_order(N_ELEMENTS(size_order)-2) EQ  TYP_DOUBLE) THEN type = TYP_DCMPLX
   IF (size_z(N_ELEMENTS(size_z)-2) EQ  TYP_DOUBLE) THEN type = TYP_DCMPLX
   IF (size_z(N_ELEMENTS(size_z)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DCMPLX
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   sequence_cvt = IMSL_1
   IF (KEYWORD_SET(sequence)) THEN sequence_cvt = (1 > IMSL_LONG(sequence(0)))
   ;
   ; Output LONG keyword(s)
   ;
   ; Floating point arguments and keywords
   ; Result vector.
   IF (type EQ TYP_DCMPLX) THEN BEGIN
      order_cvt = double(order)
      z_cvt = imsl_cvt_arr([z], type)
      result = dblarr(2*sequence_cvt, N_ELEMENTS(z))
   END ELSE BEGIN
      order_cvt = float(order)
      z_cvt = imsl_cvt_arr([z], type)
      result = fltarr(2*sequence_cvt, N_ELEMENTS(z))
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_104, Type, err_status, order_cvt, z_cvt, sequence_cvt, $
     IMSL_N_ELEMENTS(z), fcn_name, result
;Stop
; print, type, err_status, order_cvt, sequence_cvt, imsl_n_elements(z), ' ', fcn_name,max(result)
   ;
   ; Now copy over the result to the correct storage format.
   ; If result is of length 1, return it as a scalar.
   ;
   result = imsl_cvt_arr([result], /back)
   IF (N_ELEMENTS(result) EQ 1) THEN result = result(0)
   RETURN, result
END

                   
                   
                   

  
      

  
