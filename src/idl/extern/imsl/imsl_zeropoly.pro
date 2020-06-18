; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_zeropoly.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_zeropoly, coef, $                   ;INPUT 1-d array
                companion=companion, $         ;INPUT Scalar ON/OFF flag
                jenkins_traub=jenkins_traub, $ ;INPUT Scalar ON/OFF flag
                double=double, $               ;INPUT Scalar ON/OFF flag
                complex=complex                ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;    - COEF must be a 1-D array.
   ;          - INIT_PARAMS is not allowed.
   ;          
   nargs = n_params()
   IF (nargs NE 1) THEN  message, 'Incorrect number of arguments.'
   size_coef = IMSL_SIZE(coef)
   IF (size_coef(0) NE 1) THEN $
     message, 'COEF must be a 1-D array.'
   num_coef = size_coef(1)
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_coef(N_ELEMENTS(size_coef)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_coef(N_ELEMENTS(size_coef)-2) GT TYP_DOUBLE) THEN type = TYP_COMPLEX
   IF (size_coef(N_ELEMENTS(size_coef)-2) GT TYP_COMPLEX) THEN type = TYP_DCMPLX
   IF (KEYWORD_SET(complex)) THEN BEGIN
      IF (type EQ TYP_DOUBLE) THEN type = TYP_DCMPLX ELSE type = TYP_COMPLEX
   END
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   ;
   ;
   ; Setup the parameters for the call to the system function.
   ;
   if (keyword_set(companion)) then companion_cvt = IMSL_1
   if (keyword_set(jenkins_traub)) then jenkins_cvt = IMSL_1;
   ;
   ; Floating point arguments and keywords
   ; Return type is always TYP_COMPLEX
   cmplx_scale = IMSL_2
   IF ((type eq TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      result = dblarr(cmplx_scale*(num_coef-1))
   END ELSE BEGIN
      result = fltarr(cmplx_scale*(num_coef-1))
   END
   ;
   ; Convert COEF.
   ;
   coef_cvt = imsl_cvt_arr(coef, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_197, type, err_status, coef_cvt, num_coef, $
                          companion_cvt, $
                          jenkins_cvt, $
                          result
   ;
   result = imsl_cvt_arr(result, /back)
   ; Return
   RETURN, result
END
