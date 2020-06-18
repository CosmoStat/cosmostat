; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_fftinit.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_fftinit, n, $                     ;INPUT Scalar LONG
                sine=sine, $               ;INPUT Scalar ON/OFF flag
                cosine=cosine, $           ;INPUT Scalar ON/OFF flag
                double=double, $           ;INPUT Scalar ON/OFF flag
                complex=complex            ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;          
   nargs = n_params()
   IF (nargs NE 1) THEN  message, 'Incorrect number of arguments.'
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(double)) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(complex)) THEN type = TYP_COMPLEX
   IF ((KEYWORD_SET(double) AND KEYWORD_SET(complex))) THEN type = TYP_DCMPLX
   IF ((KEYWORD_SET(cosine) AND KEYWORD_SET(sine))) THEN $
	message, 'The keywords COSINE and SINE cannot be used together.'
   IF ((KEYWORD_SET(sine) AND KEYWORD_SET(complex))) THEN $
	message, 'The keywords SINE and COMPLEX cannot be used together.'
   IF ((KEYWORD_SET(cosine) AND KEYWORD_SET(complex))) THEN $
	message, 'The keywords COSINE and COMPLEX cannot be used together.'
   cos_sine = IMSL_0
   IF (KEYWORD_SET(sine)) THEN cos_sine = IMSL_1
   IF (KEYWORD_SET(cosine)) THEN cos_sine = IMSL_2
   n_cvt = IMSL_LONG(n(0))

   cmplx_scale = IMSL_1
   IF ((type EQ TYP_COMPLEX) OR (type EQ TYP_DCMPLX)) THEN $
     cmplx_scale = IMSL_2
   length = IMSL_LONG(cmplx_scale*2*n_cvt + 15L)
   if (KEYWORD_SET(sine)) THEN length = IMSL_LONG(2.5*n_cvt + 15L)
   if (KEYWORD_SET(cosine)) THEN length = IMSL_LONG(3*n_cvt + 15L)
   IF ((type EQ TYP_COMPLEX) OR (type EQ TYP_FLOAT)) $
     THEN result = fltarr(length > 1) $
     ELSE result = dblarr(length > 1)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_141, type, err_status, n_cvt, length, cos_sine, result
   ;
   ; Return
   RETURN, result
END
