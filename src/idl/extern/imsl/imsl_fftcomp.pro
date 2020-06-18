; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_fftcomp.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_fftcomp, a, $                   ;INPUT 1-d or 2-D array
                sine=sine, $               ;INPUT Scalar ON/OFF flag
                cosine=cosine, $           ;INPUT Scalar ON/OFF flag
                backward=backward, $       ;INPUT Scalar ON/OFF flag
                init_params=init_params, $ ;INPUT 1-D array.
                double=double, $           ;INPUT Scalar ON/OFF flag
                complex=complex            ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  The two basic uses of this function are divided between the
   ;  cases when the input array is either 1-D or 2-D.
   ;
   ;  Some of the following checks take place after the precision
   ;  has been determined.
   ;
   ;  Case 1: 1-D input:
   ;          - A must be a 1-D array. 
   ;            Set dim0 and dim1  w.r.t. A
   ;          - If INIT_PARAMS is specified, then:
   ;            Must be a 1-D array of length 2*dim0+15.
   ;            Must be a compatible data type to the data type 
   ;            that will be used to call the C/Math function.
   ;            
   ;  Case 2: 2-D input:
   ;          - A must be a 2-D array. Set n & m w.r.t. A
   ;          - INIT_PARAMS is not allowed.
   ;          
   nargs = n_params()
   IF (nargs NE 1) THEN  message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   n_dim = IMSL_LONG(size_a(0))
   IF ((n_dim NE 1) AND (n_dim NE 2)) THEN $
     message, 'A must be a 1-D OR 2-D array.'
   IF (n_dim EQ 1) THEN BEGIN
      dim0 = IMSL_LONG(size_a(1))
      dim1 = IMSL_0
   END ELSE BEGIN
      dim0 = IMSL_LONG(size_a(1))
      dim1 = IMSL_LONG(size_a(2))
      IF (KEYWORD_SET(init_params)) THEN $
        message, 'INIT_PARAMS is not valid with 2-D arrays.'
   END
   ; Perform checks on COSINE and SINE keywords.
   IF ((KEYWORD_SET(cosine) AND KEYWORD_SET(sine))) THEN $
	message, 'The keywords COSINE and SINE cannot be used together.'
   IF ((KEYWORD_SET(sine) AND KEYWORD_SET(complex))) THEN $
	message, 'The keywords SINE and COMPLEX cannot be used together.'
   IF ((KEYWORD_SET(cosine) AND KEYWORD_SET(complex))) THEN $
	message, 'The keywords COSINE and COMPLEX cannot be used together.'
   cos_sine = IMSL_0
   IF (KEYWORD_SET(sine)) THEN cos_sine = IMSL_1
   IF (KEYWORD_SET(cosine)) THEN cos_sine = IMSL_2
   IF ((cos_sine GT 0) and KEYWORD_SET(backward)) THEN $
	message, 'The keywords COSINE and SINE cannot be used with BACKWARD.'
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_a(N_ELEMENTS(size_a)-2) GT TYP_DOUBLE) THEN type = TYP_COMPLEX
   IF (size_a(N_ELEMENTS(size_a)-2) GT TYP_COMPLEX) THEN type = TYP_DCMPLX
   IF (KEYWORD_SET(complex)) THEN BEGIN
      if (type EQ TYP_DOUBLE) then type = TYP_DCMPLX ELSE type = TYP_COMPLEX
   END
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   ;
   ; C/Math supports only complex and double complex 2-D FFTs.
   IF ((n_dim eq 2) and (type eq TYP_FLOAT)) then type = TYP_COMPLEX
   IF ((n_dim EQ 2) AND (type EQ TYP_DOUBLE)) THEN type = TYP_DCMPLX
   ;
   ; Now that the precision is known, check INIT_PARAMS, if present.
   IF (KEYWORD_SET(init_params)) THEN BEGIN
      size_init = IMSL_SIZE(init_params)
      IF (cos_sine EQ 0) THEN BEGIN
        IF ((type EQ TYP_FLOAT) OR (type EQ (TYP_DOUBLE))) $
          THEN length = 2*dim0+15 ELSE length = 4*dim0+15
      ENDIF
      IF (cos_sine EQ 1) THEN BEGIN
        IF ((type EQ TYP_FLOAT) OR (type EQ (TYP_DOUBLE))) $
          THEN length = IMSL_LONG(2.5*dim0+15) ELSE length = IMSL_LONG(2*2.5*dim0+15)
      ENDIF
      IF (cos_sine EQ 2) THEN BEGIN
        IF ((type EQ TYP_FLOAT) OR (type EQ (TYP_DOUBLE))) $
          THEN length = IMSL_LONG(IMSL_3*dim0+15) ELSE length = IMSL_LONG(2*IMSL_3*dim0+15)
      ENDIF
      IF ((size_init(0) NE 1) OR (N_ELEMENTS(init_params) NE length)) $
        THEN message, 'INIT_PARAMS is not the correct size.'
      init_type = size_init(N_ELEMENTS(size_init)-2)
      IF (((type EQ TYP_FLOAT) OR (type EQ TYP_COMPLEX)) $
          AND (init_type NE TYP_FLOAT)) THEN $
        message, 'INIT_PARAMS is not a compatible data type.'
      IF (((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) $
          AND (init_type NE TYP_DOUBLE)) THEN $
        message, 'INIT_PARAMS is not a compatible data type.'
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   if (keyword_set(backward)) then backward_cvt = IMSL_1;
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type eq TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      IF (n_dim EQ 1) THEN result = dblarr(cmplx_scale*dim0) $
        ELSE result = dblarr(cmplx_scale*dim1, dim0)
   END ELSE BEGIN
      IF (n_dim EQ 1) THEN result = fltarr(cmplx_scale*dim0) $
        ELSE result = fltarr(cmplx_scale*dim1, dim0)
   END
   ; Since the result of the Sine and Cosine transform is (n+1), we
   ; need to make the result array larger.
   IF (cos_sine GT 0) THEN BEGIN
	if (type eq TYP_DOUBLE) THEN result = dblarr(dim0+1) $
            ELSE result = fltarr(dim0+1)
   ENDIF
   ;
   ; Convert A.
   ;
   a_cvt = imsl_cvt_arr(a, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_140, type, err_status, a_cvt, dim0, dim1, n_dim,$
     init_params, backward_cvt, cos_sine, result
   ;
   IF (cmplx_scale EQ 2) THEN $
     result = imsl_cvt_arr(result, /back) ELSE $
     IF (n_dim EQ 2) THEN result = transpose(result)
   ; Since the result of the Sine and Cosine transform was (n+1),we
   ; remove the last element, which is actually used for workspace.
   IF (cos_sine GT 0) THEN BEGIN
	result = result(0:dim0-1)
   ENDIF
   ; Return
   RETURN, result
END
