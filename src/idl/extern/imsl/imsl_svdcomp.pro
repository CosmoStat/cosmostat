; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_svdcomp.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_svdcomp, a, $                       ;INPUT 2-D array: floating point 
                tol_rank=tol_rank, $       ;INPUT Scalar floating point
                complex=complex, $         ;INPUT Scalar ON/OFF flag
                double=double, $           ;INPUT Scalar ON/OFF flag
                inverse=inverse, $         ;OUTPUT 2-D array: floating point
                rank=rank, $               ;OUTPUT Scalar floating point
                u=u, $                     ;OUTPUT 2-D array: floating point
                v=v                        ;OUTPUT 2-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - A must be a 2-D array.
   ;          
   nargs = n_params()
   IF (nargs NE 1) THEN $
     message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D square array.'
   m = IMSL_LONG(size_a(1))
   n = IMSL_LONG(size_a(2))
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = TYP_COMPLEX
   IF (size_a(N_ELEMENTS(size_a)-2) EQ TYP_DCMPLX) THEN type = TYP_DCMPLX
   IF (KEYWORD_SET(double) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_DOUBLE
      IF (type EQ TYP_COMPLEX) THEN type = TYP_DCMPLX
   END
   IF (KEYWORD_SET(complex) EQ true) THEN begin
      IF (type EQ TYP_FLOAT) THEN type = TYP_COMPLEX
      IF (type EQ TYP_DOUBLE) THEN type = TYP_DCMPLX
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      result = dblarr(cmplx_scale*(m > n))
      ; 
      ; Input
      tmp = imsl_machine(/double)
      tol_rank_cvt = 100.*tmp.(3)
      IF (KEYWORD_SET(tol_rank)) THEN tol_rank_cvt = double(tol_rank(0))
      ; 
      ; Output
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = dblarr(cmplx_scale*m, n)
      rank_spc = IMSL_0
      IF (ARG_PRESENT(u)) THEN u_spc = dblarr(cmplx_scale*(m < n), m)
      IF (ARG_PRESENT(v)) THEN v_spc = dblarr(cmplx_scale*(m < n), n)
   END ELSE BEGIN
      result = fltarr(cmplx_scale*(m > n))
      ; 
      ; Input
      tmp = imsl_machine(/float)
      tol_rank_cvt = 100.*tmp.(3)
      IF (KEYWORD_SET(tol_rank)) THEN tol_rank_cvt = float(tol_rank(0))
      ; 
      ; Output
      IF (ARG_PRESENT(inverse)) THEN inverse_spc = fltarr(cmplx_scale*m, n)
      rank_spc = IMSL_0
      IF (ARG_PRESENT(u)) THEN u_spc = fltarr(cmplx_scale*(m < n), m)
      IF (ARG_PRESENT(v)) THEN v_spc = fltarr(cmplx_scale*(m < n), n)
   END
   ;
   ; Convert A
   ;
   a_cvt = imsl_cvt_arr(a, type)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_193, type, err_status, a_cvt, m, n, $ 
                              inverse_spc, $
                              tol_rank_cvt, $
                              rank_spc, $
                              u_spc, $
                              v_spc, $
                              result
   ;
   ; Now copy over all output keywords results.
   IF (arg_present(rank)) THEN rank = rank_spc
   IF (arg_present(inverse)) THEN $
     IF (cmplx_scale EQ 2) THEN inverse = imsl_cvt_arr(inverse_spc, /back) $
     ELSE inverse = transpose(inverse_spc)
   IF (arg_present(u)) THEN $
     IF (cmplx_scale EQ 2) THEN u = imsl_cvt_arr(u_spc, /back) $
     ELSE u = transpose(u_spc)
   IF (arg_present(v)) THEN $
     IF (cmplx_scale EQ 2) THEN v = imsl_cvt_arr(v_spc, /back) $
     ELSE v = transpose(v_spc)
   ;
   ; return
   ; Since result will *always* have a zero imaginary part,
   ; we convert to TYP_FLOAT or TYP_DOUBLE.
   ; TYP_DOUBLE. This is done a little further down in the code.
   IF (cmplx_scale EQ 2) THEN $
     IF (type EQ TYP_COMPLEX) THEN result = float(imsl_cvt_arr(result, /back)) $
     ELSE result = double(imsl_cvt_arr(result, /back))

	  outsize = (m < n)-1
	  result = result(0:outsize)
   
   RETURN, result
END
