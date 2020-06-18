; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_mvmul.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Sp_mvmul, arg1, $
                   arg2, $
                   arg3, $
                   arg4, $
                   arg5, $
                   arg6, $
                   symmetric = symmetric

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Case 1: Multiply a sparse matrix stored in coordinate
   ;         format and a dense vector.
   ;   Args:
   ;         n_rows  Number of rows.
   ;         n_cols  Number of columns.
   ;         a       Sparse matrix stored in coordinate format.
   ;         x       Vector of length ncol
   ;   Data types supported:
   ;         IMSL_D_SP_ELEM, thus double precision only.
   ;   Restrictions:
   ;   Error checks:
   ;         nrow GT 0.
   ;         ncol GT 0.
   ;         nrow EQ ncol (Currently. See 'Restrictions' below)
   ;         a is a sparse matrix stored in coordinate format.
   ;         x is 1-D array of length ncol.
   ;
   ; Case 2: Multiply a sparse matrix stored in band
   ;         format and a dense vector.
   ;   Args:
   ;         n_rows  Number of rows.
   ;         n_cols  Number of columns.
   ;         nlca   Number of lower codiagonals.
   ;         nuca   Number of upper codiagonals.
   ;         a      Sparse matrix stored in band format.
   ;         x      Vector of length ncol
   ;   Data types supported:
   ;         float/double and complex only.
   ;   Restrictions:
   ;   Error checks:
   ;         nrow GT 0.
   ;         ncol GT 0.
   ;         nlca GE 0
   ;         nuca GE 0
   ;         If SYMMETRIC is set, then
   ;         - nlca EQ nuca
   ;         a is a array of length
   ;             (nlca+nuca+1)*n
   ;             (2*nlca+1)*n if SYMMETRIC is set.
   ;         x is 1-D array of length n.
   ;
   ;
   ;
   nargs = n_params()
   IF ((nargs NE 4) AND (nargs NE 6)) THEN $
     message, "Incorrect number of arguments."

   IF (nargs EQ 4) THEN BEGIN
      ; Coordinate case
      n_rows = IMSL_LONG(arg1(0))
      n_cols = IMSL_LONG(arg2(0))
      IF (n_rows LT 1) THEN message, 'N_ROWS must be greater than 0.'
      IF (n_cols LT 1) THEN message, 'N_COLS must be greater than 0.'
      IF (n_cols NE n_rows) THEN message, 'N_ROWS and N_COLS must be equal.'
      ; Test that A is in coordinate format.
      imsl_sp_acrd, arg3, nz, type_a
      size_x = IMSL_SIZE(arg4)
      IF ((size_x(0) NE 1) OR (N_ELEMENTS(arg4) NE n_rows)) THEN  $
        message, "X must be a 1-D array of length N_ROWS."
      type  =  type_a
      x = imsl_cvt_arr(arg4, type)
      cmplx_scale = IMSL_1
      IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
      IF ((type EQ TYP_DOUBLE) OR (type EQ  TYP_DCMPLX)) THEN BEGIN
         ;
         ; Result
         ;
         result  =  DBLARR(cmplx_scale*n_rows)
      END ELSE BEGIN
         ;
         ; Result
         ;
         result  =  FLTARR(cmplx_scale*n_rows)
      END
      ;
      ; Call the system function.
      ;
      err_status = 0L
      IF (KEYWORD_SET(symmetric)) THEN BEGIN
         a = 'a'
         MATHSTAT_210, type, err_status, n_rows, n_cols, nz, arg3, x, trans, result
         IF (cmplx_scale EQ 2) THEN result_tmp = imsl_cvt_arr(result,   /back)  ELSE result_tmp = result
         trans = IMSL_1
         MATHSTAT_210, type, err_status, n_rows, n_cols, nz, arg3, x, trans, result
         IF (cmplx_scale EQ 2)  THEN result_tmp2 = imsl_cvt_arr(result,   /back)  ELSE result_tmp2 = result
         diag = imsl_sp_diag(arg3, n_rows, /symmetric)
         result = result_tmp+result_tmp2-x*diag
      END ELSE BEGIN
        MATHSTAT_210, type, err_status, n_rows, n_cols, nz, arg3, x, trans, result
        IF (cmplx_scale EQ 2) THEN  result  = imsl_cvt_arr(result,   /back)
     END
     RETURN,    result
    END ELSE BEGIN
      ; Band case
      n_rows = IMSL_LONG(arg1(0))
      n_cols = IMSL_LONG(arg2(0))
      nlca = IMSL_LONG(arg3(0))
      nuca = IMSL_LONG(arg4(0))
      IF (n_rows LT 1) THEN message, 'N_ROWS must be greater than 0.'
      IF (n_cols LT 1) THEN message, 'N_ROWS must be greater than 0.'
      IF ((nlca LT 0) OR (nlca GE n_rows)) $
        THEN message, 'NLCA must be greater than or equal to zero ' + $
                      'and less than N_ROWS.'
      IF ((nuca LT 0) OR (nuca GT n_cols)) $
        THEN message, 'NUCA must be greater than or equal to zero ' + $
                      'and less than N_COLS.'
      ;
      ; In symmetric case, NLCA must be equal to NUCA.
      IF KEYWORD_SET(symmetric) THEN $
        IF (nuca LT nlca) THEN $
        message, 'NUCA must be equal to NLCA if SYMMETRIC is set.'
      ;
      ; Check out A.
      size_a = IMSL_SIZE(arg5)
      type_a  =  size_a(N_ELEMENTS(size_a)-2)
      n = n_cols
      IF KEYWORD_SET(symmetric) THEN exp_a_lngth = (n*(nlca+1)) $
        ELSE exp_a_lngth  =  (n*(nlca+nuca+1))
      IF ((size_a(0) NE 1) OR (N_ELEMENTS(arg5) NE  exp_a_lngth)) $
        THEN message, 'A is not the correct size.'
      ;
      ; Check out X.
      size_x = IMSL_SIZE(arg6)
      type_x = size_x(N_ELEMENTS(size_x)-2)
      IF ((size_x(0) NE 1) OR (N_ELEMENTS(arg6) NE n_cols)) THEN  $
        message, 'X must be a 1-D array of length N_COLS.'
      ; Decide on the data type
      type  =  TYP_FLOAT
      IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = (type > TYP_DOUBLE)
      IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_COMPLEX) THEN type = (type > TYP_COMPLEX)
      IF (size_x(N_ELEMENTS(size_x)-2) EQ TYP_DCMPLX) THEN type = (type > TYP_DCMPLX)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = (type > TYP_DOUBLE)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN type = (type > TYP_COMPLEX)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ TYP_DCMPLX) THEN type = (type > TYP_DCMPLX)
      ;
      ; Setup the parameters for the call to the system function.
      ;
      IF (type EQ TYP_FLOAT) THEN BEGIN
         ;
         ; Using float precision
         ;
         result = fltarr(n)
         x_cvt = float(arg6)
         IF KEYWORD_SET(symmetric) THEN $
           a_cvt = fltarr((nlca+nuca+1)*n) $
           ELSE a_cvt = float(arg5)
      END
      IF (type EQ TYP_DOUBLE) THEN BEGIN
         ;
         ; Using double precision
         ;
         result = dblarr(n)
         x_cvt = double(arg6)
         IF KEYWORD_SET(symmetric) THEN $
           a_cvt = dblarr((nlca+nuca+1)*n) $
           ELSE a_cvt = double(arg5)
      END
      IF (type EQ TYP_COMPLEX) THEN  BEGIN
         ;
         ; Using complex data.
         ;
         result = complexarr(n)
         x_cvt = imsl_cvt_arr(arg6, type)
         IF KEYWORD_SET(symmetric) THEN $
           a_cvt = complexarr((nlca+nuca+1)*n) $
         ELSE $
           a_cvt = imsl_cvt_arr(arg5, type)
      END
      IF (type EQ TYP_DCMPLX) THEN  BEGIN
         ;
         ; Using dcomplex data.
         ;
         result = dcomplexarr(n)
         x_cvt = imsl_cvt_arr(arg6, type)
         IF KEYWORD_SET(symmetric) THEN $
           a_cvt = dcomplexarr((nlca+nuca+1)*n) $
         ELSE $
           a_cvt = imsl_cvt_arr(arg5, type)
      END
      ;
      ; If SYMMETRIC set, convert a to general band storage.
      ;
      IF KEYWORD_SET(symmetric) THEN BEGIN
         a_cvt(0:(nlca+1)*n - 1) = arg5
         FOR i = IMSL_0, nlca-1 DO $
           a_cvt((IMSL_2*nlca-i)*n : (IMSL_2*nlca-i+1)*n -1 - (nlca-i)) = $
           arg5(i*n+(nlca-i):((i+1)*n)-1)
      END
      ;
      ;Call thge system function.
      err_status = 0L
      MATHSTAT_209, type, err_status, n_rows, n_cols, nlca, nuca, a_cvt, $
        x_cvt,  result
   END ; End of banded case.
   ;
   ; Return.
   ;
   RETURN, result
END



