; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_matrix_norm.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_Matrix_norm,  arg1,  $
                       arg2,   $
                       arg3,   $
                       arg4,   $
                       arg5,    $
                       symmetric  =  symmetric,      $    ;INPUT Scalar ON/OFF flag
                       inf_norm=inf_norm, $               ;INPUT Scalar ON/OFF flag
                       one_norm=one_norm, $               ;INPUT Scalar ON/OFF flag
                       double=double                      ;INPUT Scalar ON/OFF flag

@imsl_init.pro ; ITTVIS
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Case 1: Norm of a dense rectangular matrix.
   ;   Args:
   ;         a       Dense matrix.
   ;   Restrictions:
   ;   Error checks:
   ;         a is a 2D array.
   ;         a is not complex.
   ;         Keyword SYMMETRIC is not allowed.
   ;   Data types supported:
   ;         float/double.
   ;
   ; Case 2: Norm of a matrix stored in band format.
   ;   Args:
   ;         n      Order of matrix.
   ;         nlca   Number of lower codiagonals.
   ;         nuca   Number of upper codiagonals.
   ;         a      Sparse matrix stored in band format.
   ;   Data types supported:
   ;         float/double
   ;   Error checks:
   ;         n GT 0
   ;         nlca GE 0 and LE n
   ;         nuca GE 0 and LE n
   ;         If SYMMETRIC is set, then
   ;         - nlca EQ nuca
   ;         a is a array of length
   ;             (nlca+nuca+1)*n
   ;             (2*nlca+1)*n if SYMMETRIC is set.
   ;
   ; Case 3: Norm of a matrix stored in coordinate format.
   ;   Args:
   ;         nrows  Number of rows.
   ;         ncols  Number of columns.
   ;         a      Sparse matrix stored in coordinate format.
   ;   Data types supported:
   ;         double for band/coordinate.
   ;   Error checks:
   ;  - A must be an array of structures of type IMSL_D_SP_ELEM.
   ;
   nargs = n_params()
   IF ((nargs NE 4) AND (nargs NE 3)AND (nargs NE 1)) THEN $
     message, "Incorrect number of arguments."

   IF (nargs EQ 1) THEN BEGIN
      size_a = IMSL_SIZE(arg1)
      IF (size_a(0) NE 2) THEN MESSAGE,  'A must be a 2-D array.'
      n_rows  =  IMSL_LONG(size_a(1))
      n_cols  =  IMSL_LONG(size_a(2))
      ;
      ; Decide on what precision to use.
      ;
      type  =  TYP_FLOAT
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type  =  TYP_DOUBLE
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN $
        MESSAGE, 'Complex data is not allowed.'
      IF (KEYWORD_SET(symmetric)) THEN $
        MESSAGE,  'The keyword SYMMETRIC is not allowed with dense matrices.'
      IF (KEYWORD_SET(DOUBLE) EQ true) THEN type   =   TYP_DOUBLE
      IF (KEYWORD_SET(DOUBLE) EQ true) THEN BEGIN
         result      =      DOUBLE(0.0)
         a_cvt  =  DOUBLE(TRANSPOSE(arg1))
      END ELSE BEGIN
         result    =    FLOAT(0.)
         a_cvt  =  FLOAT(TRANSPOSE(arg1))
      END
      ;
      ; Call the system function.
      ;
      IF (KEYWORD_SET(one_norm)) THEN one_norm_set  =  IMSL_1 ELSE one_norm_set  =  IMSL_0
      IF (KEYWORD_SET(inf_norm)) THEN inf_norm_set    =    IMSL_1 ELSE inf_norm_set    =    IMSL_0
      err_status = 0L
      MATHSTAT_215,   type,   err_status,   n_rows,   n_cols,   $
        a_cvt,   one_norm_set,  inf_norm_set, result
      ;
   END
   IF (nargs EQ 3) THEN BEGIN
      ; Coordinate case
      n_rows = IMSL_LONG(arg1(0))
      n_cols = IMSL_LONG(arg2(0))
      IF (n_rows LT 1) THEN MESSAGE,  'N_ROWS must be greater than 0.'
      IF (n_cols LT 1) THEN MESSAGE,   'N_COLS must be greater than 0.'
      ; Test that A is in coordinate format.
      imsl_sp_acrd, arg3, nz, type_a
      type  =  type_a
      IF (type GT TYP_DOUBLE) THEN MESSAGE, 'Complex data is not allowed.'
      ;
      ; Setup the parameters for the call to the system function.
      ;
      IF (type EQ TYP_DOUBLE) THEN BEGIN
         ;
         ; Using double precision
         ;
         result = DOUBLE(0.0)
      END ELSE BEGIN
         ;
         ; Using single precision
         ;
         result = FLOAT(0.0)
      END
      ;
      IF (MAX(arg3.row) GE n_rows) THEN $
        MESSAGE,  'N_ROWS is too small for the largest row index in A.'
      IF (MAX(arg3.col) GE n_cols) THEN $
        MESSAGE,  'N_COLS is too small for  the largest column index in A.'
      ;
      ; Call the system function.
      ;
      err_status = 0L
      IF (KEYWORD_SET(symmetric)) THEN symm_set  =  IMSL_1 ELSE symm_set  =  IMSL_0
      IF (KEYWORD_SET(one_norm)) THEN one_norm_set  =  IMSL_1 ELSE one_norm_set  =  IMSL_0
      IF (KEYWORD_SET(inf_norm)) THEN inf_norm_set    =    IMSL_1 ELSE inf_norm_set    =    IMSL_0
      MATHSTAT_216,      type,      err_status,      n_rows,      n_cols,      nz,      $
        arg3,   one_norm_set,  inf_norm_set,  symm_set,  result
      ;
   END
   IF (nargs EQ 4) THEN BEGIN
      ; Band case
      n = IMSL_LONG(arg1(0))
      nlca = IMSL_LONG(arg2(0))
      nuca = IMSL_LONG(arg3(0))
      IF (n LT 1) THEN message, 'N must be greater than 0.'
      IF ((nlca LT 0) OR (nlca GT n)) $
        THEN message, 'NLCA must be greater than or equal to zero ' + $
                      'and less than N.'
      IF ((nuca LT 0) OR (nuca GT n)) $
        THEN message, 'NUCA must be greater than or equal to zero ' + $
        'and less than.'
      ;
      ; In symmetric case, NLCA must be equal to NUCA.
      IF KEYWORD_SET(symmetric) THEN $
        IF (nuca LT nlca) THEN $
          MESSAGE,  'NUCA must be equal to NLCA if SYMMETRIC is set.'
      ; Check out A.
      size_a = IMSL_SIZE(arg4)
      type_a = size_a(N_ELEMENTS(size_a)-2)
      IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_COMPLEX) THEN $
        MESSAGE,  'Complex data is not allowed.'
      IF KEYWORD_SET(symmetric) THEN exp_a_lngth   =   (n*(nlca+1)) $
        ELSE exp_a_lngth = (n*(nlca+nuca+1))
      IF ((size_a(0) NE 1) OR (N_ELEMENTS(arg4) NE  exp_a_lngth)) $
        THEN message, 'A is not the correct size.'
      ; Decide on the data type
      type = TYP_FLOAT
      IF (type_a EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
      IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
      ;
      ; Setup the parameters for the call to the system function.
      ;
      IF (type EQ TYP_DOUBLE) THEN BEGIN
         ;
         ; Using double precision
         ;
         result = DOUBLE(0.0)
         IF KEYWORD_SET(symmetric) THEN $
           a_cvt = dblarr((nlca+nuca+1)*n) $
           ELSE a_cvt = double(arg4)
      END ELSE BEGIN
         ;
         ; Using single precision
         ;
         result = FLOAT(0.0)
         IF KEYWORD_SET(symmetric) THEN $
           a_cvt = fltarr((nlca+nuca+1)*n) $
           ELSE a_cvt = float(arg4)
      END
      ;
      ; If SYMMETRIC set, convert a to general band storage.
      ;
      IF KEYWORD_SET(symmetric) THEN BEGIN
         a_cvt(0:(nlca+1)*n - 1) = arg4
         FOR i = IMSL_0, nlca-1 DO $
           a_cvt((IMSL_2*nlca-i)*n : (IMSL_2*nlca-i+1)*n -1 - (nlca-i)) = $
           arg4(i*n+(nlca-i):((i+1)*n)-1)
      END
      ;
      ;Call thge system function.
      err_status = 0L
      IF (KEYWORD_SET(symmetric)) THEN symm_set  =  IMSL_1 ELSE symm_set  =  IMSL_0
      IF (KEYWORD_SET(one_norm)) THEN one_norm_set  =  IMSL_1 ELSE one_norm_set  =  IMSL_0
      IF (KEYWORD_SET(inf_norm)) THEN inf_norm_set  =  IMSL_1 ELSE inf_norm_set  =  IMSL_0
      MATHSTAT_217,  type,  err_status,  n,  nlca,  nuca,  a_cvt,  $
         one_norm_set, inf_norm_set, symm_set, result
   END ; End of banded case.
   ;
   ; Return.
   ;
   RETURN, result
END



