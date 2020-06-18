; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_pdsol.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
;
; Check_fac_pd is used by sp_pdsol to test if the argument passed
; for the LU factorization is valid.  This argument is only passed
; in the case that the system has been factored by SP_PDFAC previously.
;
Pro imsl_Check_fac_pd, factor,  nz, type_a

@imsl_init.pro

   err_str = "FACTOR must be a structure containing the sparse factorization."

   ; See if it is a structure.
   size_fac = IMSL_SIZE(factor)
   IF ((size_fac(0) NE 1) OR $
       (size_fac(N_ELEMENTS(size_fac)-2) NE 8)) THEN message, err_str

   ; Check out the tag names.
   tags = tag_names(factor)
   IF (N_ELEMENTS(tags) NE 6) THEN message, err_str
   expected_tags = ['NZSUB', 'XNZSUB', 'XLNZ', 'ALNZ', 'PERM', 'DIAG']
   FOR i = 0, 5 DO $
     IF (tags(i) NE expected_tags(i)) THEN message, err_str

   nz  =  IMSL_N_ELEMENTS(factor.alnz)+IMSL_N_ELEMENTS(factor.diag)
   size_fac_alnz = IMSL_SIZE(factor.alnz)
   type_a = (size_fac_alnz(N_ELEMENTS(size_fac_alnz)-2))
   if (type_a eq 9) then type_a = TYP_DCMPLX ;; rsi change of type
END

Function  imsl_Sp_pdsol, b, $                          ;INPUT 1-D array: floating point
                    a, $                          ;INPUT Struct
                    multifrontal=multifrontal, $  ;INPUT Scalar LONG
                    factor=factor, $              ;INPUT Struct
                    csc_col=csc_col, $            ;INPUT 1-D array: LONG
                    csc_row_ind=csc_row_ind, $    ;INPUT 1-D array: LONG
                    csc_val=csc_val, $            ;INPUT 1-D array: floating point
                    sm_diag=sm_diag, $            ;OUTPUT Scalar floating point
                    lg_diag=lg_diag, $            ;OUTPUT Scalar floating point
                    n_nonzero=n_nonzero           ;OUTPUT Scalar LONG


@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Checks on positional args:
   ; There are two different cases for the number of positional argument
   ; passed to this function.
   ; CASE 1: One (1) positional argument:
   ;  In this case, the positional argument is the RHS, b.
   ;  If CSC_* are not provided, then FACTOR
   ;  must be provided, and we are doing a solve-only call.
   ;  The following checks are made:
   ;  - B must be a 1-D array of length N_ROWS.
   ;  If CSC_* are provided:
   ;  - CSC_COL is a 1-D array of length (N_ROWS + 1)
   ;  - CSC_ROW_IND is a 1-D array of length CSC_COL(N_ROWS)
   ;  - CSC_VAL is a 1-D array of length CSC_COL(N_ROWS)
   ;  If CSC are not provided
   ;  - FACTOR must be provided.
   ;  - FACTOR is checked out to see if is the correct type of stucture.
   ;  - If NUM_FACTOR is provided, then
   ;     MULTIFRONTAL, SM_DIAG, LG_DIAG, N_NONZERO are not allowed.
   ;
   ; CASE 2: Two (2) positional arguments:
   ;  In this case, the positional arguments the RHS, b, and the sparse matrix, A,
   ;  stored as a 1-D array of type IMSL_D_SP_ELEM.
   ;  - B must be a 1-D array of length N_ROWS.
   ;  - A must be an array of structures of type IMSL_D_SP_ELEM.
   ;  - FACTOR is not allowed.
   ;  - The row/column indices are checked against teh size of the matrix
   ;    as defined by B  This check helps avoid the case when B
   ;    is input too long, which can lead to a crash within CNL.
   ;
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2)) THEN $
     message, 'Incorrect number of arguments.'
   type_a = IMSL_0
   factored = IMSL_0
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN  $
     message, 'B must be a 1-D array of length N_ROWS.'
   n_rows = IMSL_LONG(size_b(3))
   IF (nargs EQ 1) THEN BEGIN  ; Either doing solve-only, or CSC format.
      tmp = KEYWORD_SET(CSC_COL) + KEYWORD_SET(CSC_ROW_IND) + KEYWORD_SET(CSC_VAL)
      IF ((tmp NE 3) AND (tmp NE 0)) THEN BEGIN
         message, 'The keywords CSC_COL, CSC_ROW_IND, AND CSC_VAL ' + $
           'must be used together.'
      END ELSE BEGIN
         IF (tmp EQ 3) THEN BEGIN
            ; Check that A is in CSC format.
            imsl_Sp_ccsc, csc_row_ind, csc_col, csc_val, n_rows
            nz = IMSL_N_ELEMENTS(csc_row_ind)
         END ELSE BEGIN ; (Thus tmp eq 0, and we are doing a solve-only step)
            IF (NOT KEYWORD_SET(factor)) THEN $
              message, 'FACTOR is required as input if only B is supplied.'
            imsl_Check_fac_pd, factor, nz, type_a
            IF ((KEYWORD_SET(MULTIFRONTAL)+ARG_PRESENT(SM_DIAG) + $
                 ARG_PRESENT(LG_DIAG) + ARG_PRESENT(N_NONZERO)) GT 0) THEN $
              message, 'If either FACTOR is provided, then ' + $
                       'MULTIFRONTAL, SM_DIAG, LG_DIAG, N_NONZERO are not allowed.'
            factored = IMSL_1
         END
      END
   END ELSE BEGIN               ;  Two positional arguments.
      ; Test that A is in coordinate format.
      imsl_sp_acrd, a, nz, type_a
      IF ((KEYWORD_SET(CSC_COL) + KEYWORD_SET(CSC_ROW_IND) +  $
           KEYWORD_SET(CSC_VAL)) GT 0) THEN $
        message, 'The keywords CSC_COL, CSC_ROW_IND, AND CSC_VAL are not allowed ' + $
           'when two  positional argument are passed.'
      IF (KEYWORD_SET(factor)) THEN $
        message, 'FACTOR is not allowed if both B and A are present.'
      index_max   =   MAX(a.row)   >   MAX(a.col)
      IF (index_max LT N_ELEMENTS(b)-1-1) THEN $
        MESSAGE, "The largest row/column index of A not correspond to the size of B."
   END
   ;
   ; Decide on what precision to use.
   ; Note that the datatype is set based on the type of A.  This way,
   ; we will not be converting the type of A, and we will avoid creating
   ; a copy of A (which is potentially very large).
   ;
   type = TYP_FLOAT
   IF (nargs EQ 1) THEN BEGIN
      IF (factored EQ 0) THEN type = ((size_csc_val(N_ELEMENTS(size_csc_val)-2)) > TYP_FLOAT) $
                         ELSE type = type_a
   END ELSE BEGIN
      type  =  type_a
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Output LONG keyword(s)
   IF (KEYWORD_SET(n_nonzero) EQ TRUE) THEN $
     n_nonzero_spc = IMSL_0
   ;
   ; Input LONG argument(s)
   IF (KEYWORD_SET(multifrontal) EQ TRUE) THEN $
     multifrontal_cvt = IMSL_1
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      ;
      ; Result
      ;
      result  =  DBLARR(cmplx_scale*n_rows)
      IF (ARG_PRESENT(sm_diag) EQ TRUE) THEN $
        sm_diag_spc = double(0.0)
      IF (ARG_PRESENT(lg_diag) EQ TRUE) THEN $
        lg_diag_spc = double(0.0)
   END ELSE BEGIN
      ;
      ; Result
      ;
      result  =  FLTARR(cmplx_scale*n_rows)
      IF (ARG_PRESENT(sm_diag) EQ TRUE) THEN $
        sm_diag_spc = float(0.0)
      IF (ARG_PRESENT(lg_diag) EQ TRUE) THEN $
        lg_diag_spc = float(0.0)
   END

   ; Input
   IF (nargs EQ 1) THEN $
     IF (factored EQ 0) THEN csc_val_cvt  =  imsl_cvt_arr(csc_val, type)
      ;
   b_cvt    =    imsl_cvt_arr(b,   type)
   ;
   err_status = 0L
   IF (factored EQ 0) THEN BEGIN
        MATHSTAT_212, LONG(type), err_status, b_cvt, a, n_rows, nz, $
                  n_nonzero_spc, $
                  sm_diag_spc, $
                  lg_diag_spc, $
                  multifrontal_cvt, $
                  csc_col_cvt, $
                  csc_row_cvt, $
                  csc_val_cvt, $
                  factored, $
                  t, $
                  t, $
                  t, $
                  t, $
                  t, $
                  t, $
                  result
   END ELSE BEGIN
        MATHSTAT_212, LONG(type), err_status, b_cvt, a, n_rows, nz, $
                  n_nonzero_spc, $
                  sm_diag_spc, $
                  lg_diag_spc, $
                  multifrontal_cvt, $
                  csc_col_cvt, $
                  csc_row_cvt, $
                  csc_val_cvt, $
                  factored, $
                  factor.nzsub, $
                  factor.xnzsub, $
                  factor.xlnz, $
                  factor.perm, $
                  factor.alnz, $
                  factor.diag, $
                  result
   END


   IF (ARG_PRESENT(n_nonzero) EQ TRUE) THEN n_nonzero = n_nonzero_spc
   IF (ARG_PRESENT(sm_diag) EQ TRUE) THEN sm_diag = sm_diag_spc
   IF (ARG_PRESENT(lg_diag) EQ TRUE) THEN lg_diag = lg_diag_spc
   ; Return.
   ;
   IF (cmplx_scale EQ 2) $
     THEN RETURN, imsl_cvt_arr(result, /back) $
     ELSE RETURN, result


END
