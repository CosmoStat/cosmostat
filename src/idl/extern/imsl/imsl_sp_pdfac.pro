; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_pdfac.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
Function  imsl_Sp_pdfac, a, $                          ;INPUT Struct
                    n_rows, $                     ;INPUT Scalar LONG
                    multifrontal=multifrontal, $  ;INPUT Scalar LONG
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
   ; CASE 1: Zero (0) positional argument:
   ;  In this case, the coefficient matrix should be supplied by the
   ;  CSC_* keywords.
   ;  The following checks are made:
   ;  - The keywords CSC_COL CSC_ROW_IND CSC_VAL are required.
   ;  - CSC_COL is a 1-D array of length (N_ROWS + 1)
   ;  - CSC_ROW_IND is a 1-D array of length CSC_COL(N_ROWS)
   ;  - CSC_VAL is a 1-D array of length CSC_COL(N_ROWS)
   ;
   ; CASE 2: Two (2) positional arguments:
   ;  In this case, the coefficient matrix is supplied in coordinate form
   ;  by the first positional argument.  The second positional argument
   ;  is the number of rows.
   ;  The following checks are made:
   ;  - A must be an array of structures of type IMSL_D_SP_ELEM.
   ;  - The keywords CSC_COL CSC_ROW_IND CSC_VAL are not allowed.
   ;  - The row/column indices are checked against teh size of the matrix
   ;    as defined by N_ROWS.  This check helps avoid the case when the size
   ;    N_ROWS is input too large, which can lead to a crash within CNL.
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;
   nargs = n_params()
   IF ((nargs NE 0) AND (nargs NE 2)) THEN $
     message, "Incorrect number of arguments."
   type_a = IMSL_0
   size_a = IMSL_SIZE(a)
   IF (nargs EQ 0) THEN BEGIN ;The CSC_* keywords should be used.
      IF ((KEYWORD_SET(CSC_COL) + KEYWORD_SET(CSC_ROW_IND) + KEYWORD_SET(CSC_VAL)) $
          NE 3) THEN BEGIN
         message, "The keywords CSC_COL, CSC_ROW_IND, AND CSC_VAL are required " + $
           "when only one positional argument is passed."
      END ELSE BEGIN
         size_csc_col     = IMSL_SIZE(csc_col)
         sz_csc_row_ind   = IMSL_SIZE(csc_row_ind)
         size_csc_val     = IMSL_SIZE(csc_val)
         n_rows_cvt = IMSL_N_ELEMENTS(csc_col) - 1

         IF ((size_csc_col(0) NE 1) OR (N_ELEMENTS(csc_col) NE (n_rows_cvt + 1))) THEN $
           message, "CSC_COL must be a 1-D array OF length N_ELEMENTS(B)+1."
         csc_col_cvt = IMSL_LONG(csc_col)

         IF ((sz_csc_row_ind(0) NE 1) OR $
             (N_ELEMENTS(csc_row_ind) NE csc_col_cvt(n_rows_cvt))) THEN $
           message, "CSC_ROW_IND must be a 1-D array OF length CSC_COL(N_ELEMENTS(B))."
         csc_row_cvt = IMSL_LONG(csc_row_ind)

         IF ((size_csc_val(0) NE 1) OR $
             (N_ELEMENTS(csc_val) NE csc_col_cvt(n_rows_cvt))) THEN $
           message, "CSC_VAL must be a 1-D array OF length CSC_COL(N_ELEMENTS(B))."
      END
      nz = IMSL_N_ELEMENTS(csc_row_ind)
      a = 0.0 ; Just a place holder.
   END ELSE BEGIN ;  Two positional arguments ==> A supplied in coordinate form.
      n_rows_cvt = IMSL_LONG(n_rows(0))
      ; Test A to see if its in coordinate format.
      imsl_sp_acrd, a, nz, type_a
      IF ((KEYWORD_SET(CSC_COL) + KEYWORD_SET(CSC_ROW_IND) +  $
           KEYWORD_SET(CSC_VAL)) GT 0) THEN $
        message, "The keywords CSC_COL, CSC_ROW_IND, AND CSC_VAL are not allowed " + $
        "when two  positional arguments are passed."
      index_max   =   MAX(a.row)   >   MAX(a.col)
      IF (index_max LT n_rows_cvt-1) THEN $
        MESSAGE, "The largest row/column index does not correspond to the size of the matrix."
   END
   ;
   ; Decide on what precision to use.
   ; Note that the datatype is set based on the type of A.  This way,
   ; we will not be converting the type of A, and we will avoid creating
   ; a copy of A (which is potentially very large).
   IF (nargs EQ 0) THEN BEGIN
      type = (size_csc_val(N_ELEMENTS(size_csc_val)-2)) > TYP_FLOAT
   END ELSE BEGIN
      type  =  type_a
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Output LONG keyword(s)
   ; Always get n_nonzero, even if caller didn't request it since it is
   ; used inside the C interface code.
   n_nonzero_spc = IMSL_0
   ;
   ; Input LONG argument(s)
   IF (KEYWORD_SET(multifrontal) EQ TRUE) THEN $
     multifrontal_cvt = IMSL_1
   ;
   ; Double arguments and keywords
   ;
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      IF (ARG_PRESENT(sm_diag) EQ TRUE) THEN sm_diag_spc = double(0.0)
      IF (ARG_PRESENT(lg_diag) EQ TRUE) THEN lg_diag_spc  =  DOUBLE(0.0)
   END ELSE BEGIN
      IF (ARG_PRESENT(sm_diag) EQ TRUE) THEN sm_diag_spc = float(0.0)
      IF (ARG_PRESENT(lg_diag) EQ TRUE) THEN lg_diag_spc  =  float(0.0)
   END


   IF (nargs EQ 0) THEN csc_val_cvt = imsl_cvt_arr(csc_val, type)
   err_status = 0L
   MATHSTAT_211, LONG(type), err_status, a, n_rows_cvt, nz, $
                  n_nonzero_spc, $
                  xnzsub, $
                  xlnz, $
                  perm, $
                  diag, $
                  nzsub, $
                  alnz, $
                  sm_diag_spc, $
                  lg_diag_spc, $
                  multifrontal_cvt, $
                  maxsub, $
                  maxlnz, $
                  invp, $
                  multifrontal_space, $
                  csc_col_cvt, $
                  csc_row_cvt, $
                  csc_val_cvt


   IF (ARG_PRESENT(n_nonzero) EQ TRUE) THEN n_nonzero = n_nonzero_spc
   IF (ARG_PRESENT(sm_diag) EQ TRUE) THEN sm_diag = sm_diag_spc
   IF (ARG_PRESENT(lg_diag) EQ TRUE) THEN lg_diag = lg_diag_spc

   RETURN,  {nzsub:nzsub, $
              xnzsub:xnzsub, $
              xlnz:xlnz, $
              alnz:alnz, $
              perm:perm, $
              diag:diag}

END
