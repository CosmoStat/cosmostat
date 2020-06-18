; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_lufac.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_sp_lufac, a, $                           ;INPUT Struct
                   n_rows, $                      ;INPUT Scalar LONG
                   transpose=transpose, $         ;INPUT Scalar ON/OFF flag
                   pivoting=pivoting, $           ;INPUT Scalar LONG
                   n_search_rows=n_search_rows, $ ;INPUT Scalar LONG
                   iter_refine=iter_refine, $     ;INPUT Scalar ON/OFF flag
                   tol_drop=tol_drop, $           ;INPUT Scalar floating point
                   stability = stability, $       ;INPUT Scalar floating point
                   gwth_lim=gwth_lim, $           ;INPUT Scalar floating point
                   memory_block=memory_block, $   ;INPUT Scalar LONG
                   hybrid_density=hybrid_density,$;INPUT Scalar floating point
                   hybrid_order=hybrid_order, $   ;INPUT Scalar LONG
                   csc_col=csc_col, $             ;INPUT 1-D array: LONG
                   csc_row_ind=csc_row_ind, $     ;INPUT 1-D array: LONG
                   csc_val=csc_val, $             ;INPUT 1-D array: floating point
                   condition=condition, $         ;OUTPUT Scalar floating point
                   gwth_factor=gwth_factor, $     ;OUTPUT Scalar floating point
                   smallest_pvt=smallest_pvt, $   ;OUTPUT Scalar floating point
                   n_nonzero=n_nonzero            ;OUTPUT Scalar LONG

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
   ;
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;  - HYBRID_DENSITY and HYBRID_ORDER must be used together.
   ;  - STABILITY must be GE 1.0.
   ;
   nargs = n_params()
   IF ((nargs NE 0) AND (nargs NE 2)) THEN $
     message, "Incorrect number of arguments."
   type_a = 0L
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

   IF ((KEYWORD_SET(hybrid_density) + KEYWORD_SET(hybrid_order)) EQ 1) THEN $
     message, "The keywords HYBRID_DENSITY and HYBRID_ORDER must be used together."
   IF (KEYWORD_SET(stability)) THEN $
     IF (stability(0) LT 1.0) THEN $
     message, "The stability factor must be set greater than or equal to 1.0."
   ;
   ; Decide on what precision to use.
   ; Note that the datatype is set based on the type of A.  This way,
   ; we will not be converting the type of A, and we will avoid creating
   ; a copy of A (which is potentially very large).
   IF (nargs EQ 0) THEN BEGIN
      type = LONG(size_csc_val(N_ELEMENTS(size_csc_val)-2)) > TYP_FLOAT
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
   IF (KEYWORD_SET(transpose) EQ TRUE) THEN $
     transpose_cvt = IMSL_1
   IF (KEYWORD_SET(iter_refine) EQ TRUE) THEN $
     iter_refine_cvt = IMSL_1
   IF (KEYWORD_SET(pivoting) EQ TRUE) THEN BEGIN
      pivoting_cvt = IMSL_LONG(pivoting(0))
      IF ((pivoting_cvt LT 1) OR (pivoting_cvt GT 3)) THEN pivoting_cvt = 1
   END
   IF (KEYWORD_SET(n_search_rows) EQ TRUE) THEN $
     n_srch_rows_cvt = IMSL_LONG(n_search_rows(0))
   IF (KEYWORD_SET(memory_block) EQ TRUE) THEN $
     mem_block_cvt = IMSL_LONG(memory_block(0))
   IF (KEYWORD_SET(hybrid_order) EQ TRUE) THEN $
     hyb_order_cvt = IMSL_LONG(hybrid_order(0))
   ;
   ;
   ; Floating point arguments and keywords
   cmplx_scale = IMSL_1
   IF ((type EQ TYP_DCMPLX) OR (type EQ  TYP_COMPLEX)) THEN cmplx_scale = IMSL_2
   IF ((type EQ TYP_DOUBLE) OR (type EQ TYP_DCMPLX)) THEN BEGIN
      ;
      ; Input
      IF (KEYWORD_SET(tol_drop) EQ TRUE) THEN $
        tol_drop_cvt = (double(tol_drop))(0)
      IF (KEYWORD_SET(stability) EQ TRUE) THEN $
        stability_cvt = ((double(stability))(0))
      IF (KEYWORD_SET(gwth_lim) EQ TRUE) THEN $
        gwth_lim_cvt = (double(gwth_lim))(0)
      IF (KEYWORD_SET(hybrid_density) EQ TRUE) THEN $
        hyb_dens_cvt = DOUBLE(hybrid_density(0))
      ;
      ; Output
      IF (ARG_PRESENT(condition) EQ TRUE) THEN $
        condition_spc = double(0.0)
      IF (ARG_PRESENT(gwth_factor) EQ TRUE) THEN $
        gwth_factor_spc = double(0.0)
      IF (ARG_PRESENT(smallest_pvt) EQ TRUE) THEN $
        smallest_pvt_spc  =  DOUBLE(0.0)
   END ELSE BEGIN
      ;
      ; Input
      IF (KEYWORD_SET(tol_drop) EQ TRUE) THEN $
        tol_drop_cvt = (float(tol_drop))(0)
      IF (KEYWORD_SET(stability) EQ TRUE) THEN $
        stability_cvt = ((float(stability))(0))
      IF (KEYWORD_SET(gwth_lim) EQ TRUE) THEN $
        gwth_lim_cvt = (float(gwth_lim))(0)
      IF (KEYWORD_SET(hybrid_density) EQ TRUE) THEN $
        hyb_dens_cvt = FLOAT(hybrid_density(0))
      ;
      ; Output
      IF (ARG_PRESENT(condition) EQ TRUE) THEN $
        condition_spc = float(0.0)
      IF (ARG_PRESENT(gwth_factor) EQ TRUE) THEN $
        gwth_factor_spc = float(0.0)
      IF (ARG_PRESENT(smallest_pvt) EQ TRUE) THEN $
        smallest_pvt_spc  =  FLOAT(0.0)
   END

   IF (nargs EQ 0) THEN csc_val_cvt = imsl_cvt_arr(csc_val, type)
   ;
   ; Factor arguments, always present.
   fac_coord_spc = IMSL_1
   fac_coord_rpvt = IMSL_LONARR(n_rows_cvt)
   fac_coord_cpvt = IMSL_LONARR(n_rows_cvt)
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_207,  type, err_status, a, n_rows_cvt, nz, $
                   transpose_cvt, $
                   pivoting_cvt, $
                   n_srch_rows_cvt, $
                   iter_refine_cvt, $
                   tol_drop_cvt, $
                   stability_cvt, $
                   gwth_lim_cvt, $
                   memory_block_cvt, $
                   condition_spc, $
                   gwth_factor_spc, $
                   smallest_pvt_spc, $
                   n_nonzero_spc, $
                   fac_coord_spc, $
                   fac_coord_rpvt, $
                   fac_coord_cpvt, $
                   hyb_dens_cvt, $
                   hyb_order_cvt, $
                   csc_col_cvt, $
                   csc_row_cvt, $
                   csc_val_cvt

   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(condition) EQ TRUE) THEN $
     condition = condition_spc
   IF (ARG_PRESENT(gwth_factor) EQ TRUE) THEN $
     gwth_factor = gwth_factor_spc
   IF (ARG_PRESENT(smallest_pvt) EQ TRUE) THEN $
     smallest_pvt = smallest_pvt_spc
   IF (ARG_PRESENT(n_nonzero) EQ TRUE) THEN n_nonzero = n_nonzero_spc

   result  = { lu:fac_coord_spc, row_pvts:fac_coord_rpvt, col_pvts:fac_coord_cpvt}
   RETURN, result
END









