; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_lusol.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
;
; Check_fac_lu is used by sp_lusol to test if the argument passed
; for the LU factorization is valid.  This argument is only passed
; in the case that the system has been factored by SP_LUFAC previously.
;
PRO imsl_Check_fac_lu, factor_coord, type_fac

@imsl_init.pro

   err_str = "FACTOR_COORD must be a structure containing the sparse LU factorization."

   ; See if it a structures.
   size_facoord = IMSL_SIZE(factor_coord)
   size_facoord_elem = IMSL_SIZE(factor_coord(0))
   IF ((size_facoord(0) NE 1) OR $
       (size_facoord_elem(N_ELEMENTS(size_facoord_elem)-2) NE 8)) THEN message, err_str

   ; Check out the tag names.
   tags = tag_names(factor_coord)
   IF (N_ELEMENTS(tags) NE 3) THEN message, err_str
   expected_tags = ['LU', 'ROW_PVTS', 'COL_PVTS']
   FOR i = 0, 2 DO $
     IF (tags(i) NE expected_tags(i)) THEN message, err_str

   ; Now find data type of the factorization, based upon factor_coord.lu.
   size_lu = IMSL_SIZE(factor_coord.lu)
   size_lu_elem = IMSL_SIZE(factor_coord.lu(0))

   ; The .lu element of the structure must also be a structure.
   IF ((size_lu(0) NE 1) OR $
       (size_lu_elem(N_ELEMENTS(size_lu_elem)-2) NE 8)) THEN message, err_str
   size_lu_val = IMSL_SIZE(factor_coord.lu.val)
   type_fac = (size_lu_val(N_ELEMENTS(size_lu_val)-2))

END

FUNCTION imsl_Sp_lusol, b, $                           ;INPUT 1-D array: floating point
                   a, $                           ;INPUT Struct
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
                   factor_coord=factor_coord, $   ;INPUT Struct
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
   ; CASE 1: One (1) positional argument:
   ;  In this case, the positional argument is the RHS, b.
   ;  If CSC_* are not provided, then FACTOR_COORD is required as input,
   ;  and we are doing a solve-only call.
   ;  The following checks are made:
   ;  - B must be a 1-D array of length N_ROWS.
   ;  If CSC_* are provided:
   ;  - CSC_COL is a 1-D array of length (N_ROWS + 1)
   ;  - CSC_ROW_IND is a 1-D array of length CSC_COL(N_ROWS)
   ;  - CSC_VAL is a 1-D array of length CSC_COL(N_ROWS)
   ;  If CSC are not provided
   ;  - FACTOR_COORD must be an array of structures of type IMSL_D_SP_ELEM.
   ;
   ; CASE 2: Two (2) positional arguments:
   ;  In this case, the positional arguments the RHS, b, and the sparse matrix, A,
   ;  stored as a 1-D array of type IMSL_D_SP_ELEM.
   ;  - B must be a 1-D array of length N_ROWS.
   ;  - A must be an array of structures of type IMSL_D_SP_ELEM.
   ;  - The row/column indices are checked against teh size of the matrix
   ;    as defined by B  This check helps avoid the case when B
   ;    is input too long, which can lead to a crash within CNL.
   ;
   ; GENERAL CHECKS:
   ; The following checks are made regardless of the number of positional
   ; arguments.
   ;  - HYBRID_DENSITY and HYBRID_ORDER must be used together.
   ;  - STABILITY must be GE 1.0.
   ;  - the following keywords are not allowed if FACTOR is supplied:
   ;       CONDITION, GWTH_FACTOR, N_NONZERO, SMALLEST_PVT
   ;
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2)) THEN $
     message, "Incorrect number of arguments."
   type_a = IMSL_0
   factored = IMSL_0
   size_b = IMSL_SIZE(b)
   IF (size_b(0) NE 1) THEN  $
     message, "B must be a 1-D array of length N_ROWS."
   n_rows = IMSL_LONG(size_b(3))
   IF (nargs EQ 1) THEN BEGIN
      tmp = KEYWORD_SET(CSC_COL) + KEYWORD_SET(CSC_ROW_IND) + KEYWORD_SET(CSC_VAL)
      IF ((tmp NE 3) AND (tmp NE 0)) THEN BEGIN
         message, "The keywords CSC_COL, CSC_ROW_IND, AND CSC_VAL " + $
           "must be used together."
      END ELSE BEGIN
         IF (tmp EQ 3) THEN BEGIN
            size_csc_col     = IMSL_SIZE(csc_col)
            sz_csc_row_ind   = IMSL_SIZE(csc_row_ind)
            size_csc_val     = IMSL_SIZE(csc_val)

            IF ((size_csc_col(0) NE 1) OR (N_ELEMENTS(csc_col) NE (n_rows + 1))) THEN $
              message, "CSC_COL must be a 1-D array OF length N_ELEMENTS(B)+1."
            csc_col_cvt = IMSL_LONG(csc_col)

            IF ((sz_csc_row_ind(0) NE 1) OR $
                (N_ELEMENTS(csc_row_ind) NE csc_col_cvt(n_rows))) THEN $
              message, "CSC_ROW_IND must be a 1-D array OF length CSC_COL(N_ELEMENTS(B))."
            csc_row_cvt = IMSL_LONG(csc_row_ind)

            IF ((size_csc_val(0) NE 1) OR $
                (N_ELEMENTS(csc_val) NE csc_col_cvt(n_rows))) THEN $
              message, "CSC_VAL must be a 1-D array OF length CSC_COL(N_ELEMENTS(B))."
            nz = IMSL_N_ELEMENTS(csc_row_ind)
         END ELSE BEGIN ; (Thus tmp eq 0, and we are doing a solve-only step)
            IF NOT KEYWORD_SET(factor_coord) THEN $
              message, "FACTOR_COORD is required as input if only B is supplied."
            imsl_Check_fac_lu, factor_coord, type_a
            nz = IMSL_N_ELEMENTS(factor_coord.lu)
            factored = IMSL_1
            solve_only  =  IMSL_1
            IF (ARG_PRESENT(condition)) THEN $
              MESSAGE,   "Keywords CONDITION and FACTOR_COORD cannot be used together."
            IF (ARG_PRESENT(GWTH_FACTOR)) THEN $
              MESSAGE,   "Keywords GWTH_FACTOR and FACTOR_COORD cannot be used together."
            IF (ARG_PRESENT(N_NONZERO)) THEN $
              MESSAGE,   "Keywords N_NONZERO and FACTOR_COORD cannot be used together."
            IF (ARG_PRESENT(smallest_pvt)) THEN $
              MESSAGE,   "Keywords SMALLEST_PVT and FACTOR_COORD cannot be used together."
         END
      END
      a = 0.0d0                 ;  Just a place holder.
   END ELSE BEGIN               ;  Two positional arguments.
      ; Test that A is in coordinate format.
      imsl_sp_acrd, a, nz, type_a
      IF ((KEYWORD_SET(CSC_COL) + KEYWORD_SET(CSC_ROW_IND) +  $
           KEYWORD_SET(CSC_VAL)) GT 0) THEN $
        message, "The keywords CSC_COL, CSC_ROW_IND, AND CSC_VAL are not allowed " + $
           "when two  positional argument are passed."
      index_max   =   MAX(a.row)   >   MAX(a.col)
      IF (index_max LT N_ELEMENTS(b)-1-1) THEN $
        MESSAGE, "The largest row/column index of A not correspond to the size of B."
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
   ; Always get n_nonzero, even if caller didn't request it since it is
   ; used inside the C interface code.
   n_nonzero_spc = IMSL_0
   ;
   ; Input LONG argument(s)
   IF (KEYWORD_SET(transpose) EQ TRUE) THEN $
     transpose_cvt = IMSL_1
   IF (KEYWORD_SET(iter_refine) EQ TRUE) THEN $
     iter_refine_cvt = IMSL_1
   IF (KEYWORD_SET(pivoting) EQ TRUE) THEN $
     pivoting_cvt = IMSL_LONG(pivoting(0))
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
      ; Result
      ;
      result  =  DBLARR(cmplx_scale*n_rows)
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
      ; Result
      ;
      result  =  FLTARR(cmplx_scale*n_rows)
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

   IF (nargs EQ 1) THEN $
     IF (factored EQ 0) THEN csc_val_cvt   =   imsl_cvt_arr(csc_val,  type)
   b_cvt    =    imsl_cvt_arr(b,   type)

   ;
   ; Call the system function.
   ;
   err_status = 0L
   ; Have to use two seperate calls here to avoid making copies of factor_coord
   ; in the case we are doing a solve-only step.

   IF (factored EQ 1) THEN $
     MATHSTAT_208, LONG(type),   err_status,   b_cvt,   a,   n_rows,   nz,   $
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
                   csc_val_cvt, $
                   solve_only, $
                   factor_coord.lu, $
                   factor_coord.row_pvts, $
                   factor_coord.col_pvts, $
                   result $
     ELSE $
     MATHSTAT_208, LONG(type), err_status, b_cvt, a, n_rows, nz, $
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
                   csc_val_cvt, $
                   solve_only, $
                   t,t,t, $ ; place holders for unused arguments.
                   result

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
   ;
   ; Return.
   ;
   IF (cmplx_scale EQ 2) $
     THEN RETURN, imsl_cvt_arr(result, /back) $
     ELSE RETURN, result

END









