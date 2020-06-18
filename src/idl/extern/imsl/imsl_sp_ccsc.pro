; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_ccsc.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_Sp_ccsc, csc_row_ind, csc_col, csc_val, n_rows
   ;
   ; Procedure to check that the input variables make up
   ; a sparse matrix stored in CSC format.
   ;
   ;
@imsl_init.pro
   ON_ERROR, on_err_action

   size_csc_col     = IMSL_SIZE(csc_col)
   sz_csc_row_ind   = IMSL_SIZE(csc_row_ind)
   size_csc_val     = IMSL_SIZE(csc_val)

   IF ((size_csc_col(0) NE 1) OR (N_ELEMENTS(csc_col) NE (n_rows + 1))) THEN $
     message, 'CSC_COL must be a 1-D array OF length N_ELEMENTS(B)+1.'
   csc_col_cvt = IMSL_LONG(csc_col)

   IF ((sz_csc_row_ind(0) NE 1) OR $
       (N_ELEMENTS(csc_row_ind) NE csc_col_cvt(n_rows))) THEN $
     message, 'CSC_ROW_IND must be a 1-D array OF length CSC_COL(N_ELEMENTS(B)).'
   csc_row_cvt = IMSL_LONG(csc_row_ind)

   IF ((size_csc_val(0) NE 1) OR $
       (N_ELEMENTS(csc_val) NE csc_col_cvt(n_rows))) THEN $
     message, 'CSC_VAL must be a 1-D array of length CSC_COL(N_ELEMENTS(B)).'

END
