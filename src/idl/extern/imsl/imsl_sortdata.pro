; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sortdata.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_sortdata, x, $                         ;INPUT 1-D or 2-D array: floating point
                   n_keys, $                    ;INPUT Scalar LONG
                   ascending=ascending, $       ;INPUT Scalar ON/OFF flag
                   descending=descending, $     ;INPUT Scalar ON/OFF flag
                   double=double, $             ;INPUT Scalar ON/OFF flag
                   frequencies=frequencies, $   ;INPUT 1-D array: floating point
                   indices_keys=indices_keys, $ ;INPUT 1-D array: LONG
                   list_cells=list_cells, $     ;OUTPUT 2-D array: floating point
                   n_cells=n_cells, $           ;OUTPUT 1-D array: LONG
                   n_list_cells=n_list_cells, $ ;OUTPUT 2-D array: LONG
                   permutation=permutation, $   ;OUTPUT 1-D array: LONG
                   table_bal=table_bal, $       ;OUTPUT 1-D array: floating point
                   table_n=table_n, $           ;OUTPUT 1-D array: LONG
                   table_unbal=table_unbal, $   ;OUTPUT 1-D array: floating point
                   table_values=table_values    ;OUTPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - Check that the input array, X,  is 1-D or 2-D
   ;  - ASCENDING and DESCENDING are mutually exclusive
   ;  - If FREQUENCIES, it is of length x->value.arr->dim[0]
   ;  - If INDICES_KEYS, it is of length *argv[1]
   ;  - The following three keywords must appear together
   ;        TABLE_N
   ;        TABLE_VALUES
   ;        TABLE_BAL
   ;  - The following three keywords must appear together
   ;        N_LIST_CELLS
   ;        LIST_CELLS
   ;        TABLE_UNBAL
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF ((size_x(0) LT 1) or (size_x(0) GT 2)) THEN BEGIN 
      message, "The input array, X, must be 1-D or 2-D"
   END
   nobs = size_x(1) ;Define nobs
   IF (size_x(0) EQ 1) THEN nvar = IMSL_1 ELSE nvar = IMSL_LONG(size_x(2)) ; Define nvar
   n_keys_cvt = (IMSL_LONG(n_keys))(0)

   IF (KEYWORD_SET(ascending) AND KEYWORD_SET(descending)) THEN $
     message, "ASCENDING AND DESCENDING are mutually exclusive."
   
   IF (KEYWORD_SET(frequencies) EQ true) THEN BEGIN
      size_freq = IMSL_SIZE(frequencies)
      IF (size_freq(0) NE 1) THEN $
      message, "FREQUENCIES must be a 1-D array."
      IF (N_ELEMENTS(frequencies) NE nobs) THEN $
        message, "FREQUENCIES is not the correct length."
   END
   
   IF (KEYWORD_SET(indices_keys) EQ true) THEN BEGIN
      size_ind = IMSL_SIZE(indices_keys)
      IF (size_ind(0) NE 1) THEN $
      message, "INDICES_KEYS must be a 1-D array."
      IF (N_ELEMENTS(indices_keys) NE n_keys_cvt) THEN $
        message, "INDICES_KEYS is not the correct length."
   END

   l_tmp = 0
   IF (ARG_PRESENT(table_n) EQ true) THEN l_tmp = l_tmp + 1
   IF (ARG_PRESENT(table_values) EQ true) THEN l_tmp = l_tmp + 1
   IF (ARG_PRESENT(table_bal) EQ true) THEN l_tmp = l_tmp + 1
   IF ((l_tmp NE 0) AND (l_tmp NE 3)) THEN $
     message, "The keywords TABLE_N, TABLE_VALUES, TABLE_BAL must be used together."

   l_tmp = 0
   IF (ARG_PRESENT(n_list_cells) EQ true) THEN l_tmp = l_tmp + 1
   IF (ARG_PRESENT(list_cells) EQ true) THEN l_tmp = l_tmp + 1
   IF (ARG_PRESENT(table_unbal) EQ true) THEN l_tmp = l_tmp + 1
   IF ((l_tmp NE 0) AND (l_tmp NE 3)) THEN $
     message, "The keywords N_LIST_CELLS, LIST_CELLS, and TABLE_UNBAL must be used together."
   
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Here are some vaules we may need later to define array sizes.
   ;
   max_n_cells = IMSL_LONG(nobs^n_keys_cvt)
   max_value = n_keys_cvt*nobs
   max_table =  IMSL_LONG(nobs^n_keys)
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(ascending) EQ TRUE) THEN $
     ascending_cvt = IMSL_1
   IF (KEYWORD_SET(descending) EQ TRUE) THEN $
     descending_cvt = IMSL_1
   IF (KEYWORD_SET(indices_keys) EQ TRUE) THEN $
     indices_kys_cvt = IMSL_LONG(indices_keys)
   ; 
   ; Output LONG keyword(s)
   IF ARG_PRESENT(permutation) THEN $
     permutation_spc = IMSL_LONARR(nobs)
   IF ARG_PRESENT(n_cells) THEN BEGIN
      n_cells_spc = IMSL_LONARR(nobs)
      lngth_n_cells = IMSL_0
   END
   
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      l_result = dblarr(nvar, nobs)
      ; 
      ; Input
      IF (nvar GT 1) THEN x_cvt = double(transpose(x)) ELSE x_cvt = double(x)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = double(frequencies)
      ; 
      ; Output 
      IF (ARG_PRESENT(list_cells) EQ TRUE) THEN BEGIN
         n_list_cls_spc = IMSL_0
         list_cls_spc = dblarr(n_keys_cvt, max_n_cells)
         table_unbal_spc = dblarr(max_n_cells)
      END
      IF ARG_PRESENT(table_bal) THEN BEGIN
         table_n_spc = IMSL_LONARR(n_keys_cvt)
         table_val_spc = dblarr(max_value)
         table_bal_spc = dblarr(max_table)
      END
   END ELSE BEGIN
      ; Result vector.
      l_result = fltarr(nvar, nobs)
      ; 
      ; Input 
      IF (nvar GT 1) THEN x_cvt = float(transpose(x)) ELSE x_cvt = float(x)
      IF (KEYWORD_SET(frequencies) EQ TRUE) THEN $
        frequencies_cvt = float(frequencies)
      ; 
      ; Output 
      IF (ARG_PRESENT(list_cells) EQ TRUE) THEN BEGIN
         n_list_cls_spc = IMSL_0
         list_cls_spc = fltarr(n_keys_cvt, max_n_cells)
         table_unbal_spc = fltarr(max_n_cells)
      END
      IF ARG_PRESENT(table_bal) THEN BEGIN
         table_n_spc = IMSL_LONARR(n_keys_cvt)
         table_val_spc = fltarr(max_value)
         table_bal_spc = fltarr(max_table)
      END
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_188, type, err_status, x_cvt, n_keys_cvt, nobs, nvar, $
                 ascending_cvt, $
                 descending_cvt, $
                 frequencies_cvt, $
                 indices_kys_cvt, $
                 list_cls_spc, $
                 n_list_cls_spc, $
                 table_unbal_spc, $
                 permutation_spc, $
                 table_bal_spc, $
                 table_n_spc, $
                 table_val_spc, $
                 n_cells_spc, $
                 lngth_n_cells
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(list_cells) EQ TRUE) THEN BEGIN
      n_list_cells = n_list_cls_spc
      list_cls_spc = transpose(list_cls_spc)
      list_cells = list_cls_spc(0:n_list_cells-1, 0:n_keys_cvt-1)
      table_unbal = table_unbal_spc(0:n_list_cells-1)
   IF (ARG_PRESENT(n_cells) EQ TRUE) THEN $
     n_cells = n_cells_spc(0:((lngth_n_cells-1) > 0))
   END
   IF (ARG_PRESENT(permutation) EQ TRUE) THEN $
     permutation = permutation_spc
   IF ARG_PRESENT(table_bal) THEN BEGIN
      table_n = table_n_spc
      table_values = table_val_spc(0:IMSL_LONG(total(table_n_spc))-1)
      l_tmp = IMSL_1
      FOR i = 0, (N_ELEMENTS(table_n)-1) DO l_tmp = l_tmp*table_n(i)
      table_bal = table_bal_spc(0:l_tmp-1)
   END
   ;
   ; Return.
   ;
   RETURN, transpose(x_cvt)
END
   

