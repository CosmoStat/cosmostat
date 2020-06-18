; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_acrd.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_Sp_acrd, a, nz, type_a
   ;
   ; Procedure to check that the input variable, A, is
   ; a sparse matrix stored in coordinate format.
   ;
   ; Args:
   ;  a       The sparse matrix to check out. (Input)
   ;  nz      The number of nonzeros in the sparse matrix. (output)
   ;  type_a  The data type of A. (Output)
   ;
@imsl_init.pro
   ON_ERROR, on_err_action

   size_a = IMSL_SIZE(a)
   size_a_elem = IMSL_SIZE(a(0))
   IF ((size_a(0) NE 1) OR $
       (size_a_elem(N_ELEMENTS(size_a_elem)-2) NE 8)) THEN $
     message, "A must be an array of structures for a sparse matrix stored in coordinate format."

   a_sName = tag_names(a, /structure)
   IF ((a_sName NE 'IMSL_F_SP_ELEM') AND  $
       (a_sName NE 'IMSL_D_SP_ELEM') AND  $
       (a_sName NE 'IMSL_C_SP_ELEM') AND  $
       (a_sName NE 'IMSL_Z_SP_ELEM') ) THEN $
     MESSAGE,  "A must be an array of structures for a sparse matrix stored in coordinate format."

   size_val = IMSL_SIZE(a.val)
   type_a = LONG(size_val(N_ELEMENTS(size_val)-2))
   nz = IMSL_N_ELEMENTS(a)

END


