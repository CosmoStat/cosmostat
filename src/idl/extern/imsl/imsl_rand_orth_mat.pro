; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_rand_orth_mat.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_rand_orth_mat, n, $                ;INPUT Scalar LONG
                     eigenvalues=eigenvalues, $  ;INPUT 1-D array: floating point
                     a_matrix=a_matrix, $        ;INPUT 2-D array: floating point
                     double=double               ;INPUT Scalar ON/OFF flag
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - N greater than 1.
   ;   - EIGENVALUES must be a 1_D array of length N.
   ;   - If A_MATRTIX supplied:
   ;       - EIGENVALUES must also be supplied.
   ;       - A must be N by N.
   ;
   nargs = n_params()
   IF (nargs LT 1) THEN MESSAGE,  "Incorrect number of arguments."
   n_cvt = IMSL_LONG(n(0))
   IF (n_cvt LE 1) THEN MESSAGE,  "N must be greater than 1."
   IF (KEYWORD_SET(EIGENVALUES) EQ TRUE) THEN BEGIN
     size_eigs = IMSL_LONG(size(eigenvalues))
     IF (size_eigs(0) NE 1) THEN message, 'EIGENVALUES must be a 1-D array.'
     if (N_ELEMENTS(eigenvalues) NE n_cvt) then MESSAGE, "EIGENVALUES is not the correct size."
   END
   IF (KEYWORD_SET(A_MATRIX) EQ TRUE) THEN BEGIN
     IF (NOT KEYWORD_SET(EIGENVALUES)) THEN MESSAGE, $
       "EIGENVALUES must be provided if A_MATRIX is set."  
     size_a = IMSL_LONG(size(a_matrix))
     IF (size_a(0) NE 2) THEN message, 'A_MATRIX must be a 2-D array.'
     if ((size_a(1) NE n_cvt) OR (size_a(2) NE n_cvt)) then MESSAGE, "A_MATRIX is not the correct size."
   END

   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      IF (KEYWORD_SET(EIGENVALUES) EQ TRUE) THEN eigs_cvt = DOUBLE(EIGENVALUES)
      IF (KEYWORD_SET(A_MATRIX) EQ TRUE) THEN a_cvt = DOUBLE(TRANSPOSE(a_matrix))
      result = DBLARR(n_cvt, n_cvt)
   END ELSE BEGIN
      ; Input
      ; Input
      IF (KEYWORD_SET(EIGENVALUES) EQ TRUE) THEN eigs_cvt = FLOAT(EIGENVALUES)
      IF (KEYWORD_SET(A_MATRIX) EQ TRUE) THEN a_cvt = FLOAT(TRANSPOSE(a_matrix))
      result = FLTARR(n_cvt, n_cvt)
   END
   ;
   ; Call the system function.
   ;
   ne_spc = IMSL_0
   err_status = 0L
   MATHSTAT_306,  type,  err_status,  $
                                      n_cvt,   $
                                      eigs_cvt,   $
                                      a_cvt,  $
                                      result
   ;
   ; Return.
   ;
   RETURN, TRANSPOSE(result)
END

                   
                   
                   

  
      

  
