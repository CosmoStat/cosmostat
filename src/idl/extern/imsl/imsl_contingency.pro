; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_contingency.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_contingency, x, $                             ;INPUT 2-D array: floating point
                      double=double, $                 ;INPUT Scalar ON/OFF flag
                      Chi_sq_test=chi_sq_test, $       ;OUTPUT 1-D array: floating point
                      Lrt=lrt, $                       ;OUTPUT 1-D array: floating point
                      Chi_sq_stats=chi_sq_stats, $     ;OUTPUT 1-D array: floating point
                      Expected=expected, $             ;OUTPUT 2-D array: floating point
                      Chi_sq_contrib=chi_sq_contrib, $ ;OUTPUT 2-D array: floating point
                      Table_stats=table_stats          ;OUTPUT 2-D array: floating point
                      
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; Make sure the input argument is a 2-D simple array.
   ;
   nargs = n_params()
   IF (nargs EQ 1) THEN BEGIN
      size_x = IMSL_SIZE(x)
      IF (size_x(0) NE 2) THEN BEGIN
         message, "X must be a 2-D array."
      END
   END ELSE message, "Incorrect number of arguments."
      
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   n_rows = size_x(1)
   n_columns = size_x(2)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      ; 
      ; Input 
      x_cvt = double(transpose(x))
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test_spc = dblarr(3)
      IF (ARG_PRESENT(Lrt) EQ TRUE) THEN $
        Lrt_spc = dblarr(3)
      IF (ARG_PRESENT(expected) EQ TRUE) THEN $
        expected_spc = dblarr(n_columns+1, n_rows+1)
      IF (ARG_PRESENT(chi_sq_contrib) EQ TRUE) THEN $
        chi_sq_ctrb_spc = dblarr(n_columns+1, n_rows+1)
      IF (ARG_PRESENT(chi_sq_stats) EQ TRUE) THEN $
        chi_sq_sts_spc = dblarr(5)
      IF (ARG_PRESENT(table_stats) EQ TRUE) THEN $
        table_sts_spc = dblarr(5, 23)
   END ELSE BEGIN
      result = float(0.0)
      ; 
      ; Input 
      x_cvt = float(transpose(x))
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
        Chi_sq_test_spc = fltarr(3)
      IF (ARG_PRESENT(Lrt) EQ TRUE) THEN $
        Lrt_spc = fltarr(3)
      IF (ARG_PRESENT(expected) EQ TRUE) THEN $
        expected_spc = fltarr(n_columns+1, n_rows+1)
      IF (ARG_PRESENT(chi_sq_contrib) EQ TRUE) THEN $
        chi_sq_ctrb_spc = fltarr(n_columns+1, n_rows+1)
      IF (ARG_PRESENT(chi_sq_stats) EQ TRUE) THEN $
        chi_sq_sts_spc = fltarr(5)
      IF (ARG_PRESENT(table_stats) EQ TRUE) THEN $
        table_sts_spc = fltarr(5, 23)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_125, type, err_status, x_cvt, n_rows, n_columns, $
                            chi_sq_test_spc, $
                            lrt_spc, $
                            expected_spc, $
                            chi_sq_ctrb_spc, $
                            chi_sq_sts_spc, $
                            table_sts_spc, $
                            result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(chi_sq_contrib) EQ TRUE) THEN $
     chi_sq_contrib = transpose(chi_sq_ctrb_spc)  ;NOTE: Transpose 2-D array
   IF (ARG_PRESENT(chi_sq_stats) EQ TRUE) THEN $
     chi_sq_stats = chi_sq_sts_spc
   IF (ARG_PRESENT(Chi_sq_test) EQ TRUE) THEN $
     Chi_sq_test = Chi_sq_test_spc
   IF (ARG_PRESENT(expected) EQ TRUE) THEN $
     expected = transpose(expected_spc)  ;NOTE: Transpose 2-D array
   IF (ARG_PRESENT(Lrt) EQ TRUE) THEN $
     Lrt = Lrt_spc
   IF (ARG_PRESENT(table_stats) EQ TRUE) THEN $
     table_stats = transpose(table_sts_spc)  ;NOTE: Transpose 2-D array
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
