; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_qrsol.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_qrsol, b, $                     ;INPUT 1-D array: floating point 
                a, $                       ;INPUT 2-D array: floating point 
                pivot=pivot, $             ;INPUT/OUTPUT 1-D array: LONG
                auxqr=auxqr, $             ;INPUT 1-D array: floating point
                qr=qr, $                   ;INPUT 2-D array: floating point
                tolerance=tolerance, $     ;INPUT Scalar floating point
                double=double, $           ;INPUT Scalar ON/OFF flag
                residual=residual, $       ;OUTPUT 1-D array: floating point
                basis=basis                ;OUTPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  The two basic uses of this function are divided between the
   ;  cases when there are either 1 or 2 position arguments.
   ;  Case 1: Only one positional argument.
   ;          In this case the RHS is supplied through the positional
   ;          argument, and the factored system must be supplied through
   ;          the keywords . In this case the following checks 
   ;          are performed.
   ;          - B:
   ;              - Must be a 1-D array. 
   ;              - set (m = n_ELEMENTS(B))
   ;          - AUXQR:
   ;              - Must be supplied. 
   ;              - set (n = length of this array).
   ;          - QR:
   ;              - Must be supplied. 
   ;              - Must be of size (m x n)
   ;          - PIVOT:
   ;              - Must be supplied. 
   ;              - Must be of length n
   ;          - BASIS:
   ;              - Must supplied.
   ;          
   ;  Case 2: Two positional arguments.
   ;          In this case, both the RHS and the original system have
   ;          have been supplied. In this case the following checks are 
   ;          performed.
   ;          - A must be a 2-D array of order (m x n)
   ;            m and n are set by this size of this argument.
   ;          - B must be a 1-D array of length m
   ;          - QR can't be present.
   ;          - AUXQR can't be present.
   ;          
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2))  THEN $
     message, 'Incorrect number of arguments.'
   IF (nargs EQ 1) THEN BEGIN
      size_b = IMSL_LONG(size(b))
      IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
      m = IMSL_LONG(size_b(1))
      IF (NOT KEYWORD_SET(auxqr)) THEN $
        message, 'The keyword AUXQR is required for the specified usage of QRSOL.'
      size_auxqr = IMSL_SIZE(auxqr)
      IF (size_auxqr(0) NE 1) THEN message, 'AUXQR must be a 1-D array.'
      n = IMSL_LONG(size_auxqr(1))
      IF (NOT KEYWORD_SET(qr)) THEN $
        message, 'The keyword QR is required for the specified usage of QRSOL.'
      size_qr = IMSL_SIZE(qr)
      IF (size_qr(0) NE 2) THEN message, 'AUXQR must be a 2-D array.'
      IF ((size_qr(1) NE m) OR (size_qr(2) NE n)) THEN $
        message, 'QR is NOT the correct size.'
      IF (NOT KEYWORD_SET(pivot)) THEN $
        message, 'The keyword PIVOT is required for the specified usage of QRSOL.'
      size_pivot = IMSL_SIZE(pivot)
      IF (size_pivot(0) NE 1) THEN message, 'PIVOT must be a 1-D array.'
      IF (size_pivot(1) NE n) THEN message, 'PIVOT is NOT the correct length.'
      ; BASIS = 0 is valid, so check it carefully.
      IF ((arg_present(basis) + KEYWORD_SET(basis)) EQ 0) THEN  $
        message, 'The keyword BASIS is required for the specified usage of QRSOL.'
   END ELSE BEGIN
      size_a = IMSL_SIZE(a)
      IF (size_a(0) NE 2) THEN message, 'A must be a 2-D array.'
      m = IMSL_LONG(size_a(1))
      n = IMSL_LONG(size_a(2))
      size_b = IMSL_LONG(size(b))
      IF (size_b(0) NE 1) THEN message, 'B must be a 1-D array.'
      IF (size_b(1) NE m) THEN message, 'B is NOT the correct length.'
      IF (KEYWORD_SET(qr)) THEN $
        message, 'QR is not valid for this usage of QRSOL.'
      IF (KEYWORD_SET(auxqr)) THEN $
        message, 'AUXQR is not valid for this usage of QRSOL.'
   END
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (nargs EQ 2) THEN $
     IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG
   ; Always send space for Pivot.
   IF (KEYWORD_SET(pivot)) THEN pivot_cvt = IMSL_LONG(pivot)
   ; BASIS = 0 is valid, so convert it carefully.
   CASE (arg_present(basis) + KEYWORD_SET(basis)) OF
      0:basis_cvt = IMSL_LONG((m < n))
      1:basis_cvt = IMSL_LONG((m < n))
      2:basis_cvt = IMSL_LONG(basis(0))
   END
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = dblarr(n)
      b_cvt = double(b)
      IF (nargs EQ 2) THEN a_cvt = double(transpose(a))
      IF (KEYWORD_SET(auxqr)) THEN auxqr_cvt = double(auxqr)
      IF (KEYWORD_SET(qr)) THEN qr_cvt = double(transpose(qr))
      tmp = imsl_machine(/double)
      tolerance_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
      IF (ARG_PRESENT(residual)) THEN residual_spc = dblarr(m)
   END ELSE BEGIN
      result = fltarr(n)
      b_cvt = float(b)
      IF (nargs EQ 2) THEN a_cvt = float(transpose(a))
      IF (KEYWORD_SET(auxqr)) THEN auxqr_cvt = float(auxqr)
      IF (KEYWORD_SET(qr)) THEN qr_cvt = float(transpose(qr))
      tmp = imsl_machine(/float)
      tolerance_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
      IF (ARG_PRESENT(residual)) THEN residual_spc = fltarr(m)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_177, type, err_status, b_cvt, a_cvt, m, n, $
                              pivot_cvt, $
                              auxqr_cvt, $
                              qr_cvt, $
                              tolerance_cvt, $
                              residual_spc, $
                              basis_cvt, $
                              result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(residual)) THEN residual = residual_spc
   IF (ARG_PRESENT(basis)) THEN basis = basis_cvt
   ;
   ; Return
   RETURN, result
END
