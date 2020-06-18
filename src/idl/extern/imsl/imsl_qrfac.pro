; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_qrfac.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO      imsl_qrfac, a, $                     ;INPUT 2-D array: floating point 
                pivot, $                   ;INPUT/OUTPUT 1-D array: LONG
                auxqr, $                   ;OUTPUT 1-D array: floating point
                qr, $                      ;OUTPUT 2-D array: floating point
                tolerance=tolerance, $     ;INPUT Scalar floating point
                double=double, $           ;INPUT Scalar ON/OFF flag
                q=q, $                     ;OUTPUT 2-D array: floating point
                r=r, $                     ;OUTPUT 2-D array: floating point
                ap=ap, $                   ;OUTPUT 2-D array: floating point
                basis=basis                ;OUTPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  The three basic uses of this function are divided between the
   ;  cases when there are either 1, 2, or 4  position arguments.
   ;  ALL CASES:
   ;     In all case the system to be factored is supplied
   ;     through the positional argument.
   ;     In this case the following checks are performed.
   ;     - A:
   ;         - Must be a 2-D array. 
   ;         - set (m = n_lelements(A(*, 0)))
   ;         - set (n = n_lelements(A(0, *)))
   ;  Case 1: Only one positional argument.
   ;          In this case the system to be factored is supplied
   ;          through the positional argument.
   ;          - Check A as described above.
   ;  Case 2: Two positional arguments.
   ;          In this case the system to be factored is supplied in
   ;          the first positionl argument, and the second positional
   ;          argument is an input/output vector specifying the pivot
   ;          sequence.
   ;          In this case the following checks are performed.
   ;          - Check A as described above.
   ;          - Pivot
   ;              - If it is undefined, create a 1-D array of length n,
   ;                and set all elements to 0.
   ;              - Must be a 1-D array of length n
   ;  Case 3: Four positional arguments.
   ;          In this case the system to be factored is supplied in
   ;          the first positionl argument, and the second positional
   ;          argument is an input/output vector specifying the pivot
   ;          sequence, and the last two positional arguments are output
   ;          arguments which will contain the factored system.
   ;          In this case the following checks are performed.
   ;          - Check A as described above.
   ;          - Pivot
   ;              - If it is undefined, create a 1-D array of length n,
   ;                and set all elements to 0.
   ;              - Must be a 1-D array of length n
   ;          
   nargs = n_params()
   IF (((nargs NE 1) AND (nargs NE 2)) and  (nargs NE 4)) THEN $
     message, 'Incorrect number of arguments.'
   size_a = IMSL_SIZE(a)
   IF (size_a(0) NE 2) THEN message, 'A must be a 2-D array.'
   m = IMSL_LONG(size_a(1))
   n = IMSL_LONG(size_a(2))
   IF (nargs GT 1) THEN BEGIN 
      IF (N_ELEMENTS(pivot) EQ 0) THEN pivot = IMSL_LONARR(n)
      size_pivot = IMSL_SIZE(pivot)
      IF (size_pivot(0) NE 1) THEN message, 'PIVOT must be a 1-D array.'
      IF (N_ELEMENTS(pivot) NE n) THEN message, 'PIVOT is not the correct length.'
   END
   ; It is possible to call this procedure and not expect any output.
   ; If that is the case, then just return.
   IF ((nargs EQ 1) AND $
       ((KEYWORD_SET(basis) + KEYWORD_SET(ap) + KEYWORD_SET(q) + $
         KEYWORD_SET(r)) LT 1)) THEN RETURN
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_a(N_ELEMENTS(size_a)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG
   ; Always send space for Pivot.
   IF (nargs GT 1) THEN pivot_cvt = IMSL_LONG(pivot) ELSE pivot_cvt = IMSL_LONARR(n)
   ; BASIS = 0 is valid, so convert it carefully.
   CASE (arg_present(basis) + KEYWORD_SET(basis)) OF
      0:basis_cvt = IMSL_LONG((m < n))
      1:basis_cvt = IMSL_LONG((m < n))
      2:basis_cvt = IMSL_LONG(basis(0))
   END
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      a_cvt = double(transpose(a))
      ; Always send space for Auxqr, Qr.
      auxqr_spc = dblarr(n)
      qr_spc = dblarr(n, m)
      tmp = imsl_machine(/double)
      tolerance_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = double(tolerance(0))
      IF (ARG_PRESENT(q)) THEN q_spc = dblarr(m,m)
   END ELSE BEGIN
      a_cvt = float(transpose(a))
      ; Always send space for Auxqr, Qr.
      auxqr_spc = fltarr(n)
      qr_spc = fltarr(n, m)
      tmp = imsl_machine(/float)
      tolerance_cvt = sqrt(tmp.(3))
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt = float(tolerance(0))
      IF (ARG_PRESENT(q)) THEN q_spc = fltarr(m,m)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_176, type, err_status, a_cvt, m, n, $
                              pivot_cvt, $
                              auxqr_spc, $
                              qr_spc, $
                              q_spc, $
                              tolerance_cvt, $
                              basis_cvt
   ;
   ; Now copy over all output keywords results.
   ;
   ; The keywords R and AP are produced from the
   ; the contents of A and QR
   ; Compute R if needed
   IF (arg_present(r)) THEN BEGIN
      IF (type EQ TYP_FLOAT) THEN r = fltarr(n, m) ELSE r = dblarr(n, m)
      FOR i = 0, ((n-1)<(m-1)) DO r(i:*, i) = qr_spc(i:*, i)
      r = transpose(r)
   END
   ; Compute AP if needed
   IF (arg_present(ap)) THEN BEGIN
      ap = transpose(a_cvt(pivot_cvt-1, *))
   END
   
   IF (nargs GT 1) THEN pivot = pivot_cvt
   IF (nargs EQ 4) THEN BEGIN
      auxqr = auxqr_spc
      qr = transpose(qr_spc)
   END
   IF (ARG_PRESENT(q)) THEN q = transpose(q_spc)
   IF (ARG_PRESENT(basis)) THEN basis = basis_cvt
   ;
   ; End of procedure.
END
