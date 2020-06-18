; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_kalman.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_kalman, b, $                        ;INPUT/OUTPUT 1-D array: floating point
                 covb, $                     ;INPUT/OUTPUT 2-D array: floating point
                 n, $                        ;INPUT/OUTPUT Scalar: LONG
                 ss, $                       ;INPUT/OUTPUT Scalar: floating point
                 alndet, $                   ;INPUT/OUTPUT Scalar: floating point
                 Y=y, $                      ;INPUT 1-D array: floating point
                 Z=z, $                      ;INPUT 2-D array: floating point
                 R=r, $                      ;INPUT 2-D array: floating point
                 T_matrix=t_matrix, $        ;INPUT 2-D array: floating point
                 Q_matrix=q_matrix, $        ;INPUT 2-D array: floating point
                 Tolerance=tolerance, $      ;INPUT Scalar: floating point
                 V=v, $                      ;OUTPUT 2-D array: floating point
                 Covv=covv, $                ;OUTPUT 2-D array: floating point
                 Double=double
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;  B must be a 1-D array, set NB based on B.  NB=1 if B is scalar.
   ;  COVB must be NB by NB
   ;  Keywords Y, Z, and R must be used together
   ;    If present, Y is 1-D array of lenght NY, nset NY based on Y.
   ;    If present, Z must be NY by NB.
   ;    If present, R must be NY by NY.
   ;  If present, T must be NB by NB.
   ;  If present, Q must be NB by NB.
   ;  V is only allowed if Y is present
   ;  COVV is only allowed if Y is present.
   ;
   nargs = n_params()
   IF (nargs NE 5) THEN $
         message, "Incorrect number of arguments."
   size_b = IMSL_SIZE(b)
   IF ((size_b(0) LT 0) OR (size_b(0) GT 1)) THEN BEGIN 
      message, "B must be a scalar or 1-D array"
   END
   nb = IMSL_N_ELEMENTS(b)

   size_covb = IMSL_SIZE(covb)
   IF (nb EQ 1) THEN BEGIN
      IF (N_ELEMENTS(covb) NE 1) THEN MESSAGE, "COVB is not the correct size."
   END ELSE BEGIN
      IF (size_covb(0) NE 2) THEN MESSAGE, "COVB must be a 2-D array."
      IF ((size_covb(1) NE nb) OR (size_covb(2) NE nb)) THEN MESSAGE, "COVB is not the correct size."
   END

   ; Check keywords Y, Z, and R.
   ; Need to be careful since they can be set to 0.
   IF ((ARG_PRESENT(y) EQ TRUE) AND (N_ELEMENTS(y) GT 0)) THEN y_set = TRUE ELSE y_set = FALSE
   IF ((ARG_PRESENT(z) EQ TRUE) AND (N_ELEMENTS(z) GT 0)) THEN z_set = TRUE ELSE z_set = FALSE
   IF ((ARG_PRESENT(r) EQ TRUE) AND (N_ELEMENTS(r) GT 0)) THEN r_set = TRUE ELSE r_set = FALSE
   yzr = y_set+z_set+r_set
   IF ((yzr EQ 1) OR (yzr EQ 2)) THEN MESSAGE, "Keywords Y, Z, and R must be used together."
   IF (yzr GT 0) THEN BEGIN
      size_y = IMSL_SIZE(y)
      IF (size_y(0) GT 1) THEN MESSAGE, "Y must be a scalar or 1-D array."
      ny = IMSL_N_ELEMENTS(y)

      ; z must be ny by nb.
      IF (nb EQ 1) THEN BEGIN ; z[ny]
         IF (N_ELEMENTS(z) NE ny) THEN MESSAGE, "Z is not the correct size."
         size_z = IMSL_SIZE(z)
         IF ((ny GT 1) AND (size_z(0) NE 1)) THEN MESSAGE, "Z is not the correct size."
      END ELSE IF (ny EQ 1) THEN BEGIN ; z[nb]
         IF (N_ELEMENTS(z) NE nb) THEN MESSAGE, "Z is not the correct size."
         size_z = IMSL_SIZE(z)
         IF ((nb GT 1) AND (size_z(0) NE 1)) THEN MESSAGE, "Z is not the correct size."
      END ELSE BEGIN ; z[ny][nb]
         size_z = IMSL_SIZE(z)
         IF (size_z(0) NE 2) THEN MESSAGE, "Z is not the corrrect size."
         IF ((size_z(1) NE ny) OR (size_z(2) NE nb)) THEN MESSAGE, "Z is not the correct size."
      END

      ; r must be ny by ny
      IF ((ny EQ 1) AND (N_ELEMENTS(r) NE 1)) THEN MESSAGE, "R is not the correct size."
      size_r = IMSL_SIZE(r)
      IF (ny GT 1) THEN BEGIN
         IF (size_r(0) NE 2) THEN MESSAGE, "R is not the correct size."
         IF ((size_r(1) NE ny) OR (size_r(2) NE ny)) THEN MESSAGE, "R is not the correct size."
      END
   END ; Checks on Y, Z and R.

   ; Check keyword T_matrix
   IF ((ARG_PRESENT(T_matrix) EQ TRUE) AND (N_ELEMENTS(T_matrix) GT 0)) THEN T_matrix_set = TRUE ELSE T_matrix_set = FALSE
   IF (T_matrix_set EQ TRUE)  THEN BEGIN
      IF ((nb EQ 1) AND (N_ELEMENTS(T_matrix) NE 1)) THEN MESSAGE, "T_MATRIX is not the correct size."
      size_T_matrix = IMSL_SIZE(T_matrix)
      IF (nb GT 1) THEN BEGIN
         IF (size_T_matrix(0) NE 2) THEN MESSAGE, "T_MATRIX is not the correct size."
         IF ((size_T_matrix(1) NE nb) OR (size_T_matrix(2) NE nb)) THEN MESSAGE, "T_MATRIX is not the correct size."
      END
   END
   
   ; Check keyword Q_matrix
   IF ((ARG_PRESENT(q_matrix) EQ TRUE) AND (N_ELEMENTS(q_matrix) GT 0)) THEN q_matrix_set = TRUE ELSE q_matrix_set = FALSE
   IF (q_matrix_set EQ TRUE) THEN BEGIN
      IF ((nb EQ 1) AND (N_ELEMENTS(q_matrix) NE 1)) THEN MESSAGE, "Q_MATRIX is not the correct size."
      size_q_matrix = IMSL_SIZE(q_matrix)
      IF (nb GT 1) THEN BEGIN
         IF (size_q_matrix(0) NE 2) THEN MESSAGE, "Q_MATRIX is not the correct size."
         IF ((size_q_matrix(1) NE nb) OR (size_q_matrix(2) NE nb)) THEN MESSAGE, "Q_MATRIX is not the correct size."
      END
   END

   ; Check keyword V
   IF ((ARG_PRESENT(v) EQ TRUE) AND (y_set NE TRUE)) THEN $
     MESSAGE, "Keyword V can only be used if Y is also present."
   ; Check keyword COVV
   IF ((ARG_PRESENT(covv) EQ TRUE) AND (y_set NE TRUE)) THEN $
     MESSAGE, "Keyword COVV can only be used if Y is also present."
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_b(N_ELEMENTS(size_b)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_covb(N_ELEMENTS(size_covb)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   size_ss = IMSL_SIZE(ss)
   IF (size_ss(N_ELEMENTS(size_ss)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   size_alndet = IMSL_SIZE(alndet)
   IF (size_alndet(N_ELEMENTS(size_alndet)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function
   n_cvt = IMSL_LONG(n(0))
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input 
      b_cvt = double(b)
      covb_cvt = DOUBLE(covb) ; No transpose since it is symmetric.
      n_cvt = IMSL_LONG(n(0))
      ss_cvt = DOUBLE(ss(0))
      alndet_cvt = DOUBLE(alndet(0))
      
      IF (y_set EQ TRUE) THEN BEGIN
         y_cvt = DOUBLE(y)
         IF ((ny GT 1) OR (nb GT 1)) THEN z_cvt = DOUBLE(TRANSPOSE(z)) ELSE z_cvt = DOUBLE(z)
         IF (ny GT 1)  THEN r_cvt = DOUBLE(TRANSPOSE(r)) ELSE r_cvt = DOUBLE(r)
      END
      IF (t_matrix_set EQ TRUE) THEN BEGIN
         IF (nb GT 1)  THEN t_matrix_cvt = DOUBLE(TRANSPOSE(t_matrix)) ELSE t_matrix_cvt = DOUBLE(t_matrix)
      END
      IF (q_matrix_set EQ TRUE) THEN BEGIN
         IF (nb GT 1)  THEN q_matrix_cvt = DOUBLE(TRANSPOSE(q_matrix)) ELSE q_matrix_cvt = DOUBLE(q_matrix)
      END
      mach = imsl_machine(/DOUBLE)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  DOUBLE(tolerance(0)) $
        ELSE tolerance_CVT = 100.*(mach.max_rel_space)

      ; 
      ; Output 
      IF (ARG_PRESENT(v) EQ TRUE) THEN v_spc = dblarr(ny)
      IF (ARG_PRESENT(covv) EQ TRUE) THEN covv_spc = dblarr(ny, ny)
   END ELSE BEGIN
      ; 
      ; Input 
      b_cvt = float(b)
      covb_cvt = FLOAT(covb) ; No transpose since it is symmetric.
      n_cvt = IMSL_LONG(n(0))
      ss_cvt = FLOAT(ss(0))
      alndet_cvt = FLOAT(alndet(0))
      
      IF (y_set EQ TRUE) THEN BEGIN
         y_cvt = FLOAT(y)
         IF ((ny GT 1) OR (nb GT 1)) THEN z_cvt = FLOAT(TRANSPOSE(z)) ELSE z_cvt = FLOAT(z)
         IF (ny GT 1)  THEN r_cvt = FLOAT(TRANSPOSE(r)) ELSE r_cvt = FLOAT(r)
      END
      IF (t_matrix_set EQ TRUE) THEN BEGIN
         IF (nb GT 1)  THEN t_matrix_cvt = FLOAT(TRANSPOSE(t_matrix)) ELSE t_matrix_cvt = FLOAT(t_matrix)
      END
      IF (q_matrix_set EQ TRUE) THEN BEGIN
         IF (nb GT 1)  THEN q_matrix_cvt = FLOAT(TRANSPOSE(q_matrix)) ELSE q_matrix_cvt = FLOAT(q_matrix)
      END

      mach = imsl_machine(/FLOAT)
      IF (KEYWORD_SET(tolerance)) THEN tolerance_cvt  =  FLOAT(tolerance(0)) $
        ELSE tolerance_CVT = 100.*(mach.max_rel_space)

      IF (KEYWORD_SET(tolerance) EQ TRUE) THEN tolerance_cvt = FLOAT(tolerance(0))
      ; 
      ; Output 
      IF (ARG_PRESENT(v) EQ TRUE) THEN v_spc = fltarr(ny)
      IF (ARG_PRESENT(covv) EQ TRUE) THEN covv_spc = fltarr(ny, ny)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_312, type, err_status,  $
                    nb, b_cvt, covb_cvt, n_cvt, ss_cvt, alndet_cvt, $
                    y_cvt, z_cvt, r_cvt, ny, $
                    t_matrix_cvt, $
                    q_matrix_cvt, $
                    tolerance_cvt, $
                    v_spc, $
                    covv_spc
   ;
   ; Now copy over all input/output arguments and output keywords.
   ;
   b = b_cvt
   covb = covb_cvt
   n = n_cvt
   ss = ss_cvt
   alndet = alndet_cvt
   
   IF (ARG_PRESENT(v) EQ TRUE) THEN $
     v = TRANSPOSE(v_spc)
   IF (ARG_PRESENT(covv) EQ TRUE) THEN $
     covv = TRANSPOSE(covv_spc)
   ;
   ; Return.
   ;
   RETURN
END
   

