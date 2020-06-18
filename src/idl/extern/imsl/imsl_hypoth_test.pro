; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_hypoth_test.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_hypoth_test,  info_v, $              ;INPUT structure
                       dfh, $                      ;INPUT Scalar floating point
                       scph, $                     ;INPUT 2-D array: floating point
                       Double=double, $            ;INPUT Scalar ON/OFF flag
                       u=u, $                      ;INPUT 2-D array: floating point
                       wilk_lambda=wilk_lambda, $  ;OUTPUT 1-D array: floating point
                       roy_max_root=roy_max_root, $  ;OUTPUT 1-D array: floating point
                       hotelling_trace=hotelling_trace, $  ;OUTPUT 1-D array: floating point
                       pillai_trace=pillai_trace   ;OUTPUT 1-D array: floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;  - Make some simple checks on info_v.
   ;  - Get the value of n_dependent out of the info_v structure.
   ;  - SCPH must be a 2D array, of size nu X nu. (set nu here)
   ;  - If keyword U is supplied, then it must be a 2D array
   ;    of size (n_dependent x nu)
   ;  NOTE: In the above checks on 2D arrays, some of the dimensions
   ;        may be one, which means if the dimension may go away
   ;        if it is the last dimension of the array (thanks to WAVE).
   ;        This makes the error checking much more difficult.
    ;    
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN $
         message, "Incorrect number of arguments."
   size_info_v = IMSL_SIZE(info_v)
   IF ((size_info_v(0) NE 1) OR (size_info_v(N_ELEMENTS(size_info_v)-2) NE 1)  ) THEN $
       message, "PREDICT_INFO must be a 1-D array of type BYTE"
   IF (size_info_v(1) LT 23) THEN $
       MESSAGE,  "PREDICT_INFO is not valid, it must be exactly as returned by MULTIREGRESS"
   n_dependent  =  -IMSL_1
   ; Get the value of n_dependent form the info structure.
   err_status = 0L
   info_type = IMSL_0
   MATHSTAT_272, err_status, info_v, IMSL_N_ELEMENTS(info_v), n_dependent, n_coefs, info_type
   ;
   size_scph = IMSL_SIZE(scph)
   IF ((size_scph(0) ne 1) AND (size_scph(0) ne 2)) THEN MESSAGE, "SCPH is not the correct size."
   IF (N_ELEMENTS(scph) EQ 1) THEN nu = 1 $
     ELSE BEGIN
      IF (size_scph(0) NE 2) THEN MESSAGE, "SCPH is not the correct size."
      IF (size_scph(1) NE size_scph(2)) THEN MESSAGE,   "SCPH must be square."
      nu  =  size_scph(1)
      IF (nu GT n_dependent) THEN MESSAGE, "SCPH is not the correct size."
   END
  
   IF (KEYWORD_SET(u)) THEN BEGIN
      size_u = IMSL_SIZE(u)
      IF (n_dependent EQ 1) THEN BEGIN
         ; If N_DEPENDENT is 1, then U needs to be a single value.
         IF (n_elements(u) NE 1) THEN MESSAGE,  "U is not the correct size."
      END ELSE BEGIN
         ; N_DEPENDENT greater than 1 forces U to be a 1D or 2D array.
         IF ((size_u(0) NE 1) AND (size_u(0) NE 2)) THEN MESSAGE, "U is not the correct size."
         ; If U is 1D then NU needs to be 1, and N_ELEMENTS(u) must be N_DEPENDENT.
         IF (size_u(0) EQ 1) THEN BEGIN
            IF ((nu NE 1) OR (N_ELEMENTS(u) NE n_dependent)) THEN MESSAGE,  "U is not the correct size."
         END
         ; If U is 2D, then it must be size N_DEPENDENT x NU.
         IF (size_u(0) EQ 2) THEN BEGIN
            IF ((size_u(1) NE n_dependent) OR (size_u(2) NE nu)) THEN MESSAGE,  "U is not the correct size."
         END
      END
   END
   ;
   ; Decide on what precision to use.
   ; Use the precision of PREDICT_INFO.
   ;
   type  =  TYP_FLOAT
   IF (info_type GT TYP_FLOAT) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Result array
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      l_result = DOUBLE(0)
      ; 
      ; Input 
      dfh_cvt = double(dfh(0))      
      IF (nu GT 1) THEN BEGIN
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = DOUBLE(TRANSPOSE(u))
         scph_cvt = DOUBLE(TRANSPOSE(scph))
      END ELSE BEGIN
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = DOUBLE( u )
         scph_cvt = DOUBLE( scph )
      END     
      ; 
      ; Output 
      IF (ARG_PRESENT(wilk_lambda) EQ TRUE) THEN wl_spc = dblarr(2)
      IF (ARG_PRESENT(roy_max_root) EQ TRUE) THEN rm_spc = dblarr(2)
      IF (ARG_PRESENT(hotelling_trace) EQ TRUE) THEN ht_spc = dblarr(2)
      IF (ARG_PRESENT(pillai_trace) EQ TRUE) THEN pl_spc = dblarr(2)
   END ELSE BEGIN
      l_result = float(0)
      ; 
      ; Input 
      dfh_cvt = float(dfh(0))
      IF (nu GT 1) THEN BEGIN
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = float(TRANSPOSE(u))
         scph_cvt = float(TRANSPOSE(scph))
      END ELSE BEGIN
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = float( u )
         scph_cvt = float( scph )
      END     
      ; 
      ; Output 
      IF (ARG_PRESENT(wilk_lambda) EQ TRUE) THEN wl_spc = fltarr(2)
      IF (ARG_PRESENT(roy_max_root) EQ TRUE) THEN rm_spc = fltarr(2)
      IF (ARG_PRESENT(hotelling_trace) EQ TRUE) THEN ht_spc = fltarr(2)
      IF (ARG_PRESENT(pillai_trace) EQ TRUE) THEN pl_spc = fltarr(2)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_273, type, err_status, info_v, IMSL_N_ELEMENTS(info_v), info_type,   $
                    dfh_cvt, scph_cvt,  n_dependent, nu, $
                    u_cvt, $
                    wl_spc, $
                    rm_spc, $
                    ht_spc, $
                    pl_spc, $
                    l_result
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(wilk_lambda) EQ TRUE) THEN wilk_lambda = wl_spc
      IF (ARG_PRESENT(roy_max_root) EQ TRUE) THEN roy_max_root = rm_spc
      IF (ARG_PRESENT(hotelling_trace) EQ TRUE) THEN hotelling_trace = ht_spc
      IF (ARG_PRESENT(pillai_trace) EQ TRUE) THEN pillai_trace = pl_spc
   ;
   ; Return.
   ;
   RETURN, l_result
END
   

