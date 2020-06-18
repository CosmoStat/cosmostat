; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_hypoth_partial.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_hypoth_partial,  info_v, $                ;INPUT structure
                       hp, $                       ;INPUT 2-D array: floating point
                       Double=double, $            ;INPUT Scalar ON/OFF flag
                       gp=gp, $                    ;INPUT 2-D array: floating point
                       g_matrix=g_matrix, $        ;OUTPUT 2-D array: floating point
                       h_matrix=h_matrix, $        ;OUTPUT 2-D array: floating point
                       rank_hp=rank_hp             ;OUTPUT Scalar LONG
@imsl_init.pro

   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;  - Make some simple checks on info_v.
   ;  - Get the values of n_dependent ans n_coefs out of the info_v structure.
   ;  - HP is an array of size nhp x n_coefs
   ;  - If keyword GP is supplied, then it must be an array
   ;    of size (nhp x nu), set nu here.
   ;  - If keyword U is supplied, then it must be anarray
   ;    of size (n_dependent x nu)
   ;  NOTE: In the above checks on 2D arrays, some of the dimensions
   ;        may be one, which means if the dimension may go away
   ;        if it is the last dimension of the array (thanks to WAVE).
   ;        This makes the error checking much more difficult.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_info_v = IMSL_SIZE(info_v)
   IF ((size_info_v(0) NE 1) OR (size_info_v(N_ELEMENTS(size_info_v)-2) NE 1)  ) THEN $
       message, "PREDICT_INFO must be a 1-D array of type BYTE"
   IF (size_info_v(1) LT 23) THEN $
       MESSAGE,  "PREDICT_INFO is not valid, it must be exactly as returned by MULTIREGRESS"
   n_dependent  =  -IMSL_1
   n_coefs  =  -IMSL_1
   ; Get the value of n_dependent and n_coefs from the info structure.
   err_status = 0L
   info_type  =  IMSL_0
   MATHSTAT_272, err_status, info_v, IMSL_N_ELEMENTS(info_v), n_dependent, n_coefs, info_type
   ;  Check H
   size_hp = IMSL_SIZE(hp)
   IF ((size_hp(0) ne 1) AND (size_hp(0) ne 2)) THEN MESSAGE, "HP is not the correct size."
   IF (size_hp(0) EQ 1) THEN BEGIN
      nhp   =   IMSL_N_ELEMENTS(hp)
      IF (n_coefs NE 1) THEN MESSAGE,  "HP is not the correct size."
   END
   IF (size_hp(0) EQ 2) THEN BEGIN
      nhp   =   size_hp(1)
      IF (n_coefs NE size_hp(2)) THEN MESSAGE,  "HP is not the correct size."
   END
   ; Check G, if supplied.
   nu = -1
   nu_set = 0
   IF (KEYWORD_SET(gp)) THEN BEGIN
      size_gp = IMSL_SIZE(gp)
      IF ((size_gp(0) NE 1) AND (size_gp(0) NE 2)) THEN MESSAGE,  "GP is not the correct size."
      IF (size_gp(0) EQ 1) THEN BEGIN
         IF (N_ELEMENTS(gp) NE nhp) THEN MESSAGE,  "GP is not the correct size."
         nu  =  IMSL_1
         nu_set = IMSL_1
      END
      IF (size_gp(0) EQ 2) THEN BEGIN
         IF (size_gp(1) NE nhp) THEN MESSAGE,  "GP is not the correct size."
         nu  =  size_gp(2)
         nu_set = IMSL_1
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
   result  = IMSL_LONG(0)
   rank_hp_spc  = IMSL_LONG(0)
   ; Result array
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input 
      IF (nhp GT 1) THEN BEGIN
         hp_cvt = DOUBLE(TRANSPOSE(hp))
      END ELSE BEGIN
         hp_cvt = DOUBLE( hp )
      END
      IF (nu GT 1) THEN BEGIN
         IF (KEYWORD_SET(gp) EQ TRUE) THEN gp_cvt = DOUBLE(TRANSPOSE(gp))
      END ELSE BEGIN
         IF (KEYWORD_SET(gp) EQ TRUE) THEN gp_cvt = DOUBLE( gp )
      END
      ; 
      ; Output 
      IF (ARG_PRESENT(h_matrix) EQ TRUE) THEN h_matrix_spc = DBLARR(n_coefs, nhp)
      ; Note, default value for NU is N_DEPENDENT, so G is definet w.r.t. n_dependent.
      IF (ARG_PRESENT(g_matrix) EQ TRUE) THEN g_matrix_spc = DBLARR(n_dependent)
   END ELSE BEGIN
      ; Input 
      IF (nhp GT 1) THEN BEGIN
         hp_cvt = FLOAT(TRANSPOSE(hp))
      END ELSE BEGIN
         hp_cvt = FLOAT( hp )
      END
      IF (nu GT 1) THEN BEGIN
         IF (KEYWORD_SET(gp) EQ TRUE) THEN gp_cvt = FLOAT(TRANSPOSE(gp))
      END ELSE BEGIN
         IF (KEYWORD_SET(gp) EQ TRUE) THEN gp_cvt = FLOAT( gp )
      END
      ; 
      ; Output 
      IF (ARG_PRESENT(h_matrix) EQ TRUE) THEN h_matrix_spc = FLTARR(n_coefs, nhp)
      ; Note, default value for NU is N_DEPENDENT, so G is definet w.r.t. n_dependent.
      IF (ARG_PRESENT(g_matrix) EQ TRUE) THEN g_matrix_spc = FLTARR(n_dependent)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_278, type, err_status, info_v, IMSL_N_ELEMENTS(info_v), info_type,   $
                     nhp, $
                     hp_cvt, $
                     gp_cvt, $
                     rank_hp_spc, $
                     h_matrix_spc, $
                     g_matrix_spc, $
                     result
   ;
   ; Now copy over all output keywords results.
   ;
   rank_hp=rank_hp_spc
   IF (ARG_PRESENT(h_matrix) EQ TRUE) THEN h_matrix =  TRANSPOSE(h_matrix_spc(*, 0:result-1))
   IF (ARG_PRESENT(g_matrix) EQ TRUE) THEN g_matrix =  g_matrix_spc
   ; Return.
   ;
   RETURN, result
END
   

