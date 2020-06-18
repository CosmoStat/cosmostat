; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_hypoth_scph.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_hypoth_scph,  info_v, $              ;INPUT structure
                       h, $                        ;INPUT 2-D array: floating point
                       Double=double, $            ;INPUT Scalar ON/OFF flag
                       g=g, $                      ;INPUT 2-D array: floating point
                       u=u, $                      ;INPUT 2-D array: floating point
                       dfh=dfh                     ;OUTPUT Scalar floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;  - Make some simple checks on info_v.
   ;  - Get the values of n_dependent ans n_coefs out of the info_v structure.
   ;  - H is an array of size nh x n_coefs
   ;  - If keyword G is supplied, then it must be an array
   ;    of size (nh x nu)
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
   size_h = IMSL_SIZE(h)
   IF ((size_h(0) ne 1) AND (size_h(0) ne 2)) THEN MESSAGE, "H is not the correct size."
   IF (size_h(0) EQ 1) THEN BEGIN
      nh   =   IMSL_N_ELEMENTS(h)
      IF (n_coefs NE 1) THEN MESSAGE,  "H is not the correct size."
   END
   IF (size_h(0) EQ 2) THEN BEGIN
      nh   =   size_h(1)
      IF (n_coefs NE size_h(2)) THEN MESSAGE,  "H is not the correct size."
   END
   ; Check G, if supplied.
   nu_set = 0
   IF (KEYWORD_SET(g)) THEN BEGIN
      size_g = IMSL_SIZE(g)
      IF ((size_g(0) NE 1) AND (size_g(0) NE 2)) THEN MESSAGE,  "G is not the correct size."
      IF (size_g(0) EQ 1) THEN BEGIN
         IF (N_ELEMENTS(g) NE nh) THEN MESSAGE,  "G is not the correct size."
         nu  =  IMSL_1
         nu_set = IMSL_1
      END
      IF (size_g(0) EQ 2) THEN BEGIN
         IF (size_g(1) NE nh) THEN MESSAGE,  "G is not the correct size."
         nu  =  size_g(2)
         nu_set = IMSL_1
      END
   END
   ; Check U, if supplied.
   IF (KEYWORD_SET(u)) THEN BEGIN
      size_u = IMSL_SIZE(u)
      IF ((size_u(0) NE 1) AND (size_u(0) NE 2)) THEN MESSAGE,  "U is not the correct size."
      IF (size_u(0) EQ 1) THEN BEGIN
         IF (N_ELEMENTS(u) NE n_dependent) THEN MESSAGE, "U is not the correct size."
         IF (nu_set) THEN BEGIN
            IF (nu NE 1) THEN MESSAGE, "U is not the correct size."
         END ELSE BEGIN
            nu = IMSL_1
            nu_set  =  IMSL_1
         END
      END
      IF (size_u(0) EQ 2) THEN BEGIN
         IF (size_u(1) NE n_dependent) THEN MESSAGE,  "U is not the correct size."
         IF (nu_set) THEN BEGIN
            IF (nu NE size_u(2)) THEN MESSAGE, "U is not the correct size."
         END ELSE BEGIN
            nu = size_u(2)
            nu_set  =  IMSL_1
         END
      END
   END
   IF (nu_set EQ 0) THEN nu = n_dependent
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
      l_result = DBLARR(nu, nu)
      ; Input 
      dfh_spc = double(0)
      IF (nu GT 1) THEN BEGIN
         IF (KEYWORD_SET(g) EQ TRUE) THEN g_cvt = DOUBLE(TRANSPOSE(g))
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = DOUBLE(TRANSPOSE(u))
      END ELSE BEGIN
         IF (KEYWORD_SET(g) EQ TRUE) THEN g_cvt = DOUBLE( g )
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = DOUBLE( u )
      END     
      IF (n_coefs GT 1) THEN BEGIN
         h_cvt = DOUBLE(TRANSPOSE(h))
      END ELSE BEGIN
         h_cvt = DOUBLE( h )
      END     
      ; 
      ; Output 
   END ELSE BEGIN
      l_result = FLTARR(nu, nu)
      ; Input 
      dfh_spc = FLOAT(0)
      IF (nu GT 1) THEN BEGIN
         IF (KEYWORD_SET(g) EQ TRUE) THEN g_cvt = FLOAT(TRANSPOSE(g))
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = FLOAT(TRANSPOSE(u))
      END ELSE BEGIN
         IF (KEYWORD_SET(g) EQ TRUE) THEN g_cvt = FLOAT( g )
         IF (KEYWORD_SET(u) EQ TRUE) THEN u_cvt = FLOAT( u )
      END     
      IF (n_coefs GT 1) THEN BEGIN
         h_cvt = FLOAT(TRANSPOSE(h))
      END ELSE BEGIN
         h_cvt = FLOAT( h )
      END     
      ; 
      ; Output 
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
;   INFO, type, err_status, info_v, IMSL_N_ELEMENTS(info_v), info_type,   $
;                    dfh_spc, h_cvt,  n_dependent, nu, $
;                    nh, $
;                    n_coefs, $
;                    g_cvt, $
;                    u_cvt, $
;                    l_result
   MATHSTAT_274, type, err_status, info_v, IMSL_N_ELEMENTS(info_v), info_type,   $
                    dfh_spc, h_cvt,  n_dependent, nu, $
                    nh, $
                    n_coefs, $
                    g_cvt, $
                    u_cvt, $
                    l_result
   ;
   ; Now copy over all output keywords results.
   ;
   dfh = dfh_spc
   ;
   ; Return.
   ;
   IF (nu GT 1) THEN l_result = TRANSPOSE(l_result)
   RETURN, l_result
END
   

