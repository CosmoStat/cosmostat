; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_random_npp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_random_npp, tbegin, $           ;INPUT Scalar floating point
                     tend, $                  ;INPUT Scalar floating point
                     ftheta, $                ;INPUT Scalar floating point
                     theta_min, $             ;INPUT Scalar floating point
                     theta_max, $             ;INPUT Scalar floating point
                     neub, $                  ;INPUT Scalar floating point
                     double=double            ;INPUT Scalar ON/OFF flag
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - NEUB greater than 1.
   ;   - FTHETA must be a string.
   ;
   nargs = n_params()
   IF (nargs LT 6) THEN MESSAGE,  "Incorrect number of arguments."
   neub_cvt = IMSL_LONG(neub(0))
   IF (neub_cvt LE 0) THEN MESSAGE,  "NEUB must be greater than 0."
   size_ftheta = IMSL_SIZE(ftheta)
   IF ((N_ELEMENTS(ftheta) NE 1) OR (size_ftheta(N_ELEMENTS(size_ftheta)-2) NE 7)) THEN $
     message, 'FTHETA must be a scalar string.'

   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (SIZE(tbegin(0), /Type) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (SIZE(tend(0), /Type) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (SIZE(theta_min(0), /Type) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (SIZE(theta_max(0), /Type) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      tbegin_cvt = DOUBLE(tbegin(0))
      tend_cvt = DOUBLE(tend(0))
      theta_min_cvt = DOUBLE(theta_min(0))
      theta_max_cvt = DOUBLE(theta_max(0))
      result = DBLARR(neub_cvt)
   END ELSE BEGIN
      ; Input
      tbegin_cvt = FLOAT(tbegin(0))
      tend_cvt = FLOAT(tend(0))
      theta_min_cvt = FLOAT(theta_min(0))
      theta_max_cvt = FLOAT(theta_max(0))
      result = FLTARR(neub_cvt)
   END
   ;
   ; Call the system function.
   ;
   ne_spc = IMSL_0
   err_status = 0L
   MATHSTAT_303,  type,  err_status,  $
                                      tbegin_cvt,   $
                                      tend_cvt,   $
                                      ftheta,  $
                                      theta_min_cvt,   $
                                      theta_max_cvt,  $
                                      neub_cvt,  $
                                      ne_spc,   $
                                      result
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ne_spc GT 0) THEN result = result(0:ne_spc-1)
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
