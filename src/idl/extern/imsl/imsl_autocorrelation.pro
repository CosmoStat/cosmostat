; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_autocorrelation.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_autocorrelation, x, $                  ;INPUT 1-D array: floating point
                     lagmax, $                       ;INPUT 1-D Scalar LONG
                     xmean_in=xmean_in, $            ;INPUT Scalar LONG
                     double=double, $                ;INPUT Scalar ON/OFF flag
                     acv=acv, $                      ;OUTPUT 1-D array: floating point
                     se_option=se_option, $          ;INPUT 1-D Scalar LONG
                     seac=seac, $                    ;OUTPUT 1-D array: floating point
                     xmean_out=xmean_out             ;OUTPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   X must be a 1D array
   ;   LAGMAX must be greater than 0.
   ;   SE_OPTION and SEAC must be used together.
   
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_x = IMSL_SIZE(x)
   IF (size_x(0) NE 1) THEN BEGIN
      message, "X must be a 1-D array."
   END
   nobs  =  IMSL_N_ELEMENTS(x)
   IF (lagmax LT 1) THEN MESSAGE,  "LAGMAX must be positive."
   
   IF ((KEYWORD_SET(se_option)+KEYWORD_SET(seac)) EQ 1) THEN $
   IF ((KEYWORD_SET(xmean_in)+ARG_PRESENT(xmean_out)) EQ 2) THEN $
     MESSAGE, "Keywords XMEAN_IN and XMEAN_OUT cannot be used together."
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
   ; Input LONG keyword(s)
   IF (KEYWORD_SET(se_option) EQ TRUE) THEN se_option_cvt = IMSL_LONG(se_option)
   ;
   ; Input LONG argument(s)
   lagmax_cvt = IMSL_LONG(lagmax)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      x_cvt = double(x)
      IF (KEYWORD_SET(xmean_in) EQ TRUE) THEN xmean_in_cvt = DOUBLE(xmean_in)
      ; Output
      IF (ARG_PRESENT(xmean_out) EQ TRUE) THEN xmean_out_spc = DOUBLE(0.0)
      IF (ARG_PRESENT(acv) EQ TRUE) THEN acv_spc = dblarr(lagmax_cvt+1)
      IF (ARG_PRESENT(seac) EQ TRUE) THEN seac_spc = dblarr(lagmax_cvt)
      result = dblarr(lagmax_cvt+1)
   END ELSE BEGIN
      ; Input
      x_cvt = float(x)
      IF (KEYWORD_SET(xmean_in) EQ TRUE) THEN xmean_in_cvt = FLOAT(xmean_in)
      ; Output
      IF (ARG_PRESENT(xmean_out) EQ TRUE) THEN xmean_out_spc = FLOAT(0.0)
      IF (ARG_PRESENT(acv) EQ TRUE) THEN acv_spc = fltarr(lagmax_cvt+1)
      IF (ARG_PRESENT(seac) EQ TRUE) THEN seac_spc = fltarr(lagmax_cvt)
      result = fltarr(lagmax_cvt+1)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_255,type,err_status, x_cvt, lagmax_cvt, nobs, $
                           acv_spc, $
                           seac_spc, $
                           se_option_cvt, $
                           xmean_in_cvt, $
                           xmean_out_spc, $
                           result
                           
   ;
   ; Now copy over all output keywords results.
   ;
   IF (ARG_PRESENT(xmean_out) EQ TRUE) THEN xmean_out = xmean_out_spc
   IF (ARG_PRESENT(acv) EQ TRUE) THEN acv = acv_spc
   IF (ARG_PRESENT(seac) EQ TRUE) THEN seac = seac_spc
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
