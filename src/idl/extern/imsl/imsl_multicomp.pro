; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_multicomp.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_multicomp, means, $    ;INPUT 1-D array: floating point
                    df, $            ;INPUT Scalar LONG
                    std_err, $       ;INPUT Scalar floating point
                    double=double, $ ;INPUT Scalar ON/OFF flag 
                    alpha=alpha      ;INPUT Scalar floating point
@imsl_init.pro
   ON_ERROR, on_err_action
   ; Error checking.
   ; MEANS must be a 1-D array
   ; DF is converted to a LONG
   ; STDERR is converted to a float/double after deciding on th
   ;        the precision to use..
   ;
   nargs = n_params()
   IF (nargs EQ 3) THEN BEGIN
      size_means = IMSL_SIZE(means)
      IF (size_means(0) NE 1) THEN BEGIN
         message, "MEANS must be a 1-D array."
      END
      df_cvt = (IMSL_LONG(df))(0)
      size_std_err = IMSL_SIZE(std_err)
   END ELSE message, "Incorrect number of arguments."
      
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_means(N_ELEMENTS(size_means)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (size_std_err(N_ELEMENTS(size_std_err)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   n_groups = size_means(N_ELEMENTS(size_means)-1)
   ; 
   ; Result vector
   result = IMSL_LONARR(n_groups-1)
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; 
      ; Input 
      std_err_cvt = (double(std_err))(0)
      means_cvt = double(means)
      IF (KEYWORD_SET(alpha) EQ TRUE) THEN $
        alpha_cvt = (double(alpha))(0)
   END ELSE BEGIN
      ; 
      ; Input 
      std_err_cvt = (float(std_err))(0)
      means_cvt = float(means)
      IF (KEYWORD_SET(alpha) EQ TRUE) THEN $
        alpha_cvt = (float(alpha))(0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_160, type, err_status, $
     means_cvt, n_groups, df_cvt, std_err_cvt, alpha_cvt, result
   ;
   ; Return.
   ;
   RETURN, result
END

                   
                   
                   

  
      

  
