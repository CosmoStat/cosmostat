; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_constant.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_constant, name, $               ;INPUT Scalar STRING
                units, $                   ;INPUT Scalar STRING
                double=double              ;INPUT Scalar ON/OFF flag
   
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - NAME must be a scalar string.
   ; - if UNITS provided, it  must be a scalar string.
   ;                       
   nargs = n_params()
   IF ((nargs NE 1) AND (nargs NE 2)) THEN message, 'Incorrect number of arguments.'
   size_name = IMSL_SIZE(name)
   IF ((N_ELEMENTS(name) NE 1) OR (size_name(N_ELEMENTS(size_name)-2) NE 7)) THEN $
     message, 'Name must be a scalar string.'
   IF (nargs EQ 2) THEN BEGIN
      size_units = IMSL_SIZE(units)
      IF ((N_ELEMENTS(units) NE 1) OR $
          (size_units(N_ELEMENTS(size_units)-2) NE 7)) THEN $
        message, 'Units must be a scalar string.'
   END
   ; 
   ; Decide on what precision to use.
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ; Input LONG keyword(s)
   ;
   ; Floating point arguments and keywords   
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
   END ELSE BEGIN
      result = float(0.0)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_124, type, err_status, name, units, result

   RETURN, result
END

