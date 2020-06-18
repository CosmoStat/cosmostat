; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_machine.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_machine, $ 
                double=double, $      ;INPUT Scalar ON/OFF flag
                float=float           ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;          
   nargs = n_params()
   IF (nargs NE 0) THEN $
     message, 'Incorrect number of arguments.'
   ;
   ; Decide on what precision to use.
   ;
   CASE (KEYWORD_SET(float) + KEYWORD_SET(double)) OF
      0: type = TYP_MEMINT
      1: IF (KEYWORD_SET(float)) THEN type = TYP_FLOAT ELSE type = TYP_DOUBLE
      2: type = TYP_DOUBLE
      ELSE:  type = TYP_MEMINT
   END
   ;
   ; Setup the parameters for the call to the system function.
   ;
   CASE type OF
      TYP_MEMINT:BEGIN
         result = {IMACHINE, BITS_PER_CHAR:IMSL_0, INTEGER_BASE:IMSL_0, $
                    INTEGER_DIGITS:IMSL_0, MAX_INTEGER:IMSL_0, LONG_DIGITS:IMSL_0, $
                    MAX_LONG:IMSL_0, FLOAT_BASE:IMSL_0, FLOAT_DIGITS:IMSL_0, FLOAT_MIN_EXP:IMSL_0, $
                    FLOAT_MAX_EXP:IMSL_0, DOUBLE_DIGITS:IMSL_0, DOUBLE_MIN_EXP:IMSL_0, $
                    DOUBLE_MAX_EXP:IMSL_0}
         l_result = IMSL_LONARR(13)
      END 
      TYP_FLOAT: BEGIN
         result = {FMACHINE, MIN_POS:FLOAT(0.0), MAX_POS:FLOAT(0.0), $
                    MIN_REL_SPACE:FLOAT(0.0), MAX_REL_SPACE:FLOAT(0.0), $
                    LOG_10:FLOAT(0.0), NAN:FLOAT(0.0), $
                    POS_INF:FLOAT(0.0), NEG_INF:FLOAT(0.0)}
         l_result = fltarr(8)
      END
      TYP_DOUBLE:BEGIN
         result = {DMACHINE, MIN_POS:DOUBLE(0.0), MAX_POS:DOUBLE(0.0), $
                    MIN_REL_SPACE:DOUBLE(0.0), MAX_REL_SPACE:DOUBLE(0.0), $
                    LOG_10:DOUBLE(0.0), NAN:DOUBLE(0.0), $
                    POS_INF:DOUBLE(0.0), NEG_INF:DOUBLE(0.0)}
         l_result = dblarr(8)
      END
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_159, type, err_status, l_result
   ;
   ; Now copy over the result
   FOR i = 0, N_ELEMENTS(l_result)-1 DO result.(i) = l_result(i)
   
   RETURN, result
END

                   
                   
                   

  
      

  
