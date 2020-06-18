; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_statdata.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_statdata, choice, $       ;INPUT Scalar LONG
                     double=double   ;INPUT Scalar ON/OFF flag
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - CHOICE must be in the range 1, 2,  ..., 9.
   ;
   nargs = n_params()
   IF (nargs NE 1) THEN message, "Incorrect number of arguments."
   IF ((IMSL_LONG(choice(0)) LT 1) OR (IMSL_LONG(choice(0)) GT 9)) THEN $
     message, 'The value OF CHOICE must be in the range 1,2,...,9.'
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   choice_cvt = IMSL_LONG(choice(0))
   ;
   ; Floating point arguments and keywords
   res_dims = [[7, 2, 5, 1, 5, 1,2, 4, 34], $
               [16, 176, 150, 144, 13, 197, 296, 100, 113]]
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      result = dblarr(res_dims(choice_cvt-1, 0), res_dims(choice_cvt-1, 1))
   END ELSE BEGIN
      ; Result vector.
      result = fltarr(res_dims(choice_cvt-1, 0), res_dims(choice_cvt-1, 1))
   END
   
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_191, type, err_status, choice_cvt, result
                           
   ;
   ; Transpose the result before returning.
   RETURN, transpose(result)
END

                   
                   
                   

  
      

  
