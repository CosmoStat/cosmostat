; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_datetodays.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
FUNCTION imsl_datetodays,day, $           ;INPUT Scalar LONG
                month, $               ;INPUT Scalar LONG
                year                   ;INPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;                       
   nargs = n_params()
   IF ((nargs LT 0) OR (nargs GT 3)) THEN message, 'Incorrect number of arguments.'
   ;
   IF (nargs GE 1) THEN day_cvt = IMSL_LONG(day(0)) ELSE day_cvt = IMSL_1
   IF (nargs GE 2) THEN month_cvt = IMSL_LONG(month(0)) ELSE month_cvt = IMSL_1
   IF (nargs GE 3) THEN year_cvt = IMSL_LONG(year(0)) ELSE year_cvt = IMSL_LONG(1900)
   result = IMSL_0
   ;
   ; Call the system function.
   ;
   type = TYP_MEMINT
   err_status = 0L
   MATHSTAT_130, type, err_status, day_cvt, month_cvt, year_cvt, result

   RETURN, result
END

