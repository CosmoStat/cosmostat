; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_daystodate.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
   
PRO imsl_daystodate, days, $              ;INPUT Scalar LONG
                day, $                 ;OUTPUT Scalar LONG
                month, $               ;OUTPUT Scalar LONG
                year                   ;OUTPUT Scalar LONG
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ;                       
   nargs = n_params()
   IF ((nargs LT 1) OR (nargs GT 4)) THEN message, 'Incorrect number of arguments.'
   ;
   days_cvt = IMSL_LONG(days(0))
   day_spc = IMSL_0
   month_spc = IMSL_0
   year_spc = IMSL_0
   ;
   ; Call the system function.
   ;
   type = TYP_MEMINT
   err_status = 0L
   MATHSTAT_131, type, err_status, days_cvt, day_spc, month_spc, year_spc
   day = day_spc
   month = month_spc
   year=year_spc
END

