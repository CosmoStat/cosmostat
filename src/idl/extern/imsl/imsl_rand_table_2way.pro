; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_rand_table_2way.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_rand_table_2way, row_totals, $  ;INPUT 1-D array: LONG
                   col_totals                 ;INPUT 1-D array: LONG

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking:
   ; The following checks are performed.
   ;   ROW_TOTALS must be a 1D array with atleast 2 elements
   ;   COL_TOTALS must be a 1D array with atleast 2 elements
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, "Incorrect number of arguments."
   size_rt = IMSL_SIZE(row_totals)
   IF (size_rt(0) NE 1) THEN BEGIN
      message, "ROW_TOTALS must be a 1-D array."
   END
   nrows = IMSL_N_ELEMENTS(row_totals)
   if (nrows LT 2) THEN MESSAGE, "ROW_TOTALS is not the correct length.."

   size_ct = IMSL_SIZE(col_totals)
   IF (size_ct(0) NE 1) THEN BEGIN
      message, "COL_TOTALS must be a 1-D array."
   END
   ncols = IMSL_N_ELEMENTS(col_totals)
   if (ncols LT 2) THEN MESSAGE, "COL_TOTALS is not the correct length.."

   type = TYP_MEMINT
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG keyword(s)
   ; Output LONG keyword(s)/Args
   row_totals_cvt = IMSL_LONG(row_totals)
   col_totals_cvt = IMSL_LONG(col_totals)
   result = IMSL_LONARR(ncols, nrows)
   ;
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_302,type,err_status, nrows, ncols, row_totals_cvt, col_totals_cvt, $
                           result
   ;
   ; Return.
   ;
   RETURN, transpose(result)
END

                   
                   
                   

  
      

  
