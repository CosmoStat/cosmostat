; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_random_table.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
PRO imsl_random_table, table, $          ;INPUT/OUTPUT 1-D array:
                  double=double, $       ;INPUT Scalar ON/OFF flag
                  get=get, $             ;INPUT Scalar ON/OFF flag
                  set=set, $             ;INPUT Scalar ON/OFF flag
                  gfsr=gfsr              ;INPUT Scalar ON/OFF flag
 
@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; The following checks are performed.
   ;   - Either GET or SET must be specified.
   ;   - If SET is set and GFSR is not set, then TABLE must be length 128.
   ;   - If SET and GFSR are set, then TABLE must be type long, and length 1565.
   ;
   nargs = n_params()
   IF (nargs LT 1) THEN MESSAGE,  "Incorrect number of arguments."
   IF ((KEYWORD_SET(set) + KEYWORD_SET(get)) NE 1) THEN $
      MESSAGE, "One of the keywords GET or SET must be specified."
   IF ((KEYWORD_SET(set)) AND (NOT KEYWORD_SET(gfsr))) THEN shuff_set = IMSL_1 ELSE shuff_set = IMSL_0
   IF ((KEYWORD_SET(get)) AND (NOT KEYWORD_SET(gfsr))) THEN shuff_get = IMSL_1 ELSE shuff_get = IMSL_0
   IF ((KEYWORD_SET(set)) AND (KEYWORD_SET(gfsr))) THEN gfsr_set = IMSL_1 ELSE gfsr_set = IMSL_0
   IF ((KEYWORD_SET(get)) AND (KEYWORD_SET(gfsr))) THEN gfsr_get = IMSL_1 ELSE gfsr_get = IMSL_0
  
   ; Setting shuffled table.
   IF (shuff_set EQ IMSL_1) THEN BEGIN
      size_table = IMSL_SIZE(table)
      IF (size_table(0) NE 1) THEN BEGIN
         message, "TABLE must be a 1-D array."
      END
      IF (n_elements(table) NE 128) THEN BEGIN
         message, "TABLE is not the correct size."
      END
   END
   ; Setting GFSR table.
   IF (gfsr_set EQ IMSL_1) THEN BEGIN
      size_table = IMSL_SIZE(table)
      IF (size_table(0) NE 1) THEN BEGIN
         message, "TABLE must be a 1-D array."
      END
      IF (n_elements(table) NE 1565) THEN BEGIN
         message, "TABLE is not the correct size."
      END
      table_cvt = IMSL_LONG(table)
   END
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(SET)) THEN BEGIN
     IF (SIZE(table(0), /Type) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   END
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   if (gfsr_get EQ IMSL_1) then table_cvt = IMSL_LONARR(1565)
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Input
      if (shuff_set EQ IMSL_1) then table_cvt = DOUBLE(table)
      ; Output
      if (shuff_get EQ IMSL_1) then table_cvt = DBLARR(128)
   END ELSE BEGIN
      ; Input
      if (shuff_set EQ IMSL_1) then table_cvt = FLOAT(table)
      ; Output
      if (shuff_get EQ IMSL_1) then table_cvt = FLTARR(128)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_301,  type,  err_status,  $
                                      table_cvt,   $
                                      gfsr_set,   $
                                      gfsr_get,   $
                                      shuff_set,   $
                                      shuff_get
   ;
   ; Now copy over all results.
   ;
   IF ((shuff_get EQ 1) OR (gfsr_get EQ 1)) THEN table = table_cvt
   ;
   ; Return.
   ;
   RETURN
END

                   
                   
                   

  
      

  
