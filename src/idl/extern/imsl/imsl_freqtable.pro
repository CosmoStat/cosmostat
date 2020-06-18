; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_freqtable.pro#1 $      
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_freqtable, X, $                  ;INPUT 1-D array: floating point 
           nxbins,  $                     ;INPUT Scalar LONG
           y, $                           ;INPUT 1-D array: floating point 
           nybins,  $                     ;INPUT Scalar LONG
           class_marks=class_marks, $     ;INPUT 1-D array: floating point 
           cutpoints=cutpoints, $         ;INPUT 1-D array: floating point 
           lower_bound=lower_bound, $     ;INPUT Scalar floating point
           upper_bound=upper_bound, $     ;INPUT Scalar floating point
           double=double                  ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ; If nargs EQ 2, then we are doing a oneway freqtable.
   ;  - X must be a 1-D array. (set n = length)
   ;  - nxbins must be greater than 0.
   ;  - UPPER_BOUND and LOWER_BOUND must be used together.
   ;  - If CUTPOINTS is set, then it must be a 1-D array of length (nxbins-1).
   ;  - If CUTPOINTS is set, then no other keywords are allowed.
   ;  - If CLASS_MARKS is set, then it must be a 1-D array of length (nxbins).
   ;  - If CLASS_MARKS is set, then no other keywords are allowed.
   ; If nargs EQ 4, then we are doing a oneway freqtable.
   ;  - X and Y both must be 1D arrays of same length.  (set n = length)
   ;  - nxbins must be greater than 0.
   ;  - nybins must be greater than 0.
   ;  - UPPER_BOUND and LOWER_BOUND must be used together, and both must be
   ;    two element arrays.
   ;  - If CUTPOINTS is set, then it must be a 1-D array of length (nxbins-1)+(nybins-1).
   ;  - If CUTPOINTS is set, then no other keywords are allowed.
   ;  - If CLASS_MARKS is set, then it must be a 1-D array of length (nxbins+nybins).
   ;  - If CLASS_MARKS is set, then no other keywords are allowed.
   ; 
   ;          
   nargs = n_params()
   IF ((nargs NE 2) AND (nargs NE 4))  THEN MESSAGE,  'Incorrect number of arguments.'
   size_x    =    IMSL_LONG(SIZE(x))
   IF (size_x(0) NE 1) THEN MESSAGE,   'X must be a 1-D array.'
   n    =    IMSL_LONG(size_x(1))
   IF (nargs EQ 4) THEN BEGIN
      size_y    =    IMSL_LONG(SIZE(y))
      IF (size_y(0) NE 1) THEN MESSAGE,    'Y must be a 1-D array.'
      IF (size_y(1) NE n) THEN MESSAGE,     'X and Y must be the same length.'
      nybins_cvt  =  IMSL_LONG(nybins(0))
      IF (nybins_cvt LT 1) THEN MESSAGE,   'NYBINS must be positive.'
   END
   nxbins_cvt   =   IMSL_LONG(nxbins(0))
   IF (nxbins_cvt LT 1) THEN MESSAGE,   'NXBINS must be positive.'
   IF ((KEYWORD_SET(upper_bound) + KEYWORD_SET(lower_bound)) EQ 1) THEN $
     MESSAGE,   'LOWER_BOUND AND UPPER_BOUND must be used together.'
   IF ((KEYWORD_SET(upper_bound)) AND (nargs EQ 4)) THEN BEGIN
      size_ub   =   IMSL_LONG(SIZE(upper_bound))
      IF (size_ub(0) NE 1) THEN $
        MESSAGE,     'UPPER_BOUND must be a 1-D array for a two-way frequency table.'
      IF (N_ELEMENTS(upper_bound) NE 2) THEN $
        MESSAGE,    'UPPER_BOUND must have exactly two elements for a two-way frequency table.'
      size_lb   =   IMSL_LONG(SIZE(lower_bound))
      IF (size_lb(0) NE 1) THEN $
        MESSAGE,     'LOWER_BOUND must be a 1-D array for a two-way frequency table.'
      IF (N_ELEMENTS(lower_bound) NE 2) THEN $
        MESSAGE,    'UPPER_BOUND must have exactly two elements for a two-way frequency table.'
   END
   
   IF (KEYWORD_SET(cutpoints) OR KEYWORD_SET(class_marks)) THEN BEGIN 
      IF ((KEYWORD_SET(lower_bound) + KEYWORD_SET(class_marks)) EQ 2) THEN $
        MESSAGE,   'LOWER_BOUND and CLASS_MARKS are mutually exclusive.'
      IF ((KEYWORD_SET(lower_bound) + KEYWORD_SET(cutpoints)) EQ 2) THEN $
        MESSAGE,   'LOWER_BOUND and CUTPOINTS are mutually exclusive.'
      IF ((KEYWORD_SET(class_marks) + KEYWORD_SET(cutpoints)) EQ 2) THEN $
        MESSAGE,   'CLASS_MARKS and CUTPOINTS are mutually exclusive.'
   END
   
   IF (KEYWORD_SET(cutpoints)) THEN BEGIN 
      size_ctpts   =   IMSL_LONG(SIZE(cutpoints))
      IF (size_ctpts(0) NE 1) THEN MESSAGE,    'CUTPOINTS must be a 1-D array.'
      IF (nargs EQ 2) THEN BEGIN 
         IF (size_ctpts(1) NE (nxbins_cvt-1)) THEN $
           MESSAGE,  'CUTPOINTS is not the correct size.'
      END ELSE BEGIN
         IF (size_ctpts(1) NE (nxbins_cvt-1)+(nybins_cvt-1)) THEN $
           MESSAGE,   'CUTPOINTS is not the correct size.'
      END
   END
        
   IF (KEYWORD_SET(class_marks)) THEN BEGIN 
      size_cmarks   =   IMSL_LONG(SIZE(class_marks))
      IF (size_cmarks(0) NE 1) THEN MESSAGE,   'CLASS_MARKS must be a 1-D array.'
      IF (nargs EQ 2) THEN BEGIN 
         IF (size_cmarks(1) NE (nxbins_cvt)) THEN $
           MESSAGE,    'CLASS_MARKS is not the correct size.'
      END ELSE BEGIN
         IF (size_cmarks(1) NE (nxbins_cvt+nybins_cvt)) THEN $
           MESSAGE,     'CLASS_MARKS is not the correct size.'
      END
   END         
   ;
   ; Decide on what precision to use.
   ;
   type = TYP_FLOAT
   IF (size_x(N_ELEMENTS(size_x)-2) EQ  TYP_DOUBLE) THEN type = TYP_DOUBLE
   IF (nargs EQ 4) THEN BEGIN 
      IF (size_y(N_ELEMENTS(size_y)-2) EQ  TYP_DOUBLE) THEN type   =   TYP_DOUBLE
   END
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Floating point arguments and keywords
   ;
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      IF (nargs EQ 4) THEN BEGIN
         result    =    DBLARR(nybins_cvt,    nxbins_cvt)
         y_cvt = double(y)
      END ELSE  result      =      DBLARR(nxbins_cvt)
      x_cvt = double(x)
      IF (KEYWORD_SET(lower_bound)) THEN lbound_cvt = double(lower_bound)
      IF (KEYWORD_SET(upper_bound)) THEN ubound_cvt = double(upper_bound)
      IF (KEYWORD_SET(class_marks)) THEN cmarks_cvt = double(class_marks)
      IF (KEYWORD_SET(cutpoints)) THEN cutpoints_cvt = double(cutpoints)
   END ELSE BEGIN
      IF (nargs EQ 4) THEN BEGIN
         result    =    FLTARR(nybins_cvt,    nxbins_cvt)
         y_cvt = float(y)
      END ELSE  result      =      FLTARR(nxbins_cvt)
      x_cvt = float(x)
      IF (KEYWORD_SET(lower_bound)) THEN lbound_cvt = float(lower_bound)
      IF (KEYWORD_SET(upper_bound)) THEN ubound_cvt = float(upper_bound)
      IF (KEYWORD_SET(class_marks)) THEN cmarks_cvt = float(class_marks)
      IF (KEYWORD_SET(cutpoints)) THEN cutpoints_cvt = float(cutpoints)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_144, type, err_status, X_cvt, nxbins_cvt, y_cvt, nybins_cvt,n, $
                              cmarks_cvt, $
                              cutpoints_cvt, $
                              lbound_cvt, $
                              ubound_cvt, $
                              result
   IF (nargs EQ 4) THEN result = TRANSPOSE(result)
   ; Return
   RETURN, result
END
