; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_anova1.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;

FUNCTION imsl_anova1, n, $                       ;INPUT 1-D array: LONG
                 y, $                            ;INPUT 1-D array: floating point
                 confidence=confidence, $        ;INPUT Scalar floating point
                 double=double, $                ;INPUT Scalar ON/OFF flag
                 anova_table=anova_table, $      ;OUTPUT 1-D array: floating point
                 bonferroni=bonferroni, $        ;OUTPUT 1-D array: floating point
                 dunn_sidak=dunn_sidak, $        ;OUTPUT 1-D array: floating point
                 group_counts=group_counts, $    ;OUTPUT 1-D array: LONG
                 group_means=group_means, $      ;OUTPUT 1-D array: floating point
                 group_std_dev=group_std_dev, $  ;OUTPUT 1-D array: floating point
                 one_at_a_time=one_at_a_time, $  ;OUTPUT 1-D array: floating point
                 scheffe=scheffe, $              ;OUTPUT 1-D array: floating point
                 tukey=tukey                     ;OUTPUT 1-D array: floating point

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking.
   ; - Make sure the input arguments are simple arrays.
   ;   Check the length of the second argument passed to this function.
   ;   It should be equal to the product of the elements of the first
   ;   argument.
   ;
   ; - At most one of the following keywords may be supplied:
   ;   BONFERRONI, DUNN_SIDAK, ONE_AT_A_TIME, SCHEFFE, TUKEY
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN $
         message, "Incorrect number of arguments."
   size_n = IMSL_SIZE(n)
   IF (size_n(0) NE 1) THEN BEGIN
      message, "N must be a 1-D array."
   END
   IF (size_n(1) LT 2) THEN message, "N should have at least two elements."
   size_y = IMSL_SIZE(y)
   IF (size_y(0) NE 1) THEN BEGIN
      message, "Y must be a 1-D array."
   END
   exp_lgth_y = IMSL_0
   FOR i = IMSL_0, N_elements(n)-1 DO exp_lgth_y = exp_lgth_y + (IMSL_LONG(n(i)))
   IF (N_ELEMENTS(y) NE exp_lgth_y) THEN $
     message, "The length of Y is incorrect"
   ;
   i_tmp = 0
   IF (ARG_PRESENT(tukey) EQ TRUE) then i_tmp = i_tmp + 1
   IF (ARG_PRESENT(scheffe) EQ TRUE) then i_tmp = i_tmp + 1
   IF (ARG_PRESENT(bonferroni) EQ TRUE) then i_tmp = i_tmp + 1
   IF (ARG_PRESENT(dunn_sidak) EQ TRUE) then i_tmp = i_tmp + 1
   IF (ARG_PRESENT(one_at_a_time) EQ TRUE) then i_tmp = i_tmp + 1
   IF (i_tmp GT 1) THEN $
     message, "At most one of the following keywords may be specified:" + $
     "BONFERRONI, DUNN_SIDAK, ONE_AT_A_TIME, SCHEFFE, TUKEY"
   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (size_y(N_ELEMENTS(size_y)-2) GT TYP_FLOAT) THEN type = TYP_DOUBLE
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Setup the parameters for the call to the system function.
   ;
   ; Input LONG argument(s)
   n_cvt = IMSL_LONG(N)
   n_groups = IMSL_N_ELEMENTS(n)
   ;
   ; Output LONG keywords(s)
   IF (ARG_PRESENT(group_counts) EQ TRUE) THEN $
     grp_counts_spc = IMSL_LONARR(n_groups)
   ;
   ; In case one of the above 5 keywords was set, we will need the following
   ; number to figure out the size of the returned array.
   ; This is the 'n_groups choose 2' value.
   dim_method = (n_groups*(n_groups-1))/2
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      result = double(0.0)
      ; 
      ; Input 
      y_cvt = double(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (double(confidence))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table_spc = dblarr(15)
      IF (ARG_PRESENT(bonferroni) EQ TRUE) THEN $
        bonferroni_spc = dblarr(5, dim_method)
      IF (ARG_PRESENT(dunn_sidak) EQ TRUE) THEN $
        dunn_sidak_spc = dblarr(5, dim_method)
      IF (ARG_PRESENT(one_at_a_time) EQ TRUE) THEN $
        one_a_time_spc = dblarr(5, dim_method)
      IF (ARG_PRESENT(scheffe) EQ TRUE) THEN $
        scheffe_spc = dblarr(5, dim_method)
      IF (ARG_PRESENT(tukey) EQ TRUE) THEN $
        tukey_spc = dblarr(5, dim_method)
      IF (ARG_PRESENT(group_means) EQ TRUE) THEN $
        grp_means_spc = dblarr(n_groups)
      IF (ARG_PRESENT(group_std_dev) EQ TRUE) THEN $
        grp_std_dev_spc = dblarr(n_groups)
   END ELSE BEGIN
      result = float(0.0)
      ; 
      ; Input 
      y_cvt = float(y)
      IF (KEYWORD_SET(confidence) EQ TRUE) THEN $
        confidence_cvt = (float(confidence))(0)
      ;
      ; Output keywords.
      ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table_spc = fltarr(15)
      IF (ARG_PRESENT(bonferroni) EQ TRUE) THEN $
        bonferroni_spc = fltarr(5, dim_method)
      IF (ARG_PRESENT(dunn_sidak) EQ TRUE) THEN $
        dunn_sidak_spc = fltarr(5, dim_method)
      IF (ARG_PRESENT(one_at_a_time) EQ TRUE) THEN $
        one_a_time_spc = fltarr(5, dim_method)
      IF (ARG_PRESENT(scheffe) EQ TRUE) THEN $
        scheffe_spc = fltarr(5, dim_method)
      IF (ARG_PRESENT(tukey) EQ TRUE) THEN $
        tukey_spc = fltarr(5, dim_method)
      IF (ARG_PRESENT(group_means) EQ TRUE) THEN $
        grp_means_spc = fltarr(n_groups)
      IF (ARG_PRESENT(group_std_dev) EQ TRUE) THEN $
        grp_std_dev_spc = fltarr(n_groups)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_101, type, err_status, n_cvt, y_cvt, n_groups, $
                       confidence_cvt, $
                       anova_table_spc, $
                       bonferroni_spc, $
                       dunn_sidak_spc, $
                       grp_counts_spc, $
                       grp_means_spc, $
                       grp_std_dev_spc, $
                       one_a_time_spc, $
                       scheffe_spc, $
                       tukey_spc, $
                       result
   ;
   ; Now copy over all output keywords results.
   ;
      IF (ARG_PRESENT(anova_table) EQ TRUE) THEN $
        anova_table = anova_table_spc
      IF (ARG_PRESENT(bonferroni) EQ TRUE) THEN $
        bonferroni = transpose(bonferroni_spc) ;NOTE transpose()
      IF (ARG_PRESENT(dunn_sidak) EQ TRUE) THEN $
        dunn_sidak = transpose(dunn_sidak_spc) ;NOTE transpose()
      IF (ARG_PRESENT(group_counts) EQ TRUE) THEN $
        group_counts = grp_counts_spc
      IF (ARG_PRESENT(group_means) EQ TRUE) THEN $
        group_means = grp_means_spc
      IF (ARG_PRESENT(group_std_dev) EQ TRUE) THEN $
        group_std_dev = grp_std_dev_spc
      IF (ARG_PRESENT(one_at_a_time) EQ TRUE) THEN $
        one_at_a_time = transpose(one_a_time_spc) ;NOTE transpose()
      IF (ARG_PRESENT(scheffe) EQ TRUE) THEN $
        scheffe = transpose(scheffe_spc) ;NOTE transpose()
      IF (ARG_PRESENT(tukey) EQ TRUE) THEN $
        tukey = transpose(tukey_spc) ;NOTE transpose()
   ;
   ; Return.
   ;
   RETURN, result
END
