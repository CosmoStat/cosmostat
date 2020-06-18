; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_radbe.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
function imsl_Chk_radial_struct, radial_fit, type
@imsl_init.pro
   ; The following checks are made on the argument that
   ; is supposed to be a radial structure.
   ; - Must be a structures.
   ; - Must have 7 tags
   ; - The 4-th tag determines if the precision:TYP_FLOAT or TYP_DOUBLE
   ; - Tagnames must match expected tagnames.
   ; - For each tag:
   ;       o  check the data type.
   ;
   err_str = "RADIAL_FIT must be a valid radial structure."
   size_rad_fit = IMSL_SIZE(radial_fit)
   IF ((size_rad_fit(0) NE 1) OR $
       (size_rad_fit(N_ELEMENTS(size_rad_fit)-2) NE 8)) THEN $
     message, err_str

   ntags = n_tags(radial_fit)
   IF (ntags NE 7) THEN message, err_str
   rad_fit_tagnames = tag_names(radial_fit)
   type = SIZE(radial_fit.centers, /TYPE)
   ; Ensure that the structure has the correct integer types
   ; by making a copy.
   result = {DIMENSION: IMSL_LONG(radial_fit.dimension), $
     NUM_CENTERS: IMSL_LONG(radial_fit.num_centers), $
     EXTRA_TERMS: IMSL_LONG(radial_fit.extra_terms), $
     CENTERS: radial_fit.centers, $
     COEFFICIENTS: FIX(radial_fit.coefficients, TYPE=type), $
     RADIAL_FCN: STRING(radial_fit.radial_fcn), $
     DELTA: FIX(radial_fit.delta, TYPE=type)}
   return, result
END



FUNCTION imsl_radbe, abscissa, $              ;INPUT 2-D array: floating point
                radial_fit                 ;INPUT radial struct

@imsl_init.pro
   ON_ERROR, on_err_action

   ;
   ; Error checking.
   ; - abscissa must be a 2-D array (size  ndim x ndata )
   ; - radial_fit must be a valid Radial-fit structure.
   ;
   nargs = n_params()
   IF (nargs NE 2) THEN message, 'Incorrect number of arguments.'
   size_abs = IMSL_SIZE(abscissa)
   IF (size_abs(0) NE 2) THEN message, 'ABSCISSA must be a 2-D array.'
   ndim = IMSL_LONG(size_abs(1))
   ndata = IMSL_LONG(size_abs(2))
   ; Check the contents of radial_fit.
   radial_fit_cvt = imsl_chk_radial_struct(radial_fit, type)
   ;
   ; Decide on what precision to use.
   ; The precision is dictated by the contents of radial_fit.
   ;
   ; Setup the parameters for the call to the system function.
   IF (radial_fit_cvt.RADIAL_FCN NE '') THEN $
     basis_cvt = radial_fit_cvt.RADIAL_FCN
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      abscissa_cvt = double(abscissa)
      result = dblarr(ndata)
   END ELSE BEGIN
      abscissa_cvt = float(abscissa)
      result = fltarr(ndata)
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_179, type, err_status, abscissa_cvt, radial_fit_cvt, $
                              ndim, ndata, $
                              basis_cvt, $
                              result

 RETURN, result
END

