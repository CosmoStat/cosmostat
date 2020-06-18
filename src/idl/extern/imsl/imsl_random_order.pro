; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_random_order.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
FUNCTION imsl_random_order, ifirst, $           ;INPUT Scalar LONG
                   ilast, $                     ;INPUT Scalar LONG
                   n, $                         ;INPUT Scalar LONG
                   double=double, $             ;INPUT Scalar ON/OFF flag
                   normal=normal, $             ;INPUT Scalar ON/OFF flag
                   uniform=uniform              ;INPUT Scalar ON/OFF flag

@imsl_init.pro
   ON_ERROR, on_err_action
   ;
   ; Error checking:
   ;  - Check (0 LT IFIRST LE ILAST LE N)
   ;  - NORMAL and UNIFORM are mutually exclusive.
   ;
   nargs = n_params()
   IF (nargs NE 3) THEN $
         message, "Incorrect number of arguments."
   ifirst_cvt = (IMSL_LONG(ifirst))(0)
   ilast_cvt = (IMSL_LONG(ilast))(0)
   n_cvt = (IMSL_LONG(n))(0)

   if ((ifirst_cvt LE 0 OR ifirst_cvt GT ilast_cvt) OR ilast_cvt GT n_cvt) THEN BEGIN
	message, "IFIRST, ILAST, and N must satisfy the relationship (0 LT IFIRST LE ILAST LE N)."
   END	
   IF (KEYWORD_SET(normal) AND KEYWORD_SET(uniform)) THEN $
     message, "Keywords NORMAL and UNIFORM are mutually exclusive."

   use_normal = IMSL_0
   IF (KEYWORD_SET(normal)) then use_normal = IMSL_1

   ;
   ; Decide on what precision to use.
   ; Use the highest precision of the input argument(s).
   ;
   type = TYP_FLOAT
   IF (KEYWORD_SET(double) EQ true) THEN type = TYP_DOUBLE
   ;
   ; Floating point arguments and keywords
   IF (type EQ TYP_DOUBLE) THEN BEGIN
      ; Result vector.
      l_result = dblarr(ilast_cvt-ifirst_cvt+1)
   END ELSE BEGIN
      ; Result vector.
      l_result = fltarr(ilast_cvt-ifirst_cvt+1)
      ; 
   END
   ;
   ; Call the system function.
   ;
   err_status = 0L
   MATHSTAT_300, type, err_status, ifirst_cvt, ilast_cvt, n_cvt, use_normal, l_result
   ;
   ; Return.
   ;
   RETURN, l_result
END
   

