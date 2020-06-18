function gaussian, xi, parms, pderiv
;+
; NAME:
;       GAUSSIAN
; PURPOSE:
;       Compute the 1-d Gaussian function and optionally the derivative
; EXPLANATION:
;       Compute the 1-D Gaussian function and optionally the derivative 
;       at an array of points.
;
; CALLING SEQUENCE:
;       y = gaussian( xi, parms,[ pderiv ])
;
; INPUTS:
;       xi = array, independent variable of Gaussian function.
;
;       parms = parameters of Gaussian, 2 or 3 element array:
;               parms(0) = maximum value (factor) of Gaussian,
;               parms(1) = mean value (center) of Gaussian,
;               parms(2) = standard deviation (sigma) of Gaussian.
;               (if parms has only 2 elements then sigma taken from common).
;
; OPTIONAL OUTPUT:
;       pderiv = optional output of partial derivatives,
;               computed only if parameter is present in call.
;
;               pderiv(*,i) = partial derivative at all xi absisca values
;               with respect to parms(i), i=0,1,2.
;
;       Function returns array of Gaussian evaluated at xi.
;
; EXAMPLE:
;       Evaulate a Gaussian centered at x=0, with sigma=1, and a peak value
;       of 10 at the points 0.5 and 1.5.   Also compute the derivative
;
;       IDL> f = gaussian( [0.5,1.5], [10,0,1], DERIV )
;       ==> f= [8.825,3.25].   DERIV will be a 2 x 3 array containing the
;       numerical derivative at the two points with respect to the 3 parameters.
; 
; COMMON BLOCKS:
;       common gaussian, sigma
; HISTORY:
;       Written, Frank Varosi NASA/GSFC 1992.
;       Converted to IDL V5.0   W. Landsman   September 1997
;-
  On_error,2

  if N_params() LT 2 then begin
        print,'Syntax - y = GAUSSIAN( xi, parms,[ pderiv ])'
        print,'         parms[0] = maximum value (factor) of Gaussian'
        print,'         parms[1] = mean value (center) of Gaussian'
        print,'         parms[2] = standard deviation (sigma) of Gaussian'
        return, -1
  endif

  common gaussian, sigma

        Nparmg = N_elements( parms )
        parms = float( parms )
        if (Nparmg GE 3) then sigma = parms[2]

        z = ( xi - parms[1] )/sigma
        zz = z*z
        gauss = fltarr( N_elements( zz ) )
        w = where( zz LT 172, nw )
        if (nw GT 0) then gauss[w] = exp( -zz[w] / 2 )

        if N_params() GE 3 then begin

                pderiv = fltarr( N_elements( xi ), Nparmg )
                fsig = parms[0] / sigma

                pderiv[0,0] = gauss
                pderiv[0,1] = gauss * z * fsig

                if (Nparmg GE 3) then  pderiv[0,2] = gauss * zz * fsig
           endif

return, parms[0] * gauss
end
