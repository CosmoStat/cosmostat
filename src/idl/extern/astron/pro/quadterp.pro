PRO quadterp, xtab, ytab, xint, yint, MISSING = MISSING
;+
; NAME:
;       QUADTERP     
; PURPOSE:
;       Quadratic interpolation of X,Y vectors onto a new X grid
; EXPLANATION:
;       Quadratically interpolate (3 point Lagrangian) a function Y = f(X)
;       at specified grid points.  Use LINTERP for linear interpolation
;
; CALLING SEQUENCE:
;       QUADTERP, Xtab, Ytab, Xint, Yint, [ MISSING = ]
;
; INPUT: 
;       Xtab - Vector (X TABle) containing the current independent variable 
;               Must be either monotonic increasing or decreasing
;       Ytab - Vector (Y TABle) containing the dependent variable defined
;               at each of the points of XTAB.
;       Xint - Scalar or vector giving the values of X for which interpolated 
;               Y values are sought
;
; OUTPUT: 
;       Yint - Interpolated value(s) of Y, same number of points as Xint
;
; OPTIONAL INPUT KEYWORD:
;       MISSING - Scalar specifying Yint value(s) to be assigned, when Xint
;               value(s) are outside of the range of Xtab.     Default is to
;               truncate the out of range Yint value(s) to the nearest value 
;               of Ytab.   See the help for the INTERPOLATE function.
; METHOD:
;       3-point Lagrangian interpolation.  The average of the two quadratics 
;       derived from the four nearest  points is returned in YTAB.   A single
;       quadratic is used near the end points.   The procedure TABINV is used 
;       to locate center point of the interpolation.
;
; RESTRICTIONS:
;       Unless MISSING keyword is set, points outside the range of Xtab in 
;       which valid quadratics can be computed are returned at the value 
;       of the nearest end point of Ytab (i.e. Ytab(0) and Ytab(NPTS-1) ).
;
; EXAMPLE:
;       A spectrum has been defined using a wavelength vector WAVE and a
;       flux vector FLUX.  Interpolate onto a new wavelength grid, e.g. 
;
;       IDL> wgrid = [1540.,1541.,1542.,1543.,1544.,1545.]
;       IDL> quadterp, wave, flux, wgrid, fgrid 
;     
;       FGRID will be a 5 element vector containing the quadratically
;       interpolated values of FLUX at the wavelengths given in WGRID.
;
;  EXTERNAL ROUTINES:
;       TABINV, ZPARCHECK, DATATYPE(), ISARRAY()
;  NOTES:
;       Users of IDL V5.3 can use a faster version of quadterp.pro available at
;       http://idlastro.gsfc.nasa.gov/ftp/v53/ which uses the intrinsic 
;       VALUE_LOCATE() function instead of TABINV
;  REVISION HISTORY:
;       31 October 1986 by B. Boothman, adapted from the IUE RDAF
;       12 December 1988 J. Murthy, corrected error in Xint
;       September 1992, W. Landsman, fixed problem with double precision
;       August 1993, W. Landsman, added MISSING keyword
;       June, 1995, W. Landsman, use single quadratic near end points
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Fix occasional problem with integer X table,  
;       YINT is a scalar if XINT is a scalar   W. Landsman Dec 1999
;-
 On_error,0

 if N_params() LT 4 then begin
     print,'Syntax - QUADTERP, xtab, ytab, xint, yint, [ MISSING = ]'
     return
 endif

 zparcheck,'QUADTERP',xtab,1,[1,2,3,4,5],1,'Independent (X) vector'
 zparcheck,'QUADTERP',ytab,2,[1,2,3,4,5],1,'Dependent (Y) vector'

 npts = min( [N_elements(xtab), N_elements(ytab) ] )
 m = n_elements(xint)

 if datatype(xtab) NE 'DOU' then xt = float(xtab) else xt = xtab
 if datatype(xint) NE 'DOU' then yint = fltarr(m) else yint = dblarr(m)  

 decreasing =  (xt[npts-1] - xt[0]) lt 0 

 if npts LT 3 then  $
     message,' ERROR - At least 3 points required for quadratic interpolation'

 icen = where( (xint-xt[1])*(xt[npts-2]-xint) GT 0, no )

 if no gt 0 then begin
        x = xint[icen]

; Determine index of data-points from which interpolation is made 

        tabinv, xtab, x, index 
        index = long( index < (npts-3) )

; Extract XTAB-YTAB values  

        x3 = xt[index+2]   &   y3 = ytab[index+2]
        x2 = xt[index+1]   &   y2 = ytab[index+1]
        x1 = xt[index]     &   y1 = ytab[index]
        x0 = xt[index-1]   &   y0 = ytab[index-1]

        g0 = x-x0 & g1 = x-x1 & g2 = x-x2 & g3 = x-x3
        x01 = x0-x1 & x02 = x0-x2 & x12 = x1-x2 & x23 = x2-x3 & x13 = x1-x3

; Use average of two quadratic interpolations for interior points 

        a= g1 * g2 /(x01*x02)
        b=(-g0 / x01 + g3 / x13) * g2 / x12
        c=( g0 / x02 - g3 / x23) * g1 / x12
        d= (g1*g2)/(x13*x23)

        yint[icen] = (y0*a + y1*b + y2*c + y3*d) / 2. 
 endif

; Use single quadratic interpolation near end points 

 if decreasing then begin
        ifirst = where (xint ge xtab[1],no)
 endif else begin
        ifirst = where (xint le xtab[1],no)
 endelse

 if  no gt 0 then begin
        if decreasing then xf = xint[ifirst] < xtab[0] else $
                           xf = xint[ifirst] > xtab[0]
        
        x0 = xt[0] & x1 = xt[1] & x2 = xt[2]
        y0 = ytab[0] & y1 = ytab[1] & y2 = ytab[2]
        g0 = xf-x0 & g1 = xf-x1 & g2 = xf-x2
        x01 = x0-x1 & x02 = x0-x2 & x12 = x1-x2
        
        a = g1*g2/(x01*x02)
        b = -g0*g2/(x01*x12)
        c = g0*g1/(x02*x12)
        
        yint[ifirst] = a*y0 + b*y1 + c*y2
 endif

 if decreasing then begin
        ilast = where (xint le xtab[npts-2],no)
 endif else begin
        ilast = where (xint ge xtab[npts-2],no)
 endelse

 if no gt 0 then begin
        if decreasing then xl = xint[ilast] > xtab[npts-1] else $
                           xl = xint[ilast] < xtab[npts-1]

        x0 = xt[npts-3] & x1 = xt[npts-2] & x2 = xt[npts-1]
        y0 = ytab[npts-3] & y1 = ytab[npts-2] & y2 = ytab[npts-1]
        g0 = xl-x0 & g1 = xl-x1 & g2 = xl-x2
        x01 = x0-x1 & x02 = x0-x2 & x12 = x1-x2
        
        a = g1*g2/(x01*x02)
        b = -g0*g2/(x01*x12)
        c = g0*g1/(x02*x12)
        
        yint[ilast] = a*y0 + b*y1 + c*y2
 endif

; Any points outside of table range?

 if N_elements(missing) EQ 1 then begin
        Xmin = min( [ Xtab[0],Xtab[npts-1] ], max = Xmax)
        bad = where( (Xint LT Xmin) or (Xint GT Xmax ), Nbad)
        if Nbad GT 0 then Yint[bad] = missing
 endif

 if not ISARRAY(xint) then yint = yint[0]

 return
 end
