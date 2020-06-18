pro mmm, sky_vector, skymde, sigma , skew, HIGHBAD = highbad, DEBUG = debug
;+
; NAME:
;	MMM
; PURPOSE: 
;	Estimate the sky background in a stellar contaminated field.
; EXPLANATION:  
;	MMM assumes that contaminated sky pixel values overwhelmingly display 
;	POSITIVE departures from the true value.  Adapted from DAOPHOT 
;	routine of the same name.
;
; CALLING SEQUENCE:
;	MMM, sky, [ skymde, sigma, skew, HIGHBAD = , DEBUG =  ]
;
; INPUTS:
;	SKY - Array or Vector containing sky values.  This version of
;		MMM does not require SKY to be sorted beforehand.  SKY
;		is unaltered by this program.
;
; OPTIONAL OUTPUTS:
;	SKYMDE - Scalar giving estimated mode of the sky values
;	SIGMA -  Scalar giving standard deviation of the peak in the sky
;		histogram.  If for some reason it is impossible to derive
;		SKYMDE, then SIGMA = -1.0
;	SKEW -   Scalar giving skewness of the peak in the sky histogram
;		If no output variables are supplied or if /DEBUG is set
;		then the values of SKYMDE, SIGMA and SKEW will be printed.
;
; OPTIONAL KEYWORD INPUTS:
;	HIGHBAD - scalar value of the high "bad" pixel level (e.g. cosmic rays)
;		If not supplied, then there is assumed to be no high bad
;		pixels.
;	DEBUG - If this keyword is set and non-zero, then additional information
;		is displayed at the terminal.
;
; RESTRICTIONS:
;	Program assumes that low "bad" pixels (e.g. bad CCD columns) have
;	already been deleted from the SKY vector.
;
; METHOD:
;	The algorithm used by MMM consists of roughly two parts:
;	(1) The average and sigma of the sky pixels is computed.   These values
;	are used to eliminate outliers, i.e. values with a low probability
;	given a Gaussian with specified average and sigma.   The average
;	and sigma are then recomputed and the process repeated up to 20
;	iterations:
;	(2) The amount of contamination by stars is estimated by comparing the 
;	mean and median of the remaining sky pixels.   If the mean is larger
;	than the median then the true sky value is estimated by
;	3*median - 2*mean
;         
; REVISION HISTORY:
;	Adapted to IDL from 1986 version of DAOPHOT in STSDAS, 
;	W. Landsman, STX Feb 1987
;	Adapted for IDL Version 2, J. Isensee, STX, Sept 1990
;	Added HIGHBAD keyword, W. Landsman January, 1991
;	Fixed occasional problem with integer inputs    W. Landsman  Feb, 1994
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 On_error,2               ;Return to caller
 mxiter = 30              ;Maximum number of iterations allowed
 minsky = 20              ;Minimum number of legal sky elements
 if N_params() EQ 0 then begin		
	print,"Syntax:  MMM, sky, skymde, sigma, skew, [HIGHBAD = , /DEBUG] "
	return
 endif

 nsky = N_elements( sky_vector )            ;Get number of sky elements     
 if nsky LE minsky then message, $   
    'ERROR -Input vector must contain at least '+strtrim(minsky,2)+' elements'

 nlast = nsky-1                        ;Subscript of last pixel in SKY array
 if keyword_set(DEBUG) then $
     message,'Processing '+strtrim(nsky,2) + ' element array',/INF

 sky = sky_vector[ sort( sky_vector ) ]    ;Sort SKY in ascending values

 skymid = 0.5*sky[(nsky-1)/2] + 0.5*sky[nsky/2] ;Median value of all sky values  
       
 cut1 = min( [skymid-sky[0],sky[nsky-1] - skymid] ) 
 if N_elements(highbad) EQ 1 then cut1 = min( [cut1,highbad - skymid] )
 cut2 = skymid + cut1
 cut1 = skymid - cut1
         
; Select the pixels between Cut1 and Cut2

 good = where( (sky LE cut2) and (sky GE cut1) )   
 delta = sky[good] - skymid  ;Subtract median to improve arithmetic accuracy
 sum = double(total(delta))                     
 sumsq = double(total(delta^2))

 maxmum = max( good,MIN=minmum )   ;Highest value accepted at upper end of vector
 minmum = minmum -1                ;Highest value accepted at lower end of vector

; Compute mean and sigma (from the first pass).

 skymed = 0.5*(sky[(minmum+maxmum+1)/2] + sky[(minmum+maxmum)/2  +1]) ;median 
 skymn = sum/(maxmum-minmum)	 			;mean       
 sigma = sqrt(sumsq/(maxmum-minmum)-skymn^2)             ;sigma          
 skymn = skymn + skymid 	;Add median which was subtracted off earlier 

 skymde = skymn                        

;    If mean is less than the mode, then the contamination is slight, and the
;    mean value is what we really want.
 if (skymed LT skymn) then skymde = 3.*skymed - 2.*skymn

; Rejection and recomputation loop:

 niter = 0                             
START_LOOP:
   niter = niter + 1                     
   if ( niter GT mxiter ) then begin
      sigma=-1.0 &  skew = 0.0   
      message,'ERROR - Too many ('+strtrim(mxiter,2) + ') iterations,' + $
               ' unable to compute sky',/CON
      return
   endif

   if ( maxmum-minmum LT minsky ) then begin 	;Error?	

      sigma = -1.0 &  skew = 0.0   
      message,'ERROR - Too few ('+strtrim(maxmum-minmum,2) +  $
                 ') valid sky elements, unable to compute sky'
   endif 

; Compute Chavenet rejection criterion.

    r = alog10( float( maxmum-minmum ) )      
    r = max( [ 2., ( -0.1042*r + 1.1695)*r + 0.8895 ] )

; Compute rejection limits (symmetric about the current mode).

    cut = r*sigma + 0.5*abs(skymn-skymde)   
    cut = max( [ 1.5, cut] )
    cut1 = skymde - cut   &    cut2 = skymde + cut
;       
    redo = 0B
;
    newmin = minmum		
    tst_min = sky[newmin+1] GE cut1	 ;Is MINMUM+1 above current CUT?
    done = (newmin EQ -1) and tst_min    ;Are we at first pixel of SKY?
    if not done then  $
        done =  (sky[newmin>0] LT cut1) and tst_min
    if not done then begin
	istep = 1 - 2*fix(tst_min)
	repeat begin
        	newmin = newmin + istep
	        done = (newmin EQ -1)
	        if not done then $
		    done = (sky[newmin] LE cut1) and (sky[newmin+1] GE cut1)
	endrep until done
	if tst_min then delta = sky[newmin+1:minmum] - skymid $
	           else delta = sky[minmum+1:newmin] - skymid
        sum = sum - istep*total(delta)
	sumsq = sumsq - istep*total(delta^2)
	redo = 1b
	minmum = newmin
     endif
;       
   newmax = maxmum
   tst_max = sky[maxmum] LE cut2	   ;Is current maximum below upper cut?
   done = (maxmum EQ nlast) and tst_max	   ;Are we at last pixel of SKY array?
   if not done then $ 	
       done = ( tst_max ) and (sky[(maxmum+1)<nlast] GT cut2) 
    if not done then begin		   ;Keep incrementing NEWMAX
       istep = -1 + 2*fix(tst_max)	   ;Increment up or down?
       Repeat begin
          newmax = newmax + istep
	  done = (newmax EQ nlast)
	  if not done then $
		done = ( sky[newmax] LE cut2 ) and ( sky[newmax+1] GE cut2 )
       endrep until done
       if tst_max then delta = sky[maxmum+1:newmax] - skymid $
	       else delta = sky[newmax+1:maxmum] - skymid
       sum = sum + istep*total(delta)
       sumsq = sumsq + istep*total(delta^2)
       redo = 1b
       maxmum = newmax
    endif
;       
; Compute mean and sigma (from this pass).
;
   nsky = maxmum - minmum
   skymn = sum/nsky       
   sigma = float( sqrt(sumsq/nsky - skymn^2) )
   skymn = skymn + skymid                 

;  Determine a more robust median by averaging the central 10 or 11 values

   i1 = (minmum + maxmum - 8)/2    &    i2 = (minmum +maxmum + 11)/2
   skymed = total(sky[i1:i2])/(i2-i1+1)

;  If the mean is less than the median, then the problem of contamination
;  is slight, and the mean is what we really want.

   skymde = float(skymn) 
   if ( skymed LT skymn ) then skymde = float(3.*skymed-2.*skymn)
;       
   if redo then goto, START_LOOP
;       
 skew = float( (skymn-skymde)/max([1.,sigma]) )

 if keyword_set(DEBUG) or ( N_params() EQ 1 ) then begin
        print, 'Number of unrejected sky elements: ', strtrim(nsky,2), $
              '    Number of iterations: ', strtrim(niter,2)
	print, 'MODE, SIGMA, SKEW of sky vector:', skymde, sigma, skew   
 endif
 
 return
 end
