FUNCTION POISSON_IMAGE, IMAGE, SEED, OUTPUT_KIND=OUTPUT_KIND

;+
; NAME:
;	POISSON_IMAGE
;
; PURPOSE:
;	Add Poisson noise to an array.
;
; CALLING SEQUENCE:
;	RESULT = POISSON_IMAGE ( IMAGE[, SEED] )
;
; INPUT:
;	IMAGE: A numeric array (byte, integer, long, or float) of arbitrary
;	dimensionality.  This is the array of values around which values in the
;	result will be Poisson-distributed.
;
; OPTIONAL INPUTS:
;	SEED: A longword seed for the random number generator.  If this is not
;	supplied, the value -123456789L is used for generating the first random
;	value.
;	
; KEYWORD PARAMETER:
;	OUTPUT_KIND: The data type of the output array, that is byte, integer,
;	longword, or float.  The words "byte", "int", "integer","long",
;	"longword", and "float", in upper or lower case, are accepted, as are
;	the numeric IDL values 1,2,3,4 for byte, integer, longword, and float.
;
; OUTPUT:
;	POISSON_IMAGE returns a copy of the input array, Poisson noise added.
;
; DEPENDENCIES:
;	This procedure calls a user-contributed IDL procedure:
;		STILL_BUSY
;
; RESTRICTIONS:
;	Negative input values are mapped to a result of 0.  A call is made to
;	the routine STILL_BUSY to provide messages once a minute during the
;	processing of large arrays.  You can eliminate that line of code if you
;	want to.
;
; PROCEDURE:
;	Create an image with values Poisson-distributed around the mean values
;	IMAGE, using Knuth's "Algorithm Q".  (Donald E. Knuth, The Art of
;	Computer Programming, Volume 2, "Seminumerical Algorithms", Addison-
;	Wesley (1969), page 117.  This routine IS NOT VECTORIZED, AND SHOULD RUN
;	SLOWLY.  A deft IDL'er could probably vectorize the algorithm, and
;	anyone who does so is entitled to a gold star.  Where the gaussian and
;	Poisson distributions are essentially identical (mean value > 50) a
;	normal, that is gaussian, distribution is used.
;
; EXAMPLE:
;	Here is how you can create a 100x100 array of values Poisson-distributed
;	around the mean value 5.0, and check the empirical probability against
;	the Poisson distribution:
;
;		n = 100
;		a = replicate(5.0,n,n)
;		b = poisson_image ( a, out="byte" )
;		tvscl, b
;		print, stdev(b,mean), mean, sqrt(5.0)
;		h = histogram ( b )
;		prob = float(h)/n_elements(b)
;		probi = exp(-5.0)
;		for i=0,10 do begin print, i, prob(i), probi & $
;			probi=probi*5.0/(i+1)
;
; MODIFICATION HISTORY:
; 	Written by:	James Hamill
;			Siemens Medical Systems
;			2501 N. Barrington Rd.
;			Hoffman Estates, IL  60195-7372
;			(708)304-7760
;			hamill@sgi.siemens.com
;			February, 1992
;-



on_error, 1

default_seed = -123456789L
BIG_ENOUGH = 50	; bigger than this, use normal approximation

undefined_type = 0
byte_type = 1
int_type = 2
long_type = 3
float_type = 4
string_type = 7


;-------------------------------------------------------------------------------
;  First part of the code: decide what type output to make.
;-------------------------------------------------------------------------------


sizi = size(image)
dims = sizi(1:sizi[0])
old_type = sizi(sizi[0]+1)

if n_elements(output_type) eq 0 then begin
			; if keyword doesn't specify the type ...

  case old_type of
    byte_type:	output_type = int_type
    int_type:	output_type = long_type
    else:	output_type = float_type
  endcase

endif else begin	; if keyword is set then follow the instructions


  sizok = size(output_kind)	; text or numeric parameter?
  ok_type = sizok(sizok[0]+1)


  case ok_type of

    int_type: output_type = output_kind

    long_type: output_type = output_kind

    string_type: begin

      case strupcase(output_kind) of
	"BYTE":		output_type = byte_type
	"INT":		output_type = int_type
	"INTEGER":	output_type = int_type
	"LONG":		output_type = long_type
	"LONGWORD":	output_type = long_type
	else: MESSAGE,"Invalid output type"
      endcase

    end

    else: MESSAGE, "Invalid output_kind keyword."

  endcase

  if (output_type lt byte_type) or (output_type gt float_type) then $
	message, "Invalid output_kind keyword"
		; thus we trap nonsense like string, complex, or structure


endelse

result = make_array(dimension=dims,type=output_type)		; (zeroed)


;-------------------------------------------------------------------------------
;  Check the random number generator seed, set it if necessary.
;-------------------------------------------------------------------------------


if n_params(dummy) ge 2 then begin
  sizs = size(seed)
  seed_type = sizs(sizs[0]+1)
  if seed_type eq undefined_type then begin
    seed = default_seed			; undefined parameter: set it
  endif else begin			; check that seed is a long scalar
    if (seed_type ne long_type) or (sizs[0] ne 0) then $
	MESSAGE,"Invalid seed."
  endelse
endif else seed=default_seed


;-------------------------------------------------------------------------------
;  Consider the small ones: use the Poisson distribution for these.  This is
;  Knuth's algorithm Q.
;-------------------------------------------------------------------------------


ROI = where( (image gt 0) and (image lt BIG_ENOUGH), count)
if count gt 0 then begin

  for i=0L,n_elements(ROI)-1 do begin

    j = ROI[i]
    p = exp(-double(image[j]))
    q = 1D0
    n = -1

    while q ge p do begin
      n = n + 1
      u = randomu(seed)
      q = q*u
    endwhile

    result[j] = n > 0
    if i mod 10000L eq 0 then print,$
	"Poisson distribution ... working on element #"+strtrim(i)+"."

  endfor

endif


;-------------------------------------------------------------------------------
;  Consider the region in which the normal distribution can be used.  Round to
;  the nearest whole number.
;-------------------------------------------------------------------------------


ROI = where ( image ge BIG_ENOUGH, count )

if count ne 0 then begin

  n_values = n_elements(ROI)
  result[ROI] = fix(0.5 + image[ROI] + sqrt(image[ROI])*randomn(seed,n_values))

endif



RETURN, result


END

