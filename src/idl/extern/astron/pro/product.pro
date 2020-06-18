	FUNCTION PRODUCT, ARRAY
;+
; NAME:
;	PRODUCT
; PURPOSE:
;	Calculates the product of all the elements of an array
; EXPLANATION:
;	PRODUCT() is the multiplicative equivalent of TOTAL().
; CALLING SEQUENCE:
;	Result = PRODUCT(ARRAY)
; INPUT PARAMETERS:
;	ARRAY	= Array of elements to multiply together.  For instance, ARRAY
;		  could contain the dimensions of another array--then
;		  PRODUCT(ARRAY) would be the total number of elements of that
;		  other array.
; OUTPUT:
;	The result of the function is the total product of all the elements of
;	ARRAY.
; OPTIONAL KEYWORD PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	The result will always be of at least floating point type.
; RESTRICTIONS:
;	ARRAY must be a numerical type.
; PROCEDURE:
;	Straightforward.
; MODIFICATION HISTORY:
;	William Thompson, Feb. 1992.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;
	ON_ERROR,2
;
;  Check the number of parameters.
;
	IF N_PARAMS() NE 1 THEN MESSAGE,'Syntax:  Result = PRODUCT(ARRAY)'
;
;  Check the type of ARRAY.
;
	SZ = SIZE(ARRAY)
	TYPE = SZ[SZ[0]+1]
	IF TYPE EQ 0 THEN MESSAGE,'ARRAY not defined'
	IF TYPE EQ 7 THEN MESSAGE,'Operation illegal with string arrays'
	IF TYPE EQ 8 THEN MESSAGE,'Operation illegal with structures'
;
;  Calculate the product.
;
	RESULT = 1.
	FOR I=0,N_ELEMENTS(ARRAY)-1 DO RESULT = RESULT * ARRAY[I]
;
	RETURN,RESULT
	END
