FUNCTION DCT, Data, C, INVERSE=inverse, T=T
;+
;Y=dct(X) is the discrete cosine transform of a vector
;X of length N or of a square array X of size NxN.
;
;Z=dct(Y,/INVERSE) computes the inverse
;
;If you are going to do the same transform or the inverse
;many times then it is inefficient to compute C each time.
;You can save C the first time by using keyword T.
;Y1=DCT(X1,T=C) constructs and returns the array C
;while doing the forward transform. Now, with a new
;array X2
;Y2=DCT(X2,C) uses array C without computing it. To do
;the inverse simply
;Z2=DCT(Y2,C,/INVERSE) etc.
;
;H. Rhody
;April 20, 2005
;-

;Find the array size
X=double(Data)
sx=size(X)

Norm=  1d
if keyword_set(Inverse) then X=X*Norm

;Case X is a column vector.
IF sx[0] EQ 2 AND sx[1] EQ 1 THEN BEGIN
	N=sx[2]
	SEL='COL'
ENDIF ELSE BEGIN
	IF sx[0] EQ 2 AND sx[1] NE sx[2] THEN MESSAGE,'Array must be a square array or a vector'
	N=sx[1]
	SEL='ARRAY'
ENDELSE

;Case X is a row vector.
IF sx[0] EQ 1 THEN BEGIN
	N=n_elements(X)
	SEL='ROW'
END

IF n_params() EQ 1 THEN BEGIN ;Compute the array C
;Compute the DCT matrix C
s=dindgen(n)#replicate(1d,n)
u=transpose(s)
n = double(n)
; C=cos((2d*s+1)*u*!dpi/2/(double(n)*sqrt(2d)/sqrt(double(n))
C=cos((2d*s+1)*u*!dpi/2/n)*sqrt(2d)/sqrt(n)

C[*,0]=c[*,0] / sqrt(2d)
ENDIF
T=C

CASE SEL OF
'COL':	BEGIN
	;For an inverse transform compute Transpose(C)##X
	IF KEYWORD_SET(inverse) THEN Return,Transpose(C)##X 

	;For a forward transform compute C##X
	Return,C##X  / Norm
	END

'ROW':	BEGIN
	;For an inverse transform compute X##C
	IF KEYWORD_SET(inverse) THEN Return,X##C 

	;For a forward transform compute X##Transpose(C)
	Return,X##Transpose(C) / Norm
	END

'ARRAY':	BEGIN
	;For an inverse transform compute Transpose(C)##X##C
	IF KEYWORD_SET(inverse) THEN Return,Transpose(C)##X##C

	;For a forward transform compute C##X##Transpose(C)
	Return,C##X##Transpose(C) /  Norm
	END

ELSE: MESSAGE,'Input must be a square array or a vector'

ENDCASE

END
