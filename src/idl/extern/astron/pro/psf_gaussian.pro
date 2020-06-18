function psf_gaussian, parameters, NPIXEL=npix, NDIMENSION=ndim, FWHM=fwhm,  $
					CENTROID=cntrd, ST_DEV=st_dev,  $
					XY_CORREL=xy_corr, NORMALIZE=normalize
;+
; NAME:
;	PSF_GAUSSIAN
;
; PURPOSE:
;	Create a 1-d, 2-d, or 3-d Gaussian with specified FWHM, center 
; EXPLANATION:
;	Return a point spread function having Gaussian profiles,
;	as either a 1D vector, a 2D image, or 3D volumetric-data.
;
; CALLING SEQUENCE:
;	psf = psf_Gaussian( NPIXEL=, FWHM= , [/NORMALIZE, /ST_DEV,  )
; or:
;	psf = psf_Gaussian( parameters, NPIXEL =  )
;
; REQUIRED INPUT KEYWORD:
;	NPIXEL = number pixels for each dimension, specify as an array,
;		or just one number to make all sizes equal.
;
; OPTIONAL KEYWORDS:
;
;	NDIMEN = dimension of result: 1 (vector), 2 (image), or 3 (volume),
;		default = 2 (an image result).
;
;	FWHM = the desired Full-Width Half-Max (pixels) in each dimension,
;		specify as an array, or single number to make all the same.
;
;	CENTROID = pixels numbers of PSF maximum ( 0.5 is center of a pixel ),
;		default is exact center of requested vector/image/volume.
;
;	STDEV = optional way to specify width by standard deviation param.
;
;	XY_CORREL = scalar between 0 and 1 specifying correlation coefficient
;		Use this keyword, for example, to specify an elliptical 
;		gaussian oriented at an angle to the X,Y axis
;
;	/NORMALIZE causes resulting PSF to be normalized so Total( psf ) = 1.
;
; INPUTS (optional):
;
;	parameters = an NDIMEN by 3 array giving for each dimension:
;			[ maxval, center, stdev ],  overrides other keywords.
;
; EXAMPLE:
;	Create a 31 x 31 array containing a normalized centered gaussian 
;	with an X FWHM = 4.3 and a Y FWHM = 3.6
;
;	IDL> array = PSF_GAUSSIAN( Npixel=31, FWHM=[4.3,3.6], /NORMAL
;
; EXTERNAL CALLS:
;	function Gaussian
;
; HISTORY:
;	Written, Frank Varosi NASA/GSFC 1991.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
	On_error,2

	if (N_params() LT 1 ) and not keyword_set( FWHM) then begin
		print,'Syntax - psf = PSF_GAUSSIAN( parameters, NPIXEL = )
		print, $
	'or       psf = PSF_GAUSSIAN( FWHM = ,STDEV = ,NPIXEL = ,[CENTROID = ])'
		return, -1
	endif

	sp = size( parameters )

	if (sp[0] GE 1) then begin
		ndim = sp[0]
		factor = total( parameters[*,0] )/float( ndim )
		cntrd = parameters[*,1]
		st_dev = parameters[*,2]
	   endif

	if N_elements( ndim ) NE 1 then ndim=2
	ndim = ndim>1

	if N_elements( npix ) LE 0 then begin
		message,"must specify size of result with NPIX=",/INFO
		return,(-1)
	  endif else if N_elements( npix ) LT ndim then $
			npix = replicate( npix[0], ndim )

	if (N_elements( cntrd ) LT ndim) AND (N_elements( cntrd ) GT 0) then $
			cntrd = replicate( cntrd[0], ndim )

	if N_elements( cntrd ) LE 0 then cntrd=(npix-1)/2. else cntrd=cntrd-0.5
	if N_elements( fwhm ) GT 0 then st_dev = fwhm/( 2* sqrt( 2* aLog(2) ) )

	if N_elements( st_dev ) LE 0 then begin
		message,"must specify ST_DEV= or FWHM=",/INFO
		return,(-1)
	  endif

	if N_elements( st_dev ) LT ndim then $
			st_dev = replicate( st_dev[0], ndim )
	sigfac = 1 / (2. * st_dev^2 )

	CASE ndim OF

	1: BEGIN
		x = findgen( npix[0] ) - cntrd[0]
		psf = gaussian( x, [1,0,st_dev] )
	     END

	2: BEGIN
		psf = make_array( DIM=npix[0:ndim-1], /FLOAT )
		x = findgen( npix[0] ) - cntrd[0]
		y = findgen( npix[1] ) - cntrd[1]

		if N_elements( xy_corr ) EQ 1 then begin
			y2 = sigfac[1] * y^2
			x1 = sigfac[0] * x
			yc = y * ( xy_corr/(st_dev[0]*st_dev[1]) )
			for j=0,npix[1]-1 do begin
				zz = x * (yc[j] + x1) + y2[j]
				w = where( zz LT 86, nw )
				if (nw GT 0) then psf[w,j] = exp( -zz[w] )
			  endfor
		  endif else begin
			psfx = gaussian( x, [ 1, 0, st_dev[0] ] )
			psfy = gaussian( y, [ 1, 0, st_dev[1] ] )
			for j=0,npix[1]-1 do psf[0,j] = psfx * psfy[j]
		   endelse
	     END

	3: BEGIN
		psf = make_array( DIM=npix[0:ndim-1], /FLOAT )
		x = findgen( npix[0] ) - cntrd[0]
		y = findgen( npix[1] ) - cntrd[1]
		z = findgen( npix[2] ) - cntrd[2]
		psfx = gaussian( x, [ 1, 0, st_dev[0] ] )
		psfy = gaussian( y, [ 1, 0, st_dev[1] ] )
		psfz = gaussian( z, [ 1, 0, st_dev[2] ] )
		for k=0,npix[2]-1 do begin
		    for j=0,npix[1]-1 do psf[0,j,k] = psfx * psfy[j] * psfz[k]
		 endfor
	     END

	ENDCASE

	if keyword_set( normalize ) then return, psf/total( psf )

	if N_elements( factor ) EQ 1 then begin
		if (factor NE 1) then return,factor*psf else return,psf
	   endif else return, psf
end
