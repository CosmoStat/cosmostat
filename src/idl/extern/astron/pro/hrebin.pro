 pro hrebin, oldim, oldhd, newim, newhd, newx, newy, $
            SAMPLE=sample, OUTSIZE = outsize, ERRMSG = errmsg
;+
; NAME:
;    HREBIN
; PURPOSE:
;    Expand or contract a FITS image using (F)REBIN and update the header 
; EXPLANATION:
;    If output size is a multiple of input size then REBIN is used, else
;    FREBIN is used.     User can either overwrite the input array,
;    or write to new variables.
;
; CALLING SEQUENCE:
;    HREBIN, oldhd        ;Special calling sequence to just update header
;    HREBIN, oldim, oldhd, [ newim, newhd, newx, newy, OUTSIZE = ,/SAMPLE, 
;                            ERRMSG =  ]
;
; INPUTS:
;    OLDIM - the original image array
;    OLDHD - the original image FITS header, string array
;
; OPTIONAL INPUTS:
;    NEWX - size of the new image in the X direction, integer scalar
;    NEWY - size of the new image in the Y direction, integer scalar
;            HREBIN will prompt for NEWX and NEWY if not supplied
;
; OPTIONAL OUTPUTS:
;    NEWIM - the image after expansion or contraction with REBIN
;    NEWHD - header for newim containing updated astrometry info
;            If output parameters are not supplied, the program will modify
;            the input parameters OLDIM and OLDHD to contain the new array and 
;            updated header.
;
; OPTIONAL INPUT KEYWORDS:
;    /SAMPLE - Expansion or contraction is done using REBIN which uses 
;              bilinear interpolation when magnifying and boxaveraging when 
;              minifying.   If the SAMPLE keyword is supplied and non-zero, 
;              then nearest neighbor sampling is used in both cases.   Keyword
;              has no effect when output size is not a multiple of input size.
;
;    OUTSIZE - Two element integer vector which can be used instead of the
;             NEWX and NEWY parameters to specify the output image dimensions
;
; OPTIONAL KEYWORD OUTPUT:
;       ERRMSG - If this keyword is supplied, then any error mesasges will be
;               returned to the user in this parameter rather than depending on
;               on the MESSAGE routine in IDL.   If no errors are encountered
;               then a null string is returned.               
; PROCEDURE:
;     The parameters BSCALE, NAXIS1, NAXIS2, CRPIX1, and CRPIX2 and the CD 
;     (or CDELT) parameters are updated for the new FITS header.
;
; EXAMPLE:
;     Compress a 2048 x 2048 image array IM, with header FITS HDR, to a 
;     724 x 724 array.   Overwrite the input variables with the compressed 
;     image and header.
;
;     IDL> hrebin, im, hdr, OUT = [724, 724]
;
; PROCEDURES USED:
;     CHECK_FITS, EXTAST, FREBIN, GSSS_STDAST, STRN(), SXPAR(), SXADDHIST, 
;     SXADDPAR, ZPARCHECK
;
; MODIFICATION HISTORY:
;     Written, December 1990  W. Landsman, ST System Corp.
;     Update CD1_1 keywords   W. Landsman   November 1992
;     Check for a GSSS header   W. Landsman  June 1994
;     Update BSCALE even if no astrometry present   W. Landsman  May 1997
;     Converted to IDL V5.0   W. Landsman   September 1997
;     Use FREBIN to accept sizes that are not a integer multiple of the original
;         size    W. Landsman     August 1998
;     Correct for "edge" effects when expanding with REBIN W. Landsman Apr. 1999
;     Fixed initialization of header only call broken in Apr 98 change May. 1999
;     Remove reference to obsolete !ERR  W. Landsman   February 2000
;     Use double precision formatting for CD matrix W. Landsman April 2000
;- 
 On_error,2

 npar = N_params()      ;Check # of parameters
 if (npar EQ 3) or (npar EQ 5) or (npar EQ 0) then begin
     print,'Syntax - HREBIN, oldim, oldhd,[ newim, newhd, OUTSIZE=, ' + $
                           '/SAMPLE, ERRMSG= ]'
     return
 endif

 if not keyword_set(SAMPLE) then sample = 0
 save_err = arg_present(errmsg)      ;Does user want to return error messages?

; If only 1 parameter is supplied, then assume it is a FITS header

 if ( npar EQ 1 ) then begin           

        zparcheck, 'HREBIN', oldim, 1, 7, 1, 'Image header'
        oldhd = oldim
        xsize = sxpar( oldhd,'NAXIS1' )
        ysize = sxpar( oldhd,'NAXIS2' )

 endif else begin 

     check_FITS, oldim, oldhd, dimen, /NOTYPE, ERRMSG = errmsg
     if errmsg NE '' then begin
        if not save_err then message,'ERROR - ' + errmsg,/CON
        return
     endif
     if N_elements(dimen) NE 2 then begin 
           errmsg = 'Input image array must be 2-dimensional'
           if not save_err then message,'ERROR - ' + errmsg,/CON
           return
     endif
      xsize = dimen[0]  &  ysize = dimen[1]
 endelse

 if ( npar LT 6 ) then begin

    if ( N_elements(OUTSIZE) NE 2 ) then begin
    tit = !MSG_PREFIX + 'HREBIN: '
    print, tit, 'Original array size is '+ strn(xsize) + ' by ' + strn(ysize)
    read, tit + 'Enter size of new image in the X direction: ',newx
    read, tit + 'Enter size of new image in the Y direction: ',newy
  endif else begin
     newx = outsize[0]
     newy = outsize[1]
   endelse 
 
 endif

;  If an image array supplied then apply the REBIN  or FREBIN functions
; If output size is a multiple of input size then use REBIN else use FREBIN

  exact = (((xsize mod newx) EQ 0) or ((newx mod xsize) EQ 0)) and $
          (((ysize mod newy) EQ 0) or ((newy mod ysize) EQ 0))

 if npar GT 1 then begin
 if exact then begin
   if npar GT 2 then newim = rebin( oldim, newx, newy, SAMPLE=sample) $
                else oldim = rebin( oldim, newx, newy, SAMPLE=sample)
 
 endif else begin
   if npar GT 2 then newim = frebin( oldim, newx, newy) $
                else oldim = frebin( oldim, newx, newy)
 endelse 
   endif


 if ( sample GT 0 ) then type = ' Nearest Neighbor Approximation' else begin
          if ( newx LT xsize ) then type = ' Box Averaging' else $
                                    type = ' Bilinear Interpolation'
 endelse

 newhd = oldhd
 sxaddpar, newhd, 'NAXIS1', fix(newx)
 sxaddpar, newhd, 'NAXIS2', fix(newy)
 label = 'HREBIN: '+ strmid( systime(),4,20 )
 sxaddpar,newhd,'history',label + ' Original Image Size Was '+ $
         strn(xsize) +' by ' +  strn(ysize) 
 if ( npar GT 1 ) then sxaddpar,newhd,'history',label+type

; Update astrometry info if it exists

 extast, newhd, astr, noparams
 if noparams GE 0 then begin
 if strmid(astr.ctype[0],5,3) EQ 'GSS' then begin
        gsss_stdast, newhd
        extast, newhd, astr, noparams
 endif

 xratio = float(newx) / xsize   ;Expansion or contraction in X
 yratio = float(newy) / ysize   ;Expansion or contraction in Y
 pix_ratio = xratio*yratio      ;Ratio of pixel areas

; Correct the position of the reference pixel.   Note that CRPIX values are
; given in FORTRAN (first pixel is (1,1)) convention

 crpix = astr.crpix

; When expanding with REBIN with bilinear interpolation (SAMPLE = 0), edge
; effects are introduced, which require a different calculation of the updated
; CRPIX1 and CRPIX2 values.

 if (exact) and (not keyword_set(SAMPLE)) and (xratio GT 0) then $
      crpix1 = (crpix[0]-1.0)*xratio + 1.0                  else $
      crpix1 = (crpix[0]-0.5)*xratio + 0.5

 if (exact) and (not keyword_set(SAMPLE)) and (yratio GT 0) then $
      crpix2 = (crpix[1]-1.0)*yratio + 1.0                  else $
      crpix2 = (crpix[1]-0.5)*yratio + 0.5

 sxaddpar, newhd, 'CRPIX1', crpix1, FORMAT='(G14.7)'
 sxaddpar, newhd, 'CRPIX2', crpix2, FORMAT='(G14.7)'

; Scale either the CDELT parameters or the CD1_1 parameters.

 if (noparams EQ 0) or ( noparams EQ 1) then begin 

    cdelt = astr.cdelt
    sxaddpar, newhd, 'CDELT1', CDELT[0]/xratio
    sxaddpar, newhd, 'CDELT2', CDELT[1]/yratio

 endif else if noparams EQ 2 then begin

    cd = astr.cd
    sxaddpar, newhd, 'CD1_1', cd[0,0]/xratio
    sxaddpar, newhd, 'CD1_2', cd[0,1]/yratio
    sxaddpar, newhd, 'CD2_1', cd[1,0]/xratio
    sxaddpar, newhd, 'CD2_2', cd[1,1]/yratio

 endif 
 endif

; Adjust BZERO and BSCALE for new pixel size

 bscale = sxpar( oldhd, 'BSCALE')
 if (bscale NE 0) and (bscale NE 1) then $
    sxaddpar, newhd, 'BSCALE', bscale/pix_ratio, 'Calibration Factor'
 bzero = sxpar( oldhd, 'BZERO')
    if (bzero NE 0) then sxaddpar, newhd, 'BZERO', bzero/pix_ratio, $
       ' Additive Constant for Calibration'

 pixelsiz = sxpar( oldhd,'PIXELSIZ' , Count = N_pixelsiz)
 if N_pixelsiz GT 0 then sxaddpar, newhd, 'PIXELSIZ', pixelsiz/xratio

 if npar EQ 2 then oldhd = newhd else $
    if npar EQ 1 then oldim = newhd

 return
 end
