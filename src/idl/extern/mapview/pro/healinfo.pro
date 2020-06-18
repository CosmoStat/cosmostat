PRO HEALINFO, RES=res_in, NSIDE=nside_in, NPIX=npix_in,  $
    SPACING=spacing_in, GET_RES=res_out, GET_NSIDE=nside_out,  $
    GET_NPIX=npix_out, GET_SPACING=spacing_out, HELP=help
;+
;   NAME:               
;      HEALINFO
;
;   PURPOSE:    
;        This procedure prints or returns values describing healpix formats of 
;        varying resolution.
;
;   INPUT ARGUMENTS:    
;       No positional arguments.
;
;   INPUT KEYWORDS:
;       RES:        Resolution (0 - 10)
;       NSIDE:      Number of divisions along edge of biggest
;                   healpix pixel (1, 2, ..., 1024)
;       NPIX:       Total number of pixels in full sky
;                   (12, 48, ..., 12582912)
;       SPACING:    Approximate mean linear spacing between
;                   pixel centers in degrees
;       /HELP:      Flag to display short help message.
;
;   OUTPUT KEYWORDS:
;       These mean the same thing as the corresponding input
;       keywords, except they are used to return results.
;       GET_RES
;       GET_NSIDE
;       GET_NPIX
;       GET_SPACING
;
;   USAGE:  This is intended as a general-purpose routine to get
;           healpix format info either from the command line or 
;           inside programs.  Entering healinfo with no keywords
;           prints out a table of information.  Specifying a
;           value for one of the keywords prints out one line of
;           the table, i.e., the best match to the value.
;
;           Examples:
;               healinfo
;               healinfo,res=9
;
;           You can also return values in variables, in
;           which case nothing is printed:
;
;               healinfo,res=9,get_npix=npix
;               print, npix
;                    3145728
;
;           The program accomodates a variety of inputs, by the 
;           following rules:
;               RES:        Input value must be exact.
;               NSIDE:      If not exact, next higher value
;                           is used, e.g., 257 --> 512
;               NPIX:       If not exact, nearest value is
;                           is used, e.g., 491521-->786432
;                           but 491520-->196608
;               SPACING:    Floating point and cannot be exact.
;                           the next smaller spacing down from
;                           the supplied value is used.
;           These rules aren't recommendations in any sense, just 
;           ways to get the program to be as informative as
;           possible under a variety of circumstances.
;   
;   MODIFICATION HISTORY:
;       Written by R. S. Hill, RITSS, 2 Aug 2001
;       Default behaviors specified and doc header written.
;           RSH, SSAI, 8 Nov 2001
;       Separated input and output keywords because behavior
;           was too confusing to be practical.  RSH, 9 Nov 2001
;      Fixed typo in pixarea_out  W. Landsman  Nov 2002
;-

ON_ERROR, 2

IF keyword_set(help) THEN BEGIN
    print, 'HEALINFO:  no positional arguments'
    print, 'Input keywords:  RES, NSIDE, NPIX, SPACING, HELP'
    print, 'Output keywords:  GET_RES, GET_NSIDE, GET_NPIX, GET_SPACING'
    RETURN
ENDIF

ap = bytarr(4)
inp = bytarr(4)

inp[0] = n_elements(res_in) GT 0
inp[1] = n_elements(nside_in) GT 0
inp[2] = n_elements(npix_in) GT 0
inp[3] = n_elements(spacing_in) GT 0

which = (where(inp, n_inp))[0]

IF n_inp GT 1 THEN message, 'Only one input keyword can be specified.'

ap[0] = arg_present(res_out) 
ap[1] = arg_present(nside_out) 
ap[2] = arg_present(npix_out) 
ap[3] = arg_present(spacing_out) 

w = where(ap GT 0, n_ap)

print_flag = n_ap LE 0
print_all = n_inp LE 0

res = indgen(11)
nside = 2L^res
npix = 12*nside^2
pixsterad = 4*!dpi/npix
pixsiderad = sqrt(pixsterad)
pixsidedeg = pixsiderad*180d0/!dpi

CASE which OF
    -1: BEGIN
        sel = 0
        END
    0:  BEGIN
        w = where(res_in[0] EQ res, nw)
        IF nw LT 1 THEN message, 'Invalid RES keyword'
        sel = w[0]
    END
    1: BEGIN
        w = where(nside_in[0] LE nside, nw)
        IF nw LT 1 THEN BEGIN
            IF nside_in[0] GT 2*max(nside) THEN  $
                message, 'NSIDE more the twice the largest available'
            sel = n_elements(nside) - 1
        ENDIF ELSE BEGIN
            sel = min(w)
        ENDELSE
    END
    2: BEGIN
        IF npix_in[0] GT 2*npix[n_elements(npix)-1] THEN  $
            message, 'NPIX more the twice the largest available'
        minval = min(abs(float(npix_in[0])-npix), sel)
    END
    3: BEGIN
        w = where(spacing_in[0] GT pixsidedeg, nw)
        IF nw LT 1 THEN BEGIN
            IF spacing_in[0] LT 0.5*min(pixsidedeg) THEN  $
                message, 'Spacing less than half of smallest available'
            sel = n_elements(pixsidedeg) - 1
        ENDIF ELSE BEGIN
            sel = min(w)
        ENDELSE
    END
ENDCASE
res_out = res[sel]
nside_out = nside[sel]
npix_out = npix[sel]
spacing_out = pixsidedeg[sel]
pixarea_out = pixsterad[sel]

IF print_flag THEN BEGIN
    print, '   Res     Nside        Npixels   Mean Spacing (Deg) ',  $
        + '    Area (Sterad)'
    fmt = '(I5,I10,I15,F15.4,10X,E15.7)'
    IF print_all THEN BEGIN
        FOR i=0,10 DO print, res[i],nside[i],npix[i],pixsidedeg[i],  $
            pixsterad[i], format=fmt
    ENDIF ELSE BEGIN
        print, res_out,nside_out,npix_out,spacing_out, pixarea_out, format=fmt
    ENDELSE
ENDIF

RETURN
END
