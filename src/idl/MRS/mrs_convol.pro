;
;	NAME:
;		mrs_convol
;
;	PURPOSE:
;
;		Compute the convolution of a spherical map either in Healpix NESTED format or in Glesp with a filter defined by is spectrum in harmonic space
;
;	CALLING:
;       mrs_convol, Map, Filter, Rec, nside_out = nside_out, Alm=Alm
;
;	INPUT:
;       Map -- IDL array of healpix map or GLESP structure: Input image to be transformed
;		Filter -- IDL 1D array: spectrum of the filter in harmonic space
;          
;	OUTPUTS:
;		Rec: -- IDL array of healpix map or GLESP structure: Output image filtered
;
;	KEYWORDS:
;		nside_out -- int: Nside parameter of the filtered map, by default it is the same value as the one of the input map
;		Alm -- IDL structure of the ALM of the filtered output map (see mrs_almtrans.pro)
;             
;	EXTERNAL CALLS
;			mrs_almtrans
;			mrs_almrec
;
;	EXAMPLE:
;		beam = getbeam( Fwhm=10 )
;		mrs_convol, Map, beam, Rec		Convolve the image map with a 10 arcmin beam.
;
;	HISTORY:
;       Written: Jean-Luc Starck Olivier Fourt 2009.
;       September, 2009 File creation
;
;===============================================================

PRO mrs_convol, Map, Filter, Rec, nside_out = nside_out, Alm=Alm, lmax=lmax

if N_PARAMS() LT 3  then begin 
        print, 'CALLING SEQUENCE: Map, Filter, Rec, nside_out = nside_out, Alm=OutAlm, lmax=lmax'
        goto, DONE
end
mrs_almtrans, Map, Alm, /tab, lmax=lmax
mrs_alm_convol, Alm, Filter
if keyword_set(nside_out)  then Alm.nside = nside_out
mrs_almrec, Alm, Rec

DONE:

END


