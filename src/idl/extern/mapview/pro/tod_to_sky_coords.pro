PRO TOD_to_Sky_Coords, TOD, DA_str, Res, Side=Side, MajorFrm=majfrm, Center=cent, $
    ECL=ecl, GAL=gal, CEL=cel, PIXEL=pixel
;+
;  NAME:
;     TOD_to_Sky_Coords
;  PURPOSE:
;     Extract a time series of astronomical coordinates or of A and B
;     side pixel numbers for the start of each observation in an array
;     of WMAP time-ordered data for a given DA.
; CALLING SEQUENCE:
;     TOD_to_Sky_Coords, TOD, DA_str, [Res, Side=, /MajorFrm, /Center, 
;                                      ECL=, GAL=, CEL=, PIXEL= ]
; INPUTS:
;     TOD      An array of 'N' WMAP time-ordered data records, of
;              the same structure as returned by READ_TOD_FITS.
; OPTIONAL INPUTS:
;     DA_str   Case-insensitive string containing the DA to process:
;              'K1',Ka1',Q1','Q2','V1','V2','W1','W2','W3','W4'.   If not
;              supplied or set to '', then the mean optical axis is used;
;              this can be done ONLY if the MajorFrm keyword is set.
;     Res      Integer(1-13) giving the resolution of the returned nested 
;              Healpix  pixel number.  Required if the output PIXEL 
;              keyword is supplied.
; OPTIONAL INPUT KEYWORD:
;     Side     A string containing either 'A' or 'B' or 'AB' specifying which 
;              side of the spacecraft to compute coordinates for.    Default 
;              is 'AB' to compute for both sides.
;     MajorFrm If present and non-zero the coordinates for a DA at the start
;              of each major science frame are returned.  See the comments
;              section for more details.
;     Center   If present and non-zero the coordinates returned are interpolated
;              to the center of each observation.  This keyword is ignored if
;              MajorFrm is set.
; OUTPUT KEYWORDS:
;     If only one side is specified, then the following output arrays will be
;     N x 30 x M x 2 rather than N x 30 x M x 4.  More details concerning
;     the output arrays is supplied below in COMMENTS.
;
;     ECL      Array containing ecliptic long & lat for side A,
;              ecliptic long & lat for side B.
;     GAL      Array containing galactic long & lat for side A,
;              galactic long & lat for side B.
;     CEL      Array containing RA & Dec in degrees for side A,
;              RA & Dec in degrees for side B.
;     PIXEL    N x 30 x M x 2 array containing pixel number for side A,
;              pixel number for side B.  The Side keyword is ignored.
; EXAMPLE:
;     Find the Side A pointing of the K band in Galactic coordinates for 
;     the first science frame in a TOD file
;
;     IDL> fits_read_tod, 'MAP_tod_20022162357_20022172357.fits', tod
;     IDL> tod_to_sky_coords, tod[0], 'K1', Gal=gal, /MajorFrm
;
; COMMENTS:
;     This procedure is a wrapper to the Quat_to_Sky_Coords procedure,
;     encapsulating the second two steps in the example described in that
;     procedure, and interpolating quaternions as necessary while setting
;     the call up.
;
;     There are thirty major science frames (sci elements) in each TOD
;     record (thus N x 30 x M x 4 output arrays); these are separated
;     in time by 1.536 seconds.  M is the number of observations for a DA:
;	K1 & Ka1 - 12
;	Q1 & Q2  - 15
;	V1 & V2  - 20
;	W1 - W4  - 30
;     The coordinates returned at those at the start each observation unless
;     the /Center keyword is specified; in that event the coordinates are
;     interpolated to the center of each observation.
;
;     If the /MajorFrm keyword is specified then the coordinates returned 
;     are those of the start of each major science frame in each TOD record.
;     In this event the arrays returned will be dimensioned N x 30 x 4
;     (N x 30 x 2 for the PIXEL array, and for the other arrays if only
;     one side is specified).
;
;     When computing the coordinates for each observation of a DA, this
;     procedure may take quite a while to return!
; PROCEDURES USED:
;     Quat_to_Sky_Coords, Interpolate_Quaternions     WMAP Library
; REVISION HISTORY:
;     Written by Michael R. Greason, SSAI, 16 August 2006.
;-
on_error, 2
;
;			Check arguments.
;
If (N_params() LT 1) Then Begin
  print,'Syntax - TOD_to_Sky_Coords, TOD, DA_str, [Res, Side=, /MajorFrm, /Center,'
  print,'                                          ECL=, GAL=, CEL=, PIXEL= ]'
  print,'	  Output Keywords - ECL=, GAL=, CEL=, PIXEL = '
  Return
EndIf
;
If (n_elements(DA_str) LE 0) Then DAstr = ''   $
                             Else DAstr = strupcase(strtrim(DA_str, 2))
DAS    = ['K1', 'KA1', 'Q1', 'Q2', 'V1', 'V2', 'W1', 'W2', 'W3', 'W4']
NObsDA = [ 12,   12,    15,   15,   20,   20,   30,   30,   30,   30]
Ind    = where(DAS EQ DAstr)
If ((Ind LT 0) AND (NOT keyword_set(majfrm))) Then Begin
  print,'Syntax - TOD_to_Sky_Coords, TOD, DA_str, [Res, Side=, /MajorFrm, /Center,'
  print,'                                          ECL=, GAL=, CEL=, PIXEL= ]'
  print,'	  Output Keywords - ECL=, GAL=, CEL=, PIXEL = '
  print,'  A valid DA must be identified unless /MajorFrm is supplied.'
  Return
EndIf
;
ns = 4
If (n_elements(Side) GT 0) Then Begin
  If ((Side EQ 'A') OR (Side EQ 'a') OR (Side EQ 'B') OR (Side EQ 'b')) Then $
     ns = 2 
EndIf
;
;			Compute only to the start of the major science frames.
;
nt = n_elements(tod)
If (keyword_set(majfrm)) Then Begin
  q = tod.quaternions
  q = q[*,1:30,*]
  q = reform(q, 4, (30L * nt))
  If (ARG_PRESENT(Pixel)) Then Begin
    Quat_to_Sky_Coords, q, DA_str, Res, Side=Side, $
                        ECL=ecl, GAL=gal, CEL=cel, PIXEL=pixel
  EndIf Else Begin
    Quat_to_Sky_Coords, q, DA_str, Res, Side=Side, $
                        ECL=ecl, GAL=gal, CEL=cel
  EndElse
  If (n_elements(ecl)	GT 0) Then ecl   = transpose(reform(ecl,   30, nt, ns), [1, 0, 2])
  If (n_elements(gal)	GT 0) Then gal   = transpose(reform(gal,   30, nt, ns), [1, 0, 2])
  If (n_elements(cel)	GT 0) Then cel   = transpose(reform(cel,   30, nt, ns), [1, 0, 2])
  If (n_elements(pixel) GT 0) Then pixel = transpose(reform(pixel, 30, nt, 2),  [1, 0, 2])
;
EndIf Else Begin
;
;			Compute for each observation in the DA.
;
  NObs = long(NObsDA[Ind[0]])
  q    = dblarr(4, NObs, 30, nt)
  For i = 0L, (nt-1) Do Begin
    For j = 0L, 29L Do Begin
      qt = tod[i].quaternions[*,j:(j+3)]
      For k = 0L, (NObs-1L) Do Begin
        If (keyword_set(cent)) Then offset = 1.0D0 + ((double(k) + 0.5d0) / double(NObs))  $
	                       Else offset = 1.0D0 + (double(k) / double(NObs))
	Interpolate_Quaternions, qt, offset, qout, status
	If (status GT 0) Then Begin
	  print,'tod_to_sky_coords -- Quaterion interpolation error!'
	  Return
	EndIf
	q[0,k,j,i] = qout
      EndFor
    EndFor
  EndFor
  q   = reform(q, 4, (NObs * 30L * nt))
  If (ARG_PRESENT(Pixel)) Then Begin
    Quat_to_Sky_Coords, q, DA_str, Res, Side=Side, $
                        ECL=ecl, GAL=gal, CEL=cel, PIXEL=pixel
  EndIf Else Begin
    Quat_to_Sky_Coords, q, DA_str, Res, Side=Side, $
                        ECL=ecl, GAL=gal, CEL=cel
  EndElse
  If (n_elements(ecl)	GT 0) Then ecl   = transpose(reform(ecl,   NObs, 30, nt, ns), [2, 1, 0, 3])
  If (n_elements(gal)	GT 0) Then gal   = transpose(reform(gal,   NObs, 30, nt, ns), [2, 1, 0, 3])
  If (n_elements(cel)	GT 0) Then cel   = transpose(reform(cel,   NObs, 30, nt, ns), [2, 1, 0, 3])
  If (n_elements(pixel) GT 0) Then pixel = transpose(reform(pixel, NObs, 30, nt, 2),  [2, 1, 0, 3])
;
EndElse
;
Return
End
