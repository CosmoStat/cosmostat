; $Id: read_combined_map.pro,v 1.4 2006/03/27 18:13:12 bhill Exp $
;
Pro Read_Combined_Map, FileName, Imap, Qmap, Umap, Smap,		$
         N_Obs=NObs, PolWtArr=PolWt, PriHdr=phdr, StokesHdr=shdr,		$
	 PolWtHdr=chdr, Verbose=vrb
;+
; NAME:
;	Read_Combined_Map
; PURPOSE:
;	Reads a combined temperature/polarization map FITS file.
; CALLING SEQUENCE:
;	Read_Combined_Map, FileName, Imap, Qmap, UMap [, Smap]
; INPUTS:
;	FileName - The name of the FITS file.
; OUTPUTS:
;	Imap     - The Stokes I map.
;	Qmap     - The Stokes Q map.
;	Umap     - The Stokes U map.
;	Smap     - The spurious signal map, if available.
; OUTPUT KEYWORDS:
;	N_Obs     - The N_Obs column.
;	PolWtArr  - The polarization weights matrix.  This will be
;	            a 2x2xN or 3x3xN array depending upon the contents
;	            of the weights binary table.  The 3x3xN option
;	            will be accompanied by an Smap array.
;	PriHdr    - The primary FITS header.
;	StokesHdr - The Stokes data binary table FITS header.
;	PolWtHdr  - The weights matrix binary table FITS header.
; COMMENTS:
;	'mrdfits.pro' is used to read the contents of the file.
;
;	Two spellings of the Q and U columns are supported:
;	p_POLARIZATION and p_POLARISATION
;
;	The relevant binary FITS tables must be the first two extensions
;	in the file, but the order of these two is irrelevant.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, SSAI, 06 October 2004.
;	New column names.  Covariance becomes Weight.  MRG, SSAI, 29 October 2004.
;	N_QQ, etc., changed to QQ, etc.  RSH, SSAI, 27 Mar 2006.
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 4) Then message, 				$
    'Syntax:  Read_Combined_Map, FileName, Imap, Qmap, UMap [, Smap]'
;
;			Get the primary FITS header.
;
tab1 = mrdfits(FileName, 0, phdr, /Silent)
;
;			Read the two extensions.
;
tab1 = mrdfits(FileName, 1, shdr, Silent=(NOT keyword_set(vrb)))
tab2 = mrdfits(FileName, 2, chdr, Silent=(NOT keyword_set(vrb)))
;
;			Fill the output arrays.
;
If (strpos(strupcase(sxpar(shdr, 'EXTNAME')), 'WEIGHTS') GE 0) Then Begin
;
;				PolWtiance table was the first table.
;
  hdr  = chdr
  chdr = shdr
  shdr = hdr
;
;					Get the stokes signal maps.
;
  t  = tag_names(tab2)
  qf = 0
  uf = 0
  sf = 0
  For i = 0, (n_elements(t) - 1) Do Begin
    If (strupcase(t[i]) EQ 'Q_POLARIZATION') Then qf = 1
    If (strupcase(t[i]) EQ 'U_POLARIZATION') Then uf = 1
    If (strupcase(t[i]) EQ 'SPUR_SIGNAL')    Then sf = 1
  EndFor
;
  Imap = tab2.TEMPERATURE
  NObs = tab2.N_Obs
  If (qf NE 0) Then Qmap = tab2.Q_POLARIZATION			$
               Else Qmap = tab2.Q_POLARISATION
  If (uf NE 0) Then Umap = tab2.U_POLARIZATION			$
               Else Umap = tab2.U_POLARISATION
  If (sf NE 0) Then Smap = tab2.SPUR_SIGNAL
;
;					Get the weights matrix.
;
  n  = n_elements(tab1)
  t  = tag_names(tab1)
  cf = 0
  For i = 0, (n_elements(t) - 1) Do Begin
    If ((strupcase(t[i]) EQ 'M11') AND (cf LE 0)) Then cf = -1
    If (strupcase(t[i]) EQ 'M33') Then cf = 1
  EndFor
;
  If (cf GT 0) Then Begin
    PolWt = fltarr(3,3,n)
    PolWt[0,0,*] = tab1.M11
    PolWt[1,0,*] = tab1.M12
    PolWt[2,0,*] = tab1.M13
    PolWt[0,1,*] = tab1.M12
    PolWt[1,1,*] = tab1.M22
    PolWt[2,1,*] = tab1.M23
    PolWt[0,2,*] = tab1.M13
    PolWt[1,2,*] = tab1.M23
    PolWt[2,2,*] = tab1.M33
  EndIf Else If (cf LT 0) Then Begin
    PolWt = fltarr(2,2,n)
    PolWt[0,0,*] = tab1.M11
    PolWt[1,0,*] = tab1.M12
    PolWt[0,1,*] = tab1.M12
    PolWt[1,1,*] = tab1.M22
  EndIf Else Begin
    PolWt = fltarr(2,2,n)
    PolWt[0,0,*] = tab1.QQ
    PolWt[1,0,*] = tab1.QU
    PolWt[0,1,*] = tab1.QU
    PolWt[1,1,*] = tab1.UU
  EndElse
;
EndIf Else Begin
;
;				Stokes table was the first table.
;
;					Get the stokes signal maps.
;
  t  = tag_names(tab1)
  qf = 0
  uf = 0
  sf = 0
  For i = 0, (n_elements(t) - 1) Do Begin
    If (strupcase(t[i]) EQ 'Q_POLARIZATION') Then qf = 1
    If (strupcase(t[i]) EQ 'U_POLARIZATION') Then uf = 1
    If (strupcase(t[i]) EQ 'SPUR_SIGNAL')    Then sf = 1
  EndFor
;
  Imap = tab1.TEMPERATURE
  NObs = tab1.N_Obs
  If (qf NE 0) Then Qmap = tab1.Q_POLARIZATION			$
               Else Qmap = tab1.Q_POLARISATION
  If (uf NE 0) Then Umap = tab1.U_POLARIZATION			$
               Else Umap = tab1.U_POLARISATION
  If (sf NE 0) Then Smap = tab1.SPUR_SIGNAL
;
;					Get the weights matrix.
;
  n  = n_elements(tab2)
  t  = tag_names(tab2)
  cf = 0
  For i = 0, (n_elements(t) - 1) Do Begin
    If ((strupcase(t[i]) EQ 'M11') AND (cf LE 0)) Then cf = -1
    If (strupcase(t[i]) EQ 'M33') Then cf = 1
  EndFor
;
  If (cf GT 0) Then Begin
    PolWt = fltarr(3,3,n)
    PolWt[0,0,*] = tab2.M11
    PolWt[1,0,*] = tab2.M12
    PolWt[2,0,*] = tab2.M13
    PolWt[0,1,*] = tab2.M12
    PolWt[1,1,*] = tab2.M22
    PolWt[2,1,*] = tab2.M23
    PolWt[0,2,*] = tab2.M13
    PolWt[1,2,*] = tab2.M23
    PolWt[2,2,*] = tab2.M33
  EndIf Else If (cf LT 0) Then Begin
    PolWt = fltarr(2,2,n)
    PolWt[0,0,*] = tab2.M11
    PolWt[1,0,*] = tab2.M12
    PolWt[0,1,*] = tab2.M12
    PolWt[1,1,*] = tab2.M22
  EndIf Else Begin
    PolWt = fltarr(2,2,n)
    PolWt[0,0,*] = tab2.QQ
    PolWt[1,0,*] = tab2.QU
    PolWt[0,1,*] = tab2.QU
    PolWt[1,1,*] = tab2.UU
  EndElse
;
EndElse
;
Return
End
