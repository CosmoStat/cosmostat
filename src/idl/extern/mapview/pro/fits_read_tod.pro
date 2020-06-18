;+
; NAME:
;	FITS_READ_TOD
; PURPOSE:
;	Reads a MAP time-ordered data (TOD) archive file that was written in
;	FITS format.
; CALLING SEQUENCE:
;	FITS_READ_TOD, file, arch, [/NOHK] 
; INPUTS:
;	file - The name of the FITS file to read.
; OUTPUTS:
;	arch - The TOD data structure.    See the function TOD_FORMAT()
;              for additional information about this structure
; OPTIONAL INPUT KEYWORD:
;       /NOHK - If set, Then no housekeeping information is Returned, and 
;               only the time stamp and science data are defined in the 
;               output structure
; EXAMPLE:
;       Read the time-ordered data for day 221 of year 2002 into an output
;       structure
;       IDL>  file ='MAP_tod_20022202357_20022212357.fits'
;       IDL>  FITS_READ_TOD, file, arch
; COMMENTS:
;	The sci array in the archive record contains the science data; there
;	are thirty major science frames in each record; each of the frequency
;	bands has some number of observations in each of these frames; the
;	number of observations will be mentioned below.
;
;	Further, the science data in each archive record assumes a four-channel
;	form for each observation of each differencing assembly (DA).  If the
;	file contains a single channel for each DA, each of the four channels
;	of an observation for a DA will contain the same single channel value.
;	If the file contains two channels for each DA, the first two channels
;	will contain the first channel from the file; the second two channels
;	will contain the second channel from the file.  Finally, if the file
;	contains four channels for each DA then each channel in the output
;	structure will contain a separate value.
;
;	The science data is stored in the output structure by bandpass (K, Ka,
;	Q, V, and W).  These bandpass arrays map to the DAs as follows:
;	  arch.sci.k [ 0: 3,nobs]   -- K1  (nobs=12)
;	  arch.sci.ka[ 0: 3,nobs]  -- Ka1  (nobs=12)
;	  arch.sci.q [ 0: 3,nobs]   -- Q1  (nobs=15)
;	  arch.sci.q [ 4: 7,nobs]   -- Q2  (nobs=15)
;	  arch.sci.v [ 0: 3,nobs]   -- V1  (nobs=20)
;	  arch.sci.v [ 4: 7,nobs]   -- V2  (nobs=20)
;	  arch.sci.w [ 0: 3,nobs]   -- W1  (nobs=30)
;	  arch.sci.w [ 4: 7,nobs]   -- W2  (nobs=30)
;	  arch.sci.w [ 8:11,nobs]   -- W3  (nobs=30)
;	  arch.sci.w [12:15,nobs]   -- W4  (nobs=30)
; PROCEDURES USED:
;       tod_format(), jul2ts                              MAP Library
;       fits_open,fits_read,fits_close, tbinfo,tbget()    IDLAstro Library
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, SSAI, 27 September 2002.
;       Additional Vectorization, NoHK keyword W. Landsman SSAI 21 November 2002
;	Modified to automatically handle 1, 2, and 4 channel data.
;	   MRG, SSAI, 15 August 2006.
;-
; =================================================================================
Function FITS_Read_TOD_Sci_Id, tb_str
;
on_error, 2
;
w = where(strpos(tb_str.ttype, 'K1') GE 0)
If (w[0] LT 0) Then return, 0
;
Return, n_elements(w)
End
; =================================================================================
Pro FITS_Read_TOD_Sci_1ch, tab, tb_str, arch
;
on_error, 2
;
dimen = size(arch.sci.k,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'K1'),dimen[1:3] )
For i=0,3 Do t[0,0,0,i] = d
arch.sci.k = transpose(t,[3,0,1,2])
;
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'KA1'),dimen[1:3] )
For i=0,3 Do t[0,0,0,i] = d
arch.sci.ka = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.q,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'Q1'),dimen[1:3] )
For i=0,3 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'Q2'),dimen[1:3] )
For i=4,7 Do t[0,0,0,i] = d
arch.sci.q = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.v,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'V1'),dimen[1:3] )
For i=0,3 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'V2'),dimen[1:3] )
For i=4,7 Do t[0,0,0,i] = d
arch.sci.v = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.w,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'W1'),dimen[1:3] )
For i=0,3   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W2'),dimen[1:3] )
For i=4,7   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W3'),dimen[1:3] )
For i=8,11  Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W4'),dimen[1:3] )
For i=12,15 Do t[0,0,0,i] = d
arch.sci.w = transpose(t,[3,0,1,2])
;
Return
End
; =================================================================================
Pro FITS_Read_TOD_Sci_2ch, tab, tb_str, arch
;
on_error, 2
;
dimen = size(arch.sci.k,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'K11'),dimen[1:3] )
For i=0,1 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'K12'),dimen[1:3] )
For i=2,3 Do t[0,0,0,i] = d
arch.sci.k = transpose(t,[3,0,1,2])
;
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'KA11'),dimen[1:3] )
For i=0,1 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'KA12'),dimen[1:3] )
For i=2,3 Do t[0,0,0,i] = d
arch.sci.ka = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.q,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'Q11'),dimen[1:3] )
For i=0,1 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'Q12'),dimen[1:3] )
For i=2,3 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'Q21'),dimen[1:3] )
For i=4,5 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'Q22'),dimen[1:3] )
For i=6,7 Do t[0,0,0,i] = d
arch.sci.q = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.v,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'V11'),dimen[1:3] )
For i=0,1 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'V12'),dimen[1:3] )
For i=2,3 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'V21'),dimen[1:3] )
For i=4,5 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'V22'),dimen[1:3] )
For i=6,7 Do t[0,0,0,i] = d
arch.sci.v = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.w,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
d = reform( tbget(tb_str,tab,'W11'),dimen[1:3] )
For i=0,1   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W12'),dimen[1:3] )
For i=2,3   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W21'),dimen[1:3] )
For i=4,5   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W22'),dimen[1:3] )
For i=6,7   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W31'),dimen[1:3] )
For i=8,9   Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W32'),dimen[1:3] )
For i=10,11 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W41'),dimen[1:3] )
For i=12,13 Do t[0,0,0,i] = d
d = reform( tbget(tb_str,tab,'W42'),dimen[1:3] )
For i=14,15 Do t[0,0,0,i] = d
arch.sci.w = transpose(t,[3,0,1,2])
;
Return
End
; =================================================================================
Pro FITS_Read_TOD_Sci_4ch, tab, tb_str, arch
;
on_error, 2
;
dimen = size(arch.sci.k,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
t[0,0,0,0] = reform( tbget(tb_str,tab,'K113'),dimen[1:3] )
t[0,0,0,1] = reform( tbget(tb_str,tab,'K114'),dimen[1:3] )
t[0,0,0,2] = reform( tbget(tb_str,tab,'K123'),dimen[1:3] )
t[0,0,0,3] = reform( tbget(tb_str,tab,'K124'),dimen[1:3] )
arch.sci.k = transpose(t,[3,0,1,2])
;
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
t[0,0,0,0] = reform( tbget(tb_str,tab,'KA113'),dimen[1:3] )
t[0,0,0,1] = reform( tbget(tb_str,tab,'KA114'),dimen[1:3] )
t[0,0,0,2] = reform( tbget(tb_str,tab,'KA123'),dimen[1:3] )
t[0,0,0,3] = reform( tbget(tb_str,tab,'KA124'),dimen[1:3] )
arch.sci.ka = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.q,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
t[0,0,0,0] = reform( tbget(tb_str,tab,'Q113'),dimen[1:3] )
t[0,0,0,1] = reform( tbget(tb_str,tab,'Q114'),dimen[1:3] )
t[0,0,0,2] = reform( tbget(tb_str,tab,'Q123'),dimen[1:3] )
t[0,0,0,3] = reform( tbget(tb_str,tab,'Q124'),dimen[1:3] )
t[0,0,0,4] = reform( tbget(tb_str,tab,'Q213'),dimen[1:3] )
t[0,0,0,5] = reform( tbget(tb_str,tab,'Q214'),dimen[1:3] )
t[0,0,0,6] = reform( tbget(tb_str,tab,'Q223'),dimen[1:3] )
t[0,0,0,7] = reform( tbget(tb_str,tab,'Q224'),dimen[1:3] )
arch.sci.q = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.v,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
t[0,0,0,0] = reform( tbget(tb_str,tab,'V113'),dimen[1:3] )
t[0,0,0,1] = reform( tbget(tb_str,tab,'V114'),dimen[1:3] )
t[0,0,0,2] = reform( tbget(tb_str,tab,'V123'),dimen[1:3] )
t[0,0,0,3] = reform( tbget(tb_str,tab,'V124'),dimen[1:3] )
t[0,0,0,4] = reform( tbget(tb_str,tab,'V213'),dimen[1:3] )
t[0,0,0,5] = reform( tbget(tb_str,tab,'V214'),dimen[1:3] )
t[0,0,0,6] = reform( tbget(tb_str,tab,'V223'),dimen[1:3] )
t[0,0,0,7] = reform( tbget(tb_str,tab,'V224'),dimen[1:3] )
arch.sci.v = transpose(t,[3,0,1,2])
;
dimen = size(arch.sci.w,/dimen)
t = make_array(dimen=[dimen[1:3],dimen[0]],/float,/nozero)
t[0,0,0, 0] = reform( tbget(tb_str,tab,'W113'),dimen[1:3] )
t[0,0,0, 1] = reform( tbget(tb_str,tab,'W114'),dimen[1:3] )
t[0,0,0, 2] = reform( tbget(tb_str,tab,'W123'),dimen[1:3] )
t[0,0,0, 3] = reform( tbget(tb_str,tab,'W124'),dimen[1:3] )
t[0,0,0, 4] = reform( tbget(tb_str,tab,'W213'),dimen[1:3] )
t[0,0,0, 5] = reform( tbget(tb_str,tab,'W214'),dimen[1:3] )
t[0,0,0, 6] = reform( tbget(tb_str,tab,'W223'),dimen[1:3] )
t[0,0,0, 7] = reform( tbget(tb_str,tab,'W224'),dimen[1:3] )
t[0,0,0, 8] = reform( tbget(tb_str,tab,'W313'),dimen[1:3] )
t[0,0,0, 9] = reform( tbget(tb_str,tab,'W314'),dimen[1:3] )
t[0,0,0,10] = reform( tbget(tb_str,tab,'W323'),dimen[1:3] )
t[0,0,0,11] = reform( tbget(tb_str,tab,'W324'),dimen[1:3] )
t[0,0,0,12] = reform( tbget(tb_str,tab,'W413'),dimen[1:3] )
t[0,0,0,13] = reform( tbget(tb_str,tab,'W414'),dimen[1:3] )
t[0,0,0,14] = reform( tbget(tb_str,tab,'W423'),dimen[1:3] )
t[0,0,0,15] = reform( tbget(tb_str,tab,'W424'),dimen[1:3] )
arch.sci.w = transpose(t,[3,0,1,2])
;
Return
End
; =================================================================================
Pro FITS_Read_TOD, file, arch, NoHK=nohk
;
on_error, 2
;
; Check arguments.
;
If (n_params() LT 2) Then Begin 
  print, 'Syntax - FITS_READ_TOD, file, arch, [/NoHK]'
  Return
EndIf
;
; Open the FITS file and find the number of extensions.
;
fits_open, file, fcb		;Extract the file control block
If (fcb.nextend LT 5) Then message, 'Invalid MAP TOD FITS file!'
;
; Extract the number of records in the file and allocate space.
;
narch = sxpar(fcb.hmain, 'NUMREC')
If (narch LE 0) Then message, 'Invalid FITS file -- No data'
;
arch = TOD_Format(narch, NoHK=nohk)         ;Return empty output structure
;
; Extension 2: Fill the array from the science data table.
;
fits_read,fcb,tab,eh,extname = 'Science Data Table',message=msg,/No_Abort
If (msg NE '') Then message, 'Invalid FITS file -- Science'
;
tbinfo,eh,tb_str
tstab = Jul2TS(tbget(tb_str,tab, 'Time'))
dimen = size(arch.sci,/dimen)
arch.sci.Day  = reform(tstab[0,*], dimen)
arch.sci.Time = reform(tstab[1,*], dimen)
arch.sci.GenFlags = reform( tbget(tb_str,tab,'GenFlags'), dimen)
dimen = size(arch.sci.daflags,/dimen)
arch.sci.DAFlags  = reform ( tbget(tb_str,tab,'DAFlags'), dimen)
If (NOT keyword_set(nohk)) Then Begin
  dimen = size(arch.sci.error1,/dimen)
  arch.sci.Error1   = reform( tbget(tb_str,tab,'Error1'), dimen)
  arch.sci.Error2   = reform( tbget(tb_str,tab,'Error2'), dimen)
EndIf
;
Case (FITS_Read_TOD_Sci_Id(tb_str)) Of
  4L   : FITS_Read_TOD_Sci_4ch, tab, tb_str, arch
  2L   : FITS_Read_TOD_Sci_2ch, tab, tb_str, arch
  1L   : FITS_Read_TOD_Sci_1ch, tab, tb_str, arch
  Else : message, 'Invalid FITS file -- Science'
EndCase
;
If (NOT keyword_set(nohk)) Then Begin
;
; Extension 1:  Fill the array from the Meta data table.
;
  arch.calflag = 1
  fits_read,fcb,tab,eh,extname = 'Meta Data Table',message=msg,/No_Abort
  If (msg NE '') Then message, 'Invalid FITS file -- Meta data'

  If ((size(tab,/dimen))[1] GT narch) Then tab = tab[*,0:narch-1]
  tbinfo,eh,tb_str
;
  s = size(arch[0].quaternions,/dimen)
  arch.Quaternions = reform(tbget(tb_str,tab,'QUATERN'), s[0],s[1], narch)
;
  arch.Position    = tbget(tb_str,tab,'Position') 
  arch.Velocity    = tbget(tb_str,tab,'Velocity') 
  arch.GeoPosition = tbget(tb_str,tab, 'GeoPos') 
  arch.GeoVelocity = tbget(tb_str,tab, 'GeoVel') 
;
  arch.nsci        =  tbget(tb_str,tab,'NSCI')
  arch.ndeu        =  tbget(tb_str,tab,'NDIHK')
;
; Extension 3: Fill the array from the Analog Instrument HouseKeeping (AIHK) 
; data table.
;
  fits_read,fcb,tab,eh,extname = 'AIHK Data Table',message=msg,/No_Abort
  If (msg NE '') Then message, 'Invalid FITS file -- AIHK Data Table'
;
  tbinfo,eh,tb_str
  sweep = tbget(tb_str,tab,'Sweep')
  pdu = tbget(tb_str,tab,'PDU')
  aeu1 = tbget(tb_str,tab,'AEU1')
  aeu2 = tbget(tb_str,tab,'AEU2')
  awin1 = tbget(tb_str,tab,'AWIN1')
  awin2 = tbget(tb_str,tab,'AWIN2')
  counters = tbget(tb_str, tab, 'COUNTERS')
;
  tstab = Jul2TS(tbget(tb_str,tab,'Time'))
  g = where(sweep EQ 1, Ng)
  If Ng GT 0 Then Begin
    arch.aeu.Day = reform(tstab[0,g])
    arch.aeu.Time = reform(tstab[1,g])
    arch.aeu.SWEEP1.PDU_Analog = pdu[*,g]
    arch.aeu.Sweep1.AEU_Analog1 = aeu1[*,g] 
    arch.aeu.Sweep1.AEU_Analog2 = aeu2[*,g]
    arch.aeu.Sweep1.AEU_Window1 = awin1[*,g]
    arch.aeu.Sweep1.AEU_Window2 = awin2[*,g]
    arch.aeu.Counters = counters[*,g]
  EndIf

  g = where(sweep NE 1, Ng)
  If Ng GT 0 Then Begin 
    arch.aeu.Sweep2.PDU_Analog  = pdu[*,g]
    arch.aeu.Sweep2.AEU_Analog1 = aeu1[*,g]
    arch.aeu.Sweep2.AEU_Analog2 = aeu2[*,g]
    arch.aeu.Sweep2.AEU_Window1 = awin1[*,g]
    arch.aeu.Sweep2.AEU_Window2 = awin2[*,g]
  EndIf 
;
; Extension 4: Fill the array from the Digital Instrument HouseKeeping (DIHK)
; data table.
;
  fits_read,fcb,tab,eh,extname = 'DIHK Data Table',message=msg,/No_Abort
  If (msg NE '') Then message, 'Invalid FITS file -- DIHK Data Table'
;
  tbinfo,eh,tb_str
  f_ind = tbget(tb_str,tab,'FrmInd') - 1
  d_ind = tbget(tb_str,tab,'DihkInd') - 1
  data  = tbget(tb_str,tab,'Data')
  N     = sxpar(eh, 'NAXIS2') - 1
  tstab = Jul2TS(tbget(tb_str,tab,'Time'))
 
  For i = 0, N Do Begin
    arch[f_ind[i]].deu[d_ind[i]].Day  = tstab[0,i]
    arch[f_ind[i]].deu[d_ind[i]].Time = tstab[1,i]
    arch[f_ind[i]].deu[d_ind[i]].Data = data[*,i]
  EndFor
EndIf
;
; Close the file and return.
;
fits_close,fcb
Return
End
