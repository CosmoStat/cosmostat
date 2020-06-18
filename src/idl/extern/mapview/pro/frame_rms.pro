
FUNCTION Frame_RMS, Data, nframes=nframes
;+
; NAME:
;    FRAME_RMS()
; PURPOSE:
;    Compute the RMS of a major frame (1.536 s) or multiple major frames of
;     science data (if nframes is specified and >1)
; CALLING SEQUENCE:
;     avg = FRAME_AVG( data, Nframes = )
; INPUTS:
;     data  -- the input array Data to be averaged: must be two- or 
;         three-dimensional and of the form  minor frame x major frame 
;         [ x uberframe], eg. given as REFORM(Data.K(0,*)).
;  OUTPUT: 
;     result - the data appropiately averaged over the number of frames
;              the first dimension of the input data array will be removed
;  OPTIONAL INPUT:
;      Nframes - Integer giving number of frames to average over, default is 1
; REVISION HISTORY:
;   initial version, G. Hinshaw, April 1999?
;   added nframes keyword, J. Weiland, 8 June 1999
;   fixed logic for no keyword use, J. Weiland 11 June 1999
;   added test for npts lt nframes
;   Support for a 3-D array.  MRG, SSAI, 28 February 2002.
;   Use structure definition of SIZE(), STDDEV()   W. Landsman    January 2003
;- 

 sz = size(data,/str)
 dimen = sz.dimensions
 Case sz.n_dimensions Of
;
   2: locdata = data
    
   3: Begin
      locdata = reform(data, dimen[0], (dimen[1] * dimen[2]))
      sz = size(locdata,/str)
      End
 ;
Else: message,' ERROR The input data must be a 2-d or 3-d array'
;
EndCase

Npts = sz.dimensions[1]

if not keyword_set(nframes) then nframes = 1

if ((nframes gt 1) and (npts ge nframes)) then begin
  nframes = long(nframes)
  nparts = Npts/nframes
  istart = lindgen(nparts)*nframes
  istop   = istart + (nframes-1L)
  rms = fltarr(nparts)
  for i = 0L,nparts-1L do rms[i] = stddev(locdata[*,istart[i]:istop[i]])
  RETURN, RMS
endif 
;
;  no keyword set or keyword set to 1 or npts lt nframes
;
 RMS = FLTARR(Npts)
 FOR I = 0L,Npts-1L DO RMS[I] = STDDEV(locdata[*,I])
 RETURN, RMS

END
