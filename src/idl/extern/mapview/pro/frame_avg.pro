FUNCTION Frame_Avg, Data, nframes = nframes
;+
; NAME:
;    FRAME_AVG()
; PURPOSE:
;    Compute the average of a major frame (1.536 s) or multiple major frames of
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
;  REVISION HISTORY:
;   initial version, J. Weiland, 8 June 1999 (based on frame_rms.pro)
;   Support for a 3-D array.  MRG, SSAI, 28 February 2002.
;   Use structure definition of SIZE()   W. Landsman    January 2003
;-
 sz = size(data,/str)
 dimen = sz.dimensions
 Case sz.n_dimensions Of
;
   2: locdata = data
;
   3: Begin
      locdata = reform(data, dimen[0], (dimen[1] * dimen[2]))
      sz = size(locdata,/str)
      End
;
Else: message,' ERROR - The input data must be a 2-d or 3-d array'
;
EndCase

Npts = sz.dimensions[1]
dim1 = sz.dimensions[0]

 if N_elements(nframes) EQ 0 then nframes = 1
 if (nframes gt 1) then begin
  nframes = long(nframes)
  Nav = Npts/nframes
  istop = (nav*nframes)-1L
  new = reform(locdata[*,0:istop],nframes*dim1,nav)
  Avg = total(new,1)/float(dim1*nframes)
 endif else begin
  Avg = TOTAL(locdata,1)/float(dim1)
endelse

RETURN, AVG
END
