FUNCTION dbindex_blk, unit, nb, bsz, ofb, dtype
;+
; NAME:
;	DBINDEX_BLK
; PURPOSE:
;	Subroutine of DBINDEX to create associated variable of correct datatype
; EXPLANATION:
;	DBINDEX_BLK will offset into the file by a specified amount in 
;	preparation for writing to the file.
;
; CALLING SEQUENCE:
;	res = dbindex_blk(unit, nb, bsz, ofb, dtype)
;
; INPUTS:
;	unit   The unit number assigned to the file.
;	nb     The number of blocks to offset into the file.
;	bsz    The size of each block, in bytes, to offset into the file.
;	ofb    The offset into the block, in bytes.
;	dtype  The IDL datatype as defined in the SIZE function
;
; OUTPUTS:
;	res    The returned variable.  This is an associated variable.
;
; RESTRICTIONS:
;	The file must have been previously opened.
;
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, STX, 14 June 1990.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
offset = long(nb) * long(bsz) + long(ofb)
case dtype of
	7: datarec=assoc(unit,bytarr(1),offset)		; string
	1: datarec=assoc(unit,bytarr(1),offset)		; byte
	2: datarec=assoc(unit,intarr(1),offset)		; integer
	4: datarec=assoc(unit,fltarr(1),offset)		; floating point
        3: datarec=assoc(unit,lonarr(1),offset)		; longword
        5: datarec=assoc(unit,dblarr(1),offset)		; double
        6: datarec=assoc(unit,complexarr(1),offset)	; complex
endcase
;
RETURN, datarec
END
