pro ftcreate, MAXCOLS,MAXROWS,H,TAB
;+
; NAME:
;	FTCREATE
; PURPOSE:
;	Create a new (blank) FITS ASCII table and header with specified size.
;
; CALLING SEQUENCE:
;	ftcreate, maxcols, maxrows, h, tab
;
; INPUTS:
;	maxcols - number of character columns allocated, integer scalar
;	maxrows - maximum number of rows allocated, integer scalar
;
; OUTPUTS:
;	h - FITS header, string array
;	tab - empty table, byte array 
; HISTORY:
;	version 1  D. Lindler   July. 87
;
;  21-Sep-88:  Because the degenerative dimension is deleted in Sun IDL,
;	       this procedure has been modified to create table with at
;	       least two rows.  If this isn't done, the other FT routines
;	       choke on a table of one row.
;
;  24-Oct-88:  Changed length of header strings from 81 to 80.  This conforms
;	       to the latest format for FITS header strings.  The 81 character
;	       format was dropped  due to problems it caused when data was
;		transferred back to the VAX.
;
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;----------------------------------------------------------------------
 On_error,2

 if n_params() lt 4 then begin
      print,'Syntax - ftcreate, maxcols, maxrows, h, tab'
      return
 endif

; create blank table

 if maxrows eq 1 then begin
	tmp_rows = 2
 endif else begin
	tmp_rows = maxrows
 endelse

 tab = replicate(32B, maxcols, tmp_rows)

; create header

 h = string(replicate(32b,80,25))
 tmp_head = h[0]
 strput, tmp_head, 'END', 0
 h[0] = tmp_head 
;
; add keywords
;
 sxaddpar, h, 'XTENSION', 'TABLE   '
 sxaddpar, h, 'BITPIX', 8
 sxaddpar, h, 'NAXIS', 2
 sxaddpar, h, 'NAXIS1', 0
 sxaddpar, h, 'NAXIS2', 0
 sxaddpar, h, 'PCOUNT', 0
 sxaddpar, h, 'GCOUNT', 1
 sxaddpar, h, 'TFIELDS', 0

 return
 end
