pro dbput,item,val,entry
;+
; NAME:
;	DBPUT
; PURPOSE:
;	Procedure to place a new value for a specified item into
;	a data base file entry.
;
; CALLING SEQUENCE:	
;	dbput, item, val, entry
;
; INPUTS:
;	item - item name or number
;	val - item value(s)
;
; INPUT/OUTPUT:
;	entry - entry (byte array) or scalar entry number.
;	        if entry is a scalar entry number then the data
;	        base file will be updated.  Otherwise the change
;	        will be only made to the entry array which must
;	        be written latter using DBWRT.
;
; OPERATIONAL NOTES:
;	If entry is a scalar entry number or the input file name
;	is supplied, the entry in the data base will be updated
;	instead of a supplied entry variable.  In this case, !priv
;	must be greater than 1.
; HISTORY:
;	version 2  D. Lindler  Feb 1988 (new db formats)
;	modified to convert blanks into zeros correctly D. Neill Jan 1991
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
;-----------------------------------------------------------------------
;
; get item number
;
 db_item, item, inum, ivalnum, dtype, sbyte, numvals, nbytes
;   
; convert val to correct type and check size
;
 s = size(val)
 if (dtype[0] NE 7) and ( s[s[0]+1] EQ 7) then val = strtrim(val)
 case dtype[0] of
	1: v = byte(fix(val))
	2: v = fix(val)
	3: v = long(val)
	4: v = float(val)
	5: v = double(val)
	7: v = string(val)
 endcase
;
 if N_elements(v) NE numvals[0] then begin
	print,'DBPUT - Invalid number of data values'
	print,'Item '+item+' requires ',strtrim(numvals[0],2),' values'
	print,'DBPUT aborting'
	retall
 endif
;
; determine if entry number supplied
;
 s = size(entry)
 if s[0] EQ 0 then begin			;scalar entry number supplied
	dbrd,entry,e
	dbxput,v,e,dtype[0],sbyte[0],nbytes[0]*numvals[0] ;update entry
	dbwrt,e					;update file
  end else begin				;array supplied, just update it
	dbxput,v,entry,dtype[0],sbyte[0],nbytes[0]*numvals[0]
 end

 return
 end
