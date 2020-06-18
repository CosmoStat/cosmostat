function dbxval,entry,idltype,nvalues,sbyte,nbytes
;+
; NAME: 
;       DBXVAL
;
; PURPOSE:      
;       Quickly return a value of the specified item number     
; EXPLANATION:
;       Procedure to quickly return a value of the specified item number
;       from the entry.
;
; CALLING SEQUENCE:     
;       result = dbxval( entry, idltype, nvalues, sbyte, nbytes )
;
; INPUTS        
;       entry - entry or entries from data base (bytarr) 
;       idltype - idl data type (obtained with db_item_info)
;       nvalues - number of values to return (obtained with db_item)
;       sbyte - starting byte in the entry (obtained with db_item)
;       nbytes - number of bytes (needed only for string type)
;                       (obtained with db_item)
;
; OUTPUTS:      
;       function value is value of the specified item in entry
;
; RESTRICTIONS: 
;       To increase speed the routine assumes that entry and item are
;       valid and that the data base is already opened using dbopen.
;
; REVISION HISTORY:     
;       version 0  D. Lindler Nov. 1987  (for new db format)
;       Version 1, William Thompson, GSFC, 28 March 1994.
;                       Incorporated into CDS library.
;       Version 2, Richard Schwartz, GSFC/SDAC, 23 August 1996
;                       Allowed Entry to have 2 dimensions
;       Version 2.1, 22 Feb 1997, JK Feggans, 
;                               avoid reform for strings arrays.
;       Version 2.2     Use overwrite with REFORM(),  W. Landsman,  May 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Work for multiple-valued strings   W. Landsman   October 2000
;-
;----------------------------------------------------------------
;
;
nentry = n_elements(entry[0,*])

case idltype of                 ;case of data type
  1: val = byte(entry[sbyte:sbyte+nvalues-1,*],0,nvalues,nentry)
  2: val = fix(entry[sbyte:sbyte+nvalues*2-1,*],0,nvalues,nentry)
  3: val = long(entry[sbyte:sbyte+nvalues*4-1,*],0,nvalues,nentry)
  4: val = float(entry[sbyte:sbyte+nvalues*4-1,*],0,nvalues,nentry)
  5: val = double(entry[sbyte:sbyte+nvalues*8-1,*],0,nvalues,nentry)
  7: val = string( reform( entry[sbyte:sbyte+nbytes-1,*], nbytes/nvalues, $
                   nvalues, nentry))
endcase
;

if ( nvalues EQ 1 and nentry EQ 1) then return,val[0] else $
        if idltype eq 7 then return,val else return,reform(val,/overwrite)
end
