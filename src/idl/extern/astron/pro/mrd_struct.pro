;+
; NAME:
;       MRD_STRUCT
; PURPOSE:
;       Return a structure as defined in the names and values data.
; CALLING SEQUENCE:
;       struct = MRD_STRUCT(NAMES, VALUES, NROW,                $
;                   STRUCTYP=structyp,                          $
;                   TEMPDIR=tempdir, /OLD_STRUCT)
; INPUT PARAMETERS:
;       NAMES   = A string array of names of structure fields.
;       VALUES  = A string array giving the values of the structure
;                 fields.  See examples below.
;       NROW    = The number of elements in the structure array.
;       
; RETURNS:
;       A structure as described in the parameters or 0 if an error
;       is detected.
;
; OPTIONAL KEYWORD PARAMETERS:
;       STRUCTYP = The structure type.  Since IDL does not allow the
;                  redefinition of a named structure it is an error
;                  to call MRD_STRUCT with different parameters but
;                  the same STRUCTYP in the same session.  If this
;                  keyword is not set an anonymous structure is created.
;       TEMPDIR  = If the structure is more than modestly complex a
;                  temporary file is created.  This file will be
;                  created in the current directory unless the TEMPDIR
;                  keyword is specified.  Note that the temporary directory
;                  should also be in the IDL search path.
;       OLD_STRUCT=Use old format structures.
; COMMON BLOCKS:
;       MRD_COMMON
; SIDE EFFECTS:                                                            
;       May create a temporary file if the structure definition is too long 
;       for the EXECUTE function and using old style structures
;
; RESTRICTIONS:
;       By default this program uses a series of execute
;       commands and create_struct's to create the structure.
;       If the old_struct keyword is set, then a program may
;       be dynamically compiled.  The nominal maximum length
;       of the execute string is 131 characters, but many systems
;       seem to allow longer values.  This code may execute more
;       efficiently with a longer execute buffer.
; PROCEDURE:
;       A structure definition is created using the parameter values.
;       MRD_NSTRUCT is called if the OLD_STRUCT keyword is not specified
;       and generates the structure in pieces using the
;       execute and create_struct keywords.
;
;       If the old_struct flag is set, then the program tries to compile
;       the structure with a single execute command.  If the structure
;       definition is too long  MRD_FSTRUCT is called to write, compile and
;       execute a function which will define the structure.
; EXAMPLES:
;       str = mrd_struct(['fld1', 'fld2'], ['0','dblarr(10,10)'],3)
;       print, str(0).fld2(3,3)
;
;       str = mrd_struct(['a','b','c','d'],['1', '1.', '1.d0', "'1'"],1)
;               ; returns a structure with integer, float, double and string
;               ; fields.
; PROCEDURE CALLS:
;       CONCAT_DIR - Used to concatenate temporary directory with filename
; MODIFICATION HISTORY:
;       Created by T. McGlynn October, 1994.
;       Modified by T. McGlynn September, 1995.
;          Added capability to create substructures so that structure
;          may contain up to 4096 distinct elements.  [This can be
;          increased by futher iteration of the process used if needed.]
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Removed V4.0 reference to common block  October 1997
;       Allowed unlimited number of structure elements if the version
;       is greater than 5.0.  Put back in code to handle prior versions.
;       The [] will need to be translated back to () for this to
;       work.  T. McGlynn December 15 1998.
;       Add MRD_NSTRUCT since IDL has mysterious problems compiling
;       very large structures.
;-
function mrd_fstruct, names, values, nrow, $
    structyp=structyp, tempdir=tempdir, silent=silent, old_struct=old_struct
;
; Create a structure with column names given in names, a string 
; giving the value of the columns in values and optionally the
; structure name in structyp.  This function differs from mrd_struct
; in that it creates a temporary procedure.  It is used when the
; string length for the execute function exceeds 131 characters.
;
; 
; Inputs:
;       Names: a string array giving the names of the structure
;             columns.
;       Values: a string array giving the prototype value for 
;               the corresponding column.
;       Nrec: The number of instances of the structure
;             to be returned.
;       Structyp=  An optional name if the structure is to be a named
;                  structure.
;       Tempdir=   Optional directory where temporary files are to be
;                  created.
;       /Silent    If set, the !quiet system variable is set to
;                  suppress compilation messages.
;
;  Modified September 7, 1995 to handle more than structures
;  with more than 127 elements (the IDL limit).  This uses
;  structures of structures and may fail on earlier versions
;  of IDL.              TAM.
;  Modified October 16, 1995 to add SILENT keyword.
;  Modified April 25, 1996 by TAM to use RESOLVE_ROUTINE function
;    when !version.release > 4.0
;  Use /usr/bin/rm rather than rm to delete temporary file, W. Landsman 6/1997
;  Use IDL features to delete rather than operating system dependent
;  capabilities. T. McGlynn 12/98
;  Use concat_dir to concatenate directory name and file name W. Landsman 12/00

common mrd_common, usage

rel_no = float(!version.release)

if (rel_no lt 4) then begin
	if n_elements(usage) eq 0 then usage = 0 else usage = usage+1
	proname = 'mrdt_'+strcompress(string(usage),/remo)
endif else begin
        proname = 'mrd_structtemp'
endelse

 if keyword_set(tempdir) then $
      filename = concat_dir(tempdir, proname+'.pro') $
 else $
     filename = proname+'.pro'

openw, lun, filename, /get_lun
        
printf,lun, 'function '+proname
printf, lun, 'return, $'

; If we are going to return an array of structures, then we use
; the replicate function.
if nrow gt 1 then begin
        printf, lun, "replicate({ $"
endif else begin
        printf, lun, "{$"
endelse

; Put in the structure name if the user has specified one.
if keyword_set(structyp) then begin
        printf, lun,  structyp + ",$"
endif

; Now for each element put in a name/value pair.
nel = n_elements(names)

comma = ' '
exten = ''

; Do the output for the structure.

; If Oldstruct is set and the  structure contains more than
; 127 elements, then the first 64 elements will be processed normally, and the
; remaining elements will be as substructures of up to 64 elements in the last
; 63 elements of the primary structure.  The structure elements will have the
; name ss## where ## ranges from 1 to 63.  This allows a total of 64^2 or 4096
; elements to be defined in a structure.


substruct = nel gt 127  and  keyword_set(old_struct)

for i=0, nel-1 do begin

        if (i mod 64) eq 0 and i ne 0 and substruct then begin
                if i ne 64 then begin
                    ; Close the previous substructure.
                    printf, lun, "} $"
                endif
                ; Create substructure to contain next 64 columns.
                exten = "ss"+strcompress(string(i/64),/remo)
                printf, lun, comma+exten+':{ $'
                comma = ' '  ; Don't put a comma before this element
        endif
        
        printf, lun, comma+names[i] + ':'+ values[i] + '$'
        comma = ','
endfor

; Close last substructure.
if substruct then begin
    printf, lun, "} $"
endif
        
if nrow gt 1 then begin
        printf, lun, "}$" 
        printf, lun,  "," + strtrim(long(nrow),2)+ ")"
endif else begin
        printf, lun, "}"
endelse

printf, lun, "end"

free_lun, lun

a = 0
if keyword_set(silent) then begin
    old_quiet = !quiet
    !quiet = 1
endif

if rel_no ge 4.0 then begin
    resolve_routine, proname,  /is_function
endif
a = call_function(proname)
        
if keyword_set(silent) then !quiet = old_quiet

; Open the filename with the Temporary qualifier and
; then close it.  This will delete the file.
openr, lun, filename, /delete
free_lun, lun

return, a
end


function mrd_nstruct, names, values, nrow, structyp=structyp

; Build up the structure use a combination of execute and
; create_struct calls.  Basically we build as many rows as
; will fit in an execute call and create that structure.  Then
; we append that structure to whatever we've done before using
; create_struct


; Check that the number of names is the same as the number of values.
nel = n_elements(names)

a = 0
cind = 0
; Start formatting the string.
strng = "a={"

comma = ' '
for i=0,nel-1 do  begin
  
    ; Now for each element put in a name/value pair.
    tstrng = strng + comma+names[i] + ':' + values[i]
    
    ; The nominal max length of the execute is 131
    ; We need one chacacter for the "}"
    if strlen(tstrng) gt 130 then begin
        strng = strng + "}"
        res = execute(strng)
	if  res eq 0 then return, 0
	if n_elements(bigstr) eq 0 then begin
	    bigstr = a
	endif else begin
	    bigstr = create_struct(bigstr, a)
	endelse
	strng = "a={" + names[i] + ":" + values[i]
	
    endif else begin
        strng = tstrng
    endelse
    comma = ","
	 
endfor

if strlen(strng) gt 3 then begin
    strng = strng + "}"
    res = execute(strng)
    if  res eq 0 then return, 0
    if n_elements(bigstr) eq 0 then begin
	bigstr = a
    endif else begin
	bigstr = create_struct(bigstr, a)
    endelse
  
endif

if keyword_set(structyp) then begin
     bigstr = create_struct(bigstr, name=structyp)
endif
  
if nrow le 1 then return, bigstr

return, replicate(bigstr, nrow)
end

function mrd_struct, names, values, nrow, $
    structyp=structyp, tempdir=tempdir, silent=silent, old_struct=old_struct


; Create an instance of A, since an execute function cannot
; create an instance of an undefined variable.

a = 0

; Check that the number of names is the same as the number of values.
nel = n_elements(names)
if nel ne n_elements(values) then return, 0

if not keyword_set(old_struct) then begin
    return, mrd_nstruct(names, values, nrow, structyp=structyp)
endif

; Start formatting the string.
strng = "a="

; If we are going to return an array of structures, then we use
; the replicate function.
if nrow gt 1 then begin
        strng = strng + "replicate({"
endif else begin
        strng = strng + "{"
endelse

; Put in the structure name if the user has specified one.
if keyword_set(structyp) then begin
        strng = strng + structyp + ","
endif

; Now for each element put in a name/value pair.
for i=0, nel-1 do begin
        if i ne 0 then strng = strng + ','
        strng = strng + names[i] + ':'
        strng = strng + values[i]
endfor

strng = strng + "}"

; Put in the second argument to the REPLICATE function if
; needed.
if nrow gt 1 then begin
        strng = strng + "," + strtrim(long(nrow),2)+ ")"
endif

; The IDL documentation implies that 131 is the maximum length
; for EXECUTE although many implementations seem to support longer
; strings.  We'll use this value though to be safe.

if strlen(strng) gt 131 then begin
        return, mrd_fstruct(names, values, nrow, structyp=structyp, $
           tempdir=tempdir, silent=silent, old_struct=old_struct)
endif else begin
        ; Execute the string.  RES should be 1 if the execution was successful.
        res = execute(strng)

        if res eq 0 then return, 0 else return, a
endelse

end
