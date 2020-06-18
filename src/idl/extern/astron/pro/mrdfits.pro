;+
; NAME:
;     MRDFITS
;
; PURPOSE:
;     Read all standard FITS data types into arrays or structures.
;
; EXPLANATION:
;      Further information on MRDFITS is available at
;      http://idlastro.gsfc.nasa.gov/mrdfits.html 
;
; CALLING SEQUENCE:
;      Result = MRDFITS( Filename/FileUnit,[Extension, Header],
;                       /FSCALE , /DSCALE , /USE_COLNUM, /NO_TDIM, /OLD_STRUCT,
;                       RANGE=[a,b], COLUMNS=[a,b,...]), ERROR_ACTION=x,
;                       COMPRESS=comp_prog, STATUS=status
;
; INPUTS:
;      Filename = String containing the name of the file to be read or
;                 file number of an open unit.  If a unit is specified
;                 if will be left open positioned to read the next HDU.
;                 Note that the file name may be of the form
;                 name.gz or name.Z on UNIX systems.  If so
;                 the file will be dynamically decompressed.
;      FiluUnit = An integer file unit which has already been
;                 opened for input.  Data will be read from this
;                 unit and the unit will be left pointing immediately
;                 after the HDU that is read.  Thus to read a compressed
;                 file with many HDU's a user might do something like:
;                      lun=fxposit(filename, 3)  ; Skip the first three HDU's
;                      repeat begin
;                          thisHDU = mrdfits(lun, 0, hdr, status=status)
;                          ... process the HDU ...
;                      endrep until status lt 0
;
;      Extension= Extension number to be read, 0 for primary array.
;                 Assumed 0 if not specified.
;                 If a unit rather than a filename
;                 is specified in the first argument, this is
;                 the number of HDU's to skip from the current position.
;
; OUTPUTS:
;      Result = FITS data array or structure constructed from
;               the designated extension.  The format of result depends
;               upon the type of FITS data read.
;             Non-group primary array or IMAGE extension:
;               A simple multidimensional array is returned with the
;               dimensions given in the NAXISn keywords.
;             Grouped image data with PCOUNT=0.
;               As above but with GCOUNT treated as NAXIS(n+1).
;             Grouped image data with PCOUNT>0.
;               The data is returned as an array of structures.  Each
;               structure has two elements.  The first is a one-dimensional
;               array of the group parameters, the second is a multidimensional
;               array as given by the NAXIS2-n keywords.
;             ASCII and BINARY tables.
;               The data is returned as a structure with one column for
;               each field in the table.  The names of the columns are
;               normally taken from the TTYPE keywords (but see USE_COLNUM).
;               Bit field columns
;               are stored in byte arrays of the minimum necessary
;               length.  Column names are truncated to 15 characters
;               if longer, spaces are removed, and invalid characters
;               are replaced by underscores.
;               Columns specified as variable length columns are stored
;               with a dimension equal to the largest actual dimension
;               used.  Extra values in rows are filled with 0's or blanks.
;               If the size of the variable length column is not
;               a constant, then an additional column is created
;               giving the size used in the current row.  If the length
;               of each element of a variable length column is 0 then
;               the column is deleted.
;
;               Prior to V5.0, IDL structures were limited to 128 tags.
;               If the version is before V5.0, or the /OLD_STRUCT is set, then
;               for FITS files with more than 127 columns, data in the first
;               64 elements of the structure are stored in the primary
;               structure, the next 64 as a substructure of the 65th
;               element, the next 64 as a substructure of the 66th element
;               and so forth.
;
; OPTIONAL OUTPUT:
;       Header = String array containing the header from the FITS extenion.
;
; OPTIONAL INPUT KEYWORDS:
;       FSCALE - If present and non-zero then scale data to float
;                numbers for arrays and columns which have either
;                non-zero offset or non-unity scale.
;                If scaling parameters are applied, then the corresponding
;                FITS scaling keywords will be modified.
;       DSCALE - As with FSCALE except that the resulting data is
;                stored in doubles.
;       /SILENT - Suppress informative messages.
;       RANGE  - A scalar or two element vector giving the start
;                and end rows to be retrieved.  For ASCII and BINARY
;                tables this specifies the row number.  For GROUPed data
;                this will specify the groups.  For array images, this
;                refers to the last non-unity index in the array.  E.g.,
;                for a 3 D image with NAXIS* values = [100,100,1], the
;                range may be specified as 0:99, since the last axis
;                is suppressed.  Note that the range uses IDL indexing
;                So that the first row is row 0.
;                If only a single value, x, is given in the range,
;                the range is assumed to be [0,x-1].
;       USE_COLNUM - When creating column names for binary and ASCII tables
;                MRDFITS attempts to use the appropriate TTYPE keyword
;                values.  If USE_COLNUM is specified and non-zero then
;                column names will be generated as 'C1, C2, ... 'Cn'
;                for the number of columns in the table.
;       STRUCTYP - The structyp keyword specifies the name to be used
;                for the structure defined when reading ASCII or binary
;                tables.  Generally users will not be able to conveniently
;                combine data from multiple files unless the STRUCTYP
;                parameter is specified.  An error will occur if the
;                user specifies the same value for the STRUCTYP keyword
;                in calls to MRDFITS in the same IDL session for extensions
;                which have different structures.
;       NO_TDIM  - Disable processing of TDIM keywords.  If NO_TDIM
;                is specified MRDFITS will ignore TDIM keywords in
;                binary tables.
;       OLD_STRUCT- Use the recursive structures formats required
;                prior to IDL 5.0 for tables with more than 127 columns.
;       TEMPDIR - The tempdir keyword allows the user to specify
;                the directory where temporary files may be created.
;                This directory should be both in the IDL path
;                and writable by the user.    Generally only needed for IDL
;          
;       COLUMNS - This keyword allows the user to specify that only a
;                subset of columns is to be returned.  The columns
;                may be specified either as number 1,... n or by
;                name or some combination of these two.
;                If USE_COLNUM is specified names should be C1,...Cn.
;                The use of this keyword will not save time or internal
;                memory since the extraction of specified columns
;                is done after all columns have been retrieved from the
;                FITS file.
;       COMPRESS - This keyword allows the user to specify a
;                decompression program to use to decompress a file that
;                will not be automatically recognized based upon
;                the file name.
;                
;       ERROR_ACTION - Set the on_error action to this value (defaults
;                to 2).
; OPTIONAL OUTPUT KEYWORDS:
;       STATUS - A integer status indicating success or failure of
;                the request.  A status of >=0 indicates a successful read.
;                Currently
;                    0 -> successful completion
;                   -1 -> error
;                   -2 -> end of file
;
; EXAMPLE:
;       Read a FITS primary array:
;               a = mrdfits('TEST.FITS')    or
;               a = mrdfits('TEST.FITS', 0, header)
;       The second example also retrieves header information.
;
;       Read rows 10-100 of the second extension of a FITS file.
;               a = mrdfits('TEST.FITS', 2, header, range=[10,100])
;
;       Read a table and ask that any scalings be applied and the
;       scaled data be converted to doubles.  Use simple column names,
;       suppress outputs.
;               a = mrdfits('TEST.FITS', 1, /dscale, /use_colnum, /silent)
;
; RESTRICTIONS:
;       (1)     Cannot handle data in non-standard FITS formats.
;       (2)     Doesn't do anything with BLANK or NULL values or
;               NaN's.  They are just read in.  They may be scaled
;               if scaling is applied.
; NOTES:
;       This multiple format FITS reader is designed to provide a
;       single, simple interface to reading all common types of FITS data.
;       MRDFITS DOES NOT scale data by default.  The FSCALE or DSCALE
;       parameters must be used.
;
; PROCEDURES USED:
;       The following procedures are contained in the main MRDFITS program.
;           MRD_IMAGE           -- Generate array/structure for images.
;           MRD_READ_IMAGE      -- Read image data.
;           MRD_ASCII           -- Generate structure for ASCII tables.
;           MRD_READ_ASCII      -- Read an ASCII table.
;           MRD_TABLE           -- Generate structure for Binary tables.
;           MRD_READ_TABLE      -- Read binary table info.
;           MRD_READ_HEAP       -- Read variable length record info.
;           MRD_SCALE           -- Apply scaling to data.
;           MRD_COLUMNS         -- Extract columns.
;
;        Other ASTRON Library routines used
;           FXPAR(), FXADDPAR, IEEE_TO_HOST, FXPOSIT, FXMOVE(), IS_IEEE_BIG()
;           MRD_STRUCT(), MRD_SKIP
;
; MODIfICATION HISTORY:
;       V1.0 November 9, 1994 ----  Initial release.
;          Creator: Thomas A. McGlynn
;       V1.1 January 20, 1995 T.A. McGlynn
;          Fixed bug in variable length records.
;          Added TDIM support -- new routine mrd_tdim in MRD_TABLE.
;       V1.2
;          Added support for dynamic decompression of files.
;          Fixed further bugs in variable length record handling.
;       V1.2a
;          Added NO_TDIM keyword to turn off TDIM processing for
;          those who don't want it.
;          Bug fixes: Handle one row tables correctly, use BZERO rather than
;               BOFFSET.     Fix error in scaling of images.  
;       V1.2b 
;          Changed MRD_HREAD to handle null characters in headers.
;       V2.0 April 1, 1996
;          -Handles FITS tables with an arbitrary number of columns.
;          -Substantial changes to MRD_STRUCT to allow the use of
;          substructures when more than 127 columns are desired.
;          -All references to table columns are now made through the
;          functions MRD_GETC and MRD_PUTC.  See description above.
;          -Use of SILENT will now eliminate compilation messages for
;          temporary functions.
;          -Bugs in handling of variable length columns with either
;          a single row in the table or a maximum of a single element
;          in the column fixed.
;          -Added support for DCOMPLEX numbers in binary tables (M formats) for
;          IDL versions above 4.0.  
;          -Created regression test procedure to check in new versions.
;          -Added error_action parameter to allow user to specify
;          on_error action.  This should allow better interaction with
;          new CHECK facility.  ON_ERROR statements deleted from
;          most called routines.
;          - Modified MRDFITS to read in headers containing null characters
;          with a warning message printed.
;       V2.0a April 16, 1996
;          - Added IS_IEEE_BIG() checks (and routine) so that we don't
;          worry about IEEE to host conversions if the machine's native
;          format is IEEE Big-endian.
;       V2.1 August 24, 1996
;          - Use resolve_routine for dynamically defined functions
;          for versions > 4.0
;          - Fix some processing in random groups format.
;          - Handle cases where the data segment is--legally--null.
;          In this case MRDFITS returns a scalar 0.
;          - Fix bugs with the values for BSCALE and BZERO (and PSCAL and
;          PZERO) parameters set by MRDFITS.
;       V2.1a April 24, 1997  Handle binary tables with zero length columns
;       V2.1b May 13,1997 Remove whitespace from replicate structure definition
;       V2.1c May 28,1997 Less strict parsing of XTENSION keyword
;       V2.1d June 16, 1997 Fixed problem for >32767 entries introduced 24-Apr
;       V2.1e Aug 12, 1997 Fixed problem handling double complex arrays
;       V2.1f Oct 22, 1997 IDL reserved words can't be structure tag names
;       Converted to IDL V5.0   W. Landsman  2-Nov-1997
;       V2.1g Nov 24, 1997 Handle XTENSION keywords with extra blanks.
;       V2.1h Jul 26, 1998 More flexible parsing of TFORM characters
;       V2.2 Dec 14, 1998 Allow fields with longer names for
;                        later versions of IDL.
;                        Fix handling of arrays in scaling routines.
;                        Allow >128 fields in structures for IDL >4.0
;                        Use more efficient structure copying for
;                        IDL>5.0
;       V2.2b June 17, 1999 Fix bug in handling case where
;                           all variable length columns are deleted
;                           because they are empty.
;       V2.3 March 7, 2000 Allow user to supply file handle rather
;                          than file name.
;                          Added status field.
;                          Now needs FXMOVE routine
;       V2.3b April 4, 2000
;                          Added compress option (from D. Palmer)
;       V2.4  July 4, 2000 Added STATUS=-1 for "File access error" (Zarro/GSFC)
;
;       Note to users of IDL prior to V5.0:  This version is compiled
;       with the [] array syntax.  To convert this version to run under
;       the earlier syntax you use W. Landsman's IDL5to4 routine
;       at http://idlastro.gsfc.nasa.gov/ftp/contrib/landsman/v5.
;       You need to change all array subscripts from the []
;       to () but this cannot be done as a global replace since
;       array initializations of the form x=[1,2] must be left unchanged.
;       Subroutines called by this routine may need similar modification.
;-
; Utility Functions ==========================================================
;=============================================================================

;=============================================================================
; Substructure handling:
;      These functions are used for versions of IDL prior to 5.0
;      or if the user has specified the old_struct keyword.
;
;      mrd_substruct - Is this column in a substructure
;      mrd_getc ------ Get the value of a column
;      mrd_putc ------ Put the value into a column
;      mrd_sizcol ---- Get the size of a column
;=============================================================================
function mrd_substruct, table, column
;
; This function returns a value indicating whether the specified
; column will be found in a substructure or in the main structure.

if column lt 64 then return, 0

tagn = tag_names(table)
sz = size(table[0].(64))
if sz[sz[0]+1] eq 8 then return, 1 else return, 0

end


function mrd_getc, table, column

; This function returns the n'th column of a table.  It allows the table
; to contain more than 127 columns, the IDL limit for the number of fields
; in a structure.
if not mrd_substruct(table, column) then begin
    return, table.(column)  
endif else begin
    return, table.(63+column/64).(column mod 64)
endelse

end

function mrd_sizcol, table, column

; This function returns the size function
; of the n'th column of a table.  It allows the table
; to contain more than 127 columns, the IDL limit for the number of fields
; in a structure.
if not mrd_substruct(table, column) then begin
    return, size(table[0].(column))
endif else begin
    return, size(table[0].(63+column/64).(column mod 64))
endelse
end

pro mrd_putc, table, column, value
common mrd_versions, idl_version, mrd_version

; This function updates a column in a table.

if not mrd_substruct(table, column) then begin
   table.(column) = value
endif else begin
   table.(63+column/64).(column mod 64) = value
endelse
end



function  mrd_dofn, name, index, use_colnum 
 
; Convert the string name to a valid variable name.  If name 
; is not defined generate the string Cnn when nn is the index 
; number. 

 reserved_names = ['AND','BEGIN','CASE','COMMON','DO','ELSE','END','ENDCASE', $
 'ENDELSE','ENDFOR','ENDIF','ENDREP','ENDWHILE','EQ','FOR','FUNCTION','GE', $
 'GOTO','GT','IF','INHERITS','LE','LT','MOD','NE','NOT','OF','ON_IOERROR', $
 'OR','PRO','REPEAT','THEN','UNTIL','WHILE','XOR' ]

table = 0
sz = size(name) 
nsz = n_elements(sz) 
if not use_colnum and (sz[nsz-2] ne 0) then begin 
        if sz[nsz-2] eq 7 then begin 
                str = name[0] 
        endif else begin 
                str = 'C'+string(index) 
        endelse 
endif else begin 
        str = 'C'+string(index) 
endelse 
 
str = strcompress(str, /remove_all)  
reserved = where(strupcase(str) EQ reserved_names, Nreserv)
if Nreserv GT 0 then str = str + '_'

len = strlen(str) 
c = strmid(str, 0, 1) 
 
; First character must be alphabetic. 
 
if  not (('a' le c  and 'z' ge c) or ('A' le c  and 'Z' ge c)) then begin 
        str = 'X' + str 
endif 
 
; 
; Replace invalid characters with underscores. 
; We assume ASCII ordering.
;
for i=1, len-1 do begin 
        c = strmid(str, i, 1) 
        if not (('a' le c  and 'z' ge c)  or  $ 
                ('A' le c  and 'Z' ge c)  or  $ 
                ('0' le c  and '9' ge c)  or  $ 
                (c eq '_') or (c eq '$')      $ 
               ) then strput, str, '_', i 
endfor 
 
return, str 
end 

;***************************************************************



pro mrd_doff, form, dim, type 
 
; 
; Parse the TFORM keyword and return the type and dimension of the 
; data. 

; Find the first non-numeric character. 
 
len = strlen(form) 
 
if len le 0 then return 
 
for i=0, len-1 do begin 
 
        c = strmid(form, i, 1) 
        if c lt '0'  or c gt '9' then goto, not_number 
         
endfor 
 
not_number: 

if i ge len then return              ;Modified from len-1 on 26-Jul-1998
 
if i gt 0 then begin 
        dim = long(strmid(form, 0, i)) 
        if dim EQ 0l then dim = -1l
endif else begin 
        dim = 0 
endelse 
 
type = strmid(form, i, 1) 
end 



;*********************************************************************

function mrd_chkfn, name, namelist, index 
 
; 
;  Check that this name is unique with regard to other column names.
; 

if float(!version.release) ge 5 then maxlen = 127 else maxlen=15

if strlen(name) gt maxlen then name = strmid(name, 0, maxlen) 
w = where(name eq strmid(namelist, 0, maxlen) )
if w[0] ne -1 then begin
        ; We have found a name conflict. 
        ; 
                name = 'gen$name_'+strcompress(string(index+1),/remove_all) 
endif 
 
return, name 
 
end

;=====================================================================
; END OF GENERAL UTILITY FUNCTIONS ===================================
;=====================================================================

; Utility functions for ASCII tables.

pro mrd_atype, form, type, slen

;
; Parse the TFORM keyword and return the type and dimension of the
; data.

; Find the first non-numeric character.


; Get rid of blanks.
form = strcompress(form,/remove_all)
len = strlen(form)
if len le 0 then return

type = strmid(form, 0,1)
length = strmid(form,1,len-1)
;
; Ignore the number of decimal places.  We assume that there
; is a decimal point.
;
p = strpos(length, '.')
if p gt 0 then length = strmid(length,0,p)

if strlen(length) gt 0 then slen = fix(length) else slen = 1

end



pro mrd_read_ascii, unit, range, nbytes, nrows, nfld, typarr, posarr,   $
     lenarr, nullarr, table, old_struct=old_struct

; Read in the table information.
;
; Unit          Unit to read data from.
; Range         Range of rows to be read
; Nbytes        Number of bytes per row.
; Nrows         Number of rows.
; Nfld          Number of fields in structure.
; Typarr        Array indicating type of variable.
; Posarr        Starting position of fields (first char at 0)
; Lenarr        Length of fields
; Nullarr       Array of null values
; Table         Table to read information into.
; Old_struct    Should recursive structure format be used?

bigstr = bytarr(nbytes, range[1]-range[0]+1)

if range[0] gt 0 then mrd_skip, unit, nbytes*range[0]
readu,unit, bigstr

; Skip to the end of the data area.

nSkipRow = nrows - range[1] - 1
nskipB = 2880 - (nbytes*nrows) mod 2880
if nskipB eq 2880 then nskipB = 0

mrd_skip, unit, nskipRow*nbytes+nskipB

s1 = posarr-1
s2 = s1 + lenarr - 1
for i=0, nfld-1 do begin

        flds = strtrim( bigstr[s1[i]:s2[i],* ] )
        
        if strtrim(nullarr[i]) ne '' then begin
                ; next statement used to be outside this If block
                ; but curr_col is only used when we've got to check
                ; fo nulls.
                ;    (TAM 5/6/96)
                curr_col = mrd_getc(table, i)
                w = where(flds ne strtrim(nullarr[i]))
                if w[0] ne -1 then begin
                        if N_elements(w) EQ 1 then w = w[0]
                        if typarr[i] eq 'I' then begin
                                curr_col[w] = long(flds[w])
                        endif else if typarr[i] eq 'E' or typarr[i] eq 'F' then begin
                                curr_col[w] = float(flds[w])
                        endif else if typarr[i] eq 'D' then begin
                                curr_col[w] = double(flds[w])
                        endif else if typarr[i] eq 'A' then begin
                                curr_col[w] = flds[w]
                        endif
                endif

                if (not keyword_set(old_struct)) then begin
                    table.(i) = curr_col
                endif else begin
                    mrd_putc, table, i, curr_col
                endelse
                
        endif else begin
                
            if typarr[i] eq 'I' then begin
                if (not keyword_set(old_struct)) then begin
                        table.(i) = long(flds)
                endif else begin
                        mrd_putc, table, i, long(flds)
                endelse
            endif else if typarr[i] eq 'E' or typarr[i] eq 'F' then begin
                if (not keyword_set(old_struct)) then begin
                        table.(i) = float(flds)
                endif else begin
                        mrd_putc, table, i, float(flds)
                endelse
            endif else if typarr[i] eq 'D' then begin
                if (not keyword_set(old_struct)) then begin
                        table.(i) = double(flds)
                endif else begin
                        mrd_putc, table, i, double(flds)
                endelse
            endif else if typarr[i] eq 'A' then begin
                if (not keyword_set(old_struct)) then begin
                        table.(i) = flds
                endif else begin
                        mrd_putc, table, i, flds
                endelse
            endif
        endelse
endfor

end



pro mrd_ascii, header, structyp, use_colnum,   $
    range, table, $
    nbytes, nrows, nfld, typarr, posarr, lenarr, nullarr, $
    fnames, fvalues, scales, offsets, scaling, status, $
    silent=silent, tempdir=tempdir, columns=columns, old_struct=old_struct

; Define a structure to hold a FITS ASCII table.
;
; Header                FITS header for table.
; Structyp              IDL structure type to be used for
;                       structure.
; Use_colnum            Use column numbers not names.
; Range                 Range of rows of interest
; Table                 Structure to be defined.
; Nbytes                Bytes per row
; Nrows                 Number of rows in table
; Nfld                  Number of fields
; Typarr                Array of field types
; Posarr                Array of field offsets
; Lenarr                Array of field lengths
; Nullarr               Array of field null values
; Fname                 Column names
; Fvalues               Formats for columns
; Scales/offsets        Scaling factors for columns
; Scaling               Do we need to scale?
; Status                Return status.

table = 0

types = ['I', 'E', 'F', 'D', 'A']
sclstr = ['0l', '0.0', '0.0', '0.0d0', ' ']
status = 0

if strmid(fxpar(header, 'XTENSION'),0,8) ne 'TABLE   ' then begin
        print, 'MRDFITS: Header is not from ASCII table.'
        status = -1;
        return
endif

nfld = fxpar(header, 'TFIELDS')
nrows = fxpar(header, 'NAXIS2')
nbytes = fxpar(header, 'NAXIS1')

if range[0] ge 0 then begin
        range[0] = range[0] < (nrows-1)
        range[1] = range[1] < (nrows-1)
endif else begin
        range[0] = 0
        range[1] = nrows-1
endelse

nrows = range[1] - range[0] + 1

if nrows le 0 then begin
        if not keyword_set(silent) then begin
                print,'MRDFITS: ASCII table.  ',strcompress(string(nfld)),  $
                 ' columns, no rows'
        endif
        return
endif

;
;  Loop over the columns

typarr = strarr(nfld)
lenarr = intarr(nfld)
posarr = intarr(nfld)
nullarr = strarr(nfld)
fnames = strarr(nfld)
fvalues = strarr(nfld)
scales = dblarr(nfld)
offsets = dblarr(nfld)


for i=0, nfld-1 do begin
        suffix = strcompress(string(i+1), /remove_all)
        fname = fxpar(header, 'TTYPE' + suffix)
        fform = fxpar(header, 'TFORM' + suffix)
        fpos = fxpar(header, 'TBCOL' + suffix)
        fnull = fxpar(header, 'TNULL' + suffix)
        scales[i] = fxpar(header, 'TSCAL' + suffix)
        if scales[i] eq 0.0d0 then scales[i] = 1.0d0
        offsets[i] = fxpar(header, 'TZERO'+suffix)
        
        fname = mrd_dofn(fname,i+1, use_colnum)
        fnames[i] = fname

        fname = mrd_chkfn(fname, fnames, i)
        
        mrd_atype, fform, ftype, flen
        typarr[i] = ftype
        lenarr[i] = flen
        posarr[i] = fpos
        nullarr[i] = fnull
        
        for j=0, n_elements(types) - 1 do begin
                if ftype eq types[j] then begin
                        if ftype ne 'A' then begin
                                val = sclstr[j]
                        endif else begin
                                val = 'string(replicate(32b,'+strtrim(flen,2)+'))'
                        endelse
                        
                        fvalues[i] = val
                        
                        goto, next_col
                endif
        endfor
        
        print, 'MRDFITS: Invalid format code:',ftype, ' for column ', i+1
        status = -1
        return
next_col:
endfor

if scaling then begin
        w = where(scales ne 1.0d0 or offsets ne 0.0d0)
        if w[0] eq -1 then scaling = 0
endif
if not scaling and not keyword_set(columns) then begin
        table = mrd_struct(fnames, fvalues, nrows, structyp=structyp, $
          tempdir=tempdir, silent=silent, old_struct=old_struct)
endif else begin
        table = mrd_struct(fnames, fvalues, nrows, tempdir=tempdir, $
          silent=silent, old_struct=old_struct)
endelse

if not keyword_set(silent) then begin
        print,'MRDFITS: ASCII table.  ',strcompress(string(nfld)),  $
         ' columns by ',strcompress(string(nrows)), ' rows.'
endif
status = 0
return

end
pro  mrd_columns, table, columns, fnames, fvalues, $
    vcls, vtpes, scales,  offsets, scaling,        $
    structyp=structyp, tempdir=tempdir, silent=silent, $
    old_struct=old_struct

; Eliminate columns from the table that do not match the
; user specification.


sz = size(columns)

type = sz[sz[0]+1]
nele = sz[sz[0]+2]
if type eq 8 or type eq 6 or type eq 0 then return  ; Can't use structs
                                                    ; or complex.

if type eq 4 or type eq 5 then tcols = fix(columns)
if type eq 1 or type eq 2 or type eq 3 then tcols = columns

; Convert strings to uppercase and compare with column names.

if type eq 7 then begin
        for i=0, nele-1 do begin
                cname = strupcase(columns[i])
                w = where(cname eq strupcase(fnames))
                if w[0] ne -1 then begin
                        if n_elements(tcols) eq 0 then begin
                                tcols = w[0]+1
                        endif else begin
                                tcols = [tcols, w[0]+1]
                        endelse
                endif
        endfor
endif

; Subtract one from column indices and check that all indices >= 0.
if n_elements(tcols) gt 0 then begin
        tcols = tcols-1
        w = where(tcols ge 0)
        if w[0] eq -1 then begin
                dummy = temporary(tcols)
        endif
endif

if n_elements(tcols) le 0 then begin
        print, 'MRDFITS:  No columns match'
        
        ; Undefine variables.  First ensure they are defined, then
        ; use temporary() to undefine them.
        table = 0
        fnames = 0
        fvalues = 0
        vcls = 0
        vtpes = 0
        scales = 0
        offsets = 0
        dummy = temporary(fnames)
        dummy = temporary(fvalues)
        dummy = temporary(vcls)
        dummy = temporary(vtpes)
        dummy = temporary(scales)
        dummy = temporary(offsets)
        scaling = 0
        
endif else begin

        ; Replace arrays with only desired columns.
        
        fnames = fnames[tcols]
        fvalues = fvalues[tcols]
        
        ; Check if there are still variable length columns.
        if n_elements(vcls) gt 0 then begin
                vcls = vcls[tcols]
                vtpes = vtpes[tcols]
                w = where(vcls eq 1)
                if w[0] eq -1 then begin
                        dummy = temporary(vcls)
                        dummy = temporary(vtpes)
                endif
        endif
        
        ; Check if there are still columns that need scaling.
        if n_elements(scales) gt 0 then begin
                scales = scales[tcols]
                offsets = offsets[tcols]
                w = where(scales ne 1.0d0  or offsets ne 0.0d0)
                if w[0] eq -1 then scaling = 0
        endif
        

        ndim = n_elements(table)
        
        if scaling or n_elements(vcls) gt 0 then begin
                tabx = mrd_struct(fnames, fvalues, ndim, tempdir=tempdir, $
                  silent=silent, old_struct=old_struct)
        endif else begin
                tabx = mrd_struct(fnames, fvalues, ndim, structyp=structyp, $
                  tempdir=tempdir, silent=silent, old_struct=old_struct)
        endelse
        
        for i=0, n_elements(tcols)-1 do begin
            if  not keyword_set(old_struct) then begin
                tabx.(i) = table.(tcols[i]);
            endif else begin
                mrd_putc, tabx, i, mrd_getc(table,tcols[i])
            endelse
        endfor

        table = tabx
endelse
end
pro mrd_read_image, unit, range, maxd, rsize, table 
 
; Read in the table information. 
; 
; Unit          Unit to read data from. 
; Table         Table/array to read information into. 
; 

; If necessary skip to beginning of desired data. 


if range[0] gt 0 then mrd_skip, unit, range[0]*rsize
 
if rsize eq 0 then return 
 
readu, unit, table

; Skip to the end of the data

skipB = 2880 - (maxd*rsize) mod 2880
if skipB eq 2880 then skipB = 0

if range[1] lt maxd-1 then begin
    skipB = skipB + (maxd-range[1]-1)*rsize
endif

mrd_skip, unit, skipB

if not is_ieee_big() then ieee_to_host, table 
end 

pro mrd_axes_trunc,naxis, dims, silent

mysilent = silent

for i=naxis-1,1,-1 do begin 

    if dims[i] eq 1 then begin
        if not mysilent then begin
            print, 'MRDFITS: Truncating unused dimensions'
            mysilent = 1
        endif
        dims = dims[0:i-1] 
        naxis = naxis - 1 
        
    endif else begin 
        return
    endelse 
    
endfor 
 
return
end

pro mrd_image, header, range, maxd, rsize, table, scales, offsets, scaling, $
  status, silent=silent, tempdir=tempdir
 
; Define structure/array to hold a FITS image. 
; 
; Header                FITS header for table. 
; Range                 Range of data to be retrieved. 
; Rsize                 Size of a row or group. 
; Table                 Structure to be defined. 
; Status                ; Silent=silent         Suppress info messages?
 
table = 0

lens = [1,2,4,4,8] 
typstrs=['Byte', 'Int*2', 'Int*4', 'Real*4', 'Real*8']
 
status = 0 
 
 
naxis = fxpar(header, 'NAXIS') 
bitpix= fxpar(header, 'BITPIX') 
if naxis gt 0 then dims = fxpar(header, 'NAXIS*') else dims = 0 
gcount = fxpar(header, 'GCOUNT') 
pcount = fxpar(header, 'PCOUNT')
isgroup = fxpar(header, 'GROUPS')
gcount = long(gcount)
 
if bitpix eq 8 then type = 0            $ 
else if bitpix eq 16 then type = 1      $ 
else if bitpix eq 32 then type = 2      $ 
else if bitpix eq -32 then type = 3     $ 
else if bitpix eq -64 then type = 4 
 
; Note that for random groups data we must ignore the first NAXISn keyword. 
if isgroup GT 0  then begin 


        range[0] = range[0] > 0
        if (range[1] eq -1) then begin
                range[1] = gcount-1
        endif else begin
                range[1] = range[1] < gcount - 1
        endelse
	maxd = gcount
        
        if (n_elements(dims) gt 1) then begin
            dims = dims[1:*]
            naxis = naxis-1
        endif else begin
            print, 'Warning: No data specified for group data.'
            dims = [0]
            naxis = 0
        endelse
        
        ; The last entry is the scaling for the sample data.
        
        if (pcount gt 0) then begin
            scales = dblarr(pcount+1)
            offsets = dblarr(pcount+1)
        endif
        
        values = strarr(2)
        typarr=["bytarr", "intarr", "lonarr", "fltarr", "dblarr"] 
        
        
        mrd_axes_trunc, naxis, dims, keyword_set(silent)
        
        values[0] = typarr[type] + "("+string(pcount)+")" 
        rsize = dims[0] 
        sarr = "(" + strcompress(string(dims[0]), /remo )
         
        for i=1, naxis-1 do begin 
                sarr = sarr + "," + strcompress(string(dims[i]),/remo)
                rsize = rsize*dims[i] 
        endfor 
         
        sarr = sarr + ")"

        if not keyword_set(silent) then print,'MRDFITS--Image with groups:', $
          ' Ngroup=',strcompress(string(gcount)),' Npar=',                   $
          strcompress(string(pcount),/remo), ' Group=', sarr, '  Type=',typstrs[type]

        sarr = typarr[type] + sarr
        values[1] = sarr 
        rsize = (rsize + pcount)*lens[type] 
         
        table = mrd_struct(['params','array'], values, range[1]-range[0]+1, $
          tempdir=tempdir, silent=silent)
           

        for i=0, pcount-1 do begin
                istr = strcompress(string(i+1),/remo)
                scales[i] = fxpar(header, 'PSCAL'+istr)
                if scales[i] eq 0.0d0 then scales[i] =1.0d0
                offsets[i] = fxpar(header, 'PZERO'+istr)
                scales[pcount] = fxpar(header, 'BSCALE')
                if scales[pcount] eq 0.0d0 then scales[pcount] = 1.0d0
                offsets[pcount] = fxpar(header, 'BZERO')
        endfor

        if scaling then begin
                w = where(scales ne 1.0d0 or offsets ne 0.0d0)
                if w[0] eq -1 then scaling = 0
        endif
        
endif else begin 
 
        if naxis eq 0 then begin 
                rsize = 0 
                table = 0
                if not keyword_set(silent) then begin
                        print, 'MRDFITS: Null image, NAXIS=0'
                endif
                return 
        endif 
         
        if gcount gt 1 then begin 
                dims = [dims, gcount] 
                naxis = naxis + 1 
        endif 
         
        mrd_axes_trunc, naxis, dims, keyword_set(silent)

        if not keyword_set(silent) then begin
                str = '('
                for i=0, naxis-1 do begin
                        if i ne 0 then str = str + ','
                        str = str + strcompress(string(dims[i]),/remo)
                endfor
                str = str+')'
                print, 'MRDFITS: Image array ',str, '  Type=', typstrs[type]
        endif
                
        maxd = dims[naxis-1] 
         
        if range[0] ne -1 then begin 
                range[0] = range[0]<(maxd-1) 
                range[1] = range[1]<(maxd-1) 
        endif else begin 
                range[0] = 0 
                range[1] = maxd - 1 
        endelse 
         
        dims[naxis-1] = range[1]-range[0]+1 
         
        rsize = 1 
        if naxis gt 1 then for i=0, naxis - 2 do rsize=rsize*dims[i] 
        rsize = rsize*lens[type] 
         
        sz = lonarr(naxis+3) 
        sz[0] = naxis 
        sz[1:naxis] = dims 
        nele = 1l 
         
        for i=0, naxis-1 do begin 
                nele = nele*dims[i] 
        endfor 
         
        sz[naxis+1] = type + 1   
        sz[naxis+2] = nele 
         
        if nele gt 0 then  begin 
                table = make_array(size=sz) 
        endif else begin 
                table = 0 
        endelse 
        
        scales = dblarr(1)
        offsets = dblarr(1)
        scales[0] = fxpar(header, 'BSCALE')
        offsets[0] = fxpar(header, 'BZERO')
        if scales[0] eq 0.0d0 then scales[0] = 1.0d0
        if scaling and scales[0] eq 1.0d0 and offsets[0] eq 0.0d0 then $
            scaling = 0
         
endelse 
         
status = 0 
return 
 
end 
pro mrd_scale, type, scales, offsets, table, header,  $
               fnames, fvalues, nrec, dscale = dscale, structyp=structyp, $
               tempdir=tempdir, silent=silent, old_struct=old_struct
;
; Scale a FITS array or table.
;
; Type:         FITS file type, 0=image/primary array
;                               1=ASCII table
;                               2=Binary table
;
; scales:       An array of scaling info
; offsets:      An array of offset information
; table:        The FITS data.
; header:       The FITS header.
; dscale:       Should data be scaled to R*8?
; fnames:       Names of table columns.
; fvalues:      Values of table columns.
; nrec:         Number of records used.
; structyp:     Structure name.
; old_struct:   Use old style structures
;
; Modified: 1-Aug-1995 to fix typo.

w = where(scales ne 1.d0  or offsets ne 0.d0)
if w[0] eq -1 then return
ww = where(scales eq 1.d0 and offsets eq 0.d0)

; First do ASCII and Binary tables.
if type ne 0 then begin
        
        if type eq 1 then begin
                if keyword_set(dscale) then begin
                        fvalues[w] = '0.0d0'
                endif else begin
                        fvalues[w] = '0.0'
                endelse
        endif else if type eq 2 then begin

                if keyword_set(dscale) then begin
                        sclr = '0.d0'
                        vc = 'dblarr'
                endif else begin
                        sclr = '0.0'
                        vc = 'fltarr'
                endelse
                
                for i=0, n_elements(w)-1 do begin
                        col = w[i]
                        if not keyword_set(old_struct) then begin
                            sz = size(table[0].(col))
                        endif else begin
                            sz = mrd_sizcol(table,col)
                        endelse

                        if sz[0] eq 0 then begin
                                fvalues[col] = sclr
                        endif else begin
                                str = vc + '('
                                for j=0, sz[0]-1 do begin
                                        if j ne 0 then str = str + ','
                                        str = str + string(sz[j+1])
                                endfor
                                str = str + ')'
                                fvalues[col] = str
                        endelse
                endfor
        endif

        tabx = mrd_struct(fnames, fvalues, nrec, structyp=structyp, $
          tempdir=tempdir, silent=silent, old_struct=old_struct)

        if ww[0] ne -1 then begin
                for i=0, n_elements(ww)-1 do begin
                    if (not keyword_set(old_struct)) then begin
                        tabx.(ww[i]) = table.(ww[i])
                    endif else begin
                        mrd_putc, tabx, ww[i], mrd_getc(table,ww[i])
                    endelse
                endfor
        endif
        
        for i=0, n_elements(w)-1 do begin
                if not keyword_set(old_struct) then begin
                    tabx.(w[i]) = table.(w[i])*scales[w[i]] + offsets[w[i]]
                endif else begin
                    mrd_putc, tabx, w[i], mrd_getc(table,w[i])*scales[w[i]] + offsets[w[i]]
                endelse
                istr = strcompress(string(w[i]+1), /remo)
                fxaddpar, header, 'TSCAL'+istr, 1.0, 'Set by MRD_SCALE'
                fxaddpar, header, 'TZERO'+istr, 0.0, 'Set by MRD_SCALE'
        endfor

        table = temporary(tabx)

endif else begin
; Now process images and random groups.

        sz = size(table[0])
        if sz[sz[0]+1] ne 8 then begin
                ; Not a structure so we just have an array of data.
                if keyword_set(dscale) then begin
                        table = table*scales[0]+offsets[0]
                endif else begin
                        table = table*float(scales[0]) + float(offsets[0])
                endelse
                fxaddpar, header, 'BSCALE', 1.0, 'Set by MRD_SCALE'
                fxaddpar, header, 'BZERO', 0.0, 'Set by MRD_SCALE'

        endif else begin
                ; Random groups.  Get the number of parameters by looking
                ; at the first element in the table.
                nparam = n_elements(table[0].(0))
                if keyword_set(dscale) then typ = 'dbl' else typ='flt'
                s1 = typ+'arr('+string(nparam)+')'
                ngr = n_elements(table)
                sz = size(table[0].(1))
                if sz[0] eq 0 then dims = [1] else dims=sz[1:sz[0]]
                s2 = typ + 'arr('
                for i=0, n_elements(dims)-1 do begin 
                        if i ne 0 then s2 = s2+ ','
                        s2 = s2+string(dims[i])
                endfor
                s2 = s2+')'
                tabx = mrd_struct(['params', 'array'],[s1,s2],ngr, $
                  tempdir=tempdir, silent=silent)
                for i=0, nparam-1 do begin
                        istr = strcompress(string(i+1),/remo)
                        fxaddpar, header, 'PSCAL'+istr, 1.0, 'Added by MRD_SCALE'
                        fxaddpar, header, 'PZERO'+istr, 0.0, 'Added by MRD_SCALE'
                        tabx.(0)[i] = table.(0)[i]*scales[i]+offsets[i]
                endfor
                tabx.(1) = table.(1)*scales[nparam] + offsets[nparam]
                fxaddpar, header, 'BSCALE', 1.0, 'Added by MRD_SCALE'
                fxaddpar, header, 'BZERO', 0.0, 'Added by MRD_SCALE'
                table = temporary(tabx)
        endelse
endelse

end
pro mrd_read_heap, unit, header, range, fnames, fvalues, vcls, vtpes, table, $ 
   structyp, scaling, scales, offsets, status, silent=silent,         $
   tempdir=tempdir, columns=columns, no_fix_varcols=no_fix_varcols,   $
   old_struct=old_struct

; This program reads the heap area to get the actual values of variable 
; length arrays. 
; 
; Unit:         FITS unit number. 
; header:       FITS header. 
; fnames:       Column names. 
; fvalues:      Column values. 
; vcols:        Column numbers of variable length columns. 
; vtypes:       Actual types of variable length columns 
; table:        Table of data from standard data area, on output 
;               contains the variable length data. 
; structyp:     Structure name. 
; scaling:      Is there going to be scaling of the data?
; status:       Set to -1 if an error occurs.
; old_struct:   Use old structures?
;
; Modified 17-Feb-1995 by TAM to fix bug when there were
; multiple actual varying columns.
; Modified 2-June-1992 by TAM to fix bug where there are multiple actual
; varying columns and columns deleted.

typstr = 'LXBIJAEDCM' 
prefix = ['bytarr(', 'bytarr(', 'bytarr(', 'intarr(',     $ 
          'lonarr(', 'string(bytarr(', 'fltarr(',         $ 
          'dblarr(', 'cmplxarr(', 'dblarr(2,' ] 
status = 0 

; Convert from a list of indicators of whether a column is variable
; length to pointers to only the variable columns.

vcols = where(vcls eq 1)
vtypes = vtpes[vcols]

nv = n_elements(vcols) 
 
; Find the beginning of the heap area. 
 
heapoff = fxpar(header, 'THEAP') 
sz = fxpar(header, 'NAXIS1')*fxpar(header, 'NAXIS2') 
if heapoff ne 0 and heapoff lt sz then begin 
        print, 'MRDFITS: ERROR Heap begins within data area' 
        status = -1 
        return 
endif

; Skip to beginning.
if (heapoff > sz) then begin
    mrd_skip, unit, heapoff-sz
endif
 
; Get the size of the heap. 
pc = fxpar(header, 'PCOUNT') 
if heapoff eq 0 then heapoff = sz 
hpsiz = pc - (heapoff-sz) 
 
if (hpsiz gt 0) then heap = bytarr(hpsiz) 
 
 
; Read in the heap 
readu, unit, heap

; Skip to the end of the data area.
skipB = 2880 - (sz+pc) mod 2880
if skipB ne 2880 then begin
    mrd_skip, unit, skipB
endif
 
; Find the maximum dimensions of the arrays. 
; 
; Note that the variable length column currently has fields which 
; are I*4 2-element arrays where the first element is the 
; length of the field on the current row and the second is the 
; offset into the heap. 
 
vdims = lonarr(nv)
for i=0, nv-1 do begin 
    col = vcols[i]
    curr_col = mrd_getc(table, col)
    vdims[i] = max(curr_col[0,*])
    w = where(curr_col[0,*] ne vdims[i])
    if w[0] ne -1 then begin
        if n_elements(lencols) eq 0 then begin
                lencols = [col]
        endif else begin
                lencols=[lencols,col]
        endelse
    endif

    if vtypes[i] eq 'X' then vdims[i]=(vdims[i]+7)/8 
    ind = strpos(typstr, vtypes[i])
    
    ; Note in the following that we ensure that the array is
    ; at least one element long.
    
    fvalues[col] = prefix[ind] + string((vdims[i] > 1)) + ')'     
    if vtypes[i] eq 'A' then fvalues[col] = fvalues[col] + ')' 
    
endfor 
 
nfld = n_elements(fnames) 
 
; Get rid of columns which have no actual data. 
w= intarr(nfld) 
w[*] = 1
corres = indgen(nfld)

if not keyword_set(no_fix_varcols) then begin
  
    ww = where(vdims eq 0) 
    if ww[0] ne -1 then  begin
        w[vcols[ww]] = 0
        if not keyword_set(silent) then begin
                print, 'MRDFITS: ', strcompress(string(n_elements(ww))),  $
                  ' unused variable length columns deleted'
        endif
    endif

    ; Check if all columns have been deleted...
    wx = where(w gt 0)
    if (wx[0] eq -1) then begin
        if not keyword_set(silent) then begin
                print, 'MRDFITS: All columns have been deleted'
        endif
	table = 0
	return
    endif
    

    ; Get rid of unused columns.
    corres = corres[wx]
    fnames = fnames[wx] 
    fvalues = fvalues[wx]
    scales = scales[wx]
    offsets = offsets[wx]

    wx = where(vdims gt 0)
    
    if (wx[0] eq -1) then begin
        vcols=[-9999]
	x=temporary(vtypes)
	x=temporary(vdims)
    endif else begin 
        vcols = vcols[wx]
        vtypes = vtypes[wx]
        vdims = vdims[wx]
    endelse

    ; Now add columns for lengths of truly variable length records.
    if n_elements(lencols) gt 0 then begin
        if not keyword_set(silent) then begin
                print, 'MRDFITS: ', strcompress(string(n_elements(lencols))), $
                  ' length column[s] added'
        endif
        

        for i=0, n_elements(lencols)-1 do begin
                col = lencols[i]
                w = where(col eq corres)
                ww = where(col eq vcols)
                w = w[0]
                ww = ww[0]
                fvstr = '0l'
                fnstr = 'L'+strcompress(string(col),/remo)+'_'+fnames[w]
                nf = n_elements(fnames)
                
                ; Note that lencols and col refer to the index of the
                ; column before we started adding in the length
                ; columns.
                
                if w eq nf-1 then begin
                        ; Subtract -1 for the length columns so 0 -> -1 and
                        ; we can distinguish this column.

                        corres = [corres, -col-1 ]
                        fnames = [fnames, fnstr ]
                        fvalues = [fvalues, fvstr ]
                        scales = [scales, 1.0d0 ]
                        offsets = [offsets, 0.0d0 ]
                endif else begin
                        corres = [corres[0:w],-col-1,corres[w+1:nf-1] ]
                        fnames = [fnames[0:w],fnstr,fnames[w+1:nf-1] ]
                        fvalues = [fvalues[0:w],fvstr,fvalues[w+1:nf-1] ]
                        scales = [scales[0:w], 1.0d0, scales[w+1:nf-1] ]
                        offsets = [offsets[0:w],0.0d0, offsets[w+1:nf-1] ]
                endelse
        endfor
    endif
endif

; Generate a new table with the appropriate structure definitions 
if not scaling and not keyword_set(columns) then begin
  tablex = mrd_struct(fnames, fvalues, n_elements(table), structyp=structyp, $
     tempdir=tempdir, silent=silent, old_struct=old_struct)
endif else begin
  tablex = mrd_struct(fnames, fvalues, n_elements(table), tempdir=tempdir, $
     silent=silent, old_struct=old_struct)
endelse


nrow = range[1]-range[0]+1

; I loops over the new table columns, col loops over the old table.
; When col is negative, it is a length column.
for i=0, n_elements(fnames)-1 do begin
        
        col = corres[i]
                
        if col ge 0 then begin
                w = where(vcols eq col)
                
                ; First handle the case of a column that is not
                ; variable length -- just copy the column.
                
                if w[0] eq -1 then begin
                        if  not keyword_set(old_struct) then begin
                            tablex.(i) = table.(col)
                        endif else begin
                            mrd_putc, tablex, i, mrd_getc(table,col)
                        endelse
                        
                 
                endif else begin 
                        vc = w[0]
                        ; Now handle the variable length columns
                        
                        ; If only one row in table, then
                        ; IDL will return curr_col as one-dimensional.
                        ; Since this is a variable length pointer column we
                        ; know that the dimension of the column is 2.
                        if not keyword_set(old_struct) then begin
                            curr_col = table.(col)
                        endif else begin
                            curr_col = mrd_getc(table, col)
                        endelse
                        if (nrow eq 1) then curr_col = reform(curr_col,2,1)
                        
                        siz = curr_col[0,*] 
                        off = curr_col[1,*] 
                        w = where(siz gt 0) 
                        
                        ; Don't process rows where the length is 0.
                        off = off[w] 
                        siz = siz[w] 
                        nw = n_elements(w)-1 
                        
                        ; Now process each type.
                        if not keyword_set(old_struct) then begin
                            curr_colx = tablex.(i)
                            sz = size(curr_colx)
                            if (sz[0] lt 2) then begin
                                 curr_colx = reform(curr_colx, 1, n_elements(curr_colx), /overwrite)
                            endif
                        endif else begin
                            curr_colx = mrd_getc(tablex, i)
                        endelse
                                
                        
                        ; As above we have to worry about IDL truncating
                        ; dimensions.  This can happen if either
                        ; nrow=1 or the max dimension of the column is 1.
                        
                        is_ieee = is_ieee_big()

                        if not keyword_set(old_struct) then begin
                            sz = size(tablex.(i))
                        endif else begin
                            sz = mrd_sizcol(tablex,i)
                        endelse
                            
                        nel = sz[sz[0]+2]
                        if nrow eq 1 and nel eq 1 then begin
                            curr_colx = make_array(1,1,value=curr_colx)
                        endif else if nrow eq 1 then begin
                            curr_colx = reform(curr_colx,[nel, 1], /overwrite)
                        endif else if nel eq 1 then begin
                            curr_colx = reform(curr_colx,[1, nrow], /overwrite)
                        endif
                
                        case vtypes[vc] of 
                          'L':begin 
                                for j=0, nw do begin 
                                   curr_colx[0:siz[j]-1,w[j]] = byte(heap,off[j],siz[j])  
                                endfor 
                              end 
                          'X':begin 
                                siz = 1+(siz-1)/8 
                                for j=0, nw do begin 
                                   curr_colx[0:siz[j]-1,w[j]] = byte(heap,off[j],siz[j])
                                endfor 
                              end 
                          'B':begin 
                                for j=0, nw do begin 
                                   curr_colx[0:siz[j]-1,w[j]] = byte(heap,off[j],siz[j]) 
                                endfor 
                              end 
                          'I':begin 
                                for j=0, nw do begin 
                                 curr_colx[0:siz[j]-1,w[j]] = fix(heap, off[j], siz[j]) 
                                endfor
                                if  not is_ieee then byteorder, curr_colx, /ntohs 
                              end 
                          'J':begin 
                                for j=0, nw do begin 
                                 curr_colx[0:siz[j]-1,w[j]] = long(heap, off[j], siz[j]) 
                                endfor 
                                if  not is_ieee then byteorder, curr_colx, /ntohl 
                              end 
                          'E':begin 
                                for j=0, nw do begin 
                                  curr_colx[0:siz[j]-1,w[j]] = float(heap, off[j], siz[j]) 
                                endfor 
                                if  not is_ieee then byteorder, curr_colx, /xdrtof 
                              end 
                           'D':begin 
                                for j=0, nw do begin 
                                  curr_colx[0:siz[j]-1,w[j]] = double(heap, off[j], siz[j]) 
                                endfor 
                                if  not is_ieee then ieee_to_host, curr_colx
                               end 
                           'C':begin 
                                for j=0, nw do begin 
                                  curr_colx[0:siz[j]-1,w[j]] = complex(heap, off[j], siz[j]) 
                                endfor 
                                if  not is_ieee then byteorder, curr_colx, /xdrtof 
                               end 
                           'M':begin 
                                for j=0, nw do begin 
                                
                              curr_colx[0:siz[j]-1,w[j]] = dcomplex(heap, off[j], siz[j]) 
                                  
                                endfor 
                                if  not is_ieee then ieee_to_host, curr_colx 
                               end 
                           'A':begin 
                                curr_colx[w] = string(heap[off:off+siz]) 
                               end 
                        endcase
                
                        if nel eq 1 and nrow eq 1 then begin
                            curr_colx = curr_colx[0]
                        endif else if nrow eq 1 then begin
                            curr_colx = reform(curr_colx, nel, /overwrite)
                        endif else if nel eq 1 then begin
                            curr_colx = reform(curr_colx, nrow, /overwrite)
                        endif

                        if not keyword_set(old_struct) then begin
                            sz = size(curr_colx)
                            if sz[1] eq 1 then begin
                                tablex.(i) = reform(curr_colx,n_elements(curr_colx))
                            endif else begin
                                tablex.(i) = curr_colx
                            endelse
                        endif else begin
                            mrd_putc, tablex,i, curr_colx
                        endelse
                
                endelse
                
          endif else begin
                ; Now handle the added columns which hold the lengths
                ; of the variable length columns.
                
                ncol = -col - 1 ; Remember we subtracted an extra one.
                if not keyword_set(old_struct) then begin
                    xx = table.(ncol)
                    tablex.(i) = reform(xx[0,*])
                endif else begin
                    xx = mrd_getc(table, ncol)
                    mrd_putc, tablex, i, reform(xx[0,*])
                endelse
          endelse
endfor 
 
; Finally get rid of the initial table and return the table with the 
; variable arrays read in. 
; 
table = temporary(tablex) 
return 
end 

pro mrd_read_table, unit, range, rsize, structyp, nrows, nfld, typarr, table, $
                                      old_struct=old_struct
 
; 
; Read in the binary table information. 
; 
; Unit          Unit to read data from. 
; Range         Desired range 
; Rsize         Size of row. 
; structyp      Structure type. 
; Nfld          Number of fields in structure. 
; Typarr        Field types 
; Table         Table to read information into.
; Old_struct    Use old structures?
; 

if range[0] gt 0 then mrd_skip, unit, rsize*range[0]
readu,unit, table

; Move to the beginning of the heap -- we may have only read some rows of
; the data.
if range[1] lt nrows-1 then begin
    skip_dist = (nrows-range[1]-1)*rsize
    mrd_skip, unit, skip_dist
endif


; If necessary then convert to native format.
if not is_ieee_big() then begin
    for i=0, nfld-1 do begin 
 
        typ = typarr[i] 
        if typ eq 'B' or typ eq 'A'  or typ eq 'X' or typ eq 'L' $ 
           then goto,nxtfld 
        fld = mrd_getc(table,i)
        if typ eq 'I' then byteorder, fld, /htons 
        if typ eq 'J' or typ eq 'P' then byteorder, fld, /htonl 
        if typ eq 'E' or typarr[i] eq 'C' then byteorder, fld, /xdrtof 
        if typ eq 'D' or typarr[i] eq 'M' then begin 
           ieee_to_host, fld 
        endif
        if n_elements(fld) gt 1 then begin
            if not keyword_set(old_struct) then begin
                table.(i) = fld
            endif else begin
                mrd_putc, table, i, fld
            endelse
          endif else begin
            if not keyword_set(old_struct) then begin
                table.(i) = fld[0]
            endif else begin
                mrd_putc, table, i, fld[0]
            endelse
        endelse
         
  nxtfld: 
    endfor 
endif 
end 
; Check the values of TDIM keywords to see that they have valid
; dimensionalities.  If the TDIM keyword is not present or valid
; then the a one-dimensional array with a size given in the TFORM
; keyword is used.

pro mrd_tdim, header, index, flen, arrstr, no_tdim=no_tdim

; HEADER        Current header array.
; Index         Index of current parameter
; flen          Len given in TFORM keyword
; arrstr        String returned to be included within paren's in definition.
; no_tdim       Disable TDIM processing

arrstr = strcompress(string(flen),/remo)

if keyword_set(no_tdim) then return

tdstr = fxpar(header, 'TDIM'+strcompress(string(index),/remo))
if tdstr eq '' then return

;
; Parse the string.  It should be of the form '(n1,n2,...nx)' where
; all of the n's are positive integers and the product equals flen.
;
tdstr = strcompress(tdstr,/remo)
len = strlen(tdstr)
if strmid(tdstr,0,1) ne '(' and strmid(tdstr,len-1,1) ne ')' or len lt 3 then begin
        print, 'MRDFITS: Error: invalid TDIM for column', index
        return
endif

; Get rid of parens.
tdstr = strmid(tdstr,1,len-2)
len = len-2

nind = 0
cnum = 0

for nchr=0, len-1 do begin
        c = strmid(tdstr,nchr, 1)
        
        if c ge '0' and c le '9' then begin
                cnum = 10*cnum + long(c)
                
        endif else if c eq ',' then begin
        
                if cnum le 0 then begin
                        print,'MRDFITS: Error: invalid TDIM for column', index
                        return
                endif
                
                if n_elements(numbs) eq 0 then  $
                        numbs = cnum $
                else    numbs = [numbs,cnum]
                
                cnum = 0
                
       endif else begin
       
                print,'MRDFITS: Error: invalid TDIM for column', index
                return
                
       endelse

endfor

; Handle the last number.
if cnum le 0 then begin
        print,'MRDFITS: Error: invalid TDIM for column', index
        return
endif

if n_elements(numbs) eq 0 then numbs = cnum else numbs = [numbs,cnum]

prod = 1

for i=0, n_elements(numbs)-1 do prod = prod*numbs[i]

if prod ne flen then begin
        print,'MRDFITS: Error: TDIM/TFORM dimension mismatch'
        return
endif

arrstr = tdstr
end
 
pro mrd_table, header, structyp, use_colnum,           $ 
    range, rsize, table, nrows, nfld, typarr, fnames, fvalues,   $ 
    vcls, vtpes, scales, offsets, scaling, status, $
    silent=silent, tempdir=tempdir, columns=columns, no_tdim=no_tdim, $
    old_struct=old_struct
 
; Define a structure to hold a FITS binary table. 
; 
; Header                FITS header for table. 
; Structyp              IDL structure type to be used for 
;                       structure. 
; N_call                Number of times this routine has been called. 
; Table                 Structure to be defined. 
; Status                Return status.
; No_tdim               Disable TDIM processing.

table = 0

types = ['L', 'X', 'B', 'I', 'J', 'A', 'E', 'D', 'C', 'M', 'P']
arrstr = ['bytarr(', 'bytarr(', 'bytarr(', 'intarr(', 'lonarr(',      $ 
          'string(replicate(32b,', 'fltarr(', 'dblarr(', 'complexarr(',    $ 
          'dcomplexarr(', 'lonarr(2*'] 

sclstr = ["'T'", '0B', '0B', '0', '0L', '" "', '0.', '0.d0', 'complex(0.,0.)', $ 
          'dcomplex(0.d0,0.d0)', 'lonarr(2)'] 
 

status = 0 

xten = strmid(fxpar(header, 'XTENSION'),0,8)
if xten ne 'BINTABLE' and xten ne 'A3DTABLE' then begin 
        print, 'MRDFITS: ERROR - Header is not from binary table.' 
        nfld = 0 & status = -1 
        return 
endif 
 
nfld = fxpar(header, 'TFIELDS') 
nrow = fxpar(header, 'NAXIS2')
nrows = nrow
 
if range[0] ge 0 then begin 
        range[0] = range[0] < (nrow-1) 
        range[1] = range[1] < (nrow-1) 
endif else begin 
        range[0] = 0 
        range[1] = nrow - 1 
endelse 
nrow = range[1] - range[0] + 1 

if nrow le 0 then begin
        if not keyword_set(silent) then begin
                print, 'MRDFITS: Binary table. ', $
                 strcompress(string(nfld)), ' columns, no rows.'
        endif
        return
endif

rsize = fxpar(header, 'NAXIS1') 
 
; 
;  Loop over the columns           
 
typarr = strarr(nfld) 
fnames = strarr(nfld) 
fvalues = strarr(nfld) 
dimfld = strarr(nfld)
scales = dblarr(nfld)
offsets = dblarr(nfld)
vcls = intarr(nfld)
vtpes = strarr(nfld)
 
for i=0, nfld-1 do begin 
        istr = strcompress(string(i+1), /remo)
        fname = fxpar(header, 'TTYPE' +  istr) 
        fform = fxpar(header, 'TFORM' +   istr)
        scales[i] = fxpar(header, 'TSCAL'+istr)
        if scales[i] eq 0.d0 then scales[i] = 1.d0
        offsets[i] = fxpar(header, 'TZERO'+istr)
        fname = mrd_dofn(fname,i+1, use_colnum) 
        fname = mrd_chkfn(fname, fnames, i) 
        fnames[i] = fname 
        mrd_doff, fform, dim, ftype
        
        ; Treat arrays of length 1 as scalars.
        
        if dim eq 1 then begin
                dim = 0
        endif else if dim EQ -1 then begin 
                dimfld[i] = -1
        endif else begin
                mrd_tdim, header, i+1, dim, str, no_tdim=no_tdim
                dimfld[i] = str
        endelse
                
        typarr[i] = ftype 
        
        
        ; Find the number of bytes in a bit array. 
 
        if ftype eq 'X' and dim gt 0 then begin
            dim = (dim+7)/8 
            dimfld[i] = strtrim(string(dim),2)
        endif
         
        ; Add in the structure label. 
        ; 
         
        ; Handle variable length columns. 
        if ftype eq 'P' then begin 
 
                if dim ne 0  and dim ne 1 then begin 
                    print, 'MRDFITS: Invalid dimension for variable array column '+string(i+1) 
                    status = -1 
                    return 
                endif 
                ppos = strpos(fform, 'P') 
                vf = strmid(fform, ppos+1, 1); 
                if strpos('LXBIJAEDCM', vf) eq -1 then begin 
                    print, 'MRDFITS: Invalid type for variable array column '+string(i+1) 
                    status = -1 
                    return 
                endif 

                vcls[i] = 1
                vtpes[i] = vf
                dim = 0
                         
        endif 
         
         
        for j=0, n_elements(types) - 1 do begin 
                if ftype eq types[j] then begin 
                 
                        if dim eq 0 then begin 
                                fvalues[i] = sclstr[j] 
                        endif else begin 
                                line = arrstr[j] + dimfld[i] + ')'
                                if ftype eq 'A' then line = line + ')' 
                                fvalues[i] = line 
                        endelse 
                        goto, next_col 
                endif 
        endfor 
         
        print, 'MRDFITS: Invalid format code:',ftype, ' for column ', i+1 
        status = -1 
        return 
next_col: 
endfor 

; Check if there are any variable length columns.  If not then
; undefine vcls and vtpes
w = where(vcls eq 1)
if w[0] eq -1 then begin
        dummy = temporary(vcls)
        dummy = temporary(vtpes)
        dummy = 0
endif

if scaling then begin 
        w = where(scales ne 1.0d0 or offsets ne 0.0d0)
        if w[0] eq -1 then scaling = 0
endif

zero = where(long(dimfld) LT 0L, N_zero)
if N_zero GT 0 then begin
        if N_zero Eq nfld then begin
                print,'MRDFITS: Error - All fields have zero length'
                return
        endif
        for i=0, N_zero-1 do print, $
                'MRDFITS: Table column ' + fnames[zero[i]] + ' has zero length'
        nfld = nfld - N_zero
        good = where(dimfld GE 0)
        fnames = fnames[good]
        fvalues = fvalues[good]
endif

if n_elements(vcls) eq 0  and  (not scaling) and not keyword_set(columns) $
  then begin
        table = mrd_struct(fnames, fvalues, nrow, structyp=structyp, $
          tempdir=tempdir, silent=silent, old_struct=old_struct)
endif else begin
        table = mrd_struct(fnames, fvalues, nrow, tempdir=tempdir, $
          silent=silent, old_struct=old_struct)
endelse

if not keyword_set(silent) then begin
        print, 'MRDFITS: Binary table. ',strcompress(string(nfld)), ' columns by ',  $
          strcompress(string(nrow)), ' rows.'
        if n_elements(vcls) gt 0 then begin
                print, 'MRDFITS: Uses variable length arrays'
        endif
endif
        
status = 0 
return 
 
end 

function mrdfits, file, extension, header,      $
        structyp = structyp,                    $
        use_colnum = use_colnum,                $
        range = range,                          $
        dscale = dscale, fscale=fscale,         $
        silent = silent,                        $
        columns = columns,                      $
        no_tdim = no_tdim,                      $
        tempdir = tempdir,                      $
        error_action = error_action,            $
        old_struct=old_struct,                  $
	compress=compress,                      $
        status=status

  
;
;  Can't use keyword_set since default is 2, not 0.

   if n_elements(error_action) eq 0 then begin
        error_action = 2
   endif
   on_error, error_action

   if float(!version.release lt 5) then old_struct=1
   
; Check positional arguments.

   if n_params() le 0  or n_params() gt 3 then begin
        print, 'MRDFITS: Usage'
        print, '   a=mrdfits(file/unit, [extension, header],         $'
        print, '       /fscale, /dscale, /use_colnum, /silent    $'
        print, '       range=, tempdir=, structyp=, columns=, error_action=, status= )'
        return, 0
   endif
   
   if n_params() eq 1 then extension = 0
   
; Check optional arguments.
;
;  *** Structure name ***

if keyword_set(structyp) then begin
        sz = size(structyp)
        if sz[0] ne 0 then begin
                ; Use first element of array
                structyp = structyp[0]
                sz = size(structyp[0])
        endif
        if sz[1] ne 7 then begin
                print, 'MRDFITS: stucture type must be a string'
                return, 0
        endif
endif

;  *** Use column numbers not names?
if not keyword_set(use_colnum) then use_colnum = 0

;  *** Get only a part of the FITS file.
if keyword_set(range) then begin
        if n_elements(range) eq 2 then arange = range $
        else if n_elements(range) eq 1 then arange = [0,range[0]-1] $
        else if n_elements(range) gt 2 then arange = range[0:1] $
        else if n_elements(range) eq 0 then arange = [-1,-1]
endif else arange = [-1,-1]

arange = long(arange)

; Open the file and position to the appropriate extension then read
; the header.

sz = size(file)
if (sz[0] ne 0) then begin
    print, 'MRDFITS: Vector input not supported'
    return, 0
endif

inputUnit = 0
if sz[1] gt 0 and sz[1] lt 4 then begin
    inputUnit = 1
    unit = file
    
    if fxmove(unit,extension) lt 0 then begin
        return, -1
    endif
    
endif else begin 
    unit = fxposit(file, extension, compress=compress, /readonly)

    if unit lt 0 then begin
        print, 'MRDFITS: File access error'
        status = -1
        return, 0
    endif
endelse

if eof(unit) then begin
        print,'MRDFITS: Extension past EOF'
	if inputUnit eq 0 then begin
            free_lun, unit
	endif
	
	status = -2
	
        return, 0
endif

mrd_hread, unit, header, status
if status lt 0 then begin
        print, 'MRDFITS: Unable to read header for extension'
	if inputUnit eq 0 then begin
            free_lun, unit
	endif
        return, 0
endif

; If this is primary array then XTENSION will have value
; 0 which will be converted by strtrim to '0'

xten = strtrim( fxpar(header,'XTENSION'), 2)
if xten eq '0' or xten eq 'IMAGE' then begin
        type = 0
endif else if xten eq 'TABLE' then begin
	type = 1
endif else if xten eq 'BINTABLE'  or xten eq 'A3DTABLE' then begin
	type = 2
endif else begin 
        print, 'MRDFITS: Unable to process extension type:', xten
	if inputUnit eq 0 then begin
	        free_lun,unit
	endif
	status = -1
        return, 0
endelse

scaling = keyword_set(fscale) or keyword_set(dscale)

if type eq 0 then begin

        ;*** Images/arrays
        
        mrd_image, header, arange, maxd, rsize, table, scales, offsets, status, $
          silent=silent, tempdir=tempdir
        if status ge 0 and rsize gt 0 then mrd_read_image, unit, arange, maxd, rsize, table
        size = rsize
        
endif else if type eq 1 then begin

        ;*** ASCII tables.
        
        mrd_ascii, header, structyp, use_colnum,   $
            arange, table, nbytes, nrows, nfld,  $
            typarr, posarr, lenarr, nullarr, fnames, fvalues, $
            scales, offsets, scaling, status, silent=silent, tempdir=tempdir, $
            columns=columns,old_struct=old_struct
        size = nbytes*nrows
        
        if status ge 0   and  size gt 0  then begin
        
                ;*** Read data.
                mrd_read_ascii, unit,  arange, nbytes, nrows, $
                  nfld, typarr, posarr, lenarr, nullarr, table, $
                  old_struct=old_struct
              
                ;*** Extract desired columns.
                if status ge 0 and keyword_set(columns) then $
                        mrd_columns, table, columns, fnames, fvalues, vcls, vtps, $
                          scales, offsets, $
                          scaling, structyp=structyp, tempdir=tempdir, $
                          silent=silent, old_struct=old_struct
        endif
        
endif else begin

     ; *** Binary tables.

     mrd_table, header, structyp, use_colnum,   $
         arange, rsize, table, nrows, nfld, typarr,        $ 
         fnames, fvalues, vcls, vtpes, scales, offsets, scaling, status,   $
         silent=silent, tempdir=tempdir, columns=columns, no_tdim=no_tdim, $
         old_struct=old_struct


     size = nfld*(arange[1] - arange[0] + 1)
     if status ge 0  and  size gt 0  then begin
     
        ;*** Read data.
        mrd_read_table, unit, arange, rsize,  $
          structyp, nrows, nfld, typarr, table, old_struct=old_struct

        if status ge 0 and keyword_set(columns) then $
        
           ;*** Extract desired columns.
           mrd_columns, table, columns, fnames, fvalues, $
                vcls, vtpes, scales, offsets, scaling, structyp=structyp, $
                tempdir=tempdir, silent=silent, old_struct=old_struct
         
     
        if status ge 0 and n_elements(vcls) gt 0 then begin 
        
                ;*** Get variable length columns
                mrd_read_heap, unit, header, arange, fnames, fvalues, $
                 vcls, vtpes, table, structyp, scaling, scales, offsets, status, $
                 silent=silent, tempdir=tempdir, old_struct=old_struct
		
	endif else begin

	        ; Skip remainder of last data block
	        sz = fxpar(header, 'NAXIS1')*fxpar(header,'NAXIS2') +  $
		     fxpar(header, 'PCOUNT')
		skipB = 2880 - sz mod 2880
		if (skipB ne 2880) then mrd_skip, unit, skipB
        endelse
		     
      endif

endelse


if unit gt 0 and inputUnit eq 0 then begin
        free_lun, unit
endif

if  status ge 0  and  scaling  and  size gt 0  then begin
        w = where(scales ne 1.d0  or  offsets ne 0.0d0)
        
        ;*** Apply scalings.
        if w[0] ne -1 then mrd_scale, type, scales, offsets, table, header,$
            fnames, fvalues, 1+arange[1]-arange[0], structyp=structyp, $
            dscale=dscale, tempdir=tempdir, silent=silent, old_struct=old_struct
endif

                              
if status ge 0 then return, table else return, 0

end

