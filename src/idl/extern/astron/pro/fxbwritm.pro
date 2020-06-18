	PRO FXBWRITM, UNIT, COL, $
                      D0,  D1,  D2,  D3,  D4,  D5,  D6,  D7,  D8,  D9, $
                      D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, $
                      D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, $
                      D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, $
                      D40, D41, D42, D43, D44, D45, D46, D47, D48, D49, $
                      NOIEEE=NOIEEE, $
                      ROW=ROW, NANVALUE=NANVALUE, BUFFERSIZE=BUFFERSIZE, $
                      ERRMSG=ERRMSG, WARNMSG=WARNMSG, STATUS=OUTSTATUS
;+
; NAME: 
;	FXBWRITM
; Purpose     : 
;       Write multiple columns/rows to a disk FITS binary table file.
; Explanation : 
;       A call to FXBWRITM will write multiple rows and multiple
;       columns to a binary table in a single procedure call.  Up to
;       fifty columns may be read in a single pass.  The file should
;       have already been opened with FXBOPEN (with write access) or
;       FXBCREATE.  FXBWRITM optimizes writing multiple columns by
;       first writing a large chunk of data to the FITS file all at
;       once.  FXBWRITM cannot write variable-length arrays; use
;       FXBWRITE instead.
; Use         : 
;	FXBWRITM, UNIT, COL, D0, D1, D2, ..., ROW=ROW
; Inputs      : 
;	UNIT	= Logical unit number corresponding to the file containing the
;		  binary table.
;	D0,...	= An IDL data array to be written to the file, one for
;                 each column.
;	COL	= Column in the binary table to place data in, starting from
;		  column one.
; Opt. Inputs : 
;	None.
; Outputs     : 
;	None.
; Opt. Outputs: 
;	None.
; Keywords    : 
;	ROW	= Either row number in the binary table to writedata to,
;		  starting from row one, or a two element array containing a
;		  range of row numbers to write.  If not passed, then
;		  the entire column is written.
;	NANVALUE= Value signalling data dropout.  All points corresponding to
;		  this value are set to be IEEE NaN (not-a-number).  Ignored
;		  unless DATA is of type float, double-precision or complex.
;	ERRMSG	= If defined and passed, then any error messages will be
;		  returned to the user in this parameter rather than
;		  depending on the MESSAGE routine in IDL.  If no errors are
;		  encountered, then a null string is returned.  In order to
;		  use this feature, ERRMSG must be defined first, e.g.
;
;			ERRMSG = ''
;			FXBWRITE, ERRMSG=ERRMSG, ...
;			IF ERRMSG NE '' THEN ...
;       WARNMSG = Messages which are considered to be non-fatal
;                 "warnings" are returned in this  output string.
;       BUFFERSIZE = Data are transferred in chunks to conserve
;                 memory.  This is the size in bytes of each chunk.
;                 If a value of zero is given, then all of the data
;                 are transferred in one pass.  Default is 32768 (32
;                 kB).
;       STATUS  = An output array containing the status for each
;                 read, 1 meaning success and 0 meaning failure.
;
; Calls       : 
;	HOST_TO_IEEE
; Common      : 
;	Uses common block FXBINTABLE--see "fxbintable.pro" for more
;	information.
; Restrictions: 
;	The binary table file must have been opened with FXBCREATE or
;       FXBOPEN (with write access).
;
;	The data must be consistent with the column definition in the binary
;	table header.
;
;	The row number must be consistent with the number of rows stored in the
;	binary table header.
;
; Side effects: 
;	None.
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;       C. Markwardt, based on FXBWRITE and FXBREADM (ver 1), Jan 1999
; Written     : 
;	Craig Markwardt, GSFC, January 1999.
; Modified    :
;       Version 1, Craig Markwardt, GSFC 18 January 1999.
;               Documented this routine, 18 January 1999. 
; Version     :
;       Version 1, 18 January 1999.
;-
;
@fxbintable
	ON_ERROR, 2
;
;  Check the number of parameters.
;
        IF N_PARAMS() LT 2 THEN BEGIN
		MESSAGE = 'Syntax:  FXBWRITM, UNIT, COL, DATA1, DATA2, ' $
                  +' ..., ROW=ROW'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF
        IF N_ELEMENTS(BUFFERSIZE) EQ 0 THEN BUFFERSIZE = 32768L

;
;  COL may be one of several descriptors:
;     * a list of column numbers, beginning with 1
;     * a list of column names
;
        MYCOL = [ COL ]    ; Make sure it is an array

        SC = SIZE(MYCOL)
        NUMCOLS = N_ELEMENTS(MYCOL)
        OUTSTATUS = LONARR(NUMCOLS)

;
;  Find the logical unit number in the FXBINTABLE common block.
;
	ILUN = WHERE(LUN EQ UNIT,NLUN)
	ILUN = ILUN[0]
	IF NLUN EQ 0 THEN BEGIN
		MESSAGE = 'Unit ' + STRTRIM(UNIT,2) +	$
			' not opened properly'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF

;
;  Make sure the file was opened for write access.
;
	IF STATE[ILUN] NE 2 THEN BEGIN
		MESSAGE = 'Unit ' + STRTRIM(UNIT,2) +	$
			' not opened for write access'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF

;
;  Check the number of columns.  It should be fewer than 50
;
        IF NUMCOLS GT 50 THEN BEGIN
            MESSAGE = 'Maximum of 50 columns exceeded'
            IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                ERRMSG = MESSAGE
                RETURN
            END ELSE MESSAGE, MESSAGE
        ENDIF
;        Commented out because too much data is not a problem
;        IF NUMCOLS LT N_PARAMS()-2 THEN BEGIN
;            MESSAGE = 'ERROR: too few data parameters passed'
;            IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
;                ERRMSG = MESSAGE
;                RETURN
;            END ELSE MESSAGE, MESSAGE
;        ENDIF

        ICOL    = LONARR(NUMCOLS)
        FOUND   = BYTARR(NUMCOLS)
        NOTFOUND = ''
        NNOTFOUND = 0L
        IF N_ELEMENTS(WARNMSG) NE 0 THEN WARNMSG = ''

;
;  If COL is of type string, then search for a column with that label.
;
        IF SC[SC[0]+1] EQ 7 THEN BEGIN
            MYCOL = STRUPCASE(STRTRIM(MYCOL,2))
            FOR I = 0, NUMCOLS-1 DO BEGIN
		XCOL = WHERE(TTYPE[*,ILUN] EQ MYCOL[I], NCOL)
		ICOL[I] = XCOL[0]
                IF NCOL GT 0 THEN FOUND[I] = 1
                IF NOT FOUND[I] THEN BEGIN
                    IF NOTFOUND EQ '' THEN NOTFOUND = MYCOL[I] $
                    ELSE NOTFOUND = NOTFOUND +', ' + MYCOL[I]
                    NNOTFOUND = NNOTFOUND + 1
                ENDIF
            ENDFOR

            IF NNOTFOUND EQ NUMCOLS THEN BEGIN
                MESSAGE = 'ERROR: None of the requested columns were found'
                IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                    ERRMSG = MESSAGE
                    RETURN
                END ELSE MESSAGE, MESSAGE
            ENDIF ELSE IF NNOTFOUND GT 0 THEN BEGIN
                MESSAGE = 'WARNING: Columns ' + NOTFOUND + ' were not found'
                IF N_ELEMENTS(WARNMSG) NE 0 THEN WARNMSG = MESSAGE $
                ELSE MESSAGE, MESSAGE, /INFO
            ENDIF
;
;  Otherwise, a numerical column was passed.  Check its value.
;
	ENDIF ELSE BEGIN
            ICOL[*] = LONG(MYCOL) - 1
        ENDELSE
                

;
;  Step through each column index, and check for validity
;
        MESSAGE = ''
        FOR I = 0, NUMCOLS-1 DO BEGIN
            IF NOT FOUND[I] THEN GOTO, LOOP_END_COLCHECK

            IF (ICOL[I] LT 0) OR (ICOL[I] GE TFIELDS[ILUN]) THEN BEGIN
                MESSAGE = 'COL "'+STRTRIM(MYCOL[I],2)+$
                  '" must be between 1 and ' +	$
                  STRTRIM(TFIELDS[ILUN],2)
                FOUND[I] = 0
            ENDIF
;
;  If there are no elements in the array, then set !ERR to -1.
;
            IF FOUND[I] AND N_ELEM[ICOL[I],ILUN] EQ 0 THEN BEGIN
                FOUND[I] = 0
                MESSAGE = MESSAGE + '; Number of elements to write in "'+$
                  STRTRIM(MYCOL[I],2)+'" should be zero'
            ENDIF

;
;  Do not permit variable-length columns
;
            IF MAXVAL[ICOL[I],ILUN] GT 0 THEN BEGIN
                MESSAGE = MESSAGE + 'FXBWRITM cannot write ' +	$
                  'variable-length column "'+STRTRIM(MYCOL[I],2)+'"'
                FOUND[I] = 0
            ENDIF

            LOOP_END_COLCHECK:

        ENDFOR
;
;  If ROW was not passed, then set it equal to the entire range.  Otherwise,
;  extract the range.
;
        IF N_ELEMENTS(ROW) EQ 0 THEN BEGIN
            ROW = [1L, NAXIS2[ILUN]]
        ENDIF
	CASE N_ELEMENTS(ROW) OF
		1:  ROW2 = LONG(ROW[0])
		2:  ROW2 = LONG(ROW[1])
		ELSE:  BEGIN
			MESSAGE = 'ROW must have one or two elements'
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
			END
	ENDCASE
	ROW1 = LONG(ROW[0])

;
;  If ROW represents a range, then make sure that the row range is legal, and
;  that reading row ranges is allowed (i.e., the column is not variable length.
;
	IF ROW1 NE ROW2 THEN BEGIN
		MAXROW = NAXIS2[ILUN]
		IF (ROW1 LT 1) OR (ROW1 GT MAXROW) THEN BEGIN
			MESSAGE = 'ROW[0] must be between 1 and ' +	$
				STRTRIM(MAXROW,2)
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		END ELSE IF (ROW2 LT ROW1) OR (ROW2 GT MAXROW) THEN BEGIN
			MESSAGE = 'ROW[1] must be between ' +	$
				STRTRIM(ROW1,2) + ' and ' + STRTRIM(MAXROW,2)
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		ENDIF
;
;  Otherwise, if ROW is a single number, then just make sure it's valid.
;
	END ELSE BEGIN
		IF (ROW1 LT 1) OR (ROW1 GT NAXIS2[ILUN]) THEN BEGIN
			MESSAGE = 'ROW must be between 1 and ' +	$
				STRTRIM(NAXIS2[ILUN],2)
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		ENDIF
	ENDELSE

;
;  Check the type of the data against that defined for this column.
;
        COLNDIM = LONARR(NUMCOLS)
        COLDIM  = LONARR(NUMCOLS, 20) ;; Maximum of 20 dimensions in output
        COLTYPE = LONARR(NUMCOLS)
        BOFF1   = LONARR(NUMCOLS)
        BOFF2   = LONARR(NUMCOLS)
        NOUTPUT = LONARR(NUMCOLS)
        NROWS = ROW2-ROW1+1
        MESSAGE = ''
        DTYPENAMES = [ 'BAD TYPE', 'BYTE', 'FIX', 'LONG', $
                       'FLOAT', 'DOUBLE', 'COMPLEX', 'STRING', $
                       'BAD TYPE', 'DCOMPLEX' ]
        FOR I = 0L, NUMCOLS-1 DO BEGIN

            IF NOT FOUND[I] THEN GOTO, LOOP_END_DIMS
            ;;  Data type of the input.
            COLTYPE[I] = IDLTYPE[ICOL[I],ILUN]

            SZ = 0
            RESULT = EXECUTE('SZ = SIZE(D'+STRTRIM(I,2)+')')
            IF RESULT EQ 0 THEN BEGIN
                MESSAGE = MESSAGE + '; Could not extract type info (column '+$
                  STRTRIM(MYCOL[I],2)+')'
                FOUND[I] = 0
            ENDIF

            TYPE = SZ[SZ[0]+1]
            IF TYPE NE COLTYPE[I] THEN BEGIN
                CASE COLTYPE[I] OF
                    1: STYPE = 'byte'
                    2: STYPE = 'short integer'
                    3: STYPE = 'long integer'
                    4: STYPE = 'floating point'
                    5: STYPE = 'double precison'
                    6: STYPE = 'complex'
                    7: STYPE = 'string'
                ENDCASE
                FOUND[I] = 0
                MESSAGE = '; Data type (column '+STRTRIM(MYCOL[I],2)+$
                  ') should be ' + STYPE
            ENDIF

            DIMS = N_DIMS[*,ICOL[I],ILUN]
            NDIMS = DIMS[0]
            DIMS  = DIMS[1:NDIMS]

            IF NDIMS EQ 1 AND DIMS[0] EQ 1 THEN BEGIN

                ;; Case of only one output element, try to return a
                ;; scalar.  Otherwise, it is a vector equal to the
                ;; number of rows to be read

                COLNDIM[I] = 1L
                COLDIM[I,0] = NROWS
            ENDIF ELSE BEGIN

                COLNDIM[I] = NDIMS
                COLDIM[I,0:(NDIMS-1)] = DIMS
                IF NROWS GT 1 THEN BEGIN
                    COLDIM[I,NDIMS] = NROWS
                    COLNDIM[I] = COLNDIM[I]+1
                ENDIF

            ENDELSE

;
;  Check the number of elements in the input
;
            NOUTP = ROUND(PRODUCT(COLDIM[I,0:COLNDIM[I]-1]))
            IF SZ[SZ[0]+1] EQ 7 THEN BEGIN
                NOUTP = NOUTP / COLDIM[I,0]
                IF NOUTP NE SZ[SZ[0]+2] THEN GOTO, ERR_NELEM
                NOUTPUT[I] = NOUTP
            ENDIF ELSE IF SZ[SZ[0]+2] NE NOUTP THEN BEGIN
                ERR_NELEM:
                MESSAGE = MESSAGE+'; Data array (column '+STRTRIM(MYCOL[I],2)+$
                  ') should have ' + STRTRIM(LONG(NOUTP),2) + ' elements'
                FOUND[I] = 0
            ENDIF ELSE NOUTPUT[I] = NOUTP

            ;; Byte offsets
            BOFF1[I] = BYTOFF[ICOL[I],ILUN]
            IF ICOL[I] EQ TFIELDS[ILUN]-1 THEN BOFF2[I] = NAXIS1[ILUN]-1 $
            ELSE BOFF2[I] = BYTOFF[ICOL[I]+1,ILUN]-1

            LOOP_END_DIMS:

        ENDFOR

;
;  Check to be sure that there are columns to be written
;
        W = WHERE(FOUND EQ 1, COUNT)
        IF COUNT EQ 0 THEN BEGIN
            STRPUT, MESSAGE, ':', 0
            MESSAGE = 'ERROR: No requested columns could be written'+MESSAGE
            IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                ERRMSG = MESSAGE
                RETURN
            END ELSE MESSAGE, MESSAGE
        ENDIF ELSE IF MESSAGE NE '' THEN BEGIN
            STRPUT, MESSAGE, ':', 0
            MESSAGE = 'WARNING: Some columns could not be written'+MESSAGE
            IF N_ELEMENTS(WARNMSG) NE 0 THEN WARNMSG = MESSAGE $
            ELSE MESSAGE, MESSAGE, /INFO
        ENDIF

        ;; I construct a list of unique column names here.  Why?
        ;; Because if *all* the columns are named, then there is no
        ;; need to read the data from disk first.  Since columns can
        ;; be given more than once in MYCOL, we need to uniq-ify it.
        CC = MYCOL[UNIQ(MYCOL, SORT(MYCOL))]
        NC = N_ELEMENTS(CC)

;
;  Find the position of the first byte of the data array in the file.
;
	OFFSET0 = NHEADER[ILUN] + NAXIS1[ILUN]*(ROW1-1)

        POS = 0L
        NROWS0 = NROWS
        J = 0L
        ;; Here, we constrain the buffer to be at least 16 rows long.
        ;; If we fill up 32 kB with fewer than 16 rows, then there
        ;; must be a lot of (big) columns in this table.  It's
        ;; probably a candidate for using FXBREAD instead.
        BUFFROWS = LONG((BUFFERSIZE/NAXIS1[ILUN]) > 16L)
        IF BUFFERSIZE LE 0 THEN BUFFROWS = NROWS0

;
;  Loop through the data in chunks
;
        WHILE NROWS GT 0 DO BEGIN
        J = J + 1
        NR  = NROWS < BUFFROWS
        OFFSET1 = NAXIS1[ILUN]*POS
;
;  Proceed by reading a byte array from the input data file
;  FXBREADM reads all columns from the specified rows, and
;  sorts out the details of which bytes belong to which columns
;  in the next FOR loop.
;
        BB = BYTARR(NAXIS1[ILUN], NR)
;  If *all* columns are being filled, then there is no reason to 
;  read from the file
        
        IF NC LT TFIELDS[ILUN] THEN BEGIN
            POINT_LUN,UNIT,OFFSET0+OFFSET1
            READU, UNIT, BB
        ENDIF

;
;  Now select out the desired columns to write
;
        FOR I = 0, NUMCOLS-1 DO BEGIN
            IF NOT FOUND[I] THEN GOTO, LOOP_END_WRITE

            ;; Copy data into DD
            DD = 0
            RESULT = EXECUTE('DD = D'+STRTRIM(I,2))
            IF RESULT EQ 0 THEN BEGIN
                GOTO, LOOP_END_WRITE
;                MESSAGE = 'ERROR: Could not get data (column '+$
;                  STRTRIM(MYCOL(I),2)+')'
;                IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
;                    ERRMSG = MESSAGE
;                    RETURN
;                ENDIF ELSE MESSAGE, MESSAGE
            ENDIF
            IF N_ELEMENTS(DD) EQ 1 THEN DD = [DD]
            DD = REFORM(DD, NOUTPUT[I]/NROWS0, NROWS0)
            IF POS GT 0 OR NR LT NROWS0 THEN $
              DD = DD[*,POS:(POS+NR-1)]

            ;; Now any conversions to FITS format must be done
            COUNT = 0L
            CT = COLTYPE[I]
            CASE 1 OF
                ;; Integer types
                (CT EQ 1): CT = CT
                (CT EQ 2 OR CT EQ 3): BEGIN
                    IF NOT KEYWORD_SET(NOIEEE) THEN HOST_TO_IEEE, DD 
                END

                ;; Floating and complex types
                (CT GE 4 AND CT LE 6 OR CT EQ 9): BEGIN
                    IF NOT KEYWORD_SET(NOIEEE) THEN BEGIN
                        IF N_ELEMENTS(NANVALUE) EQ 1 THEN BEGIN
                            W=WHERE(DD EQ NANVALUE,COUNT)
                            NAN = REPLICATE('FF'XB,16)
                            NAN = CALL_FUNCTION(DTYPENAMES,NAN,0,1)
                        ENDIF
                        HOST_TO_IEEE, DD
                        IF COUNT GT 0 THEN DD[W] = NAN
                    ENDIF
                END

                ;; String type, needs to be padded with spaces
                (CT EQ 7): BEGIN
                    N_CHAR = N_DIMS[1,ICOL[I],ILUN]
                    ;; Largest string determines size of array
                    MAXLEN = MAX(STRLEN(DD))  
                    ;; Convert to bytes
                    DD = BYTE(TEMPORARY(DD))
                    IF N_ELEMENTS(DD) EQ 1 THEN DD = [DD]
                    DD = REFORM(DD, MAXLEN, NR, /OVERWRITE)

                    ;; Put it into the output array
                    IF MAXLEN GT N_CHAR THEN BEGIN
                        DD = DD[0:(N_CHAR-1),*]
                    ENDIF ELSE BEGIN
                        DB = BYTARR(N_CHAR, NR)
                        DB[0:(MAXLEN-1),*] = TEMPORARY(DD)
                        DD = TEMPORARY(DB)
                    ENDELSE

                    ;; Pad any zeroes with spaces
                    WB = WHERE(DD EQ 0b, WCOUNT)
                    IF WCOUNT GT 0 THEN DD[WB] = 32B
                    
                    ;; Pretend that it is a byte array
                    CT = 1
                END
            ENDCASE
            IF CT NE 1 THEN $
              DD = BYTE(TEMPORARY(DD),0,(BOFF2[I]-BOFF1[I]+1),NR)
            IF N_ELEMENTS(DD) EQ 1 THEN DD = [DD]
            DD = REFORM(DD, BOFF2[I]-BOFF1[I]+1, NR, /OVERWRITE)
            
            ;; Now place the data into the byte array
            BB[BOFF1[I],0] = DD
            
            OUTSTATUS[I] = 1
            LOOP_END_WRITE:
        END

        ;; Finally, write byte array to output file
        POINT_LUN, UNIT, OFFSET0+OFFSET1
        BB = REFORM(BB, N_ELEMENTS(BB), /OVERWRITE)
        WRITEU, UNIT, BB

        NROWS = NROWS - NR
        POS   = POS + NR
        ENDWHILE

;
	IF N_ELEMENTS(ERRMSG) NE 0 THEN ERRMSG = ''
	RETURN
	END
