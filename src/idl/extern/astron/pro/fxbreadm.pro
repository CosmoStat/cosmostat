PRO FXBREADM, UNIT, COL, $
              D0,  D1,  D2,  D3,  D4,  D5,  D6,  D7,  D8,  D9, $
              D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, $
              D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, $
              D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, $
              D40, D41, D42, D43, D44, D45, D46, D47, D48, D49, $
              ROW=ROW, VIRTUAL=VIR, DIMENSIONS=DIM, $
              NOSCALE=NOSCALE, NOIEEE=NOIEEE, $
              NANVALUE=NANVALUE, BUFFERSIZE=BUFFERSIZE, $
              ERRMSG=ERRMSG, WARNMSG=WARNMSG, STATUS=OUTSTATUS

;+
; NAME: 
;	FXBREADM
; Purpose     : 
;	Read multiple columns/rows from a disk FITS binary table file.
; Explanation : 
;       A call to FXBREADM will read data from multiple rows and
;       multiple columns in a single procedure call.  Up to fifty
;       columns may be read in a single pass; the number of rows is
;       limited essentially by available memory.  The file should have
;       already been opened with FXBOPEN.  FXBREADM optimizes reading
;       multiple columns by first reading a large chunk of data from
;       the FITS file directly, and then slicing the data into columns
;       within memory.  FXBREADM cannot read variable-length arrays;
;       use FXBREAD instead.
; Use         : 
;	FXBREADM, UNIT, COL, DATA1, DATA2, ... [, ROW=ROW ]
; Inputs      : 
;	UNIT	= Logical unit number corresponding to the file containing the
;		  binary table.
;       COL     = An array of columns in the binary table to read data
;                 from, either as character strings containing column
;                 labels (TTYPE), or as numerical column indices
;                 starting from column one.
; Opt. Inputs : 
;       None.
; Outputs     : 
;	D0, ... = A named variable to accept the data values, one for
;                 each column.  The columns are stored in order of the
;                 list in COL.  If the read operation fails for a
;                 particular column, then the corresponding output Dn
;                 variable is not altered.  See the STATUS keyword.
;
; Opt. Outputs: 
;	None.
; Keywords    : 
;	ROW	= Either row number in the binary table to read data from,
;		  starting from row one, or a two element array containing a
;		  range of row numbers to read.  If not passed, then the entire
;		  column is read in.
;	NOSCALE	= If set, then the ouput data will not be scaled using the
;		  optional TSCAL and TZERO keywords in the FITS header.
;		  Default is to scale.
;	VIRTUAL	= If set, and COL is passed as a name rather than a number,
;		  then if the program can't find a column with that name, it
;		  will then look for a keyword with that name in the header.
;		  Such a keyword would then act as a "virtual column", with the
;		  same value for every row.
;	DIMENSIONS = FXBREADM ignores this keyword.  It is here for
;	          compatibility only.
;	NANVALUE= Value signalling data dropout.  All points corresponding to
;		  IEEE NaN (not-a-number) are converted to this number.
;		  Ignored unless DATA is of type float, double-precision or
;		  complex.
;	ERRMSG	= If defined and passed, then any error messages will be
;		  returned to the user in this parameter rather than
;		  depending on the MESSAGE routine in IDL.  If no errors are
;		  encountered, then a null string is returned.  In order to
;		  use this feature, ERRMSG must be defined first, e.g.
;
;			ERRMSG = ''
;			FXBREAD, ERRMSG=ERRMSG, ...
;			IF ERRMSG NE '' THEN ...
;       WARNMSG = Messages which are considered to be non-fatal
;                 "warnings" are returned in this  output string.
;       BUFFERSIZE = Raw data are transferred from the file in chunks
;                 to conserve memory.  This is the size in bytes of
;                 each chunk.  If a value of zero is given, then all
;                 of the data are transferred in one pass.  Default is
;                 32768 (32 kB).
;       STATUS  = An output array containing the status for each
;                 column read, 1 meaning success and 0 meaning failure.
;
; Calls       : 
;	IEEE_TO_HOST, FXPAR, WHERENAN
; Common      : 
;	Uses common block FXBINTABLE--see "fxbintable.pro" for more
;	information.
; Restrictions: 
;	The binary table file must have been opened with FXBOPEN.
;
;	The data must be consistent with the column definition in the binary
;	table header.
;
;	The row number must be consistent with the number of rows stored in the
;	binary table header.
;
;       No variable-length columns may be read with FXBREADM.
;
;       Generaly speaking, FXBREADM will be faster than iterative
;       calls to FXBREAD when (a) a large number of columns is to be
;       read or (b) the size in bytes of each cell is small, so that
;       the overhead of the FOR loop in FXBREAD becomes significant.
;
; Side effects: 
;	If there are no elements to read in (the number of elements is zero),
;	then the program sets !ERR to -1, and DATA is unmodified.
;
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;       C. Markwardt, January 1999, based in concept on FXBREAD version 12 from
;                              IDLASTRO, but with significant and
;                              major changes to accomodate the
;                              multiple row/column technique.  Mostly
;                              the parameter checking and general data
;                              flow remain.
;
; Written     : 
;       Craig Markwardt, GSFC, January 1999.
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
		MESSAGE = 'Syntax:  FXBREADM, UNIT, COL, D0, D1, ... [, ROW= ]'
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
;  Check the number of columns.  It should be fewer than 50
;
        IF NUMCOLS GT 50 THEN BEGIN
            MESSAGE = 'Maximum of 50 columns exceeded'
            IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                ERRMSG = MESSAGE
                RETURN
            END ELSE MESSAGE, MESSAGE
        ENDIF
        IF NUMCOLS LT N_PARAMS()-2 AND N_ELEMENTS(ERRMSG) EQ 0 THEN BEGIN
            MESSAGE, 'WARNING: number of data parameters less than columns', $
              /INFO
        ENDIF
            
        ICOL    = LONARR(NUMCOLS)
        VIRTUAL = BYTARR(NUMCOLS)
        VIRTYPE = LONARR(NUMCOLS)
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
;
;  If the column was not found, and VIRTUAL was set, then search for a keyword
;  by that name.
;
                IF NCOL GT 0 THEN FOUND[I] = 1
                IF NOT FOUND[I] AND KEYWORD_SET(VIR) THEN BEGIN
                    HEADER = HEAD[*,ILUN]
                    VALUE = FXPAR(HEADER,MYCOL[I])
                    IF !ERR GE 0 THEN BEGIN
                        RESULT = EXECUTE('D'+STRTRIM(I,2)+$
                                         ' = VALUE')
                        SV = SIZE(VALUE)
                        VIRTYPE[I] = SV[SV[0]+1]
                        VIRTUAL[I] = 1
                        FOUND[I] = 1
                    ENDIF
                ENDIF ELSE IF NOT FOUND[I] THEN BEGIN
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

;  Step through each column index
        MESSAGE = ''
        FOR I = 0, NUMCOLS-1 DO BEGIN
            IF NOT FOUND[I] THEN GOTO, LOOP_END_COLCHECK
            IF VIRTUAL[I] THEN GOTO, LOOP_END_COLCHECK

            IF (ICOL[I] LT 0) OR (ICOL[I] GE TFIELDS[ILUN]) THEN BEGIN
                MESSAGE = MESSAGE + '; COL "'+STRTRIM(MYCOL(I),2)+$
                  '" must be between 1 and ' +	$
                  STRTRIM(TFIELDS[ILUN],2)
                FOUND[I] = 0
;                IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
;                    ERRMSG = MESSAGE
;                    RETURN
;                END ELSE MESSAGE, MESSAGE
            ENDIF
;
;  If there are no elements in the array, then set !ERR to -1.
;
            IF FOUND[I] AND N_ELEM[ICOL[I],ILUN] EQ 0 THEN BEGIN
                FOUND[I] = 0
                MESSAGE = MESSAGE + '; Number of elements to read in "'+$
                  STRTRIM(MYCOL[I],2)+'" is zero'
;                !ERR = -1
;                RETURN
            ENDIF

;
;  Do not permit variable-length columns
;
            IF MAXVAL[ICOL[I],ILUN] GT 0 THEN BEGIN
                MESSAGE = MESSAGE + '; FXBREADM cannot read ' + $
                  'variable-length column "'+STRTRIM(MYCOL[I],2)+'"'
                FOUND[I] = 0
;                IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
;                    ERRMSG = MESSAGE
;                    RETURN
;                END ELSE MESSAGE, MESSAGE
            ENDIF

            LOOP_END_COLCHECK:

        ENDFOR

;
;  Check to be sure that there are columns to be read
;
        W = WHERE(FOUND EQ 1, COUNT)
        IF COUNT EQ 0 THEN BEGIN
            STRPUT, MESSAGE, ':', 0
            MESSAGE = 'ERROR: No requested columns could be read'+MESSAGE
            IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                ERRMSG = MESSAGE
                RETURN
            END ELSE MESSAGE, MESSAGE
        ENDIF ELSE IF MESSAGE NE '' THEN BEGIN
            STRPUT, MESSAGE, ':', 0
            MESSAGE = 'WARNING: Some columns could not be read'+MESSAGE
            IF N_ELEMENTS(WARNMSG) NE 0 THEN WARNMSG = MESSAGE $
            ELSE MESSAGE, MESSAGE, /INFO
        ENDIF
            
;
;  If ROW was not passed, then set it equal to the entire range.  Otherwise,
;  extract the range.
;
        IF N_ELEMENTS(ROW) EQ 0 THEN ROW = [1L, NAXIS2[ILUN]]
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
;  Compose information about the output
;
        COLNDIM = LONARR(NUMCOLS)
        COLDIM  = LONARR(NUMCOLS, 20) ;; Maximum of 20 dimensions in output
        COLTYPE = LONARR(NUMCOLS)
        BOFF1   = LONARR(NUMCOLS)
        BOFF2   = LONARR(NUMCOLS)
        NROWS = ROW2-ROW1+1
        DTYPENAMES = [ 'BAD TYPE', 'BYTE', 'FIX', 'LONG', $
                       'FLOAT', 'DOUBLE', 'COMPLEX', 'STRING', $
                       'BAD TYPE', 'DCOMPLEX' ]
        FOR I = 0L, NUMCOLS-1 DO BEGIN

            IF NOT FOUND[I] THEN GOTO, LOOP_END_DIMS
            ;;  Data type of the input.
            IF VIRTUAL[I] THEN BEGIN
                COLTYPE[I] = VIRTYPE[I] 
                GOTO, LOOP_END_DIMS
            ENDIF ELSE $
              COLTYPE[I] = IDLTYPE[ICOL[I],ILUN]

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
            
            ;; For strings, the number of characters is the first
            ;; dimension.  This information is useless to us now,
            ;; since the STRING() type cast which will appear below
            ;; handles the array conversion automatically.
            IF COLTYPE[I] EQ 7 THEN BEGIN
                IF COLNDIM[I] GT 1 THEN BEGIN
                    COLDIM[I,0:COLNDIM[I]-2] = COLDIM[I,1:COLNDIM[I]-1]
                    COLDIM[I,COLNDIM[I]-1]   = 0
                    COLNDIM[I] = COLNDIM[I] - 1
                ENDIF ELSE BEGIN  ;; Case of a single row
                    COLNDIM[I] = 1L
                    COLDIM[I,0] = NROWS
                ENDELSE
            ENDIF

            ;; Byte offsets
            BOFF1[I] = BYTOFF[ICOL[I],ILUN]
            IF ICOL[I] EQ TFIELDS[ILUN]-1 THEN BOFF2[I] = NAXIS1[ILUN]-1 $
            ELSE BOFF2[I] = BYTOFF[ICOL[I]+1,ILUN]-1

            LOOP_END_DIMS:

        ENDFOR

;
;  Construct any virtual columns first
;
        WC = WHERE(FOUND EQ 1 AND VIRTUAL EQ 1, WCCOUNT)
        IF WCCOUNT GT 0 THEN BEGIN
            FOR I = 0, WCCOUNT-1 DO BEGIN
                ;; If it's virtual, then the value only needs to be
                ;; replicated
                EXTCMD = 'D'+STRTRIM(WC[I],2)+$
                  ' = REPLICATE(D'+STRTRIM(WC[I],2)+',NROWS)'
                ;; Run the command that selects the data
                RESULT = EXECUTE(EXTCMD)
                IF RESULT EQ 0 THEN BEGIN
                    MESSAGE = 'ERROR: Could not extract data (column '+$
                      STRTRIM(MYCOL[WC[I]],2)+')'
                    IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                        ERRMSG = MESSAGE
                        RETURN
                    ENDIF ELSE MESSAGE, MESSAGE
                ENDIF
                OUTSTATUS[I] = 1
            ENDFOR
        ENDIF


;  Skip to the end if all columns are virtual
        IF WCCOUNT EQ NUMCOLS THEN GOTO, PROC_CLEANUP

        IF N_ELEMENTS(NANVALUE) GE NUMCOLS THEN BEGIN
            NANVALUES = NANVALUE[0:NUMCOLS-1]
        ENDIF ELSE IF N_ELEMENTS(NANVALUE) GT 0 THEN BEGIN
            NANVALUES = REPLICATE(NANVALUE[0], NUMCOLS)
            NANVALUES[0] = NANVALUE
            I = N_ELEMENTS(NANVALUE)
            IF I LT NUMCOLS THEN $
              NANVALUES[I:*] = NANVALUE[0]
        ENDIF

;
;  Find the position of the first byte of the data array in the file.
;
	OFFSET0 = NHEADER[ILUN] + NAXIS1[ILUN]*(ROW1-1)
        POS = 0L
        NROWS0 = NROWS
        J = 0L
        FIRST = 1
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
        POINT_LUN, UNIT, OFFSET0+OFFSET1
        READU, UNIT, BB

;
;  Now select out the desired columns
;
        FOR I = 0, NUMCOLS-1 DO BEGIN
           
            ;; Extract the proper rows and columns
            IF NOT FOUND[I] THEN GOTO, LOOP_END_STORE
            IF VIRTUAL[I] THEN GOTO, LOOP_END_STORE

            ;; Extract the data from the byte array and convert it
            ;; The inner CALL_FUNCTION is to one of the coercion
            ;; functions, such as FIX(), DOUBLE(), STRING(), etc.,
            ;; which is called with an offset to force a conversion
            ;; from bytes to the data type.
            ;; The outer CALL_FUNCTION is to REFORM(), which makes
            ;; sure that the data structure is correct.
            ;;
            DIMS = COLDIM[I,0:COLNDIM[I]-1]
            PERROW = ROUND(PRODUCT(DIMS)/NROWS0)

            IF COLTYPE[I] EQ 7 THEN BEGIN
                DD = STRING(BB[BOFF1[I]:BOFF2[I], *])
            ENDIF ELSE BEGIN
                DD = CALL_FUNCTION(DTYPENAMES[COLTYPE[I]], $
                                   BB[BOFF1[I]:BOFF2[I],*],0, PERROW*NR)
            ENDELSE
            IF N_ELEMENTS(DD) EQ 1 THEN DD = [DD]
            DD = REFORM(DD, PERROW, NR, /OVERWRITE)

            ;; Now perform any type-specific conversions, etc.
            COUNT = 0L
            CT = COLTYPE[I]
            CASE 1 OF
                ;; Integer types
                (CT EQ 2 OR CT EQ 3): BEGIN
                    IF NOT KEYWORD_SET(NOIEEE) THEN $
                      IEEE_TO_HOST, DD 
                END

                ;; Floating and complex types
                (CT GE 4 OR CT LE 6 OR CT EQ 9): BEGIN
                    IF NOT KEYWORD_SET(NOIEEE) THEN BEGIN
                        IF N_ELEMENTS(NANVALUES) GT 0 THEN W=WHERENAN(DD,COUNT)
                        IEEE_TO_HOST, DD
                    ENDIF
                END

                ;; String types (CT EQ 7) have already been converted
                ;; in the above CALL_FUNCTION.  No further conversion
                ;; is necessary here.
            ENDCASE

;
;  If the parameters TZERO and TSCAL are non-trivial, then adjust the array by
;  these values.
;
            IF NOT KEYWORD_SET(NOIEEE) AND NOT KEYWORD_SET(NOSCALE) THEN BEGIN
                BZERO  = TZERO[ICOL[I],ILUN]
                BSCALE = TSCAL[ICOL[I],ILUN]
                IF (BSCALE NE 0) AND (BSCALE NE 1) THEN DD = BSCALE*DD
		IF BZERO NE 0 THEN DD = DD + BZERO
            ENDIF

;
;  Store NANVALUE everywhere where the data corresponded to IEEE NaN.
;
            IF COUNT GT 0 THEN DD[W] = NANVALUES[I]

            ;; Initialize the output variable on the first chunk
            IF FIRST THEN BEGIN
                RESULT = EXECUTE('D'+STRTRIM(I,2)+' = 0')
                DA = MAKE_ARRAY(PERROW, NROWS0, TYPE=COLTYPE[I])
                RESULT = EXECUTE('D'+STRTRIM(I,2)+' = '+$
                                 'REFORM(DA,PERROW, NROWS0,/OVERWRITE)')
            ENDIF

            ;; Finally, store this in the output variable
            RESULT = EXECUTE('D'+STRTRIM(I,2)+'(0,POS) = DD')
            DD = 0
            IF RESULT EQ 0 THEN BEGIN
                MESSAGE = 'ERROR: Could not compose output data D'+$
                  STRTRIM(I,2)
                IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                    ERRMSG = MESSAGE
                    RETURN
                ENDIF ELSE MESSAGE, MESSAGE
            ENDIF

            OUTSTATUS[I] = 1

            LOOP_END_STORE:
        ENDFOR

        FIRST = 0
        NROWS = NROWS - NR
        POS   = POS + NR
        ENDWHILE

        FOR I = 0, NUMCOLS-1 DO BEGIN
            IF OUTSTATUS[I] NE 1 THEN GOTO, LOOP_END_FINAL
            DIMS = COLDIM[I,0:COLNDIM[I]-1]
            NEL  = PRODUCT(DIMS)
            IF NEL GT 1 THEN $
              RESULT = EXECUTE('D'+STRTRIM(I,2)+' = '+$
                               'REFORM(D'+STRTRIM(I,2)+',DIMS,/OVERWRITE)') $
            ELSE $
              RESULT = EXECUTE('D'+STRTRIM(I,2)+' = D'+STRTRIM(I,2)+'(0)')
            LOOP_END_FINAL:
        ENDFOR

        PROC_CLEANUP:
;
        IF N_ELEMENTS(ERRMSG) NE 0 THEN ERRMSG = ''
	RETURN
	END
