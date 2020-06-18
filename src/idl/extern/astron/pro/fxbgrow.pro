	PRO FXBGROW, UNIT, HEADER, NROWS, ERRMSG=ERRMSG, NOZERO=NOZERO
;+
; NAME: 
;	FXBGROW
; Purpose     : 
;	Increase the number of rows in a binary table.
; Explanation : 
;       Call FXBGROW to increase the size of an already-existing FITS
;       binary table.  The number of rows increases to NROWS (or does
;       not change if NROWS is less than the number of rows already
;       existing).  WARNING:  the table to be grown must be the *last*
;       extension in the FITS file.  FXBGROW does *not* preserve any
;       following extensions.  This procedure is useful when a table
;       with an unknown number of rows must be created.  The caller
;       would then call FXBCREATE to construct a table of some base
;       size, and follow with calls to FXBGROW to lengthen the table
;       as needed.
;
; Use         : 
;	FXBGROW, UNIT, HEADER, NROWS[, ERRMSG=ERRMSG, NOZERO=NOZERO]
; Inputs      : 
;	UNIT     = Logical unit number of an already-opened file.
;	HEADER	 = String array containing the FITS binary table extension
;		   header.  The header is modified in place.
;       NROWS    = New number of rows, always more than the previous
;                  number.
; Opt. Inputs : 
;	None.
; Outputs     : 
;	None.
; Opt. Outputs: 
;	None.
; Keywords    : 
;       NOZERO   = when set, FXBGROW will not zero-pad the new data if
;                  it doesn't have to.
;	ERRMSG	  = If defined and passed, then any error messages will be
;		    returned to the user in this parameter rather than
;		    depending on the MESSAGE routine in IDL.  If no errors are
;		    encountered, then a null string is returned.  In order to
;		    use this feature, ERRMSG must be defined first, e.g.
;
;			ERRMSG = ''
;			FXBGROW, ERRMSG=ERRMSG, ...
;			IF ERRMSG NE '' THEN ...
;
; Calls       : 
;	FXADDPAR, FXHREAD
; Common      : 
;	Uses common block FXBINTABLE--see "fxbintable.pro" for more
;	information.
; Restrictions: 
;       The file must be open with write permission.
;
;       The binary table extension in question must already by written
;       to the file (using FXBCREATE), and must be the last extension
;       in the file.
;
;       A table can never shrink via this operation.
;
;       This operation is not well optimized for tables with large
;       heap usage, such as large variable-length columns.  Since the
;       procedure must move the entire heap upon every call, it could
;       be (1) memory intensive and (2) I/O intensive.
;
; Side effects: 
;	The FITS file will grow in size, and heap areas are
;	preserved by moving them to the end of the file.
;
;       The header is modified to reflect the new number of rows.
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;       Initially written, C. Markwardt, GSFC, Nov 1998
; Written     : 
;	Craig Markwardt, GSFC, Nov 1998
; Version     :
;       Version 1, 17 Nov 1998
;-
;
@fxbintable
	ON_ERROR, 0
;
;  Check the number of parameters.
;
	IF N_PARAMS() NE 3 THEN BEGIN
		MESSAGE = 'Syntax:  FXBGROW, UNIT, HEADER, NROWS'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF

;
;  Find the index of the file.
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
;  Don't shrink the file.
;
        IF NAXIS2[ILUN] GE NROWS THEN GOTO, FINISH
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
;  Ugghhh.  Have to read the HEAP and DHEAP portions, especially to
;  preserve variable length columns.  Seek to the end of the data...
;
        DBYTES = NAXIS1[ILUN]*NAXIS2[ILUN]
        HBYTES = HEAP[ILUN] + DHEAP[ILUN] - DBYTES
        IF HBYTES GT 0 THEN BEGIN
            POINT_LUN, UNIT, NHEADER[ILUN] + DBYTES
            HBUFFER = BYTARR(HBYTES)
            READU, UNIT, HBUFFER
        ENDIF
;
;  Grow the array by padding with zeroes.
;
        IF HBYTES GT 0 OR NOT KEYWORD_SET(NOZERO) THEN BEGIN
            POINT_LUN, UNIT, NHEADER[ILUN] + DBYTES
            BUFFER = BYTARR(NAXIS1[ILUN])
            FOR I = NAXIS2[ILUN]+1,NROWS  DO WRITEU,UNIT,BUFFER
        ENDIF
;
;  Write the heap/whatever data back to disk, in its new position
        IF HBYTES GT 0 THEN BEGIN
            POINT_LUN, UNIT, NHEADER[ILUN] + NAXIS1[ILUN]*NROWS
            WRITEU, UNIT, HBUFFER
            HBUFFER = 0
        ENDIF
;
;  Update the internal state.
;
        HEAP[ILUN] = HEAP[ILUN] + (NROWS-NAXIS2[ILUN])*NAXIS1[ILUN]
        NAXIS2[ILUN] = NROWS

;  Modify the header in place and on disk.
        IF N_ELEMENTS(HEADER) GT 0 THEN $
          FXADDPAR, HEADER, 'NAXIS2', LONG(NROWS), 'Number of rows (grown)'
        POINT_LUN, UNIT, MHEADER[ILUN]
        FXHREAD, UNIT, DHEADER, STATUS
        IF STATUS NE 0 THEN BEGIN
            MESSAGE = 'Could not load header from file'
            IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
                ERRMSG = MESSAGE
                RETURN
            END ELSE MESSAGE, MESSAGE
        ENDIF
        FXADDPAR, DHEADER, 'NAXIS2', LONG(NROWS), 'Number of rows (grown)'
        ;; Don't worry about the header increasing in size, since
        ;; every binary table has to have NAXIS2 already.
        SLEN = STRLEN(DHEADER[0])
        FULL = STRING(REPLICATE(32B, 80))
        IF SLEN LT 80 THEN DHEADER[0] = DHEADER[0] + STRMID(FULL,0,80-SLEN)
        BHDR = BYTE(DHEADER)
        BHDR = BHDR[0:79,*]
        POINT_LUN, UNIT, MHEADER[ILUN]
        WRITEU, UNIT, BHDR

FINISH:
	IF N_ELEMENTS(ERRMSG) NE 0 THEN ERRMSG = ''
	RETURN
	END
