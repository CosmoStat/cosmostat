; $Id: //depot/idl/IDL_70/idl_src/libs/imsl/imsl_6.0/lib/imsl_sp_readhb.pro#1 $
;
; Copyright (c) 1970-2006, VISUAL NUMERICS Inc. All Rights Reserved.
; This software is confidential information which is proprietary to and a
; trade secret of Visual Numerics Inc.  Use, duplication or disclosure is
; subject to the terms of an appropriate license agreement.
;
;
; Read in a Harwell-Boeing data file containing a sparse matrix.
; File formats are from User's Guide for the Harwell-Boeing
; Sparse Matrix Colloection (Release I)
;
PRO imsl_Sp_readhb, filename, header, pointr, rowind, $
               Values=values, $
               Rhsval=rhsval, $
               Rhsptr=rhsptr, $
               Rhsind=rhsind, $
               Sguess=sguess, $
               Xexact=xexact, $
               Header_only=header_only
;+
; Name:
;       SP_READHB
; Purpose:
;       Reads a Harwell-Boeing data file containing a sparse matrix.
; Category:
;       Sparse Matrices
; Calling Sequence:
;       Sp_readhb, filename, header [, pointr, rowind]
; Inputs:
;       Filename = Name of the file to be read.
; Outputs:
;       Header = Structure containing the information from the
;                file header.
;       Pointr = One dimensional array containing the column start
;                pointers.
;       Rowind = One dimensional array containing the row indices.
; Keyword Parameters:
;       Header_only = If present and nonzero, then return after reading
;                only the header. (Input)
;       Values = One dimensional array containing the numerical
;                values, if present in the file. (Output)
;       Rhsval = One dimensional array containing the RHS
;                values, if present in the file. (Output)
;       Rhsptr = One dimensional array containing the RHS pointers.
;                This keyword is only valid if the RHS(s) are supplied
;                in sparse or elemental form. (Output)
;       Rhsind = One dimensional array containing the RHS indices.
;                This keyword is only valid if the RHS(s) are supplied
;                in sparse or elemental form. (Output)
;       Sguess = One dimensional array containing the Starting guesses.
;                (Output)
;       Xexact = One dimensional array containing the exact solutions.
;                (Output)
; Side effects:
;       None known
; Restrictions:
;       None known
; Procedure:
;       File formats are from User's Guide for the Harwell-Boeing
;       Sparse Matrix Collection (Release I)
; Reference:
;        Duff, Iian S., Grimes, Roger G., Lewis, John G.,
;        User's Guide for the Harwell-Boeing Sparse Matrix
;        Collection (Release I), October 1992.
;
;
;
;-
@imsl_init.pro
  ON_ERROR, on_err_action

   ; Declare some variables that will be used when reading in the
   ; Harwell Boeing data.
   title  = STRING(1)
   key    = STRING(1)
   totcrd = IMSL_0
   ptrcrd = IMSL_0
   indcrd = IMSL_0
   valcrd = IMSL_0
   rhscrd = IMSL_0
   mxtype = STRING(1)
   nrow   = IMSL_0
   ncol   = IMSL_0
   nnzero = IMSL_0
   neltvl = IMSL_0
   ptrfmt = STRING(1)
   indfmt = STRING(1)
   valfmt = STRING(1)
   rhsfmt = STRING(1)
   rhstyp = STRING(1)
   nrhs   = IMSL_0
   nrhsix = IMSL_0
   ;
   OPENR, unit, filename, /GET_LUN
   ;
   ; Read in the first four lines of the header all at once.
   header_format = '(a72, a8, /, 5I14, /, a3, 11x, 4I14, /, 2a16, 2a20)'
   READF, unit, title, key,  $                   ;LINE 1
     totcrd, ptrcrd, indcrd, valcrd, rhscrd,  $  ;LINE 2
     mxtype, nrow, ncol, nnzero, neltvl, $       ;Line 3
     ptrfmt, indfmt, valfmt, rhsfmt, $           ;LINE 4
     Format = header_format
   ;
   ; If rhscrd is positive, then there is at least one RHS, so we must
   ; read in a fifth line in the header.
   IF (rhscrd GT 0) THEN READF, unit, rhstyp, nrhs, nrhsix, $
     Format = '(A3, 11x, 2I14)'                  ;LINE 5
   ;
   ; Define the structure that contains the header information.
   header = { title  : title,  $
              key    : key,    $
              totcrd : totcrd, $
              ptrcrd : ptrcrd, $
              indcrd : indcrd, $
              valcrd : valcrd, $
              rhscrd : rhscrd, $
              mxtype : mxtype, $
              nrow   : nrow,   $
              ncol   : ncol,   $
              nnzero : nnzero, $
              neltvl : neltvl, $
              ptrfmt : ptrfmt, $
              indfmt : indfmt, $
              valfmt : valfmt, $
              rhsfmt : rhsfmt, $
              rhstyp : rhstyp, $
              nrhs   : nrhs,   $
              hrhsix : nrhsix}
   ;
   ; If we are only supposed to read the header, return.
   IF (KEYWORD_SET(header_only)) THEN RETURN
   ;
   ; Start by reading in the matrix structure.
   ;
   ; 1st data block: (NROW+1) Column Pointers
   IF (ARG_PRESENT(pointr)) THEN BEGIN
      pointr = IMSL_LONARR(nrow+1)
      READF, unit, pointr, Format = ptrfmt
   END
   ;
   ; 2nd data block: NNZERO Row Indices
   IF (ARG_PRESENT(pointr)) THEN BEGIN
      rowind = IMSL_LONARR(nnzero)
      READF, unit, rowind, FORMAT = indfmt
   END
   ;
   ; 3rd data block (optional): Numerical values.
   IF (ARG_PRESENT(values)) THEN BEGIN
      IF (strmid(mxtype, 0, 1) NE 'P') THEN BEGIN ;Any numerical values?
         CASE strmid(mxtype, 2, 1) OF
            ;
            ;Assembled
            'A': BEGIN
               values = DBLARR(nnzero)
               READF, unit, values, Format = valfmt
            END
            ;
            ;Elemental
            'E': BEGIN
               values = DBLARR(neltvl)
               READF, unit, values, Format = valfmt

            END
            ;
            ; Incorrect matrix type.
            ELSE: MESSAGE, 'Incorrect third character in matrix type.'
         END
      END
   END
   ;
   ; Start on RHS's if present.
   If (nrhs GT 0) THEN BEGIN
      IF ((ARG_PRESENT(rhsptr)) OR (ARG_PRESENT(rhsptr)) OR $
          (ARG_PRESENT(rhsptr))) THEN BEGIN
         IF (strmid(rhstyp, 0, 1) EQ 'F') THEN BEGIN
            ; Read dense RHS's.
            rhsval = DBLARR(nrow*nrhs) ;Store all RHS's in a 1-D array.
            READF, unit, rhsval, Format = rhsfmt
         END ELSE BEGIN
            ; Read sparse of elemental RHS's.
            IF (strmid(mxtype, 2, 1) EQ 'A') THEN BEGIN
               ; Sparse RHS's
               ;
               ; Pointer array
               rhsptr = IMSL_LONARR(nrhs+1)
               READF, unit, rhsptr, Format = ptrfmt
               ;
               ; Sparsity pattern
               rhsind = IMSL_LONARR(nrhsix)
               READF, unit, rhsind, Format = indfmt
               ;
               ; RHS values
               rhsval = dblarr(nrhsix)
               READF, unit, rhsval, Format = rshfmt
            END ELSE BEGIN
               ; Elemental RSH's
               rhsval = dblarr(nnzero*nrhs) ;Store all RHS's in a 1-D array.
               readf, unit, rhsval, Format = rhsfmt
            END
         END
      END
      ;
      ; Read in guesses if they are present.
      IF (ARG_PRESENT(sguess) AND (strmid(rhstyp, 1, 1) EQ 'G')) THEN BEGIN
         sguess = dblarr(nrow*nrhs) ;Store all guesses in a 1-D array.
         READF, unit, sguess, Format = rhsfmt
      END
      ;
      ; Read in solutions if they are present.
      IF (ARG_PRESENT(xexact) AND (strmid(rhstyp, 2, 2) EQ 'X')) THEN BEGIN
         xexact = dblarr(nrow*nrhs) ;Store all solutions in a 1-D array.
         READF, unit, xexact, Format = rhsfmt
      END
   END
END
