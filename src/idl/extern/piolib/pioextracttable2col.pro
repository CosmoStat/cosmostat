
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOExtractTable2col

FUNCTION PIOEXTRACTTABLE2COL,X_Value,ERR_X_Value,Y_Value,ERR_Y_Value,Type,Index,MyTable2col
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 CASE (STRING(TYPE)) OF
"PIOLONG": X_Value = LONG64(X_Value)
"PIOFLAG": X_Value = BYTE(X_Value)
"PIOBYTE": X_Value = BYTE(X_Value)
"PIOINT": X_Value = LONG(X_Value)
"PIOSTRING":BEGIN
    X_Value=PIOSTRING2C(X_Value)
     NBDATA=N_ELEMENTS(X_Value)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(X_Value(I))-1)=BYTE(X_Value(I))
   X_Value=RES
END 
"PIOSHORT": X_Value = FIX(X_Value)
 "PIOFLOAT": X_Value = FLOAT(X_Value)
 "PIODOUBLE": X_Value = DOUBLE(X_Value)
 "PIOCOMPLEX": X_Value = COMPLEX(X_Value)
 "PIODBLCOMPLEX": X_Value = DBLCOMPLEX(X_Value)
ELSE: RETURN,-1006 ; wrong type definition
END
 CASE (STRING(TYPE)) OF
"PIOLONG": ERR_X_Value = LONG64(ERR_X_Value)
"PIOFLAG": ERR_X_Value = BYTE(ERR_X_Value)
"PIOBYTE": ERR_X_Value = BYTE(ERR_X_Value)
"PIOINT": ERR_X_Value = LONG(ERR_X_Value)
"PIOSTRING":BEGIN
    ERR_X_Value=PIOSTRING2C(ERR_X_Value)
     NBDATA=N_ELEMENTS(ERR_X_Value)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(ERR_X_Value(I))-1)=BYTE(ERR_X_Value(I))
   ERR_X_Value=RES
END 
"PIOSHORT": ERR_X_Value = FIX(ERR_X_Value)
 "PIOFLOAT": ERR_X_Value = FLOAT(ERR_X_Value)
 "PIODOUBLE": ERR_X_Value = DOUBLE(ERR_X_Value)
 "PIOCOMPLEX": ERR_X_Value = COMPLEX(ERR_X_Value)
 "PIODBLCOMPLEX": ERR_X_Value = DBLCOMPLEX(ERR_X_Value)
ELSE: RETURN,-1006 ; wrong type definition
END
 CASE (STRING(TYPE)) OF
"PIOLONG": Y_Value = LONG64(Y_Value)
"PIOFLAG": Y_Value = BYTE(Y_Value)
"PIOBYTE": Y_Value = BYTE(Y_Value)
"PIOINT": Y_Value = LONG(Y_Value)
"PIOSTRING":BEGIN
    Y_Value=PIOSTRING2C(Y_Value)
     NBDATA=N_ELEMENTS(Y_Value)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(Y_Value(I))-1)=BYTE(Y_Value(I))
   Y_Value=RES
END 
"PIOSHORT": Y_Value = FIX(Y_Value)
 "PIOFLOAT": Y_Value = FLOAT(Y_Value)
 "PIODOUBLE": Y_Value = DOUBLE(Y_Value)
 "PIOCOMPLEX": Y_Value = COMPLEX(Y_Value)
 "PIODBLCOMPLEX": Y_Value = DBLCOMPLEX(Y_Value)
ELSE: RETURN,-1006 ; wrong type definition
END
 CASE (STRING(TYPE)) OF
"PIOLONG": ERR_Y_Value = LONG64(ERR_Y_Value)
"PIOFLAG": ERR_Y_Value = BYTE(ERR_Y_Value)
"PIOBYTE": ERR_Y_Value = BYTE(ERR_Y_Value)
"PIOINT": ERR_Y_Value = LONG(ERR_Y_Value)
"PIOSTRING":BEGIN
    ERR_Y_Value=PIOSTRING2C(ERR_Y_Value)
     NBDATA=N_ELEMENTS(ERR_Y_Value)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(ERR_Y_Value(I))-1)=BYTE(ERR_Y_Value(I))
   ERR_Y_Value=RES
END 
"PIOSHORT": ERR_Y_Value = FIX(ERR_Y_Value)
 "PIOFLOAT": ERR_Y_Value = FLOAT(ERR_Y_Value)
 "PIODOUBLE": ERR_Y_Value = DOUBLE(ERR_Y_Value)
 "PIOCOMPLEX": ERR_Y_Value = COMPLEX(ERR_Y_Value)
 "PIODBLCOMPLEX": ERR_Y_Value = DBLCOMPLEX(ERR_Y_Value)
ELSE: RETURN,-1006 ; wrong type definition
END
    Type_TMP=BYTARR(128)
    Type_TMP(*)=0
    if (N_ELEMENTS(Type) GT 0) THEN if (STRLEN(Type) GT 0) THEN Type_TMP(0:STRLEN(Type)-1)=BYTE(Type)
    if (N_ELEMENTS(Index) EQ 0) THEN     Index=LONG64(0) ELSE    Index=LONG64(Index)
  IF (N_ELEMENTS(MyTable2col) EQ 0) THEN MyTable2col=LONG64(0) ELSE MyTable2col=LONG64(MyTable2col)

 PIOEXTRACTTABLE2COL=call_external(PIOLibIDLSO,'pioextracttable2col_tempoidl', $
        X_Value, $
        ERR_X_Value, $
        Y_Value, $
        ERR_Y_Value, $
        Type_TMP, $
        Index, $
        MyTable2col, $
               /L64_VALUE) 

 RETURN,PIOEXTRACTTABLE2COL
END
