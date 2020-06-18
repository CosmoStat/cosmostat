
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOExtractTable

FUNCTION PIOEXTRACTTABLE,Value,ERR_Value,Type,Index,MyTable
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 CASE (STRING(TYPE)) OF
"PIOLONG": Value = LONG64(Value)
"PIOFLAG": Value = BYTE(Value)
"PIOBYTE": Value = BYTE(Value)
"PIOINT": Value = LONG(Value)
"PIOSTRING":BEGIN
    Value=PIOSTRING2C(Value)
     NBDATA=N_ELEMENTS(Value)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(Value(I))-1)=BYTE(Value(I))
   Value=RES
END 
"PIOSHORT": Value = FIX(Value)
 "PIOFLOAT": Value = FLOAT(Value)
 "PIODOUBLE": Value = DOUBLE(Value)
 "PIOCOMPLEX": Value = COMPLEX(Value)
 "PIODBLCOMPLEX": Value = DBLCOMPLEX(Value)
ELSE: RETURN,-1006 ; wrong type definition
END
 CASE (STRING(TYPE)) OF
"PIOLONG": ERR_Value = LONG64(ERR_Value)
"PIOFLAG": ERR_Value = BYTE(ERR_Value)
"PIOBYTE": ERR_Value = BYTE(ERR_Value)
"PIOINT": ERR_Value = LONG(ERR_Value)
"PIOSTRING":BEGIN
    ERR_Value=PIOSTRING2C(ERR_Value)
     NBDATA=N_ELEMENTS(ERR_Value)
    RES = BYTARR(NBDATA*128L)
    RES(*)= 32B
    FOR I=0l,NBDATA-1 DO RES(128*I:128*I+STRLEN(ERR_Value(I))-1)=BYTE(ERR_Value(I))
   ERR_Value=RES
END 
"PIOSHORT": ERR_Value = FIX(ERR_Value)
 "PIOFLOAT": ERR_Value = FLOAT(ERR_Value)
 "PIODOUBLE": ERR_Value = DOUBLE(ERR_Value)
 "PIOCOMPLEX": ERR_Value = COMPLEX(ERR_Value)
 "PIODBLCOMPLEX": ERR_Value = DBLCOMPLEX(ERR_Value)
ELSE: RETURN,-1006 ; wrong type definition
END
    Type_TMP=BYTARR(128)
    Type_TMP(*)=0
    if (N_ELEMENTS(Type) GT 0) THEN if (STRLEN(Type) GT 0) THEN Type_TMP(0:STRLEN(Type)-1)=BYTE(Type)
    if (N_ELEMENTS(Index) EQ 0) THEN     Index=LONG64(0) ELSE    Index=LONG64(Index)
  IF (N_ELEMENTS(MyTable) EQ 0) THEN MyTable=LONG64(0) ELSE MyTable=LONG64(MyTable)

 PIOEXTRACTTABLE=call_external(PIOLibIDLSO,'pioextracttable_tempoidl', $
        Value, $
        ERR_Value, $
        Type_TMP, $
        Index, $
        MyTable, $
               /L64_VALUE) 

 RETURN,PIOEXTRACTTABLE
END
