
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetList

FUNCTION PIOGETLIST,Value,Type,ListName,Index,MyPar
     ;ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
Value = 0
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
   Value=RES
END 
"PIOSHORT": Value = FIX(Value)
 "PIOFLOAT": Value = FLOAT(Value)
 "PIODOUBLE": Value = DOUBLE(Value)
 "PIOCOMPLEX": Value = COMPLEX(Value)
 "PIODBLCOMPLEX": Value = DBLCOMPLEX(Value)
ELSE: RETURN,-1006 ; wrong type definition
END
    Type_TMP=BYTARR(128)
    Type_TMP(*)=0
    if (N_ELEMENTS(Type) GT 0) THEN if (STRLEN(Type) GT 0) THEN Type_TMP(0:STRLEN(Type)-1)=BYTE(Type)
    ListName_TMP=BYTARR(128)
    ListName_TMP(*)=0
    if (N_ELEMENTS(ListName) GT 0) THEN if (STRLEN(ListName) GT 0) THEN ListName_TMP(0:STRLEN(ListName)-1)=BYTE(ListName)
    if (N_ELEMENTS(Index) EQ 0) THEN     Index=LONG64(0) ELSE    Index=LONG64(Index)
  IF (N_ELEMENTS(MyPar) EQ 0) THEN MyPar=LONG64(0) ELSE MyPar=LONG64(MyPar)

 PIOGETLIST=call_external(PIOLibIDLSO,'piogetlist_tempoidl', $
        Value, $
        Type_TMP, $
        ListName_TMP, $
        Index, $
        MyPar, $
               /L64_VALUE) 

 if (STRING(TYPE) EQ "PIOSTRING") then Value=string(Value)
 RETURN,PIOGETLIST
END
