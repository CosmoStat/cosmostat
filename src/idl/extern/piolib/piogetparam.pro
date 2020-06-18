
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetParam

FUNCTION PIOGETPARAM,Value,Type,ParamName,MyPar
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
 Value=0
CASE (STRING(TYPE)) OF
"PIOLONG": Value = LONG64(Value)
"PIOFLAG": Value = BYTE(Value)
"PIOBYTE": Value = BYTE(Value)
"PIOINT": Value = LONG(Value)
"PIOSTRING":BEGIN
    RES = BYTARR(128L)
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
    ParamName_TMP=BYTARR(128)
    ParamName_TMP(*)=0
    if (N_ELEMENTS(ParamName) GT 0) THEN if (STRLEN(ParamName) GT 0) THEN ParamName_TMP(0:STRLEN(ParamName)-1)=BYTE(ParamName)
 IF (N_ELEMENTS(MyPar) EQ 0) THEN MyPar=LONG64(0) ELSE MyPar=LONG64(MyPar)

 PIOGETPARAM=call_external(PIOLibIDLSO,'piogetparam_tempoidl', $
        Value, $
        Type_TMP, $
        ParamName_TMP, $
        MyPar, $
               /L64_VALUE) 
 IF (STRING(TYPE) EQ "PIOSTRING") THEN Value=STRING(Value)

 RETURN,PIOGETPARAM
END
