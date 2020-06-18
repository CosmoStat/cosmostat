
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOWriteKeywordGrp

FUNCTION PIOWRITEKEYWORDGRP,Value,Comment,InKeyword,type,MyGroup
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
    Comment_TMP=BYTARR(128)
    Comment_TMP(*)=0
    if (N_ELEMENTS(Comment) GT 0) THEN if (STRLEN(Comment) GT 0) THEN Comment_TMP(0:STRLEN(Comment)-1)=BYTE(Comment)
    InKeyword_TMP=BYTARR(128)
    InKeyword_TMP(*)=0
    if (N_ELEMENTS(InKeyword) GT 0) THEN if (STRLEN(InKeyword) GT 0) THEN InKeyword_TMP(0:STRLEN(InKeyword)-1)=BYTE(InKeyword)
    type_TMP=BYTARR(128)
    type_TMP(*)=0
    if (N_ELEMENTS(type) GT 0) THEN if (STRLEN(type) GT 0) THEN type_TMP(0:STRLEN(type)-1)=BYTE(type)

 PIOWRITEKEYWORDGRP=call_external(PIOLibIDLSO,'piowritekeywordgrp_tempoidl', $
        Value, $
        Comment_TMP, $
        InKeyword_TMP, $
        type_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOWRITEKEYWORDGRP
END
