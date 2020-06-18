
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOWriteKeywordObject

FUNCTION PIOWRITEKEYWORDOBJECT,Value,Comment,InKeyword,Type,Object,MyGroup,INDEX=INDEX,MASK=MASK
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    IF (KEYWORD_SET(INDEX) AND KEYWORD_SET(MASK)) THEN BEGIN
        RETURN,PIOWRITEKeywordOBJECTIDXMASK(DATA,INDEX,MASK,N_ELEMENTS(INDEX),OBJECTNAME,TYPE,COMMAND,MYGROUP)
    END
    IF (KEYWORD_SET(MASK)) THEN BEGIN 
        RETURN,PIOWRITEKeywordOBJECTMASK(DATA,MASK,OBJECTNAME,TYPE,COMMAND,MYGROUP)
    END
    IF (KEYWORD_SET(INDEX)) THEN BEGIN 
        RETURN,PIOWRITEKeywordOBJECTIDX(DATA,INDEX,N_ELEMENTS(INDEX),OBJECTNAME,TYPE,COMMAND,MYGROUP)
    END
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
    Type_TMP=BYTARR(128)
    Type_TMP(*)=0
    if (N_ELEMENTS(Type) GT 0) THEN if (STRLEN(Type) GT 0) THEN Type_TMP(0:STRLEN(Type)-1)=BYTE(Type)
    Object_TMP=BYTARR(128)
    Object_TMP(*)=0
    if (N_ELEMENTS(Object) GT 0) THEN if (STRLEN(Object) GT 0) THEN Object_TMP(0:STRLEN(Object)-1)=BYTE(Object)
IF (N_PARAMS() LT 6) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOWRITEKEYWORDOBJECT=call_external(PIOLibIDLSO,'piowritekeywordobject_tempoidl', $
        Value, $
        Comment_TMP, $
        InKeyword_TMP, $
        Type_TMP, $
        Object_TMP, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOWRITEKEYWORDOBJECT
END
