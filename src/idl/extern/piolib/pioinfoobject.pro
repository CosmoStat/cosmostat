
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOInfoObject

FUNCTION PIOINFOOBJECT,TOItype,Datatype,BeginIndx,EndIndx,Author,Date,MyObjectName,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    TOItype_TMP=BYTARR(128)
    TOItype_TMP(*)=0
    if (N_ELEMENTS(TOItype) GT 0) THEN if (STRLEN(TOItype) GT 0) THEN TOItype_TMP(0:STRLEN(TOItype)-1)=BYTE(TOItype)
    Datatype_TMP=BYTARR(128)
    Datatype_TMP(*)=0
    if (N_ELEMENTS(Datatype) GT 0) THEN if (STRLEN(Datatype) GT 0) THEN Datatype_TMP(0:STRLEN(Datatype)-1)=BYTE(Datatype)
    if (N_ELEMENTS(BeginIndx) EQ 0) THEN     BeginIndx=LONG64(0) ELSE    BeginIndx=LONG64(BeginIndx)
     if (N_ELEMENTS(EndIndx) EQ 0) THEN     EndIndx=LONG64(0) ELSE    EndIndx=LONG64(EndIndx)
     Author_TMP=BYTARR(128)
    Author_TMP(*)=0
    if (N_ELEMENTS(Author) GT 0) THEN if (STRLEN(Author) GT 0) THEN Author_TMP(0:STRLEN(Author)-1)=BYTE(Author)
    Date_TMP=BYTARR(128)
    Date_TMP(*)=0
    if (N_ELEMENTS(Date) GT 0) THEN if (STRLEN(Date) GT 0) THEN Date_TMP(0:STRLEN(Date)-1)=BYTE(Date)
    MyObjectName_TMP=BYTARR(128)
    MyObjectName_TMP(*)=0
    if (N_ELEMENTS(MyObjectName) GT 0) THEN if (STRLEN(MyObjectName) GT 0) THEN MyObjectName_TMP(0:STRLEN(MyObjectName)-1)=BYTE(MyObjectName)
IF (N_PARAMS() LT 8) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOINFOOBJECT=call_external(PIOLibIDLSO,'pioinfoobject_tempoidl', $
        TOItype_TMP, $
        Datatype_TMP, $
        BeginIndx, $
        EndIndx, $
        Author_TMP, $
        Date_TMP, $
        MyObjectName_TMP, $
        MyGroup, $
               /L64_VALUE) 
  TOItype=STRING(TOItype_TMP)
  Datatype=STRING(Datatype_TMP)
  Author=STRING(Author_TMP)
  Date=STRING(Date_TMP)
  MyObjectName=STRING(MyObjectName_TMP)

 RETURN,PIOINFOOBJECT
END
