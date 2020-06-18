
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOKeywordListObject

FUNCTION PIOKEYWORDLISTOBJECT,Keyword,Type,Object,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Keyword=long64(0)
    Type=long64(0)
    Object_TMP=BYTARR(128)
    Object_TMP(*)=0
    if (N_ELEMENTS(Object) GT 0) THEN if (STRLEN(Object) GT 0) THEN Object_TMP(0:STRLEN(Object)-1)=BYTE(Object)
IF (N_PARAMS() LT 4) THEN BEGIN
    MyGroup=long64(0)
END ELSE BEGIN
    MyGroup=long64(MyGroup)
END

 PIOKEYWORDLISTOBJECT=call_external(PIOLibIDLSO,'piokeywordlistobject_tempoidl', $
        Keyword, $
        Type, $
        Object_TMP, $
        MyGroup, $
               /L64_VALUE) 
  IF (PIOKEYWORDLISTOBJECT GT 0) THEN BEGIN
     TP=Keyword
     TMPNAME=bytarr(128)
     Keyword=STRARR(PIOKEYWORDLISTOBJECT)
        FOR I=0L,PIOKEYWORDLISTOBJECT-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            Keyword(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (PIOKEYWORDLISTOBJECT GT 0) THEN BEGIN
     TP=Type
     TMPNAME=bytarr(128)
     Type=STRARR(PIOKEYWORDLISTOBJECT)
        FOR I=0L,PIOKEYWORDLISTOBJECT-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            Type(I)=strmid(string(TMPNAME),0,A)
            END
    END

 RETURN,PIOKEYWORDLISTOBJECT
END
