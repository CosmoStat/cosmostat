
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetObjectList

FUNCTION PIOGETOBJECTLIST,ObjList,ObjType,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    ObjList=long64(0)
    ObjType=long64(0)

 PIOGETOBJECTLIST=call_external(PIOLibIDLSO,'piogetobjectlist_tempoidl', $
        ObjList, $
        ObjType, $
        MyGroup, $
               /L64_VALUE) 
  IF (PIOGETOBJECTLIST GT 0) THEN BEGIN
     TP=ObjList
     TMPNAME=bytarr(128)
     ObjList=STRARR(PIOGETOBJECTLIST)
        FOR I=0L,PIOGETOBJECTLIST-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            ObjList(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (PIOGETOBJECTLIST GT 0) THEN BEGIN
     TP=ObjType
     TMPNAME=bytarr(128)
     ObjType=STRARR(PIOGETOBJECTLIST)
        FOR I=0L,PIOGETOBJECTLIST-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            ObjType(I)=strmid(string(TMPNAME),0,A)
            END
    END

 RETURN,PIOGETOBJECTLIST
END
