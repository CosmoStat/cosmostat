
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetGrpList

FUNCTION PIOGETGRPLIST,GrpList,GrpType,Database
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    GrpList=long64(0)
    GrpType=long64(0)
    Database_TMP=BYTARR(128)
    Database_TMP(*)=0
    if (N_ELEMENTS(Database) GT 0) THEN if (STRLEN(Database) GT 0) THEN Database_TMP(0:STRLEN(Database)-1)=BYTE(Database)

 PIOGETGRPLIST=call_external(PIOLibIDLSO,'piogetgrplist_tempoidl', $
        GrpList, $
        GrpType, $
        Database_TMP, $
               /L64_VALUE) 
  IF (PIOGETGRPLIST GT 0) THEN BEGIN
     TP=GrpList
     TMPNAME=bytarr(128)
     GrpList=STRARR(PIOGETGRPLIST)
        FOR I=0L,PIOGETGRPLIST-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            GrpList(I)=strmid(string(TMPNAME),0,A)
            END
    END
  IF (PIOGETGRPLIST GT 0) THEN BEGIN
     TP=GrpType
     TMPNAME=bytarr(128)
     GrpType=STRARR(PIOGETGRPLIST)
        FOR I=0L,PIOGETGRPLIST-1L DO BEGIN
            A=call_external(PIOLibIDLSO,'piogetstringidl',TP+I*128,TMPNAME)
            GrpType(I)=strmid(string(TMPNAME),0,A)
            END
    END

 RETURN,PIOGETGRPLIST
END
