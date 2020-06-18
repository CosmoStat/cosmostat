
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetGroupType

FUNCTION PIOGETGROUPTYPE,Typein,GroupName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Typein_TMP=BYTARR(128)
    Typein_TMP(*)=0
    if (N_ELEMENTS(Typein) GT 0) THEN if (STRLEN(Typein) GT 0) THEN Typein_TMP(0:STRLEN(Typein)-1)=BYTE(Typein)
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)

 PIOGETGROUPTYPE=call_external(PIOLibIDLSO,'piogetgrouptype_tempoidl', $
        Typein_TMP, $
        GroupName_TMP, $
               /L64_VALUE) 
  Typein=STRING(Typein_TMP)
  GroupName=STRING(GroupName_TMP)

 RETURN,PIOGETGROUPTYPE
END
