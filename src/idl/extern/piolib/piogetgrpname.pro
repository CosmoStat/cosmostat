
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetGrpName

FUNCTION PIOGETGRPNAME,GroupName,ObjectName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)
    ObjectName_TMP=BYTARR(128)
    ObjectName_TMP(*)=0
    if (N_ELEMENTS(ObjectName) GT 0) THEN if (STRLEN(ObjectName) GT 0) THEN ObjectName_TMP(0:STRLEN(ObjectName)-1)=BYTE(ObjectName)

 PIOGETGRPNAME=call_external(PIOLibIDLSO,'piogetgrpname_tempoidl', $
        GroupName_TMP, $
        ObjectName_TMP, $
               /L64_VALUE) 
  GroupName=STRING(GroupName_TMP)
  ObjectName=STRING(ObjectName_TMP)

 RETURN,PIOGETGRPNAME
END
