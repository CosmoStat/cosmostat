
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetLinkBase

FUNCTION PIOGETLINKBASE,DataBaseName,GroupName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    DataBaseName_TMP=BYTARR(128)
    DataBaseName_TMP(*)=0
    if (N_ELEMENTS(DataBaseName) GT 0) THEN if (STRLEN(DataBaseName) GT 0) THEN DataBaseName_TMP(0:STRLEN(DataBaseName)-1)=BYTE(DataBaseName)
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)

 PIOGETLINKBASE=call_external(PIOLibIDLSO,'piogetlinkbase_tempoidl', $
        DataBaseName_TMP, $
        GroupName_TMP, $
               /L64_VALUE) 
  DataBaseName=STRING(DataBaseName_TMP)
  GroupName=STRING(GroupName_TMP)

 RETURN,PIOGETLINKBASE
END
