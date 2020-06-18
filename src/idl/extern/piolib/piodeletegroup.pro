
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIODeleteGroup

FUNCTION PIODELETEGROUP,GroupName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)

 PIODELETEGROUP=call_external(PIOLibIDLSO,'piodeletegroup_tempoidl', $
        GroupName_TMP, $
               /L64_VALUE) 

 RETURN,PIODELETEGROUP
END
