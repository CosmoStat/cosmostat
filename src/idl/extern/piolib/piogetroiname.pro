
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetROIName

FUNCTION PIOGETROINAME,RoiName,GroupName
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    RoiName_TMP=BYTARR(128)
    RoiName_TMP(*)=0
    if (N_ELEMENTS(RoiName) GT 0) THEN if (STRLEN(RoiName) GT 0) THEN RoiName_TMP(0:STRLEN(RoiName)-1)=BYTE(RoiName)
    GroupName_TMP=BYTARR(128)
    GroupName_TMP(*)=0
    if (N_ELEMENTS(GroupName) GT 0) THEN if (STRLEN(GroupName) GT 0) THEN GroupName_TMP(0:STRLEN(GroupName)-1)=BYTE(GroupName)

 PIOGETROINAME=call_external(PIOLibIDLSO,'piogetroiname_tempoidl', $
        RoiName_TMP, $
        GroupName_TMP, $
               /L64_VALUE) 
  RoiName=STRING(RoiName_TMP)
  GroupName=STRING(GroupName_TMP)

 RETURN,PIOGETROINAME
END
