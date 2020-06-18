
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCreatePBRGrp

FUNCTION PIOCREATEPBRGRP,Groupname,ROIname,RingSize,Init_Size
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    Groupname_TMP=BYTARR(128)
    Groupname_TMP(*)=0
    if (N_ELEMENTS(Groupname) GT 0) THEN if (STRLEN(Groupname) GT 0) THEN Groupname_TMP(0:STRLEN(Groupname)-1)=BYTE(Groupname)
    ROIname_TMP=BYTARR(128)
    ROIname_TMP(*)=0
    if (N_ELEMENTS(ROIname) GT 0) THEN if (STRLEN(ROIname) GT 0) THEN ROIname_TMP(0:STRLEN(ROIname)-1)=BYTE(ROIname)
    if (N_ELEMENTS(RingSize) EQ 0) THEN     RingSize=LONG64(0) ELSE    RingSize=LONG64(RingSize)
 IF (N_PARAMS() LT 4) THEN BEGIN
    Init_Size=long64(0)
END ELSE BEGIN
    Init_Size=long64(Init_Size)
END

 PIOCREATEPBRGRP=call_external(PIOLibIDLSO,'piocreatepbrgrp_tempoidl', $
        Groupname_TMP, $
        ROIname_TMP, $
        RingSize, $
        Init_Size, $
               /L64_VALUE) 

 RETURN,PIOCREATEPBRGRP
END
