
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseTAB3DGrp

FUNCTION PIOCLOSETAB3DGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyGroup=long64(MyGroup)

 PIOCLOSETAB3DGRP=call_external(PIOLibIDLSO,'pioclosetab3dgrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCLOSETAB3DGRP
END
