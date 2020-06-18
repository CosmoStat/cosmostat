
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseTAB2DGrp

FUNCTION PIOCLOSETAB2DGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyGroup=long64(MyGroup)

 PIOCLOSETAB2DGRP=call_external(PIOLibIDLSO,'pioclosetab2dgrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCLOSETAB2DGRP
END
