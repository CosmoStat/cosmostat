
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseIMG2DGrp

FUNCTION PIOCLOSEIMG2DGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyGroup=long64(MyGroup)

 PIOCLOSEIMG2DGRP=call_external(PIOLibIDLSO,'piocloseimg2dgrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCLOSEIMG2DGRP
END
