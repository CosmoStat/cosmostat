
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOClosePBRGrp

FUNCTION PIOCLOSEPBRGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyGroup=long64(MyGroup)

 PIOCLOSEPBRGRP=call_external(PIOLibIDLSO,'pioclosepbrgrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCLOSEPBRGRP
END
