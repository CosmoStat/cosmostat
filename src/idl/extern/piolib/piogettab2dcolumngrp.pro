
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetTAB2DColumnGrp

FUNCTION PIOGETTAB2DCOLUMNGRP,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')

 PIOGETTAB2DCOLUMNGRP=call_external(PIOLibIDLSO,'piogettab2dcolumngrp_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETTAB2DCOLUMNGRP
END
