
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseIMO

FUNCTION PIOCLOSEIMO,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyGroup=long64(MyGroup)

 PIOCLOSEIMO=call_external(PIOLibIDLSO,'piocloseimo_tempoidl', $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOCLOSEIMO
END
