
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseVECTObject

FUNCTION PIOCLOSEVECTOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyObject=long64(MyObject)

 PIOCLOSEVECTOBJECT=call_external(PIOLibIDLSO,'pioclosevectobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOCLOSEVECTOBJECT
END
