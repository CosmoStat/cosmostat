
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOClosePOIObject

FUNCTION PIOCLOSEPOIOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyObject=long64(MyObject)

 PIOCLOSEPOIOBJECT=call_external(PIOLibIDLSO,'pioclosepoiobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOCLOSEPOIOBJECT
END
