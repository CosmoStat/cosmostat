
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOCloseTOIObject

FUNCTION PIOCLOSETOIOBJECT,MyObject
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    MyObject=long64(MyObject)

 PIOCLOSETOIOBJECT=call_external(PIOLibIDLSO,'pioclosetoiobject_tempoidl', $
        MyObject, $
               /L64_VALUE) 

 RETURN,PIOCLOSETOIOBJECT
END
