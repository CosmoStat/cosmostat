
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOSetLonpoleGrp

FUNCTION PIOSETLONPOLEGRP,LonPole,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(LonPole) EQ 0) THEN     LonPole=DOUBLE(0) ELSE    LonPole=DOUBLE(LonPole)
 
 PIOSETLONPOLEGRP=call_external(PIOLibIDLSO,'piosetlonpolegrp_tempoidl', $
        LonPole, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOSETLONPOLEGRP
END
