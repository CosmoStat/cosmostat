
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOGetLonpoleGrp

FUNCTION PIOGETLONPOLEGRP,LonPole,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(LonPole) EQ 0) THEN     LonPole=DOUBLE(0) ELSE    LonPole=DOUBLE(LonPole)
 
 PIOGETLONPOLEGRP=call_external(PIOLibIDLSO,'piogetlonpolegrp_tempoidl', $
        LonPole, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOGETLONPOLEGRP
END
