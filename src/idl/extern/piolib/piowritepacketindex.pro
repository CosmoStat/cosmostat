
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOWritePacketIndex

FUNCTION PIOWRITEPACKETINDEX,data_timesec,data_period,data_ndata,Ring_num,NbNum,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    if (N_ELEMENTS(data_timesec) EQ 0) THEN     data_timesec=DOUBLE(0) ELSE    data_timesec=DOUBLE(data_timesec)
     if (N_ELEMENTS(data_period) EQ 0) THEN     data_period=LONG64(0) ELSE    data_period=LONG64(data_period)
     if (N_ELEMENTS(data_ndata) EQ 0) THEN     data_ndata=LONG64(0) ELSE    data_ndata=LONG64(data_ndata)
     if (N_ELEMENTS(Ring_num) EQ 0) THEN     Ring_num=LONG64(0) ELSE    Ring_num=LONG64(Ring_num)
     if (N_ELEMENTS(NbNum) EQ 0) THEN     NbNum=LONG64(0) ELSE    NbNum=LONG64(NbNum)
 
 PIOWRITEPACKETINDEX=call_external(PIOLibIDLSO,'piowritepacketindex_tempoidl', $
        data_timesec, $
        data_period, $
        data_ndata, $
        Ring_num, $
        NbNum, $
        MyGroup, $
               /L64_VALUE) 

 RETURN,PIOWRITEPACKETINDEX
END
