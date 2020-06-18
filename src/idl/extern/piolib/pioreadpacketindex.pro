
;--------------------------------------------------------------#
;
; AUTOMATICALLY GENERATED DO NOT MODIFY 

;--------------------------------------------------------------#




;	Wrapping PIOReadPacketIndex

FUNCTION PIOREADPACKETINDEX,data_timesec,data_period,data_ndata,MyGroup
     ON_ERROR,1
     PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
    data_timesec=long64(0)
    data_period=long64(0)
    data_ndata=long64(0)

 PIOREADPACKETINDEX=call_external(PIOLibIDLSO,'pioreadpacketindex_tempoidl', $
        data_timesec, $
        data_period, $
        data_ndata, $
        MyGroup, $
               /L64_VALUE) 
  IF (PIOREADPACKETINDEX GT 0) THEN BEGIN
     TP=data_timesec
     data_timesec=DBLARR(PIOREADPACKETINDEX)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,data_timesec, LONG64(PIOREADPACKETINDEX*8))
    END
  IF (PIOREADPACKETINDEX GT 0) THEN BEGIN
     TP=data_period
     data_period=LON64ARR(PIOREADPACKETINDEX)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,data_period, LONG64(PIOREADPACKETINDEX*8))
    END
  IF (PIOREADPACKETINDEX GT 0) THEN BEGIN
     TP=data_ndata
     data_ndata=LON64ARR(PIOREADPACKETINDEX)
     A=call_external(PIOLibIDLSO,'piomemcpyidl',TP,data_ndata, LONG64(PIOREADPACKETINDEX*8))
    END

 RETURN,PIOREADPACKETINDEX
END
