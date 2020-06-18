

PRO PIOERRMESS,ERRUNIT

PIOLibIDLSO=shared_lib_path('PIOLibIDL.so')
size=long64(0)
mess=bytarr(128)

a=call_external(PIOLibIDLSO,'pioerrmessidl',long64(ERRUNIT),long64(size),mess,/L64_VALUE)

mess=strmid(string(mess),0,size)
if (ERRUNIT NE 0) THEN BEGIN
    MESSAGE,mess,/continue
END
END
