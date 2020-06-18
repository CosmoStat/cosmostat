function piostring2c,stin
;================================================================
;
; This function ensure the 128 character size of the input string
;
;================================================================

st=string(stin)
a=strlen(st)
tab=bytarr(129)
tab(127:128)=0

if (a gt 0) then begin
    if (a gt 126) then begin
        tt=byte(stin)
        tab(0:126)=tt(0:126)
    end ELSE begin
        tab(0:a-1)=byte(stin)
        tab(a:*)=0
    end
end else begin
    tab(*)=0
end
return,tab
end
