function master_wigner3j2,il1,il2,il3,lnwpro
l1=long(il1)
l2=long(il2)
l3=long(il3)
L=l1+l2+l3
L_2=L/2
min=abs(l1-l2)
max=l1+l2
c=l3*0d
w=where((L_2*2-L) eq 0 and l3 ge min and l3 le max)
if w(0) ne -1 then begin
    lnw1=lnwpro(L_2(w)-l1)
    lnw2=lnwpro(L_2(w)-l2)
    lnw3=lnwpro(L_2(w)-l3(w))
    lnwl=lnwpro(L_2(w))
    lnc=-alog(L(w)+1d)-lnwl+lnw1+lnw2+lnw3
    c(w)=exp(lnc)
endif
return,c
end
