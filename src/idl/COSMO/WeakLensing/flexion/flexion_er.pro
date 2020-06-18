pro flexion_er,img,flexion,xc,yc

nx=(size(img))[1]

i0=lonarr(nx*nx)
j0=lonarr(nx*nx)

for i=0L,nx*nx-1 do begin
    i0[i]=i mod nx
    j0[i]=i/nx
endfor

flux=total(img)
xc=total(i0*img)/flux
yc=total(j0*img)/flux

Q11=total((i0-xc)^2*img)/flux
Q12=total((i0-xc)*(j0-yc)*img)/flux
Q22=total((j0-yc)^2*img)/flux

Q111=total((i0-xc)^3*img)/flux
Q112=total((i0-xc)^2*(j0-yc)*img)/flux
Q122=total((i0-xc)*(j0-yc)^2*img)/flux
Q222=total((j0-yc)^3*img)/flux

Q1111=total((i0-xc)^4*img)/flux
Q1112=total((i0-xc)^3*(j0-yc)*img)/flux
Q1122=total((i0-xc)^2*(j0-yc)^2*img)/flux
Q1222=total((i0-xc)*(j0-yc)^3*img)/flux
Q2222=total((j0-yc)^4*img)/flux

q2=complex(Q11-Q22,2*Q12)
q0=Q11+Q22

T3=complex(Q111-3*Q122,3*Q112-Q222)
T1=complex(Q111+Q122,Q112+Q222)

F4=complex(Q1111-6*Q1122+Q2222,4*(Q1112-Q1222))
F2=complex(Q1111-Q2222,2*(Q1112+Q1222))
F0=Q1111+2*Q1122+Q2222

flexion=4./(9*F0+12*q0^2)*T1

return
end
