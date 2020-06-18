FUNCTION bwt01_lift_1D,z,INVERSE=inverse

common cf01,a,b,c,d,e

w=float(z)
sw=size(w)

mn1= sw(1) mod 2
n1=sw(1)/2+mn1

ie=indgen(n1)*2
io=ie+1

if(not keyword_set(inverse)) then begin

    ze=fltarr(n1+4)
    zo=fltarr(n1+4)
    ze(2:n1+1)=w(ie) 
    zo(2:n1+1-mn1)=w(io(0:n1-1-mn1))

    ze(1)=ze(3) & ze(0)=ze(4) 
    zo(1)=zo(2) & zo(0)=zo(3) 
    if ( mn1 eq 0) then begin
        ze(n1+2)=ze(n1+1) & ze(n1+3)=ze(n1)
        zo(n1+2)=zo(n1) & zo(n1+3)=zo(n1-1)
    end else begin
        ze(n1+2)=ze(n1) & ze(n1+3)=ze(n1-1)
        zo(n1+1)=zo(n1) & zo(n1+2)=zo(n1-1) & zo(n1+3)=zo(n1-2)
    end

    zo=zo+a*(ze+shift(ze,-1))
    ze=ze+b*(zo+shift(zo,1))
    zo=zo+c*(ze+shift(ze,-1))
    ze=ze+d*(zo+shift(zo,1))
    ze=ze*e
    zo=zo/e

    w(0:n1-1)=ze(2:n1+1)
    w(n1:2*n1-1-mn1)=zo(2:n1+1-mn1)     
endif

if (keyword_set(inverse)) then begin

    ze=fltarr(n1+4)
    ze(2:n1+1)=w(0:n1-1) 
    zo=fltarr(n1+4)
    zo(2:n1+1-mn1)=w(n1:2*n1-1-mn1)

    ze(1)=ze(3) & ze(0)=ze(4) 
    zo(1)=zo(2) & zo(0)=zo(3) 
    if ( mn1 eq 0) then begin
        ze(n1+2)=ze(n1+1) & ze(n1+3)=ze(n1)
        zo(n1+2)=zo(n1) & zo(n1+3)=zo(n1-1)
    end else begin
        ze(n1+2)=ze(n1) & ze(n1+3)=ze(n1-1)
        zo(n1+1)=zo(n1) & zo(n1+2)=zo(n1-1) & zo(n1+3)=zo(n1-2)
    end
    ze=ze/e
    zo=zo*e
    ze=ze-(zo+shift(zo,1))*d
    zo=zo-(ze+shift(ze,-1))*c
    ze=ze-(zo+shift(zo,1))*b
    zo=zo-(ze+shift(ze,-1))*a

    w(ie)=ze(2:n1+1)
    w(io(0:n1-1-mn1))=zo(2:n1+1-mn1)       
endif

return,w
end


FUNCTION bwt01_lift,z,LV,INVERSE=inverse

common cf01,a,b,c,d,e

w=float(z)
sw0=size(z)

a=-1.586134342
b=-.05298011854
c=0.8829110762
d=0.4435068511
e=1.149604398

case sw0(0) of
    1: begin
        nx=sw0(1)
        ny=1
;        PYRAMIDAL=1
    end
    2: begin
        if(sw0(1) eq 1) then begin 
;            PYRAMIDAL=1
            w=transpose(w)
            nx=sw0(2)
            ny=1
        end else begin
            nx=sw0(1)
            ny=sw0(2)
        end
    end

    else : begin    print,'Le tableau doit etre uni ou bidimensionnel'
        return,-1
    end

end

n1=nx
n2=ny

n1l=lonarr(lv)
for l=0,lv-1 do begin 
    n1l(l)=n1
    mn1= n1 mod 2
    n1=n1/2+mn1    
end   

 

if(ny gt 1) then begin
    n2l=lonarr(lv)
    for l=0,lv-1 do begin 
        n2l(l)=n2
        mn2= n2 mod 2
        n2=n2/2+mn2
    end        
end

mn1=min(n1l) & imn1=lv & nnv=0
if (mn1 le 4) then begin
    im=where(n1l gt 4,imn1)
    lvn=imn1
    n1l=n1l(0:lvn-1)
    nnv=1
end
if(ny gt 1) then begin
    mn2=min(n2l) & imn2=lv 
    if (mn2 le 4) then begin
        im=where(n2l gt 4,imn2)
        nnv=1
        lvn=min([imn1,imn2])
        n1l=n1l(0:lvn-1)
        n2l=n2l(0:lvn-1)
    end
end 
if(nnv eq 1) then begin
    lv=lvn
    print,'Nombre de niveaux trop grand'
    print,format='("Le nombre de niveaux est fixé à",i3)',lvn
end

if (not keyword_set(inverse)) then begin 
    if(ny eq 1) then begin       
        for l=0,lv-1 do begin    
            w(0:n1l(l)-1)=bwt01_lift_1d(w(0:n1l(l)-1))
        end
    end else begin
        for l=0,lv-1 do begin    
            for i=0,n2l(l)-1 do w(0:n1l(l)-1,i)=bwt01_lift_1d(w(0:n1l(l)-1,i))
            for i=0,n1l(l)-1 do w(i,0:n2l(l)-1)=bwt01_lift_1d(transpose(w(i,0:n2l(l)-1)))
        end
    endelse
end

if (keyword_set(inverse)) then begin 
    if(ny eq 1) then begin
        for l=lv-1,0,-1 do begin    
            w(0:n1l(l)-1)=bwt01_lift_1d(w(0:n1l(l)-1),/inv)
        end
    end else begin
        for l=lv-1,0,-1 do begin    
            for i=0,n2l(l)-1 do w(0:n1l(l)-1,i)=bwt01_lift_1d(w(0:n1l(l)-1,i),/inv)
            for i=0,n1l(l)-1 do w(i,0:n2l(l)-1)=bwt01_lift_1d(transpose(w(i,0:n2l(l)-1)),/inv)
        end
    endelse
end

return,w

end


