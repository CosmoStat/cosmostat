;+
; NAME:
;    BWT01_INVERSE
;
; PURPOSE:
;    Computes the inverse wavelet transform of bidimensional wavelet 
;    coefficient array. 
;
; CATEGORY:
;    Statistical tools for map analysis.
;
; CALLING SEQUENCE:
;     W = BWT01_INVERSE( ARRAY, N_LEVELS)
; 
; INPUTS:
;     ARRAY:       Input wavelet coefficient map as a bidimensional array
;     N_LEVELS:    Number of decomposition levels as an integer 
;
; OPTIONAL INPUTS:
;     
;
;	
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;     W:          Map, bidimensional array, same size as the wavelet 
;                 coefficient map (NX, NY).     
;
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;     If too many decomposition levels are requested: 
;          A warning is issued giving the maximum allowed N_LEVELS 
;          The wavelet decomposition is performed with this new value
;
; RESTRICTIONS:
;     Procedure is stopped if: The input array is not bidimensional
;
; PROCEDURE: BWT01_INV_STEP
;
; EXAMPLE:
;      map = RANDOMN(seed, 512, 512)
;      W = BWT01_DIRECT( map, 4)
;      map_inverse = BWT01_INVERSE( W, 4)
;
; MODIFICATION HISTORY:
;      15/05/2002      Olivier FORNI (IAS) For Planck L2-DPC 
;-

FUNCTION bwt01_inv_step,z,INVERSE=inverse

common cf01,a,b,c,d,e

w=float(z)
sw=size(w)

mn1= sw(1) mod 2
n1=sw(1)/2+mn1

ie=indgen(n1)*2
io=ie+1

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


return,w
end


FUNCTION bwt01_inverse,z,LV

common cf01,a,b,c,d,e

w=float(z)
sw0=size(z)

a=-1.586134342
b=-.05298011854
c=0.8829110762
d=0.4435068511
e=1.149604398

case sw0(0) of
   2: begin
            nx=sw0(1)
            ny=sw0(2)
    end

    else : begin    print,'The input array must be bidimensional'
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

n2l=lonarr(lv)
for l=0,lv-1 do begin 
    n2l(l)=n2
    mn2= n2 mod 2
    n2=n2/2+mn2
end        

mn1=min(n1l) & imn1=lv & nnv=0
if (mn1 le 4) then begin
    im=where(n1l gt 4,imn1)
    lvn=imn1
    n1l=n1l(0:lvn-1)
    nnv=1
end
mn2=min(n2l) & imn2=lv 
if (mn2 le 4) then begin
    im=where(n2l gt 4,imn2)
    nnv=1
    lvn=min([imn1,imn2])
    n1l=n1l(0:lvn-1)
    n2l=n2l(0:lvn-1)
end

if(nnv eq 1) then begin
    lv=lvn
    print,'Too many decomposition levels'
    print,format='("The number of levels N_LEVELS is set to",i3)',lvn
end



for l=lv-1,0,-1 do begin    
    for i=0,n2l(l)-1 do w(0:n1l(l)-1,i)=bwt01_inv_step(w(0:n1l(l)-1,i))
    for i=0,n1l(l)-1 do w(i,0:n2l(l)-1)=bwt01_inv_step(transpose(w(i,0:n2l(l)-1)))
end

return,w

end

 
