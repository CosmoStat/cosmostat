;+
; NAME:
;    BWT01_DIRECT
;
; PURPOSE:
;    Computes the wavelet transform of bidimensional array. 
;    Returns the first four moments of the wavelet coefficients.
;
; CATEGORY:
;    Statistical tools for map analysis.
;
; CALLING SEQUENCE:
;     W = BWT01_DIRECT( ARRAY, N_LEVELS, MOMENT)
; 
; INPUTS:
;     ARRAY:       Input map as a bidimensional array (NX, NY)
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
;     W:          Wavelet coefficient map, bidimensional array, same
;                 size as input map (NX, NY).     
;     MOMENT:     A named variable containing the first, second, third
;                 and fourth moment of the distribution of wavelet 
;                 coefficients.
;                 Array (4, 3, N_LEVELS). The first dimension contains
;                 the moments ordered as mentioned above. The second 
;                 dimension corresponds to the horizontal, vertical and
;                 cross coefficients. The third dimension represents
;                 the decomposition level.
;
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;     If too many decomposition levels are requested: 
;          A warning is issued giving the maximum allowed N_LEVELS 
;          The wavelet decomposition is performed with this new value
;
; RESTRICTIONS:
;    If the input array is not bidimensional, the procedure stops and returns
;    an expression = -1.
;                          
; PROCEDURE: BWT01_DIR_STEP, CSCAL
;
; EXAMPLE:
;      map = RANDOMN(seed, 512, 512)
;      W = BWT01_DIRECT( map, 4, moment)
;      print, moment(*, 0, 2)
;      ; Prints the four moments of the horizontal wavelet
;      ; coefficients at the third decomposition level
;
; MODIFICATION HISTORY:
;      15/05/2002      Olivier FORNI (IAS) For Planck L2-DPC 
;-
function cscal,w,sc,nc

b=size(w)
n1=b(1)
n2=b(2)
nt=min([n1,n2])
lv=fix(alog(nt)/alog(2))
 
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

i2=n1l(sc)
i1=n1l(sc-1)-1
j2=n2l(sc)
j1=n2l(sc-1)-1


case nc of
    1:wsc=w(i2:i1,0:j2-1)
    2:wsc=w(0:i2-1,j2:j1)
    3:wsc=w(i2:i1,j2:j1)
endcase    


return,wsc
end

FUNCTION bwt01_dir_step,z

common cf01,a,b,c,d,e

w=float(z)
sw=size(w)

mn1= sw(1) mod 2
n1=sw(1)/2+mn1

ie=indgen(n1)*2
io=ie+1
    
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

return,w
end


FUNCTION bwt01_direct,z,LV,mom

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


for l=0,lv-1 do begin    
    for i=0,n2l(l)-1 do w(0:n1l(l)-1,i)=bwt01_dir_step(w(0:n1l(l)-1,i))
    for i=0,n1l(l)-1 do w(i,0:n2l(l)-1)=bwt01_dir_step(transpose(w(i,0:n2l(l)-1)))
end

mom=fltarr(4,3,lv)

for l=0,lv-1 do begin
    for i=0,2 do begin
        w0=cscal(w,l+1,i+1)
        mom(*,i,l)=moment(w0)
    end
end

return,w

end

 
