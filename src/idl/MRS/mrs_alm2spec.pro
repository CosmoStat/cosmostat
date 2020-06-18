;+
; NAME:
;        mrs_alm2spec
;
; PURPOSE:
;
;   Compute the power spectrum C(l) from the A_{l,m} coefficients.
;              C(l) = 1/(2l+1)  sum_{m=-l}^{m=l} | A_{l,m} |^2   
;
; CALLING:
;     P = mrs_alm2spec(ALM, StdPS=StdPS)
;
; INPUTS:
;     ALM -- IDL  structure: Alm coefficient (see mrs_almtransform.pro). 
;    
; OUTPUTS:
;     P -- 1D IDL fltarr: Power Spectrum (C(l)).
;
; OUTPUT KEYWORDS:
;    StdPS  -- 1D IDL fltarr:standard deviation of Cl:   StdPS[l] = stddev( PowSpecIma[l, -l:l ]
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the power spectrum of an image. 
;               P = mrs_powspec(Imag) 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2006
;	January, 2006 File creation
;--------------------------------------------------------------------------------------------------------

;=====================================================================

function mrs_alm2spec, ALM, StdPS=StdPS

NbrAlm = ALM.LMAX+1
TabM = ALM.TABNBRM
;d = ALM.ALM

if ALM.COMPLEX_ALM EQ 1 then begin

	UScomplex = 1
	if ALM.TAB EQ 0 then begin
		alm2tab_complex_in, ALM.ALM, TabALM, complex=UScomplex, TabNbrM=TabM
	end else begin
		TabALM = ALM.ALM
	end
	
end else begin

	UScomplex=0
	if ALM.TAB EQ 0 then begin
		alm2tab, ALM.ALM, TabALM, complex=UScomplex, TabNbrM=TabM
	end else begin
		TabALM = ALM.ALM
	end
	
end

if ALM.complex_alm  EQ 1 then begin
	TabALM = TabALM * conj(TabALM)
end else begin
       if ALM.complex_alm EQ 0 then begin
       		TabALM = TabALM[*,*,0]^2 + TabALM[*,*,1]^2
       end else TabALM = TabALM[*,*,0]
end

PowSpecIma = TabALM
Ret = fltarr(NbrAlm)
StdPS = fltarr(NbrAlm)
for i=0, NbrAlm-1 do begin
   allalm = 2*TabM[i]-1
   if allalm GT 1 then begin
      P = fltarr(allalm)
      P[0:TabM[i]-2] = reverse( PowSpecIma[i, 1:TabM[i]-1])
      P[TabM[i]-1:*] = PowSpecIma[i, 0:TabM[i]-1]
      Ret[i] = mean( P)      
      StdPS[i] = sigma(P)
   end else begin
       Ret[i] = PowSpecIma[i, 0]
       StdPS[i] = 0.
   end
end

return, Ret
   
end

;=====================================================================
 
function mrs_powalm, ALM, nalm=nalm

NbrAlm = ALM.LMAX+1
TabM = ALM.TABNBRM
d = ALM.ALM

if ALM.COMPLEX_ALM EQ 1 then UScomplex = 1 else UScomplex=0


if ALM.TAB EQ 0 then  alm2tab, ALM.ALM, TabALM, complex=UScomplex, TabNbrM=TabM $
else TabALM = ALM.ALM

if ALM.complex_alm  EQ 1 then  TabALM = TabALM * conj(TabALM) $
else begin
       if ALM.complex_alm  EQ 0 then TabALM = TabALM[*,*,0]^2 + TabALM[*,*,1]^2 $
       else TabALM = TabALM[*,*,0]
end

PowSpecIma = TabALM
Ret = fltarr(NbrAlm)
nalm=0
for i=0, NbrAlm-1 do begin
   allalm = 2*TabM[i]-1
   nalm=nalm+allalm
   if allalm GT 1 then begin
      P = fltarr(allalm)
      P[0:TabM[i]-2] = reverse( PowSpecIma[i, 1:TabM[i]-1])
      P[TabM[i]-1:*] = PowSpecIma[i, 0:TabM[i]-1]
      Ret[i] = total( P)
   end else Ret[i] = PowSpecIma[i, 0]
end
TotPow = total(Ret)
return, TotPow
   
end
