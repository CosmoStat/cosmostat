;+
; NAME:
;        mrsp_alm2spec
;
; PURPOSE:
;
;   Computes the power spectrums and cross spectrums of a polarized map from the A_{l,m} coefficients.  
;
; CALLING:
;     P = mrsp_alm2spec(ALM, StdPS=StdPS)
;
; INPUTS:
;     ALM -- IDL  structure: Alm coefficients (see mrsp_almtrans.pro). 
;    
; OUTPUTS:
;      P -- IDL fltarr: the TT, EE, BB, TE, TB, EB spectrums. P[k,i] = Mean(  SPECTRUM[*,l,i]  )	i=0...5
;
; OUTPUT KEYWORDS:
;    StdPS  -- 2D IDL fltarr:standard deviation of spectrums:   StdPS[l,i] = stddev( PowSpecIma[l, -l:l ] )
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the spectrums of ALM from a polarized image. 
;               P = mrsp_alm2spec( ALM ) 
;         
; HISTORY:
;	Written: Olivier Fourt, 2009
;	July, 2009 File creation
;--------------------------------------------------------------------------------------------------------

;=====================================================================

function mrsp_alm2spec, ALM, StdPS=StdPS

COMMON MR1ENV
if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: mrsp_alm2spec, ALM, StdPS=StdPS'
	Ret=-1
        goto, DONE
        end

NbrAlm = ALM.LMAX+1
TabM = ALM.TABNBRM

if ALM.COMPLEX_ALM EQ 1 then begin

	UScomplex = 1
	if ALM.TAB EQ 0 then begin
		alm_pola2tab_complex_in, ALM.ALM, TabALM, complex=UScomplex, TabNbrM=TabM
	end else begin
		TabALM = ALM.ALM
	end
	
end else begin

	UScomplex=0
	if ALM.TAB EQ 0 then begin
		alm_pola2tab, ALM.ALM, TabALM, complex=UScomplex, TabNbrM=TabM
	end else begin
		TabALM = ALM.ALM
	end
	
end

if ALM.complex_alm EQ 1 then begin
	
	alm_t = TabALM[*,*,0]
	alm_e = TabALM[*,*,1]
	alm_b = TabALM[*,*,2]
	
	TabALM_TT = re(alm_t[*,*])*re(alm_t[*,*]) + im(alm_t[*,*])*im(alm_t[*,*])
	TabALM_EE = re(alm_e[*,*])*re(alm_e[*,*]) + im(alm_e[*,*])*im(alm_e[*,*])
	TabALM_BB = re(alm_b[*,*])*re(alm_b[*,*]) + im(alm_b[*,*])*im(alm_b[*,*])
	TabALM_TE = re(alm_t[*,*])*re(alm_e[*,*]) + im(alm_t[*,*])*im(alm_e[*,*])
	TabALM_TB = re(alm_t[*,*])*re(alm_b[*,*]) + im(alm_t[*,*])*im(alm_b[*,*])
	TabALM_EB = re(alm_e[*,*])*re(alm_b[*,*]) + im(alm_e[*,*])*im(alm_b[*,*])

end else begin

	alm_t = TabALM[*,*,*,0]
	alm_e = TabALM[*,*,*,1]
	alm_b = TabALM[*,*,*,2]
	
	TabALM_TT = alm_t[*,*,0]*alm_t[*,*,0] + alm_t[*,*,1]*alm_t[*,*,1]
	TabALM_EE = alm_e[*,*,0]*alm_e[*,*,0] + alm_e[*,*,1]*alm_e[*,*,1]
	TabALM_BB = alm_b[*,*,0]*alm_b[*,*,0] + alm_b[*,*,1]*alm_b[*,*,1]
	TabALM_TE = alm_t[*,*,0]*alm_e[*,*,0] + alm_t[*,*,1]*alm_e[*,*,1]
	TabALM_TB = alm_t[*,*,0]*alm_b[*,*,0] + alm_t[*,*,1]*alm_b[*,*,1]
	TabALM_EB = alm_e[*,*,0]*alm_b[*,*,0] + alm_e[*,*,1]*alm_b[*,*,1]

end

Ret = fltarr(NbrAlm,6)
StdPS = fltarr(NbrAlm,6)

for l=0, NbrAlm-1 do  begin

	allalm = 2*TabM[l]-1
	
	if allalm GT 1 then begin
	
		P = fltarr(allalm)
		
    	P[0:TabM[l]-2] = reverse( TabALM_TT[l, 1:TabM[l]-1])
    	P[TabM[l]-1:*] = TabALM_TT[l, 0:TabM[l]-1]
    	Ret[l,0] = mean(P);	Spec TT
    	StdPS[l,0] = sigma(P)
    	
    	P[0:TabM[l]-2] = reverse( TabALM_EE[l, 1:TabM[l]-1])
    	P[TabM[l]-1:*] = TabALM_EE[l, 0:TabM[l]-1]
    	Ret[l,1] = mean(P);	Spec EE
    	StdPS[l,1] = sigma(P)
    	
    	P[0:TabM[l]-2] = reverse( TabALM_BB[l, 1:TabM[l]-1])
    	P[TabM[l]-1:*] = TabALM_BB[l, 0:TabM[l]-1]
    	Ret[l,2] = mean(P);	Spec BB
    	StdPS[l,2] = sigma(P)
    	
    	P[0:TabM[l]-2] = reverse( TabALM_TE[l, 1:TabM[l]-1])
    	P[TabM[l]-1:*] = TabALM_TE[l, 0:TabM[l]-1]
    	Ret[l,3] = mean(P);	Spec TE
    	StdPS[l,3] = sigma(P)
    	
    	P[0:TabM[l]-2] = reverse( TabALM_TB[l, 1:TabM[l]-1])
    	P[TabM[l]-1:*] = TabALM_TB[l, 0:TabM[l]-1]
    	Ret[l,4] = mean(P);	Spec TB
    	StdPS[l,4] = sigma(P)
    	
    	P[0:TabM[l]-2] = reverse( TabALM_EB[l, 1:TabM[l]-1])
    	P[TabM[l]-1:*] = TabALM_EB[l, 0:TabM[l]-1]
    	Ret[l,5] = mean(P);	Spec EB
    	StdPS[l,5] = sigma(P)
    	
    end else begin
    
    	Ret[l,0] = TabALM_TT[l,0]
    	Ret[l,1] = TabALM_EE[l,0]
    	Ret[l,2] = TabALM_BB[l,0]
    	Ret[l,3] = TabALM_TE[l,0]
    	Ret[l,4] = TabALM_TB[l,0]
    	Ret[l,5] = TabALM_EB[l,0]
    	
    	StdPS[l,0] = 0.0
    	StdPS[l,1] = 0.0
    	StdPS[l,2] = 0.0
    	StdPS[l,3] = 0.0
    	StdPS[l,4] = 0.0
    	StdPS[l,5] = 0.0
    	
	end
end

DONE:
return, Ret
end

