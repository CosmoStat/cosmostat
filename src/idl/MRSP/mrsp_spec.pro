;+
; NAME:
;        mrsp_spec
;
; PURPOSE:
;   Computes the power spectrums and cross spectrums of a polarized map, 
;   using the HEALPix representation (nested data representation by default)
;
; CALLING:
;     P = mrsp_spec( Imag, nonorm=nonorm, teb=teb, NormVal=NormVal, StdPS=StdPS, lmax=lmax )
;
; INPUTS:
;     Imag -- IDL array of healpix polarized map (ie3 maps) in TQU or TEB (default TQU): Input image to be transformed 
;    
; OUTPUTS:
;     P -- IDL fltarr: the TT, EE, BB, TE, TB, EB spectrums. P[k,i] = Mean(  SPECTRUM[*,l,i]  )	i=0...5
;
; INPUT KEYWORDS:
;		teb -- map is in TEB
;		nonorm -- no normalization of alm
;		Lmax  -- int: Number of spherical harmonics computed in the decomposition and size of the computed spectrum (Lmax+1)
;					  (HEALPIX==> default is 3*nside, should be between 2*nside and 4*nside)
;
;	OUTPUT KEYWORDS:
;		StdPS  -- 2D IDL fltarr:standard deviation of spectrums:   StdPS[l,i] = stddev( PowSpecIma[l, -l:l ] )
;		NormVal -- float: normalization value of the alm coefficients
;
;
; EXTERNAL CALLS:
;       mrs_almtrans.pro
;		mrsp_almtrans.pro
;
; EXAMPLE:
;       Compute the power spectrum of an image. 
;               P = mrsp_spec(Imag) 
;         
; HISTORY:
;	Written: Olivier Fourt, 2009
;	June, 2009 File creation
;--------------------------------------------------------------------------------------------------------

;=====================================================================
 
function mrsp_spec, Imag, nonorm=nonorm, teb=teb, NormVal=NormVal, StdPS=StdPS, lmax=lmax
COMMON MR1ENV

if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: mrsp_spec, Imag, nonorm=nonorm, teb=teb, NormVal=NormVal, StdPS=StdPS, lmax=lmax'
	Ret=-1
        goto, DONE
        end

L2PowSpecNorm=0	

if DEF_NORM_POWSPEC EQ 1 then L2PowSpecNorm=1

if keyword_set(nonorm) then L2PowSpecNorm=0

if keyword_set(teb) then begin

	mrs_almtrans, Imag[*,0], ALM_T, lmax=lmax, /tab, norm=L2PowSpecNorm
	mrs_almtrans, Imag[*,1], ALM_E, lmax=lmax, /tab, norm=L2PowSpecNorm
	mrs_almtrans, Imag[*,2], ALM_B, lmax=lmax, /tab, norm=L2PowSpecNorm
	
	NormVal = ALM_T.NormVal
	NbrAlm = ALM_T.LMAX+1
	
	TabM = ALM_T.TABNBRM
	alm_t = ALM_T.ALM
	alm_e = ALM_E.ALM
	alm_b = ALM_B.ALM
	
end else begin

	mrsp_almtrans, Imag, ALM, lmax=lmax, /tab, norm=L2PowSpecNorm
	
	NormVal = ALM.NormVal
	NbrAlm = ALM.LMAX+1
	
	TabM = ALM.TABNBRM
	alm_t = ALM.ALM[*,*,*,0]
	alm_e = ALM.ALM[*,*,*,1]
	alm_b = ALM.ALM[*,*,*,2]
	
end

Ret = dblarr(NbrAlm,6)
StdPS = dblarr(NbrAlm,6)

TabALM_TT = alm_t[*,*,0]*alm_t[*,*,0] + alm_t[*,*,1]*alm_t[*,*,1]
TabALM_EE = alm_e[*,*,0]*alm_e[*,*,0] + alm_e[*,*,1]*alm_e[*,*,1]
TabALM_BB = alm_b[*,*,0]*alm_b[*,*,0] + alm_b[*,*,1]*alm_b[*,*,1]
TabALM_TE = alm_t[*,*,0]*alm_e[*,*,0] + alm_t[*,*,1]*alm_e[*,*,1]
TabALM_TB = alm_t[*,*,0]*alm_b[*,*,0] + alm_t[*,*,1]*alm_b[*,*,1]
TabALM_EB = alm_e[*,*,0]*alm_b[*,*,0] + alm_e[*,*,1]*alm_b[*,*,1]

for l=0, NbrAlm-1 do  begin

	allalm = 2*TabM[l]-1
	
	if allalm GT 1 then begin
	
		P = dblarr(allalm)
		
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
END
