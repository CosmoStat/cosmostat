;+
; NAME:
;        MR_MODPHASE_DWT_TRANS
;
; PURPOSE:
;	Compute the decimated wavelet transform of a vector field (ex Q-U CMB data) image.
;   It computes first the modulus and the phase, and run a decimated WT on each.
;   The output is a IDL structure.
;
; CALLING:
;      MR_MODPHASE_DWT_TRANS, Imag, Trans, NbrScale=NbrScale
;
; INPUTS:
;     Imag -- IDL array of a vector field image fltarr[*,*,2] : Input image be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;         ModCoeff      -- 2D IDL array [*,*] : Wavelet coefficients of the modulus
;         AngCoeff      -- 2D IDL array [*,*] : Wavelet coefficients of the phase
;
; KEYWORDS:
;         NBRSCALE  -- LONG: Number of scales of the wavelet transform
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the orthogonal wavelet transform of a vector field I with five scales
;       The result is stored in WT
;               mr_modphase_dwt_trans, Imag, WT, NbrScale=5
;               tvscl, WT.coef[*,*,0] ; plot the wavelet transform (f = 0..11) of the modulus
;         
; HISTORY:
;	Written:  Jerome Bobin, May 2008
;-
;-----------------------------------------------------------------


;######################  ShapeQU ##############################

pro ShapeQU,Q,U,S1_QU,M_QU

M_QU = sqrt(Q^2. + U^2.)
S1_QU =  complex(Q/(M_QU),U/(M_QU))

end

;###############################################################

pro extract_scale,Xin,J,T1,T2,T3

nx = size(Xin);
N = nx(1)

Jmax = floor(alog(double(N))/alog(2.));

Nl = 2^(Jmax-J+1)-1;
Nd = 2^(Jmax-J);

T1 = Xin[Nd:Nl,Nd:Nl]
T2 = Xin[0:Nd-1,Nd:Nl]
T3 = Xin[Nd:Nl,0:Nd-1]

end


;###############################################################

pro ComplexMult1,Cp1,Cp2,Cp

Th1 = imaginary(alog(Cp1))
Th2 = imaginary(alog(Cp2))

Cp = complex(cos(Th1+Th2),sin(Th1+Th2))

end

;######################  gene_vec ##############################

function GeneVec,Vmin,Vmax,Np

Vec = (Vmax - Vmin)*indgen(Np)/(Np-1) + Vmin;

return,Vec

end

;###############################################################

pro Lift1Sphere,V,FiltCoeff,ManiWC,WaveCoeff,P_Ce,C_o

;--- Sur la 1-Sphere, avec des complexes

L = N_ELEMENTS(V);
K = floor(double(L)/2.);
FiltCoeff = 0.*DCOMPLEXARR(K);
WaveCoeff = 0.*dblarr(K);

;--- Splitting

Eind = 2*indgen(K);
Oind = Eind + 1;

C_e = V[Eind];
C_o = V[Oind];

;--- Prediction

P_Ce = C_e;

Wco = 0.*dblarr(K);

Wco = imaginary(alog(P_Ce^(-1)*C_o))

WaveCoeff = Wco;

;---Update

FiltCoeff = C_e*complex(cos(0.5*Wco),sin(0.5*Wco))

ManiWC = complex(cos(Wco),sin(Wco))

end

;###############################################################

pro FoWT_Ang,Xin,J,wc

;--- Xin contains the manifold-based values for each pixel


;--- First Normalize the input :

nx = size(Xin);
N = nx(1)
wc = Xin;

;--- WT for each column/line:

;for each scale

Jmax = floor(alog(double(N))/alog(2.));

for s=0,J-1 do begin

for ll=0,2^(Jmax-s+1)-1 do begin

            ;-- Extract the line
                                    
            indi = indgen(2^(Jmax-s+1));
            
            Vline = wc[ll,indi];
            
            Lift1Sphere,Vline,FiltCoeff,ManiWC,WaveCoeff;
                        
            indi1 = indgen(2^(Jmax-s));
            
            wc[ll,indi1] = FiltCoeff;
                        
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
            
            wc[ll,indi2] = ManiWC;

           ; print,imaginary(alog(ManiWC))

endfor

for mm = 0,2^(Jmax-s+1)-1 do begin
        
            ;-- Extract the column
            
            indi = indgen(2^(Jmax-s+1));
            
            Vcol = wc[indi,mm];
            
            Lift1Sphere,Vcol,FiltCoeff,ManiWC,WaveCoeff;
            
            indi1 = indgen(2^(Jmax-s));
            
            wc[indi1,mm] = FiltCoeff;
            
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
   
            wc[indi2,mm] = ManiWC;
            
endfor

endfor

Trans = wc

end

;###############################################################

pro Lift1Mod,V,FiltCoeff,ManiWC,WaveCoeff,P_Ce,C_o

;--- Sur la 1-Sphere, avec des complexes

L = N_ELEMENTS(V);
K = floor(double(L)/2.);
FiltCoeff = 0.*DBLARR(K);
WaveCoeff = 0.*dblarr(K);

;--- Splitting

Eind = 2*indgen(K);
Oind = Eind + 1;

C_e = V[Eind];
C_o = V[Oind];

;--- Prediction

P_Ce = C_e

Wco = C_o - P_Ce

WaveCoeff = Wco;

;---Update

FiltCoeff = C_e + Wco/2.

ManiWC = Wco

end

;###############################################################

pro FoWT_Mod,Xin,J,wc

;--- Xin contains the manifold-based values for each pixel


;--- First Normalize the input :

nx = size(Xin);
N = nx(1)
wc = Xin;

;--- WT for each column/line:

;for each scale

Jmax = floor(alog(double(N))/alog(2.));

for s=0,J-1 do begin

for ll=0,2^(Jmax-s+1)-1 do begin

            ;-- Extract the line
                                    
            indi = indgen(2^(Jmax-s+1));
            
            Vline = wc[ll,indi];
            
            Lift1Mod,Vline,FiltCoeff,ManiWC,WaveCoeff;
                        
            indi1 = indgen(2^(Jmax-s));
            
            wc[ll,indi1] = FiltCoeff;
                        
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
            
            wc[ll,indi2] = ManiWC;

           ; print,imaginary(alog(ManiWC))

endfor

for mm = 0,2^(Jmax-s+1)-1 do begin
        
            ;-- Extract the column
            
            indi = indgen(2^(Jmax-s+1));
            
            Vcol = wc[indi,mm];
            
            Lift1Mod,Vcol,FiltCoeff,ManiWC,WaveCoeff;
            
            indi1 = indgen(2^(Jmax-s));
            
            wc[indi1,mm] = FiltCoeff;
            
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
   
            wc[indi2,mm] = ManiWC;
            
endfor

endfor

Trans = wc

end

;###############################################################

pro EraseOWTScale,Trans,Scale

	NbrScale = Trans.NbrScale

	N=double((size(Trans.ModCoeff))[1])

	for pp=0,NbrScale-Scale do begin ;--- for each scale

        	Nd = N/(2.^(pp+1))
        	Nf = N/(2.^pp)-1.
		
		    Trans.ModCoeff[Nd:Nf,Nd:Nf,*] = 0.
		    Trans.ModCoeff[Nd:Nf,0:Nd-1.,*] = 0.
		    Trans.ModCoeff[0:Nd-1.,Nd:Nf,*] = 0.
		    
		    Trans.AngCoeff[Nd:Nf,Nd:Nf,*] = 0.
		    Trans.AngCoeff[Nd:Nf,0:Nd-1.,*] = 0.
		    Trans.AngCoeff[0:Nd-1.,Nd:Nf,*] = 0.
				
    endfor

end

;###############################################################

pro CalcScaleRecOWT,Ima,Rec,NbrScale

mr_modphase_dwt_trans, Ima, Trans, NbrScale=NbrScale

nc = size(Ima)
Rec= 0.*dblarr(nc[1],nc[2],2,NbrScale-1.)
help,Rec

for ll=0,NbrScale-2. do begin

	TransTemp = Trans	
	EraseOWTScale,TransTemp,ll
	mr_modphase_dwt_rec,TransTemp,temp
	help,temp
	Rec[*,*,*,ll] = temp

endfor

end

;###############################################################

pro PlotScaleRec,Ima,Rec,len=len,dist=dist,wsize=wsize

nc = size(Rec)
p=max([nc[1],512.])
if keyword_set(wsize) then p = wsize

window,0,xsize=p,ysize=p
tv_vecf,Ima,len=len,dist=dist

for ll=0,nc[4]-1 do begin

	window,ll+1,xsize=p,ysize=p
	tv_vecf,Rec[*,*,*,ll],len=len,dist=dist

endfor

end

;###############################################################

pro mr_modphase_dwt_trans,Ima,Trans,NbrScale=NbrScale

if N_PARAMS() NE 2  then begin 
        print, 'CALLING SEQUENCE: mr_modphase_dwt_trans, Imag, Trans,  NbrScale=NbrScale '
        goto, DONE
        end

Ni = size(Ima)

if not keyword_set(NbrScale) then NbrScale = 4

ShapeQU,Ima[*,*,0],Ima[*,*,1],AngMap,ModMap

;--- Decomposing the modulus map

FoWT_Mod,ModMap,NbrScale,ModCoeff

;--- Decomposing the angle map

FoWT_Ang,AngMap,NbrScale,AngCoeff

;--- Declare the structure

Trans = {NbrScale : NbrScale, ModCoeff : ModCoeff, AngCoeff :AngCoeff }

DONE:

end

;###############################################################
