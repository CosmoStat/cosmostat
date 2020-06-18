;+
; NAME:
;        MR_MODPHASE_DWT_REC
;
; PURPOSE:
;	Reconstruct an vector field of an image from its decimated wavelet transform.
;
; CALLING:
;      MR_MODPHASE_DWT_REC,Trans, Rec
;
; INPUTS:
;       Trans : IDL structure; Wavelet transform structure (see MR_MODPHASE_DWT_REC) ;     
;    
; OUTPUTS:
;      Imag -- IDL array of a vector field image fltarr[*,*,2] : Output image be reconstructed 
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the orthogonal wavelet transform of a vector field I with five scales and reconstrcution
;               mr_modphase_dwt_trans, Imag, WT, NbrScale=5
;               mr_modphase_dwt_rec, WT, RecIma
;               tvscl, RecIma.coef[*,*,0] ; plot the first reconstructed component   
;               tvscl, RecIma.coef[*,*,1] ; plot the second reconstructed component   
;         
; HISTORY:
;	Written:  Jerome Bobin, May 2008
;-
;-----------------------------------------------------------------


;###############################################################

pro ILift1Sphere,V,FiltCoeff,ManiWC

;--- Sur la 1-Sphere, avec des complexes
K = N_ELEMENTS(FiltCoeff);
L = 2.*K
V = 0.*DCOMPLEXARR(L);

;--- Splitting

Eind = 2*indgen(K);
Oind = Eind + 1;

C_e = 0.*DCOMPLEXARR(K);
C_o = 0.*DCOMPLEXARR(K);

;---Update

Uwc = 0.*dblarr(K);
Wco = 0.*dblarr(K)

Wco = imaginary(alog(ManiWC));

C_e = FiltCoeff*complex(cos(-Wco/2.),sin(-Wco/2.));

;--- Prediction

P_Ce = C_e;

C_o = P_Ce*complex(cos(Wco),sin(Wco)); %--- A faire dans le plan tangent

V[Eind] = C_e;
V[Oind] = C_o;


end

;###############################################################

pro IWT_Ang,Trans,NbScale,Rec

;--- Xin contains the manifold-based values for each pixel

wc = Trans
J = NbScale

;--- First Normalize the input :

Xout = wc
nx = size(wc)
N = nx(1)

;--- WT for each column/line:

;for each scale

Jmax = floor(alog(double(N))/alog(2.));

for qq=0,J-1 do begin

s = J-1 - qq


for ll=0,2^(Jmax-s+1)-1 do begin

            ;-- Extract the line
                                    
            indi = indgen(2^(Jmax-s+1));
            
            indi1 = indgen(2^(Jmax-s));
            
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
            
            FiltCoeff = Xout[ll,indi1];
            
            ManiWC = Xout[ll,indi2];
                        
            ILift1Sphere,Vline,FiltCoeff,ManiWC;
                        
            Xout[ll,indi] = Vline  

endfor

for mm = 0,2^(Jmax-s+1)-1 do begin
        
            ;-- Extract the column
            
            indi = indgen(2^(Jmax-s+1));
            
            indi1 = indgen(2^(Jmax-s));
            
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
            
            FiltCoeff = Xout[indi1,mm];
            
            ManiWC = Xout[indi2,mm];
                        
            ILift1Sphere,Vcol,FiltCoeff,ManiWC;
            
            Xout[indi,mm] = Vcol
            
endfor

endfor

Rec = Xout

end

;###############################################################

pro ILift1Mod,V,FiltCoeff,ManiWC

;--- Sur la 1-Sphere, avec des complexes
K = N_ELEMENTS(FiltCoeff);
L = 2.*K
V = 0.*DBLARR(L);

;--- Splitting

Eind = 2*indgen(K);
Oind = Eind + 1;

C_e = 0.*DBLARR(K);
C_o = 0.*DBLARR(K);

;---Update

Uwc = 0.*dblarr(K);
Wco = 0.*dblarr(K)

Wco = ManiWC;

C_e = FiltCoeff - Wco/2.;

;--- Prediction

P_Ce = C_e; %--- Pas top

C_o = P_Ce + Wco; %--- A faire dans le plan tangent


V[Eind] = C_e;
V[Oind] = C_o;


end

;###############################################################

pro IWT_Mod,Trans,NbScale,Rec

;--- Xin contains the manifold-based values for each pixel

wc = Trans
J = NbScale

;--- First Normalize the input :

Xout = wc
nx = size(wc)
N = nx(1)

;--- WT for each column/line:

;for each scale

Jmax = floor(alog(double(N))/alog(2.));

for qq=0,J-1 do begin

s = J-1 - qq


for ll=0,2^(Jmax-s+1)-1 do begin

            ;-- Extract the line
                                    
            indi = indgen(2^(Jmax-s+1));
            
            indi1 = indgen(2^(Jmax-s));
            
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
            
            FiltCoeff = Xout[ll,indi1];
            
            ManiWC = Xout[ll,indi2];
                        
            ILift1Mod,Vline,FiltCoeff,ManiWC;
                        
            Xout[ll,indi] = Vline  

endfor

for mm = 0,2^(Jmax-s+1)-1 do begin
        
            ;-- Extract the column
            
            indi = indgen(2^(Jmax-s+1));
            
            indi1 = indgen(2^(Jmax-s));
            
            indi2 = indgen(2^(Jmax-s)) + 2^(Jmax-s)
            
            FiltCoeff = Xout[indi1,mm];
            
            ManiWC = Xout[indi2,mm];
                        
            ILift1Mod,Vcol,FiltCoeff,ManiWC;
            
            Xout[indi,mm] = Vcol
            
endfor

endfor

Rec = Xout

end
 
;###############################################################

pro mr_modphase_dwt_rec, Trans, Rec

if N_PARAMS() NE 2 then begin 
        print, 'CALL SEQUENCE:  mr_modphase_dwt_rec, WT_Struct, result'
        goto, DONE
        end
        
NbrScale = Trans.NbrScale

;--- Reconstruct the modulus map

IWT_Mod,Trans.ModCoeff,NbrScale,ModMap

;--- Reconstruct the angle map

IWT_Ang,Trans.AngCoeff,NbrScale,AngMap

Xout =ModMap * AngMap

Rec = [[[real_part(Xout)]],[[imaginary(Xout)]]]

DONE:

end

;###############################################################
 