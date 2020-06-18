;+
; NAME:
;         MR_MODPHASE_UWT_TRANS
;
; PURPOSE:
;	   Computes a multiresolution transform of a vector field.
;
; CALLING:
;
;      MR_MODPHASE_UWT_TRANS, Ima, Trans, NbrScale=NbrScale   
;
; INPUTS:
;     Ima -- [*,*,2] 3D IDL array: vector field we want transform
;    
; OUTPUTS:
;     Transf -- IDL Structure :
;			NbrScale -- scalar - number of scales	
;			ModCoeff -- 3D IDL Array : multiresolution transform of the modulus
;			AngCoeff -- 3D IDL Array : multiresolution transform of the angles 
;
; KEYWORDS:
;      NbrScale -- scalar : number of scales - default : 4
;
; EXAMPLE:
;       Compute the multiresolution of a vector field I_QU :
;              MR_MODPHASE_UWT_TRANS, I_QU, CoeffStruct
;
; HISTORY:
;	Written: Jerome BOBIN 04/2008
;-

;######################  ShapeQU ##############################

pro ShapeQU,Q,U,S1_QU,M_QU

M_QU = sqrt(Q^2. + U^2.)
S1_QU =  complex(Q/(M_QU+1e-9),U/(M_QU+1e-9))

end


;######################  gene_vec ##############################

function GeneVec,Vmin,Vmax,Np

Vec = (Vmax - Vmin)*indgen(Np)/(Np-1) + Vmin;

return,Vec

end

;######################  Complex on S1 ##############################

function expola,Vin

ExpV = complex(cos(Vin),sin(Vin))

return,ExpV

end


;######################  UPDATE COEFF ##############################

function UpCoeff_S1,Coarse,NCoarse

Mc = abs(Coarse)
MNc = abs(NCoarse)

DMod = MC - MNc

VTheta = imaginary(alog(NCoarse*Coarse^(-1)))

;W = expola(VTheta);  %--- pas necessaire aucun sens sur la manifold

W = [[[DMod]],[[VTheta]]]

help,W

return,W

end

;######################  CalcInvTheta ##############################

function CalcDiffTheta,X1,X0

P1 = X1/abs(X1 + 1e-12)
P0 = X0/abs(X0+1e-12)

V = imaginary(alog(P1*P0^(-1)))

return,V

end

;######################  CONVOL 1D ##############################

function convol1d_S1,X,F,J,MirrBord=MirrBord,ConstBord=ConstBord

N = N_ELEMENTS(X);
Nf = N_ELEMENTS(F);
Nf2 = double(floor(Nf/2));

Xcv = X;
Xcv[*] = 0.

for ll=0,N-1 do begin

        Nl = 2.^(J-1);
                            
        if (ll - Nl*Nf2) lt 0 then begin
                
                q = double(floor(ll/Nl));
                Nfc_m = max([0.,q]);
                M_Nfc = Nf2;
                                
                Np = M_Nfc + Nfc_m + 1
                indi = Nl*GeneVec(-Nfc_m,M_Nfc,Np) + ll;
                                
                indif = GeneVec(Nf2-Nfc_m,Nf-1,Np)
                
                ;--- Going to the tangent space :
                
                ;VXl = 0.*DCOMPLEXARR(Np)
                VXl = 0.*DCOMPLEXARR(2*Nf2+1)

                VXl[*] = X(ll);

                VxOthers = 0.*DCOMPLEXARR(2*Nf2+1);
                
                if keyword_set(ConstBord) then begin
                            
                            VxOthers[0:2*Nf2-Np] = X[indi[0]]
                            VxOthers[2*Nf2+1-Np:2*Nf2] = X[indi]
                                                        
                endif
                
                 if keyword_set(MirrBord) then begin
                            
                            indox = indi(bsort(indi,/REVERSE))
                            VxOthers[0:2*Nf2-Np] = X[indox[0:2*Nf2-Np]]
                            VxOthers[2*Nf2+1-Np:2*Nf2] = X[indi]
                            
                endif
                
                u1= real_part(alog(VxOthers))
                u2 = imaginary(alog(VxOthers))

               u1 = total(F*exp(u1))
               u2 = total(F*u2)
               
               Temp = complex(alog(u1),u2)           
                
                Xcv[ll] = exp(Temp)
                                             
        endif
        
        

        if ((ll + Nl*Nf2) lt N) then begin
            if (ll - Nl*Nf2) ge 0 then begin
                                            
                q = double(floor(ll/Nl));
                Nfc_m = Nf2;
                M_Nfc = Nf2;
                                
                Np = M_Nfc + Nfc_m + 1
                indi = Nl*GeneVec(-Nfc_m,M_Nfc,Np) + ll;
                                                           
                indif = GeneVec(Nf2-Nfc_m,Nf-1,Np)
                                                
                ;--- Going to the tangent space :
                
                VxOthers = X[indi]

                u1= real_part(alog(VxOthers))
                u2 = imaginary(alog(VxOthers))

               u1 = total(F*exp(u1))
               u2 = total(F*u2)
               
               Temp = complex(alog(u1),u2)           
                
                Xcv[ll] = exp(Temp)
                
             endif
        endif
                        
        if (ll + Nl*Nf2) ge N then begin
        
                Nfc_m = Nf2;
                
                q = floor((N-1-ll)/Nl);
                M_Nfc = min([Nf2,q]);
                Nfc_m = -Nf2
                Np = M_Nfc - Nfc_m +1
                
                indi = Nl*GeneVec(-Nf2,M_Nfc,Np) + ll;
                                
                indif = GeneVec(0,Np-1,Np)
                                                
                ;--- Going to the tangent space :

                VXl = 0.*DCOMPLEXARR(Np)

                VXl[*] = X(ll);
                
                VxOthers = 0.*DCOMPLEXARR(2*Nf2+1);
                
                if keyword_set(ConstBord) then begin
                            
                            VxOthers[Np:2*Nf2] = X[indi[Np-1]]
                            VxOthers[0:Np-1] = X[indi]
                                                        
                endif
                
                 if keyword_set(MirrBord) then begin
                            
                            indox = indi(bsort(indi,/REVERSE))
                            VxOthers[Np:2*Nf2] = X[indox[0:2*Nf2-Np]]
                            VxOthers[0:Np-1] = X[indi]
                            
                endif

                u1= real_part(alog(VxOthers))
                u2 = imaginary(alog(VxOthers))

               u1 = total(F*exp(u1))
               u2 = total(F*u2)
               
               Temp = complex(alog(u1),u2)           
                
                Xcv[ll] = exp(Temp) 

        endif

end

return,Xcv

end


 
;######################  FORWARD TRANSFORM ##############################

pro FoUWT_Angle,Xin,J,wc

nx = size(Xin);
N = nx(1);
M = nx(2);

wc = 0.*dblarr(N,M,J);

Coarse = Xin;
NCoarse = Coarse;

h = 1/16.*[1,4,6,4,1];

for ll=0,J-2 do begin

        ;--- For each scale - compute the coarse image
                
        for cc = 0,N-1 do begin
               
               Vo = convol1d_S1(NCoarse[cc,*],h,ll,/ConstBord);
               NCoarse[cc,*] = Vo;
                              
        end
                
        for cc = 0,M-1 do begin
        
               V = convol1d_S1(transpose(NCoarse[*,cc]),h,ll,/ConstBord);
               NCoarse[*,cc] = transpose(V);
        
        end
        
        ;--- And the wavelet coeffs

        wc[*,*,ll] = imaginary(alog(Coarse) - alog(NCoarse));
        Coarse = NCoarse;

endfor

wc[*,*,J-1] = imaginary(alog(NCoarse));

end

;######################  CONVOL 1D MODULUS ##############################

function convol1d_Mod,X,F,J,MirrBord=MirrBord,ConstBord=ConstBord

N = N_ELEMENTS(X);
Nf = N_ELEMENTS(F);
Nf2 = double(floor(Nf/2));

Xcv = X;
Xcv[*] = 0.

for ll=0,N-1 do begin

        Nl = 2.^(J-1);
                            
        if (ll - Nl*Nf2) lt 0 then begin
                
                q = double(floor(ll/Nl));
                Nfc_m = max([0.,q]);
                M_Nfc = Nf2;
                                
                Np = M_Nfc + Nfc_m + 1
                indi = Nl*GeneVec(-Nfc_m,M_Nfc,Np) + ll;
                                
                indif = GeneVec(Nf2-Nfc_m,Nf-1,Np)
                
                ;--- Going to the tangent space :
                
                ;VXl = 0.*DBLARR(Np)
                VXl = 0.*DBLARR(2*Nf2+1)

                VXl[*] = X(ll);

                VxOthers = 0.*DBLARR(2*Nf2+1);
                
                if keyword_set(ConstBord) then begin
                            
                            VxOthers[0:2*Nf2-Np] = X[indi[0]]
                            VxOthers[2*Nf2+1-Np:2*Nf2] = X[indi]
                                                        
                endif
                
                 if keyword_set(MirrBord) then begin
                            
                            indox = indi(bsort(indi,/REVERSE))
                            VxOthers[0:2*Nf2-Np] = X[indox[0:2*Nf2-Np]]
                            VxOthers[2*Nf2+1-Np:2*Nf2] = X[indi]
                            
                endif
                
               Xcv[ll] = total(F*VxOthers)
                                             
        endif
        
        

        if ((ll + Nl*Nf2) lt N) then begin
            if (ll - Nl*Nf2) ge 0 then begin
                                            
                q = double(floor(ll/Nl));
                Nfc_m = Nf2;
                M_Nfc = Nf2;
                                
                Np = M_Nfc + Nfc_m + 1
                indi = Nl*GeneVec(-Nfc_m,M_Nfc,Np) + ll;
                                                           
                indif = GeneVec(Nf2-Nfc_m,Nf-1,Np)
                                                
                ;--- Going to the tangent space :
                
                VxOthers = X[indi]

               Xcv[ll] = total(F*VxOthers)
                
                 
                
             endif
        endif
                        
        if (ll + Nl*Nf2) ge N then begin
        
                Nfc_m = Nf2;
                
                q = floor((N-1-ll)/Nl);
                M_Nfc = min([Nf2,q]);
                Nfc_m = -Nf2
                Np = M_Nfc - Nfc_m +1
                
                indi = Nl*GeneVec(-Nf2,M_Nfc,Np) + ll;
                                
                indif = GeneVec(0,Np-1,Np)
                                                
                ;--- Going to the tangent space :

                VXl = 0.*DCOMPLEXARR(Np)

                VXl[*] = X(ll);
                
                VxOthers = 0.*DCOMPLEXARR(2*Nf2+1);
                
                if keyword_set(ConstBord) then begin
                            
                            VxOthers[Np:2*Nf2] = X[indi[Np-1]]
                            VxOthers[0:Np-1] = X[indi]
                                                        
                endif
                
                 if keyword_set(MirrBord) then begin
                            
                            indox = indi(bsort(indi,/REVERSE))
                            VxOthers[Np:2*Nf2] = X[indox[0:2*Nf2-Np]]
                            VxOthers[0:Np-1] = X[indi]
                            
                endif

                Xcv[ll] = total(F*VxOthers)

        endif

end

return,Xcv

end

 

;######################  FORWARD TRANSFORM ##############################

pro FoUWT_Mod,Xin,J,wc

nx = size(Xin);
N = nx(1);
M = nx(2);

wc = 0.*DBLARR(N,M,J);

Coarse = Xin;
NCoarse = Coarse;

h = 1/16.*[1,4,6,4,1];

for ll=0,J-2 do begin

        ;--- For each scale - compute the coarse image
                
        for cc = 0,N-1 do begin
               
               Vo = convol1d_Mod(NCoarse[cc,*],h,ll,/ConstBord);
               NCoarse[cc,*] = Vo;
                              
        end
                
        for cc = 0,M-1 do begin
        
               V = convol1d_Mod(transpose(NCoarse[*,cc]),h,ll,/ConstBord);
               NCoarse[*,cc] = transpose(V);
        
        end
        
        ;--- And the wavelet coeffs

        wc[*,*,ll] = Coarse - NCoarse;
        Coarse = NCoarse;

endfor

wc[*,*,J-1] = NCoarse;

end

;###############################################################

pro EraseWTScale,Trans,Scale

Trans.ModCoeff[*,*,0:Scale] = 0.
Trans.AngCoeff[*,*,0:Scale] = 0.

end

 
 
;######################### MRP_WTTRANS #################################

pro mr_modphase_uwt_trans, Ima, Trans, NbrScale=NbrScale

if N_PARAMS() NE 2  then begin 
        print, 'CALLING SEQUENCE: mr_modphase_uwt_trans, Imag, Trans,  NbrScale=NbrScale '
        goto, DONE
        end

Ni = size(Ima)

if not keyword_set(NbrScale) then NbrScale = 4

ShapeQU,Ima[*,*,0],Ima[*,*,1],AngMap,ModMap

;--- Decomposing the modulus map

FoUWT_Mod,ModMap,NbrScale,ModCoeff

;--- Decomposing the angle map

FoUWT_Angle,AngMap,NbrScale,AngCoeff

;--- Declare the structure

Trans = {NbrScale : NbrScale, ModCoeff : ModCoeff, AngCoeff :AngCoeff }

DONE:

end