;+
; NAME:
;        MR_MODPHASE_UWT_REC
;
; PURPOSE:
;	Reconstruct an vector field of an image from its undecimated wavelet transform.
;
; CALLING:
;      MR_MODPHASE_UWT_REC, Trans, Rec
;
; INPUTS:
;       Trans : IDL structure; Wavelet transform structure (see MR_MODPHASE_UWT_REC) ;     
;    
; OUTPUTS:
;      Imag -- IDL array of a vector field image fltarr[*,*,2] : Output image be reconstructed 
;
; EXTERNAL CALLS:
;
; EXAMPLE:
;       Compute the undecimated wavelet transform of a vector field I with five scales and reconstrcution
;               mr_modphase_uwt_trans, Imag, WT, NbrScale=5
;               mr_modphase_uwt_rec, WT, RecIma
;               tvscl, RecIma.coef[*,*,0] ; plot the first reconstructed component   
;               tvscl, RecIma.coef[*,*,1] ; plot the second reconstructed component   
;         
; HISTORY:
;	Written:  Jerome Bobin, May 2008
;-
;-----------------------------------------------------------------



;###############################################################


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


;######################  INVERSE TRANSFORM ##############################

pro IUWT_Ang,wc,Xout

nw = size(wc);

J = nw(3);

wcr = 0.*wc

wc = complex(wcr,wc)

Xout = total(wc,3)

Xout = exp(Xout)

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


;######################  INVERSE TRANSFORM ##############################

pro IUWT_Mod,wc,Xout

nw = size(wc);

J = nw(3);

Xout = total(wc,3)

end

 ;###############################################################

pro EraseWTScale,Trans,Scale

Trans.ModCoeff[*,*,0:Scale] = 0.
Trans.AngCoeff[*,*,0:Scale] = 0.

end

;###############################################################

pro CalcScaleRec,Ima,Rec,NbrScale

mrp_wttrans, Ima, Trans, NbrScale=NbrScale


nc = size(Ima)
Rec= 0.*dblarr(nc[1],nc[2],2,NbrScale-1.)
help,Rec

for ll=0,NbrScale-2. do begin

	TransTemp = Trans	
	EraseWTScale,TransTemp,ll
	mrp_wtrec,TransTemp,temp
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

;######################### MRP_WTREC #################################

pro mr_modphase_uwt_rec, Trans,Rec

NbrScale = Trans.NbrScale

;--- Reconstruct the modulus map

IUWT_Mod,Trans.ModCoeff,ModMap

;--- Reconstruct the angle map

IUWT_Ang,Trans.AngCoeff,AngMap

Xout =ModMap * AngMap

Rec = [[[real_part(Xout)]],[[imaginary(Xout)]]]

end

 