;+
; NAME:
;        STAR2D
;
; PURPOSE:
;       Computes the starlet transform of an image (i.e. undecimated isotropic wavelet transform). 
;       The output is a 3D IDL array. If the keyword Gen2 is set, then it is
;       the 2nd generation starlet transform which is computed:  i.e. g = Id - h*h 
;       instead of g = Id - h.
;       
; CALLING:
;      DataTransf = STAR2D(Imag, Nscale=Nscale, Gen2=Gen2)
;       
; INPUTS:
;     Imag -- 2D IDL array: image we want to transform
;
; OUTPUTS:
;     DataTransf -- 3D IDL array: Wavelet Transform
;                                 DataTransf(*,*,i) = ith band of the 
;                                                    wavelet transform
;  
; KEYWORDS:
;      Nscale -- int: Number of scales. Default is 4.
;                     DataTransf is a cube of Nscale planes: 
;                                      DataTransf(*,*,0: Nscale-1)
;			There is no test of the validity of the number of scales with respect
;			to the size of the input images. 
; 
; EXAMPLE:
;       Compute the starlet transform of an image I with default options
;       (i.e. a trou algorithm with 4 scales).  
;               Output = STAR2D(I)
;       Reconstruction can be done by simple co-addition of all frames:
;               Rec = total(output, 3)
;
; REFERENCES:
;    [1] J.L. Starck and F. Murtagh, 
;    "Image Restoration with Noise Suppression 
;    Using the Wavelet Transform",
;    Astronomy and Astrophysics, 288, pp-343-348, 1994.
;    
;
;    For the modified STARLET transform:
;    [2] J.-L. Starck, J. Fadili and F. Murtagh, "The Undecimated Wavelet Decomposition 
;        and its Reconstruction", IEEE Transaction on Image Processing,  16,  2, pp 297--309, 2007.
; AUTHOR:
;    J.L. Starck
;    Service d'Astrophysique, Centre d'Etudes de SACLAY,
;    Orme des Merisiers, 91191 GIF-Sur-YVETTE CEDEX, France 
;    Email: jstarck@cea.fr        Tel:  +33 (0) 169 085 764
;    http://jstarck.free.fr       Fax:  +33 (0) 169 086 577
;-
;-------------------------------------------------------------------------------


function test_ind_1D, ind, N
	;this function implements mirror like limit conditions on the edges of the input signal
	;ATTENTION : the output may still be out of range ie not in [0, N-1]
	;refinements may be necessary, although not meaningful in practice
	
	ret = ind
	if ind LT 0 then ret = -ind $
	else if ind GE N then ret = 2*N-ind-2
	
	return, ret
	
end

;-------------------------------------------------------------------------------


pro b3spline_star2d, Im_in, Im_Out,  Step
   
   Nl = (size(Im_in))(1)
   Nc = (size(Im_in))(2)
   Im_Out = dblarr(Nl,Nc)

   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Buffs = fltarr(Nl,Nc)

   for i = 0, Nl-1 do begin
    for j = 0, Nc-1 do begin
      jm = test_ind(j-Step, Nc)
      jp = test_ind(j+Step, Nc)
      jm2 = test_ind(j-2*Step, Nc)
      jp2 = test_ind(j+2*Step, Nc)
      Buffs (i,j) =$ 
		C3 * Im_in(i,j) + C2 * (Im_in(i,jm) + Im_in(i,jp))$
		+ C1 * (Im_in(i,jm2) + Im_in(i,jp2)) 
	endfor
   endfor

   for i = 0, Nl-1 do begin
    for j = 0, Nc-1 do begin
      im = test_ind(i-Step, Nl)
      ip = test_ind(i+Step, Nl)
      im2 = test_ind(i-2*Step, Nl)
      ip2 = test_ind(i+2*Step, Nl)
      Im_Out (i,j) =$
			C3 * Buffs(i,j) + C2 * (Buffs(im,j) + Buffs(ip,j))$
			+ C1 * (Buffs(im2,j) + Buffs(ip2,j))
    endfor
   endfor
end



; the IDL routine convol is faster but does not allow MIRROR border which is generally the best way to handle border
pro fast_b3spline_star2d, Im_in, Im_Out,  Step
; /edge_wrap
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   KSize = 4*Step+1
   KS2 = KSize/2
   Kernel = fltarr(KSize)
   Kernel[0] = C1
   Kernel[KSize-1] = C1
   Kernel[KS2+Step] = C2
   Kernel[KS2-Step] = C2
   Kernel[KS2]=C3
   z = convol(Im_in,Kernel, /EDGE_TRUNCATE, /CENTER) 
   kernelY = transpose(Kernel)
   Im_Out = convol(z, kernelY, /EDGE_TRUNCATE, /CENTER)
end

;-------------------------------------------------------------------------------
 
function star2d, Ima, Nscale=Nscale, Gen2=Gen2 

        WT = -1
        Fast=1
	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: WT = star2d(Ima, Nscale=Nscale, Gen2=Gen2)'
        goto, DONE
	end

	if not keyword_set(Nscale) then Nscale =4

	if Nscale LT 2 then  Nscale =4
	
	vsize = size(Ima)
	if vsize[0] NE 2 then begin
        print, 'Error: First parameter is not an image ...'
        print, 'CALLING SEQUENCE:  WT = star2d(Ima, Nscale=Nscale, Gen2=Gen2)'
        goto, DONE
	end
	
	;-------------------------------------------------------------------------------
	;recovering a few parameter values, initializing the transform loop ...
	Nx = vsize[1]
	Ny = vsize[2]
	Nz = Nscale
	NStep=Nscale-1
  	WT = fltarr(Nx,Ny,Nz)
	Im_in = Ima
	Im_Out = fltarr(Nx,Ny)
	Step_trou = 1

	for i=0, NStep-1 do begin
		;B3 Spline smoothing
		if keyword_set(fast) then fast_b3spline_star2d, Im_in, Im_Out, Step_trou $
		else b3spline_star2d, Im_in, Im_Out, Step_trou 
                if keyword_set(Gen2) then begin
		   Im_Aux = Im_Out
		   if keyword_set(fast) then  fast_b3spline_star2d, Im_Out, Im_Aux, Step_trou $
		   else b3spline_star2d, Im_Out, Im_Aux, Step_trou 
		   WT(*,*,i) = Im_in - Im_Aux
		end else begin
		; Wavelet coefficient calculations
		  WT(*,*,i) = Im_in - Im_Out
		end
		Im_in = Im_Out

		; New distance between two pixels
		if i NE NStep-1 then Step_trou = Step_trou * 2
        endfor

	; Smooth array
	i =  NStep  
	WT(*,*,i) = Im_Out
	
DONE:
    return, WT
    
end

;-------------------------------------------------------------------------------
pro softth, Data, Lambda
ind = where ( Data GT 0, c)
if c GT 0 then begin
   Pos = Data[ind] - Lambda
   indN = where (Pos LT 0, c)
   if c GT 0 then Pos[indN] = 0
   Data[ind] = Pos
end 
ind = where ( Data LT 0, c)
if c GT 0 then begin
   Neg = Data[ind] + Lambda
   indP = where (Neg GT 0, c)
   if c GT 0 then Neg[indP] = 0
   Data[ind] = Neg
end 
end


function pstar2d, Ima, Nscale=Nscale, Niter=Niter, Rec=Rec, soft=soft, SigmaNoise=SigmaNoise, Nsigma=Nsigma, pos=pos, den=den, mask=Mask, TabNorm=TabNorm, $
 background=background, out_wtbackground=out_wtbackground, out_background=out_background, bgr_FirstScale=bgr_FirstScale

Gen2=1
Resi  = Ima

if keyword_set(background) then out_background = Ima*0. else out_background=0.
if not keyword_set(bgr_FirstScale) then bgr_FirstScale = Nscale-1
w  = star2d(Resi, Nscale=Nscale, Gen2=Gen2)

if keyword_set(Mask) then begin
  IndMask = where(Mask NE 0)
  Scale = w[*,*,0]
  soft=1
  if not keyword_set(SigmaNoise) then SigmaNoise = MAD(Scale[IndMask])
end else if not keyword_set(SigmaNoise) then SigmaNoise = MAD(w[*,*,0])

if not keyword_set(Niter) then Niter = 10
if not keyword_set(Nsigma) then Nsigma = 3

print, "SigmaNoise = ", SigmaNoise

MW = max(w) / 2.

if not keyword_set(TabNorm) then begin
  dirac = Ima
  dirac [*] = 0
  vs = size(dirac)
  Nx = vs[1]
  Ny = vs[2]
  dirac[Nx/2, Ny/2] = 1
  w  = star2d(dirac, Nscale=Nscale, Gen2=Gen2)
  TabNorm = fltarr(Nscale)
  for j=0, Nscale-1 do TabNorm[j] = sqrt( total(w[*,*,j]^2) )
end

for i=0,Niter+10 do begin
   Lambda = MW * ( 1. - (i+1.) / Niter)
   if Lambda LT 0 then Lambda = 0
   w  = star2d(Resi, Nscale=Nscale, Gen2=Gen2)
   if keyword_set(background) then  begin
      LastScaleResi = w[*,*,Nscale-1]
      w[*,*,Nscale-1] = 0
   end
   
   Wr = W
   if i EQ 0 then WT = w*0

   if keyword_set(Den) then begin
   for j=0, Nscale-2 do begin
     Nsig= Nsigma
     if j EQ 0 then Nsig= Nsigma+1
     Scale = WT[*,*,j]
     ScaleR = Wr[*,*,j]

     ind = where( ABS(Scale+ScaleR)  LT SigmaNoise*TabNorm[j]* Nsig, c)
     if c GT 0 then ScaleR[ind] = 0
     Wr[*,*,j] = ScaleR
   end
   end
   WT = WT + Wr
    
   if keyword_set(soft) then begin
      for j=0, Nscale-1 do begin
          Scale = WT[*,*,j]
          softth, Scale,  Lambda* TabNorm[j]
          WT[*,*,j] = Scale
      end
  end
  
  if keyword_set(background) then begin
     out_background = out_background + LastScaleResi
     softth, out_background,  Lambda* TabNorm[NScale-1]
  end
  
  if keyword_set(pos) then begin
     ind = where(WT LT 0, c)
     if c GT 0 then WT[ind] = 0
  end
  
   Rec = istar2d(WT, Gen2=Gen2)
   Resi = Ima - Rec
   if keyword_set(Mask) then Resi = Resi * Mask
   print, "Iter " , i+1, Lambda , ", Resi =  " ,  sigma(Resi)
end

info, wt
load, ima-rec
DONE:

 return, WT
 
end


;============================================================


function fistar2d, Ima, Nscale=Nscale, Niter=Niter, Rec=Rec,  SigmaNoise=SigmaNoise, Nsigma=Nsigma, pos=pos, Lambda=Lambda, TabSigmaNoise=TabSigmaNoise
Gen2=1 
Resi  = Ima
w  = star2d(Resi, Nscale=Nscale, Gen2=Gen2)
WT = w
WT[*]  = 0.
if not keyword_set(SigmaNoise) then SigmaNoise = MAD(w[*,*,0])
if not keyword_set(Niter) then Niter = 10
if not keyword_set(Nsigma) then Nsigma = 3

print, "SigmaNoise = ", SigmaNoise

   MW = max(w) / 2.
   dirac = Ima
  dirac [*] = 0
   vs = size(dirac)
   Nx = vs[1]
   Ny = vs[2]
   dirac[Nx/2, Ny/2] = 1
   w  = star2d(dirac, Nscale=Nscale, Gen2=Gen2)
   TabNorm = fltarr(Nscale)
   for j=0, Nscale-1 do TabNorm[j] = sqrt( total(w[*,*,j]^2) )
if not keyword_set(TabSigmaNoise) then begin
   TabSigmaNoise= fltarr(Nscale)
   TabSigmaNoise = TabNorm * SigmaNoise
end
AlphaPrev = 0.
T = 1.;
TPrev = 0.

for i=0,Niter do begin
   if not keyword_set(Lambda) then Lambda = MW * ( 1. - (i+1.) / Niter)
   if Lambda LT 0 then Lambda = 0
 
   TPrev = T
   T = (1. + sqrt(4.*T^2)) / 2.
   X =  WT + (TPrev-1.) / T  * (WT - AlphaPrev)

   ; Gradient Calculation
   if not keyword_set(FISTA) then X = WT
   Rec = istar2d(X, Gen2=Gen2) 
   Resi = Ima - Rec
   w  = star2d(Resi, Nscale=Nscale, Gen2=Gen2)
   X = X + w

  if i EQ Niter then begin
     for j=0, Nscale-2 do begin
        Nsig= Nsigma
        if j EQ 0 then Nsig= Nsigma+1
        Scale = X[*,*,j] 
        ind = where( ABS(Scale)  LT TabSigmaNoise[j] * Nsig, c)
        if c GT 0 then Scale[ind] = 0
        X[*,*,j] = Scale
    end
  end else begin
      for j=0, Nscale-2 do begin
          Scale = WT[*,*,j]
          softth, Scale,  Lambda* TabNorm[j]
          if keyword_set(FISTA) then   X[*,*,j] = Scale $
          else WT[*,*,j] = Scale
      end
 end
   AlphaPrev = WT
   WT = X
   
  if keyword_set(pos) then begin
     ind = where(WT LT 0, c)
     if c GT 0 then WT[ind] = 0
  end
  
       print, "Iter " , i+1, Lambda , ", Resi =  " ,  sigma(Resi)
end


 Rec = istar2d(WT, Gen2=Gen2)  
info, wt
load, ima-rec
DONE:

 return, WT
 
end




