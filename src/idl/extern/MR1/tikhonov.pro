;+
; NAME: 
;       TIKHONOV
;
; PURPOSE: 
;       Deconvolve of an image by the TIKHONOV regularization.
;       This program does not shift the maximum of the PSF at 
;       center. The minimized functional is
;         J(Result) = || Imag - Psf*Result || + RegulParam || H*Result ||
;
;       The filter H used in this program is the Laplacien.
;       If the keyword NoRegul is set, then no regularization is done, 
;       and the algorithm becomes a simple one step gradient method.
;
; CALLING:
;
;       Result = TIKHONOV(Imag, Psf, residu=residu, niter=niter, 
;                         CvgParam=CvgParam, RegulParam=RegulParam, 
;                         NoRegul=NoRegul)
;
; INPUTS:
;       Imag: IDL 2D array: image to deconvolve
;       Psf:  IDL 2D array: Point Spread Function
;     
; INPUT KEYWORD: 
;       niter -- int: number of iterations. Default is 10.
;       CvgParam -- float: Convergence parameter. Default is 0.1.
;       RegulParam -- float: Regularization parameter.
;       NoRegul -- int: if set, no regularization is performed.
;
; OUTPUT KEYWORD: 
;       Residu: IDL 2D array: residual (= Imag - Psf*Result)
;
; EXAMPLE:
;           deconvolve an image with all default options  
;                Result = TIKHONOV(Imag, Psf)
;
;          same example, but impose the number of iterations to be 30
;                Result = TIKHONOV (Imag, Psf, niter=30)
;
;          deconvolution by the one step gradient method, without 
;          any regularization, with 30 iterations
;                Result = TIKHONOV (Imag, Psf, niter=30, /NoRegul)
;
; HISTORY:
;	Written: Jean-Luc Starck 1997.
;	June, 1997 File creation
;-
 
Function TIKHONOV, ima, psf, residu=residu, niter=niter, CvgParam=CvgParam, RegulParam=RegulParam, NoRegul=NoRegul
 
 if N_PARAMS() NE 2 then BEGIN
 print, 'USAGE: res = tikhonov(ima,psf,residu=residu, niter=niter, CvgParam=CvgParam, RegulParam=RegulParam, NoRegul=NoRegul)'
 res = -1
 GOTO, CLOSING
 end
 
 if not keyword_set(niter) then niter =10
 if not keyword_set(CvgParam) then CvgParam=0.1
 if not keyword_set(RegulParam) then RegulParam = 0.2
 if not keyword_set(NoRegul) then NoRegul = 0
  
 if NoRegul EQ 1 then CvgParam = 1 
 ;-- filtre Passe-Haut pour contraintes de regularisation
 ;=======================================================
fph = fltarr(3,3)
fph[0,0:2]=[0,1,0]
fph[2,0:2]=[0,1,0]
fph[1,0:2]=[1,-4,1]

 ;-- taille de l'image 
 ;====================
sz = size(ima)
Nx = sz[1]
Ny = sz[2]

 ;--initialisation
 ;================
 res = ima
 psf_cf = dft(Psf)
 psf_cf_conj = conj(psf_cf)
 
 flux = total(ima)
 flux_psf = total(psf)


  ;-- Algorithme recursif de deconvolution au sens des moindres carres
 ;===================================================================
 FOR i=0,niter-1 do Begin
    print, '++++ Iteration ', i+1 
    ;-- Codage et convolution de l'objet 
    ;===================================
    conv_res = float ( dfti( dft(res)*Psf_cf) )
    conv_res = conv_res / flux_psf
    residu = ima - conv_res
    print, 'Residu :'
    info, residu
    
    ;-- Convolution avec le symetrique de la PSF
    ;============================================
     conv_res = float ( dfti( dft(residu)*Psf_cf_conj) )
     conv_res = conv_res / flux_psf

    ;-- Contrainte de regularite pour rendre la solution plus lisse
    ;-- on utilise une filtre passe-haut : le LAPLACIEN par exple
    ;--------------------------------------------------------------
    if (NoRegul EQ 0) then $
      res_contraint = RegulParam * convol( $
             convol(res,fph,/edge_wrap),fph, /edge_wrap) $
    else res_contraint = 0.
    
    ;-- Mise a jour de l'objet decode
    ;================================
    res = res + CvgParam*(conv_res - res_contraint)

    ind = where( res LT 0, count)
    if count GT 0 then res[ind] = 0
    
    res = res / total(res) * flux
    
    
 EndFor

CLOSING:
 return, res

end



