;+
; NAME: 
;       MRS_ALM_INPAINTING
;
; PURPOSE:
;        Apply an inpainting to a spherical map based on a sparse regularization of the spherical harmonics coefficients.
;        The iterative hard thresholding algorithm is used.
;        If a mask is not provided, all pixels with a zero value are considered as missing pixels.
;        A c++ program ($MRS/cxx/mrs_alm_inpainting) is called.
;       If file names are given for the input image and mask, these two images are not loaded into IDL.
;
; CALLING:
;     InpaintMap = mrs_alm_inpainting(Imag, Mask=Mask, Niter=Niter, OutPowSpec=OutPowSpec, lmax=lmax, gauss=gauss)
;
; INPUTS:
;     Imag -- IDL 1D array: Input Healpix image to be inpainted 
;     
; OUTPUTS:
;     InpaintMap -- IDL 1D array: Output inpainted Healpix image   
;          
; INPUT KEYWORDS:
;      niter: int: number of iterations used in the reconstruction
;      Lmax      : Number of spherical harmonics computed in the decomposition
;                                       (default is 3*nside, should be between 2*nside and 4*nside)
;      FNin: Filename containing the input image. The input image won't be read.
;      FNMask: Filename containing the mask. The mask won't be read.
;      FNOut: Filename containing the results. By default, nothing is written on the disk.
;      FNPS: Filename containing the mask. The mask won't be read.
;
; OUTPUT KEYWORDS:
;     OutPowSpec: IDL 1D array: Cl of the inpainted map. 
;       
; EXTERNAL CALLS:
;       mrs_alm_inpainting (C++ program)
;
; EXAMPLE:
;      
; HISTORY:
;       Written : Jean-Luc Starck   2009.
;-
;-----------------------------------------------------------------
; ISAP> q = mrs_alm_inpainting(i*mask, niter=200, opt='-m1 -W -v -C -N ', lmax=64)
; Option for log-normal galaxies field inpainting


function mrs_alm_inpainting,  Ima, Mask=Mask, Niter=Niter,  OutPowSpec=OutPowSpec,  OutAlm= OutAlm, lmax=lmax,  Opt=Opt,  FNMask=FNMask, FNOut=FNOut, FNIn=FNin, FNPS=FNPS, gauss=gauss, verb=verb, idl=idl, Cl=Cl, Isotropic=Isotropic, log=log, SetOutPowSpec=SetOutPowSpec, PRea=PRea,  MatMask=MatMask, PNiter=PNiter

COMMON MR1ENV
COMMON C_PLANCK

if N_PARAMS() LT 1 then begin 
        print, 'CALL SEQUENCE: InpaintMap = mrs_alm_inpainting(Imag, Mask=Mask, Niter=Niter,  OutPowSpec=OutPowSpec,  lmax=lmax,  FNMask=FNMask, FNOut=FNOut, FNIn=FNin, FNPS=FNPS, gauss=gauss'
        goto, DONE
        end

EpsLog = 1.
if  keyword_set(log) then Imag = alog(Ima+EpsLog) $
else Imag = Ima

if keyword_set(IDL) then begin 
   InpaintMap = cmb_l1_alm_inpainting(Imag, Mask=Mask,  niter=niter, FNImag=FNin, FNMask=FNMask, FNOut=FNOut, OutPowSpec=OutPowSpec, OutAlm=OutAlm)
end else BEGIN

Nx = (size(Ima))[1]
Ny = (size(Ima))[2]

if not keyword_set(Mask)  then BEGIN
     if keyword_set(FNIn) then Imag = mrs_read(FNIn)
      Mask = Imag
      Mask[*] = 0
      ind = where(Imag NE 0, c)
      if c GT 0 then Mask[ind] = 1
END

 ind = where(mask NE 1, c)
 if c GT 0 then mask[ind] = 0
 indMask = where( mask NE 0, cMask)
 MeanVal = mean(Imag[indMask])


if not keyword_set(Opt) then Opt = '  '
OptI = Opt
; if keyword_set(log) then OptI = OptI + ' -C ' +  '  '
if keyword_set(Niter) then OptI = OptI + ' -i ' + strc(Niter) +  '  '
if keyword_set(lmax) then OptI = OptI + ' -l ' + strc(lmax) +  '  '
if keyword_set(gauss) then OptI = OptI + ' -G '  +  '  '
if keyword_set(Verb) then OptI = OptI + ' -v '  +  '  '

if keyword_set(Isotropic) then begin  
    OptI = OptI + ' -I '  +  '  '
	if keyword_set(Cl) and not keyword_set(SetOutPowSpec) then begin
  	     NameCl =gettmpfilename()
  	     writefits, NameCl, double(Cl)
  	    OptI = OptI + ' -S '  +  NameCl + '  '
  	 end  ; else Cl = mrs_powspec(Imag, /nonorm)
 end

EstimReaPowSpec = 0
if keyword_set(SetOutPowSpec) then begin
  if not keyword_set(PRea) then begin
      Pi = mrs_powspec(Imag*Mask)
      P = mrs_deconv_powspec( Pi,  Mask,  MatMask=MatMask,  Niter=PNiter,  NoisePS=NoisePS, lmax=lmax, verb=verb)
      ; info, p
   end else P = Prea
   EstimReaPowSpec = P
;   if not keyword_set(Cl) then  Cli = mrs_get_cl_theo_powspec(P) else Cli= Cl
;   print, "SigmaC", sigma(Cli)
;   help, Cli
;   info, cli
   NamePrea =gettmpfilename() 
   writefits, NamePrea, P
   OptI = OptI + ' -s '  +  NamePrea + '  '
   
 ;   NameCl =gettmpfilename()
 ;   writefits, NameCl, double(Cli)
  ;  OptI = OptI + ' -S '  +  NameCl + '  '
end   

if not keyword_set(FNin)  then NameIma =gettmpfilename()    else NameIma = FNin
if not keyword_set(FNMask)  then NameMask =gettmpfilename()   else NameMask = FNMask
if not keyword_set(FNOut)  then NameResult = gettmpfilename()   else NameResult = FNOut
if not keyword_set(FNPS)  then NameResultPS = gettmpfilename()    else NameResultPS = FNPS
if not keyword_set(FNAlm)  then NameResultAlm= gettmpfilename()    else NameResultAlm = FNAlm

if not keyword_set(FNin)  then mrs_write,  NameIma,  Imag
if not keyword_set(FNMask)  then mrs_write,  NameMask,  Mask
cmd = BIN_ISAPCXX + '/mrs_alm_inpainting '+' '+ OptI + ' '+ NameIma  + '  ' +  NameMask  + '  ' +NameResult + '  ' +  NameResultPS   + ' ' + NameResultAlm

if keyword_set(Verb) then print, cmd
spawn, cmd

InpaintMap = mrs_read(NameResult)
OutPowSpec = mrdfits(NameResultPS, /silent)
fits2alm, index, Alm, NameResultAlm
alm2tab,Alm, TabALM, complex=complex, TabNbrM=TabNbrM, NbrL=lmax
ALM =  TabALM
out=0
nside = gnside(Ima)
npix = (size(imag))[1]
nel = float(N_ELEMENTS(Imag))
NormVal =  sqrt(nel / (4.*!DPI))
OutAlm = {PixelType:0, tab: 1, complex_alm: 0, nside : gnside(Ima), npix:npix, ALM : Alm, norm: 0, NormVal: NormVal,lmin:0,lmax: lmax, TabNbrM: TabNbrM, index:index , x_sky : 0 , y_sky : 0 , nx : 0 , np : 0}

if not keyword_set(FNin)  then delete, NameIma
if not keyword_set(FNMask)  then delete, NameMask
if not keyword_set(FNOut)  then delete, NameResult
if not keyword_set(FNPS)  then delete, NameResultPS
if not keyword_set(FNAlm)  then delete, NameResultAlm
if keyword_set(NamePrea) then delete, NamePrea
END 

if  keyword_set(log) then InpaintMap = exp(InpaintMap) - EpsLog 


return, InpaintMap

 DONE:
 
end

;==============================================================================================
 
