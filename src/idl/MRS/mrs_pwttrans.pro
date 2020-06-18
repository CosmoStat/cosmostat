;+
; NAME:
;        mrs_pwttrans
;
; PURPOSE:
;	Computes the pyramidal wavelet transform on the sphere, using the HEALPix representation (NESTED data
;   representation) or Glesp Data representation. The wavelet function is zonal and its spherical harmonics 
;	coefficients a_l0 follow a cubic box-spline profile.
;
;
; CALLING:
;
;     mrs_pwttrans, Imag, Trans, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave
;       
; INPUTS:
;     Imag -- IDL array of healpix map or IDL structure of a Glesp map: Input image be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;                    NbrScale : int = number of scales 
;                       nside : int = Healpix nside parameter, only present with healpix input image
;                        npix : long = Number of pixels
;                      Scale1 : finest scale (highest frequencies). A IDL array of healpix map or IDL structure of a Glesp map
;                      Scale2 : Second scale    
;                      Scalej : j th scale
;                        ...
;                      ScaleJ : with J = NbrScale, coarsest resolution
;                       lmax  : int = nlmax parameter at the first scale
;		   Tab_lmax[NbrScale] : int array Tab_lmax[j] = lmax at scale j+1, j=0...NbrScale-1
;		  Tab_nside[NbrScale] : int array Tab_nside[j] = nside parameter of the scale j+1, j=0..NbrScale-1 (Healpix input map)
;														 nx number of rings, Glesp parameter of the scale j+1, j=0..NbrScale-1 (Glesp input map)
;					 UseGLESP : int = 1 if the input image was in Glesp format or keyword Healpix_with_Glesp used, otherwise 0
;				    MeyerWave : int = 1 if the keyword MeyerWave used, otherwise 0
;					  DifInSH : int = 1 if the keyword DifInSH used, otherwise 0
;					  	   nx : int = number of rings at the first scale, Glesp parameter, only present with glesp input image
;						   np : int = max number of pixel on a ring at the first scale, Glesp parameter, only present with glesp input image
;
; KEYWORDS:
;		NbrScale -- int: Number of scale (defaut 4)
;		Lmax     -- int: Number of used spherical harmonics (defaut 3*nside, should be between 2*nside and 4*nside)
;		DifInSH		   : If set, compute be difference between two resolution in the spherical harmonic space instead of the direct space.
;		MeyerWave	   : If set, use Meyer wavelets and set the keyword DifInSH
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;   	synfast (healpix software)
;   	alm_product2 (idl)
;   	compute_g (idl)
;   	compute_h (idl)
;
; EXAMPLE:
;
;       Compute the multiresolution of an image I with default options
;        The result is stored in Output
;               mrs_pwttrans, Imag, Output, NbrScale=5
;         
; HISTORY:
;	Written: Pierrick Abrial & Jean-Luc Starck, 2005
;	December, 2004 File creation
;-------------------------------------------------------------------------------------------

pro mrs_glesp_pwttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave 

if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: mrs_pwttrans, Imag, out, NbrScale=NbrScale,lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave'
        goto, DONE
        end
	 
if  not keyword_set(NbrScale) then NbrScale = 4
 
if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end
if type_code(Imag) EQ 8 then begin
  pixel_type=1 ; glesp
  endif else pixel_type = 0 ; healix

GLESP = 1
npix = (size(imag.t_sky))[1]
nside = npix2nside(npix)
nx = imag.nx
np = imag.np
nx_orig = nx
np_orig = np
if not keyword_set(lmax)  then lmax = min([(nx-1)/2,np/4])
if not keyword_set(DifInSH)  then DifInSH = 0
if not keyword_set(MeyerWave) then MeyerWave = 0 else DifInSH = 1

ech = 1.

Hscale = imag
;alm_trans, imag, ALM_HighResolImag, nside=nside, nlmax=lmax, index=index
mrs_almtrans,imag,alm_first,lmax = lmax
ALM_HighResolImag  = alm_first.alm

ech =1
nlmax2 = lmax
nlmax = lmax
BandName = strarr(NbrScale)
Tab_lmax = intarr(NbrScale)
Tab_nside = intarr(NbrScale)

for j=0,NbrScale-2 do begin
  if keyword_set(MeyerWave) then begin
     hgmey, nlmax2, j, h, g, dif=dif 
  end else begin
     compute_g, nlmax2, ech, g
     compute_h, nlmax2, ech, h
  end
  alm_product2, ALM_HighResolImag ,h, alm_h 
  Tab_nside[j] = nx;side
;  if keyword_set(DifInSH) then begin
        alm_product2, ALM_HighResolImag, g, alm_g  
        alm_g_s = alm_first
	alm_g_s.alm = alm_g
	;alm_itrans, alm_g, WScale, nside=nside, nlmax=lmax, index=index
        mrs_almrec,alm_g_s,WScale,nx=nx,np=np  
;  endif
  ; else begin
   ;     alm_h1 = alm_h[0:(nlmax*(nlmax+1.))/8,*]
;	index1 = index[0:(nlmax*(nlmax+1.))/8]
;	nside1 = nside
;	nlmax1=nlmax
;        alm_itrans, alm_h1, LScale, nside=nside/2, nlmax=nlmax1, index=index1
;	alm_trans, LScale, ALM_LScale,nlmax=nlmax1/2, index=index2
;	alm_itrans, ALM_LScale, LRScale, nside=nside, nlmax=nlmax1/2, index=index2
;	WScale = Hscale - LRScale
;	Hscale = LScale
;  end
  Tab_lmax[j] = nlmax
  my_command = 'scale'+strcompress(string(j+1), /remove_all)
  BandName[j] = my_command 
  my_command = my_command +'=WScale'
  my_command = strcompress( my_command, /remove_all)
  ; print, 'cmd = ',  my_command
  ACK = EXECUTE( my_command) 
   
  ALM_HighResolImag = alm_h
  nlmax= nlmax/2								
  ech = ech*2
  ; if j NE NbrScale-2 then 
  np = np /2
  nx = nx/2
endfor

j=NbrScale-1
Tab_lmax[j] = nlmax
Tab_nside[j] = nx
; print,'nx ',nx
;if keyword_set(DifInSH) then alm_itrans, alm_h, LScale, nside=nside, nlmax=lmax, index=index
        alm_h_s = alm_first
	alm_h_s.alm = alm_h

;if keyword_set(DifInSH) then 
mrs_almrec, alm_h_s, LScale, nx=nx,np=np;, nlmax=lmax;, index=index


nside = Tab_nside[0]
TabNorm=[0.97,0.11,0.047,0.02246,0.011,0.0086]

my_command = 'scale'+strcompress(string(j+1), /remove_all)
BandName[j] = my_command 

my_command = my_command +'=LScale'
my_command = strcompress( my_command, /remove_all)
; print, 'cmd = ',  my_command
ACK = EXECUTE( my_command) 

  pyrtrans = 1
  my_command = 'out = { NbrScale : NbrScale,'
  my_command = my_command+'UseGLESP: GLESP, '
  my_command = my_command+'nx : nx_orig, '
  my_command = my_command+' np:np_orig, '
  my_command = my_command+' npix:npix, '
  my_command = my_command+' lmax:lmax, '
  my_command = my_command+' Tab_lmax:Tab_lmax, '
  my_command = my_command+' Tab_nside:Tab_nside, '
  my_command = my_command+' MeyerWave:MeyerWave,'
  my_command = my_command+' pyrtrans :pyrtrans,'
  my_command = my_command+' TabNorm :TabNorm,'
  my_command = my_command+' pixel_type : pixel_type,'

  for j=0, NbrScale-1 do  my_command = my_command+BandName[j]+':'+BandName[j]+','
 my_command = my_command+' DifInSH:DifInSH'
  my_command = strcompress( my_command, /remove_all)
  my_command =my_command+'}'
  ; print, my_command
  ACK = EXECUTE( my_command)
 
DONE:

END

;===================================================================================

function mrs_clean_median, Map, windowsize=windowsize, Threshold=Threshold, Nsigma=Nsigma
Npix = N_ELEMENTS(Map)
if not keyword_set(NSigma) then Nsigma=5.
 MedIma = mrs_median( Map, windowsize=WindowSize)
 Diff = Map - MedIma
if keyword_set(T) then T = Threshold else T = mad(Diff)
ind = where(ABS(Diff) LT Nsigma*T, c)
if c GT 0 then Diff[ind] = 0 else c = 0 
print,  float(c) / float(Npix) * 100.
return, MedIma+Diff
end

;===================================================================================

pro mrs_pwttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave, median=median, WindowSize=WindowSize

COMMON MR1ENV
COMMON C_PLANCK

if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: mrs_pwttrans, Imag, out, NbrScale=NbrScale,lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave'
        goto, DONE
        end
	 
if  not keyword_set(NbrScale) then NbrScale = 4

if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end

TabNorm=[0.97,0.11,0.047,0.02246,0.011,0.0086]

if type_code(Imag) EQ 8 then mrs_glesp_pwttrans, Imag, out, NbrScale=NbrScale, lmax=lmax, DifInSH=DifInSH, MeyerWave=MeyerWave  $
else BEGIN ; HEALPIX
GLESP = 0 
npix = (size(imag))[1]
nside = npix2nside(npix)
if keyword_set(lmax) eq 0 then begin
	lmax = nside *3
	if lmax GT P_Lmax then  lmax = P_Lmax
end
if not keyword_set(DifInSH)  then DifInSH = 0
if not keyword_set(MeyerWave) then MeyerWave = 0 else DifInSH = 1
if not keyword_set(windowsize)  then windowsize = 7.
ech = 1.
 if keyword_set(Median) then begin
   HealpixCXX = 1
   DifInSH=0
 end
 
Hscale = imag
MedIma = imag
 if keyword_set(Median) then  MedIma  = mrs_clean_median(MedIma, windowsize=windowsize, Threshold=Threshold)

  if keyword_set(HealpixCXX) then alm_cxxtrans, MedIma, ALM_HighResolImag, nside=nside, nlmax=lmax, index=index $
  else alm_trans, Hscale, ALM_HighResolImag, nside=nside, nlmax=lmax, index=index
       
ech =1
nlmax2 = lmax
nlmax = lmax
BandName = strarr(NbrScale)
Tab_lmax = intarr(NbrScale)
Tab_nside = intarr(NbrScale)

for j=0,NbrScale-2 do begin
  if j LT 1  and keyword_set(Median) then   begin
      MedIma  = mrs_clean_median(Hscale, windowsize=windowsize, Threshold=Threshold)
     alm_cxxtrans, MedIma, ALM_LScale, nlmax=nlmax, index=index2
  end
  
  if keyword_set(MeyerWave) then begin
     hgmey, nlmax2, j, h, g, dif=dif 
  end else begin
     compute_g, nlmax2, ech, g
     compute_h, nlmax2, ech, h
  end
  alm_product2, ALM_HighResolImag ,h, alm_h 
  Tab_nside[j] = nside
  if keyword_set(DifInSH) then begin
        alm_product2, ALM_HighResolImag, g, alm_g  
        if keyword_set(HealpixCXX) then	alm_cxxitrans, alm_g, WScale, nside=nside, nlmax=lmax, index=index $
	   else alm_itrans, alm_g, WScale, nside=nside, nlmax=lmax, index=index
  end else begin
        alm_h1 = alm_h[0:(nlmax*(nlmax+1.))/8,*]
	index1 = index[0:(nlmax*(nlmax+1.))/8]
	nside1 = nside
	nlmax1=nlmax
        if keyword_set(HealpixCXX) then	begin
	alm_cxxitrans, alm_h1, LScale, nside=nside/2, nlmax=nlmax1, index=index1
	alm_cxxtrans, LScale, ALM_LScale,nlmax=nlmax1/2, index=index2
	alm_cxxitrans, ALM_LScale, LRScale, nside=nside, nlmax=nlmax1/2, index=index2
	endif else begin
	alm_itrans, alm_h1, LScale, nside=nside/2, nlmax=nlmax1, index=index1
	alm_trans, LScale, ALM_LScale,nlmax=nlmax1/2, index=index2
	alm_itrans, ALM_LScale, LRScale, nside=nside, nlmax=nlmax1/2, index=index2
	endelse
	
	WScale = double(Hscale) - double(LRScale)
	Hscale = double(LScale)
  end
  Tab_lmax[j] = nlmax
  my_command = 'scale'+strcompress(string(j+1), /remove_all)
  BandName[j] = my_command 
  my_command = my_command +'=WScale'
  my_command = strcompress( my_command, /remove_all)
  ; print, 'cmd = ',  my_command
  ACK = EXECUTE( my_command) 
   
  ALM_HighResolImag = alm_h
  nlmax= nlmax/2								
  ech = ech*2
  ; if j NE NbrScale-2 then 
  nside = nside /2
endfor

j=NbrScale-1
Tab_lmax[j] = nlmax
Tab_nside[j] = nside
if keyword_set(DifInSH) then begin
   if keyword_set(HealpixCXX) then  alm_cxxitrans, alm_h, LScale, nside=nside, nlmax=lmax, index=index $
                               else alm_itrans, alm_h, LScale, nside=nside, nlmax=lmax, index=index
 endif
nside = Tab_nside[0]

my_command = 'scale'+strcompress(string(j+1), /remove_all)
BandName[j] = my_command 

my_command = my_command +'=LScale'
my_command = strcompress( my_command, /remove_all)
; print, 'cmd = ',  my_command
ACK = EXECUTE( my_command) 

  pyrtrans = 1
  NEEDLETWAVE=0
  B_NeedletParam=0
  Healpix_with_Glesp=0
  my_command = 'out = { NbrScale : NbrScale,'
  my_command = my_command+'UseGLESP: GLESP, '
  my_command = my_command+'nside : nside, '
  my_command = my_command+' npix:npix, '
  my_command = my_command+' lmax:lmax, '
  my_command = my_command+' Tab_lmax:Tab_lmax, '
  my_command = my_command+' Tab_nside:Tab_nside, '
  my_command = my_command+' MeyerWave:MeyerWave,'
  my_command = my_command+' NEEDLETWAVE:NEEDLETWAVE,'
  my_command = my_command+' B_NeedletParam: B_NeedletParam, '
  my_command = my_command+' Healpix_with_Glesp: Healpix_with_Glesp, '
  my_command = my_command+' pyrtrans :pyrtrans,'
  my_command = my_command+' TabNorm :TabNorm,'
  for j=0, NbrScale-1 do  my_command = my_command+BandName[j]+':'+BandName[j]+','
  my_command = my_command+' DifInSH:DifInSH'
  my_command = strcompress( my_command, /remove_all)
  my_command =my_command+'}'
  ; print, my_command
  ACK = EXECUTE( my_command)

END ; HEALPIX

DONE:

END

