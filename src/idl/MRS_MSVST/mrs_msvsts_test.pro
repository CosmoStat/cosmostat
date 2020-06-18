;+
; NAME:
;        mrs_msvsts_test
;
; PURPOSE:
;	Test the different routines of MS-VSTS with a set of FERMI simulated data.
;
; CALLING:
;
; I - Procedures for monochannel data

;    msvsts_load_data,  model, sources,  modelsources, imagepoisson
;			Loading of monochannel simulated Fermi data
;			OUTPUTS : model : simulated galactic background model
; 								sources : simulated Point sources
; 								modelsources : simulated Fermi intensity map (model + sources)
; 								 imagepoisson : simulated Fermi data with Poisson noise
;								All outputs are Healpix data with nside=128
;       
;     msvsts_ldenoisingsimple,  imagerec,  plot=plot,  updatesup=updatesup,  nbrscale=nbrscale, niter=niter
;            call msvsts_load_data to load the data and run the denoising
;			 OUTPUTS : imagerec = denoised image
;  
;    denoising_impainting,imagepoisson,imagerec,plot=plot,alm=alm,curv=curv,nbrscale=nbrscale,niter=niter, model=model, sources=sources, inputorig=inputorig, mask=mask
;			 OUTPUTS : imagepoisson = input noisy masked image
; 								imagerec : reconstructed image
; 								sources : simulated Point sources
;								model: diffuse emission
;                               inputorig: input image = sources + model
;                               mask: mask of missing data
;     
;       denoising_separation,sourcesrec,updatesup=updatesup,plot=plot,nbrscale=nbrscale,niter=niter
;
;; II - Procedures for multichannel data
;
;    mc_load_simu,simuintensity,simupoisson
 ;   multichannel_denoising,simu_denoised,nbrscale1=nbrscale1,nbrscale2=nbrscale2,niter=niter
;    multichannel_deconvolution,simu_deconv,nbrscale1=nbrscale1,nbrscale2=nbrscale2,niter=niter
;         
; HISTORY:
;	Written: Jérémy Schmitt & Jean-Luc Starck, 2012
;	February, 2009 File creation
;--------------------------------------------------------------------------------------------------------



;;;;;;;;;;;;;;;Test file

;; I - Procedures for monochannel data

pro msvsts_load_data,  model,sources,modelsources,imagepoisson
;;;;;;;;;;;;;;; Loading of monochannel simulated Fermi data
;;;;;;;;;;;;;;; OUTPUTS : model : simulated galactic background model
;;;;;;;;;;;;;;;;;;;;;;;; sources : simulated Point sources
;;;;;;;;;;;;;;;;;;;;;;;; modelsources : simulated Fermi intensity map (model + sources)
;;;;;;;;;;;;;;;;;;;;;;;; imagepoisson : simulated Fermi data with Poisson noise
;;;;;;;;;;;;;;;;;;;;;;;; All outputs are Healpix data with nside=128
FERMI_DIR = getenv("ISAP") + '/data/MSVST_FERMI_TEST_DATA/'
model=mrs_read(FERMI_DIR + 'model.expoconv_healpix.fits')
sources= mrs_read(FERMI_DIR + 'source_healpix.fits')

modelsources=model+sources
imagepoisson=poisson_image(modelsources)

end

;;;;;;;;;;;;;; Denoising of simulated Fermi data
pro msvsts_ldenoisingsimple,imagerec,plot=plot,updatesup=updatesup,nbrscale=nbrscale,niter=niter, model=model, sources=sources, inputorig=inputorig, imagepoisson=imagepoisson

if not keyword_set(nbrscale) then nbrscale=6
if not keyword_set(niter) then niter=5

msvsts_load_data, model,sources, modelsources, imagepoisson
inputorig = modelsources
if keyword_set (updatesup) then begin
 mrs_msvsts_iuwt_denoising,imagepoisson,sourcesrec,nbrscale=nbrscale,niter=niter,/hsd,/separation,back_reconstruit=back_reconstruit,/update_support,support
 imagerec=sourcesrec+back_reconstruit
endif else begin
 mrs_msvsts_iuwt_denoising,imagepoisson,imagerec,nbrscale=nbrscale,niter=niter,/hsd
endelse

if keyword_set(plot) then tvs,alog(1+imagerec)

end


;;;;;;;;;;;;;; Inpainting on simulated Fermi data
pro denoising_impainting,imagepoisson,imagerec,plot=plot,alm=alm,curv=curv,nbrscale=nbrscale,niter=niter, model=model, sources=sources, inputorig=inputorig, mask=mask

if not keyword_set(nbrscale) then nbrscale=6
if not keyword_set(niter) then niter=50

msvsts_load_data, model,sources,modelsources,imagepoisson
inputorig = modelsources

FERMI_DIR = getenv("ISAP") + '/data/MSVST_FERMI_TEST_DATA/'
mask=rims(FERMI_DIR+'mask_healpix.fits')
imagepoisson=imagepoisson*mask

if keyword_set(alm) then begin
 mrs_msvsts_iuwt_denoising,imagepoisson,imagerec,nbrscale=nbrscale,niter=niter,/hsd,mask=mask,/alm
endif else begin
 if keyword_set(curv) then begin 
  mrs_msvsts_iuwt_denoising,imagepoisson,imagerec,nbrscale=nbrscale,niter=niter,/hsd,mask=mask,/curv
 endif else begin
  mrs_msvsts_iuwt_denoising,imagepoisson,imagerec,nbrscale=nbrscale,niter=niter,/hsd,mask=mask
 endelse
endelse

if keyword_set(plot) then tvs,alog(1+imagerec)

end

;;;;;;;;;;;;;;; Background extraction on simulated Fermi data
pro denoising_separation,sourcesrec,updatesup=updatesup,plot=plot,nbrscale=nbrscale,niter=niter, model=model, sources=sources, inputorig=inputorig 

if not keyword_set(nbrscale) then nbrscale=6
if not keyword_set(niter) then niter=5

msvsts_load_data,model,sources,modelsources,imagepoisson
inputorig = modelsources

if keyword_set (updatesup) then begin
 mrs_msvsts_iuwt_denoising,imagepoisson,sourcesrec,nbrscale=nbrscale,niter=niter,/hsd,background=model,/separation,back_reconstruit=back_reconstruit,/update_support,support
endif else begin
 mrs_msvsts_iuwt_denoising,imagepoisson,sourcesrec,nbrscale=nbrscale,niter=niter,/hsd,background=model
endelse 

if keyword_set (plot) then tvs,alog(1+sourcesrec)

end


; ===================================
;; II - Procedures for multichannel data

;;;;;;;;;;;;;; Loading of multichannel simulated Fermi galactic background (healpix data of nside=512 on 14 channels)
pro mc_load_model,model,sizeimage,sizemc
FERMI_DIR = getenv("ISAP") + '/data/MSVST_FERMI_TEST_DATA/'

 temp=mrdfits(FERMI_DIR+ 'total_noIso_9_skymap_expoconv_skymap.fits',1)
 tempsp=temp.spectra
 sizedata=size(tempsp)
 sizemc=sizedata[1]
 sizeimage=sizedata[2]
 model=fltarr(sizeimage,sizemc)
 tempnest=fltarr(sizeimage)
 tempring=fltarr(sizeimage)
 for j=0,sizemc-1 do begin
  tempnest[*]=tempsp[j,*]
  ipnest=findgen(sizeimage)
  ipnest=ulong(ipnest)
  nest2ring,256,ipnest,ipring
  tempring=tempnest[ipring]
  model[*,j]=tempring
 end
end


;;;;;;;;;;;;; Loading of multichannel simulated Fermi point sources (healpix data of nside=512 on 14 channels)
pro mc_load_source,sources,sizeimage,sizemc
FERMI_DIR = getenv("ISAP") + '/data/MSVST_FERMI_TEST_DATA/'

 temp=mrdfits(FERMI_DIR+ 'src_healpix256_skymap.fits',1)
 tempsp=temp.spectra
 sizedata=size(tempsp)
 sizemc=sizedata[1]
 sizeimage=sizedata[2]
 sources=fltarr(sizeimage,sizemc)
 tempnest=fltarr(sizeimage)
 tempring=fltarr(sizeimage)
 for j=0,sizemc-1 do begin
  tempnest[*]=tempsp[j,*]
  ipnest=findgen(sizeimage)
  ipnest=ulong(ipnest)
  nest2ring,256,ipnest,ipring
  tempring=tempnest[ipring]
  sources[*,j]=tempring
 end
end

;;;;;;;;;;;;;; Loading of multichannel simulated Fermi data (healpix data of nside=512 on 14 channels)
;;;;;;;;;;;;;; OUTPUTS : simuintensity : simulated Fermi intensity map
;;;;;;;;;;;;;;;;;;;;;;;; simupoisson : simulated Fermi map with Poisson noise
pro mc_load_simu,simuintensity,simupoisson
 mc_load_model, model
 mc_load_source,sources
 simuintensity=model+sources
 simupoisson=poisson_image(simuintensity)
end

;;;;;;;;;;;;; Loading of convolution beams (4001*14 array)
pro beam_fermi,beam
 FERMI_DIR = getenv("ISAP") + '/data/MSVST_FERMI_TEST_DATA/'

 beam=fltarr(4001,14) 
 psfchan=mrdfits(FERMI_DIR+ 'psf.fits',3)
 for j=0,13 do begin
  ind=where(psfchan[179,*,j] GE max(psfchan[179,*,j])/2.)
  ;ind=where(psfchan[179,*,j] GE max(psfchan[179,*,j])/4.)
  fwhm=(max(ind)-min(ind))*0.1*60.
  beamtemp=mrs_getbeam(fwhm=fwhm)
  beam[*,j]=beamtemp
 end

end

;;;;;;;;;;;;;; Multichannel denoising of simulated Fermi data
pro multichannel_denoising,simu_denoised,nbrscale1=nbrscale1,nbrscale2=nbrscale2,niter=niter

if not keyword_set(nbrscale1) then nbrscale1=6
if not keyword_set(nbrscale2) then nbrscale2=6
if not keyword_set(niter) then niter=10

mc_load_simu,simuintensity,simupoisson

mrs_msvsts_multichannel_denoising,simupoisson,simu_denoised,NbrScale1=NbrScale1, NbrScale2=NbrScale2,niter=niter

end

;;;;;;;;;;;;;; Multichannel denconvolution of simulated Fermi data
pro multichannel_deconvolution,simu_deconv,nbrscale1=nbrscale1,nbrscale2=nbrscale2,niter=niter

if not keyword_set(nbrscale1) then nbrscale1=6
if not keyword_set(nbrscale2) then nbrscale2=6
if not keyword_set(niter) then niter=50

mc_load_simu,simuintensity,simupoisson
beam_fermi,beam
mrs_msvsts_multichannel_deconvolution,simupoisson,simu_deconv,NbrScale1=NbrScale1, NbrScale2=NbrScale2,niter=niter,beam=beam,/regularization

end