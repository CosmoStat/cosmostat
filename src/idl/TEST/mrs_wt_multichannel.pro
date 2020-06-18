;; NAME : mrs_wt_multichannel
;;
;; PURPOSE : compute the isotropic undecimated wavelet transform of multichannel data on the sphere
;;
;; CALLING : mrs_wt_multichannel,data,trans
;;
;; Input : data : multi-channel data on the sphere
;;         2D array : data[pixel,energy]
;;
;; Output : trans : multi-channel wavelet transform of data on the sphere
;;         4D array : trans[pixel,scale1,energy]    
;;
;; Keywords : NbrScale1 : number of scales along the spatial coordinates
;; Keywords : NbrScale2 : number of scales along the multichannel coordinate (time or energy)
;;
;; Written by Jérémy Schmitt, 1/04/2010

pro mrs_wt_multichannel,data,trans,NbrScale1=NbrScale1,NbrScale2=NbrScale2

if not keyword_set(nbrscale1) then nbrscale1=6
if not keyword_set(nbrscale2) then nbrscale2=6

SizeDataMC=size(data)
SizeImage=SizeDataMC[1]
SizeMC=SizeDataMC[2]
image_wt_mc=fltarr(sizeimage,nbrscale1,sizemc)
vecMC=fltarr(SizeMC)
trans=fltarr(sizeimage,nbrscale1,sizemc)

for j=0,SizeMC-1 do begin
 image=data[*,j]
 mrs_wttrans,image,image_wt,NbrScale=NbrScale1
 image_wt_mc[*,*,j]=image_wt.coef
end


for sc1=0,nbrscale1-1 do begin
 for pixel=ulong(0),SizeImage-1 do begin
  vecMC[*]=image_wt_mc[pixel,sc1,*]
  atwt1d,vecMC,vecMC1dwt,nscale=nbrscale2
  vecMC1dwt= bwt01_lift(vecMC, nbrscale2-1)
  ;help,vecMC1dwt
  ;trans[pixel,sc1,*,*]=vecMC1dwt
  trans[pixel,sc1,*]=vecMC1dwt
 end
end

end