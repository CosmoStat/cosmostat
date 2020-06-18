;+
; NAME: 
;       INVRID2D
;
; PURPOSE:
;        Reconstruct an image from its ridgelet transform.   
;
; CALLING:
;     invrid2d, RidTrans, RecImage
;
; INPUTS:
;     RidTrans -- IDL structures; Ridgelet transform structure (see RID2D)   
;     
; OUTPUTS:
;     Imag -- IDL 2D array: reconstructed   image 
;
; EXTERNAL CALLS
;           invlinogram, invpwt1d
;
; EXAMPLE:
;    Apply the ridgelet transform and the inverse ridgelet transform
;     rid2d, ima, Rid
;     invrid2d, Rid, RecIma
;
; HISTORY:
;       Written: Jean-Luc Starck & Yassir Moudden & Ludovic Poupard.
;       September, 2005 File creation
;-
;-----------------------------------------------------------------

PRO wave_recons_img_lp, trans, imrec, numblock=numblock
nligne=trans.NyRad
npoints=trans.NxRad
nbscale=trans.NBRSCALE 
typetrans=2

from=fltarr(nbscale)
to=fltarr(nbscale)

wcoef_ligne=fltarr(npoints,nbscale)
imrec=fltarr(nligne,npoints)

j=0
PosBandInStructure=11
for i=0,nbscale-1 do begin
  if keyword_set(numblock) then band = reform(trans.(PosBandInStructure+i)(j,*,numblock-1)) $
  else band = reform(trans.(PosBandInStructure+i)(j,*))
  vs = size(band)
  nx = vs[1]
  to[i] = nx
  wcoef_ligne[0:nx-1,i] = band
  trans2={ n_band : nbscale, coef : wcoef_ligne, typetrans : typetrans, from : from, to : to}  
  END
  	
for j=0, nligne-1 do begin
        wcoef_ligne[*]=0
        for i=0,nbscale-1 do  begin
	    if keyword_set(numblock) then wcoef_ligne[0:to[i]-1,i] = reform(trans.(PosBandInStructure+i)(j,*,numblock-1)) $
	    else wcoef_ligne[0:to[i]-1,i] = reform(trans.(PosBandInStructure+i)(j,*))
	end
 	trans2.coef = wcoef_ligne
	invpwt1d,trans2,recons
	imrec(j,*)=recons
endfor
END

;=====================================

pro ponderation, indl, indc, nbl, nbc, bsize, poids
poids=dblarr(bsize,bsize)
poidsc=(sin(findgen(bsize)/bsize*!Pi))^2
poidsl=(sin(findgen(bsize)/bsize*!Pi))^2
if (indl eq 0) then poidsl(0:bsize/2-1)=1.
if (indc eq 0) then poidsc(0:bsize/2-1)=1.
if (indl eq nbl-1) then poidsl(bsize/2:bsize-1)=1.
if (indc eq nbc-1) then poidsc(bsize/2:bsize-1)=1.

poids=poidsl#poidsc
end

;=====================================================================================================================

PRO invrid2d, ridtrans, imrec

if (ridtrans.LOCALRIDGELET EQ 0) $
then begin

   wave_recons_img_lp,ridtrans,radon_rec
   invlinogram,radon_rec,imrec

endif else begin 

	overlap = ridtrans.overlap
	t=size(ridtrans.band1)
	nb_block=t[3]
	nb_blockl = sqrt(nb_block)
	nb_blockc = nb_blockl
	ind_block=0
	imrec = fltarr(ridtrans.nxIma, ridtrans.nyIma)	

	for indl=0,nb_blockl-1 do begin											;Begin loop over the different blocks extracted from the input image
		for indc=0,nb_blockc-1 do begin
		
			; destructure,tabstruct,tabstruct1,ind_block
			wave_recons_img_lp, ridtrans, radon_rec, numblock=ind_block+1
			invlinogram,radon_rec,blockrec2
			
			if ind_block eq 0 then bsize=(size(blockrec2))(1) -1
						
			if (indl ne (nb_blockl - 1) ) then begin
				if (indc ne nb_blockc - 1)$
				then blockrec=blockrec2(0:bsize-1,0:bsize-1)$
				else blockrec=blockrec2(0:bsize-1,1:bsize) 
			endif else begin
				if (indc ne (nb_blockc - 1) )$
				then  blockrec=blockrec2(1:bsize,0:bsize-1)$
				else  blockrec=blockrec2(1:bsize,1:bsize)
			endelse											
											
 			if (overlap eq 1) then begin
				ponderation,indl,indc,nb_blockl,nb_blockc,bsize, poids													;maybe these weights don't need to recomputed at each step.to be modified
				imrec(indl*bsize/2:indl*bsize/2+bsize-1, indc*bsize/2:indc*bsize/2+bsize-1)  = imrec(indl*bsize/2:indl*bsize/2+bsize-1, indc*bsize/2:indc*bsize/2+bsize-1) + blockrec*poids 
			endif else imrec(indl*bsize:(indl+1)*bsize -1, indc*bsize:(indc+1)*bsize-1)=blockrec
 			
			ind_block=ind_block+1
			
		endfor
	endfor
	
endelse																			;Endloop over the different blocks extracted from the input image


DONE:
END

