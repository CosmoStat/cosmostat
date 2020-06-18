;----------------------------------------------------------------------------------------------------
; NAME: 
;       RID2D
;
; PURPOSE:
;       Compute the local ridgelet  transform of a square image. The size of the image is nside*nside and
;	nside is a power of two. The local ridgelet transform is computed on blocks of size
;	blocksize which is also a power of two, obviously smaller than nside.
;	The output is an IDL structure.
;       A band at scale j (j=0..NBRSCALE-1) can be extracted using the function
;       function rid_getband(Rid, j) (ex: Scale2 = rid_getband(RidTrans, 2))
;       and a band can be inserted in the transformation using the routine  rid_putband
;       (ex:  rid_putband, RidTrans, Scale2, 2).
;       At a given scale, the output is a 3D array where the third dimension corresponds to the block number.
;
; CALLING:
;     rid2d, image, LocalRid, NbrScale, blocksize, index_table, overlap=overlap 
;
; INPUTS:
;     Imag -- IDL 2D array: Input image to be transformed 
;     Nbrscale --  int: number of scales
;     
; OUTPUTS:
;     LocalRid -- IDL structures with the following fields:
;         NXIMA     -- Number of columns of the input image
;         NYIMA     -- Number of lines of the input image
;         NBRSCALE  -- LONG: Number of scales of the ridgelet transform (input/output parameter)
;                            
;         LOCALRIDGELET  -- Scalar: 1 for a local ridgelet transform and 0 otherwise
;         NYRID     -- LONG: Number of pixels in the y-axis direction
;         BANDj (j=1,..,NBRSCALE) -- 2D or 3D IDL array: Ridgelet coefficients
;                       For a ridgelet transform, BANDj is a 2D array
;		        For a local ridgelet transform, BANDj is a 3D array,
;		          BANDj[*,*,k] = ridgelet transform of the kth block at scale j
;         BSIZE     --  LONG: Block size used in the local ridgelet transform
;         NXB       -- LONG: Number of blocks in the x-axis direction in the local ridgelet transform
;         NYB       -- LONG: Number of blocks in the x-axis direction in the local ridgelet transform
;         OVERLAP   -- LONG: is equal to 1 if blocks are overlapping in the local ridgelet transform
;         TABNORM   -- FLOAT Array[0:NBRSCALE-1]: Normalization value for each scale
;
; KEYWORDS:
;     blocksize -- int: block size, default is 16
;     index_table -- long array, computed prior to calling rid2d, using lino_grid. 
;                   this array of indices is used in the linogram.and depends on the block size
;      Overlap -- Scalar: if set and if block size is set (i.e. local ridgelet ridgelet transform
;                         then the block are overlapping
;
; EXTERNAL CALLS
;           linogram, pwt1d, rid2d_one_block, wave_trans_img
;
; EXAMPLE:
;     rid2d, ima, LocalRid, NbrScale, blocksize, index_table
;    
; HISTORY:
;       Written: Jean-Luc Starck & Yassir Moudden & Ludovic Poupard
;       September, 2005 File creation
;
;-----------------------------------------------------------------------------------------------------

PRO rid2d, image, LocalRid, NbrScale, blocksize=blocksize, index_table=index_table, overlap=overlap   

if not keyword_set(blocksize) then blocksize = 16
if not keyword_set(index_table) then lino_grid, blocksize+1, index_table

;-----------------------------------------------------------------------------------------------------
;recovering a few parameters and defining some constants
;-----------------------------------------------------------------------------------------------------
TabNorm=[0.8546,0.215, 0.1265,0.097,0.0669,0.0483,0.0342,0.0242,0.0171,0.0121]
N=(size(image))(1)
		
bsize= blocksize 
block=fltarr(bsize+1,bsize+1)				;extracted blocks are one line and one column larger than the specified blocksize.									
											;this is required by the linogram transform which is used in the ridgelet transform.


;the number of blocks depends on whether there is or not overlapping
if not keyword_set(overlap) then begin	
	overlap=0
	nb_blockl=N/bsize
	nb_blockc=N/bsize
endif else begin
	overlap=1
	nb_blockl=2.*N/bsize-1.
	nb_blockc=2.*N/bsize-1.
endelse

nb_block=nb_blockl*nb_blockc

;---------------------------------------------------------------------------------------------------------------
; initializing the ridgelet data structure for this blocksize and computing the ridgelet transform on each block
;---------------------------------------------------------------------------------------------------------------
ind_nbloc=0L									;initial value of an index on the set of blocks
indl = 0l
indc = 0l
;block(0:bsize-1,0:bsize-1)=image(indl*bsize:(indl+1)*bsize -1, indc*bsize:(indc+1)*bsize-1)
;duplicating the last line and column to obtain a square array of odd dimensions
;block(bsize,0:bsize-1)=block(bsize-1,0:bsize-1) 
;block(0:bsize,bsize)=block(0:bsize,bsize-1)													
;since the size of the whole image is greater than the size of the blcoks in all cases considered here, 
;these above three command lines are replaced with : 

block(0:bsize,0:bsize)=image(indl*bsize:(indl+1)*bsize, indc*bsize:(indc+1)*bsize)

;for other blocks than this first one, it is necessary to check that we are not going out of the image bounds.
;if this is the case, then we decide to duplicate the last line and last column.
;better strategies would try to include lines and columns on the inside of the image but this then requires modifications of the reconstruction algorithm.
;TO BE ENHANCED IN THE FUTURE
 	
									
rid2d_one_block, block, rid, NbrScale,index_table, typetrans=2							;applying the global ridgelet transform on this first block
Tabstruct=replicate(rid, nb_block) 
		
	
for indl= 0l, nb_blockl-1 do begin				;Begin loop over the different blocks extracted from the input image
	for indc= 0l, nb_blockc-1 do begin
				
				if ind_nbloc ne 0 then begin		;BEGIN IF to avoid useless recomputation of the ridgelet transform of the first block
				
					if not keyword_set(overlap) then begin    ;if no overlap
					
						if (indl ne (nb_blockl - 1) ) then begin
							if (indc ne nb_blockc - 1)$
							then block(0:bsize,0:bsize)=image(indl*bsize:(indl+1)*bsize, indc*bsize:(indc+1)*bsize)$
							else block(0:bsize,0:bsize)=image(indl*bsize:(indl+1)*bsize, indc*bsize-1:(indc+1)*bsize-1)
						endif else begin
							if (indc ne (nb_blockc - 1) )$
							then block(0:bsize,0:bsize)=image(indl*bsize-1:(indl+1)*bsize -1 , indc*bsize:(indc+1)*bsize)$
							else block(0:bsize,0:bsize)=image(indl*bsize-1:(indl+1)*bsize -1, indc*bsize-1:(indc+1)*bsize-1)
						endelse											

					endif else begin			;if overlap
					
						if (indl ne (nb_blockl - 1) ) then begin
							if (indc ne (nb_blockc - 1) )$
							then block(0:bsize,0:bsize) = image(indl*bsize/2:indl*bsize/2+bsize, indc*bsize/2:indc*bsize/2+bsize)$
							else block(0:bsize,0:bsize) = image(indl*bsize/2:indl*bsize/2+bsize, indc*bsize/2-1:indc*bsize/2+bsize-1)
						endif else begin
							if (indc ne (nb_blockc - 1) )$
							then block(0:bsize,0:bsize) = image(indl*bsize/2-1:indl*bsize/2+bsize-1, indc*bsize/2:indc*bsize/2+bsize)$
							else block(0:bsize,0:bsize) = image(indl*bsize/2-1:indl*bsize/2+bsize-1, indc*bsize/2-1:indc*bsize/2+bsize-1) 		
						endelse					
					
					endelse
			
					;Apply the ridgelet transform on the extracted block
					rid2d_one_block, block,rid,NbrScale,index_table,typetrans=2
					Tabstruct(ind_nbloc)=rid
				
				
				endif							;END IF to avoid useless recomputation of the ridgelet transform of the first block
			
				ind_nbloc=ind_nbloc+1
			
	endfor
endfor											;Endloop over the different blocks extracted from the input image
  	
 
 
 
;building the output structure	
LocalRidgelet=1
commande="LocalRid={NxIma: N, NyIma: N, NbrScale:Tabstruct[0].n_band,LocalRidgelet:LocalRidgelet, TabNorm:TabNorm, overlap:overlap, bsize:bsize, NXB:nb_blockl, NYB:nb_blockc, NxRad: Tabstruct[0].NxRad, NyRad:Tabstruct[0].NyRad"
for i=1,NbrScale do commande=commande+",band"+ strcompress(string(i), /remove_all)+":tabstruct.band"+string(i, format="(I1)")
commande=commande+"}"
test=execute(commande)		
  
END





;=====================================================================================================
;SOME USEFUL FUNCTIONS IN THE BACKGROUND ....
;=====================================================================================================



;---------------------------------------------------------------------------------------------------------------
;ridgelet transform on one block
;---------------------------------------------------------------------------------------------------------------

pro rid2d_one_block, image, ridtrans, nbscale, tab_indices, typetrans = typetrans 

	N = (size(image))[1]
	
	no_last_ft=1		;default setting is to not compute two useless successive inverse and forward fft
	no_ft=1
	
	linogram, image, radtrans, tab_indices, no_last_ft = no_last_ft
	wave_trans_img, radtrans, ridtrans, nbscale=nbscale, typetrans=typetrans, no_ft=no_ft
			
END		
 
;---------------------------------------------------------------------------------------------------------------



;---------------------------------------------------------------------------------------------------------------
;compute a 1dwt along all the lines of a 2d array
;---------------------------------------------------------------------------------------------------------------

pro wave_trans_img, image, trans, nbscale=nbscale, typetrans=typetrans, no_ft=no_ft

taille = size(image)
nligne = taille(1)
npoints = taille(2)
wcoef=fltarr(nligne,npoints,nbscale)

for j=0, nligne-1 do begin
		ligne=reform(image(j,*))
		pwt1d, ligne , temp_trans, nbscale = nbscale, typetrans = typetrans, no_ft = no_ft
		wcoef(j,*,*) = temp_trans.coef(*,*)
endfor

n_band=temp_trans.n_band
from=temp_trans.from
to=temp_trans.to

for i=1, n_band do begin
     ma_commande = 'band'+string(i)
     IndexScale = i-1
     ma_commande = ma_commande +'=wcoef[*,from(IndexScale):to(IndexScale)-1,IndexScale]'
     ma_commande = strcompress( ma_commande, /remove_all)
     ACK = EXECUTE( MA_COMMANDE)
endfor
  
;building the output structure  
ma_commande = 'trans = { n_band: n_band,'
ma_commande = ma_commande+'typetrans: typetrans,'
ma_commande = ma_commande+'from: from,'
ma_commande = ma_commande+'to: to,'
ma_commande = ma_commande+'NxRad: npoints,'
ma_commande = ma_commande+'NyRad: nligne,'
  
for i=1, n_band do  ma_commande = ma_commande+'band'+string(i)+':band'+string(i)+','
   
ma_commande = strcompress( ma_commande, /remove_all)
ma_commande = strmid(ma_commande,0,strlen(ma_commande)-1)
ma_commande =ma_commande+'}'
 
  ACK = EXECUTE( MA_COMMANDE)

end
;---------------------------------------------------------------------------------------------------------------



;---------------------------------------------------------------------------------------------------------------
;compute a an index grid to be used for line extraction the linogram transform
;---------------------------------------------------------------------------------------------------------------

pro lino_grid, N,tab_indices

cc=(N-1)/2
tab_indices_x = lonarr(2*N-2,N)
tab_indices_y = lonarr(2*N-2,N)
tab_indices = lonarr(2*N-2,N)
tab_indices_x(0:N-1, *) = (lonarr(N,1)+1) # lindgen(1,N) 
pente =( findgen(N) - cc )/ float(cc)
tab_indices_y(0:N-1, *) = round(   ( lindgen(N) - cc )#pente )  + cc  

tab_indices_x(N:2*N-3, *) = tab_indices_y(N-2 - lindgen(N-2), *)
tab_indices_y(N:2*N-3, *) = tab_indices_x(1:N-2, *) 

tab_indices1 = tab_indices_x + tab_indices_y *N 
tab_indices2 = ( centre( lindgen(N,N) ))( tab_indices1)
tab_indices(*, 0:cc) = tab_indices2(*, cc:N-1)
tab_indices(*, cc+1:N-1) = tab_indices2(*, 0:cc-1)

end

;---------------------------------------------------------------------------------------------------------------




;---------------------------------------------------------------------------------------------------------------
;compute the linogram transform of an image which is a numerical Radon transform
;---------------------------------------------------------------------------------------------------------------

pro linogram, image, trans, tab_indices, no_last_ft = no_last_ft

TF_Image = fft(image, -1)*N_ELEMENTS(image)	

lignes = TF_Image(tab_indices)

if not keyword_set(no_last_ft) then begin
	trans = REAL_PART(  fft(lignes, +1, dimension = 2) / N_ELEMENTS(lignes(0, *)))   
endif else trans = lignes

DONE:
return


END


;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------
;---------------------------------------------------------------------------------------------------------------


