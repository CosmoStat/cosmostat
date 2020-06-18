;+
; NAME:
;      mrs_gmca
;	
; PURPOSE:
;	 Apply the blind component separation method called GMCA on spherical maps using Healpix representation (NESTED format) with the orthogonal wavelet transform on the sphere.
;		
; CALLING:
;	mrs_gmca, Data, NbrSource, RecSource, MixintMat, NbrIter, NbrScale=NbrScale, A0=A0, VARMEAN=VARMEAN, 
;			SpecConst=SpecConst, ForceGalactic=ForceGalactic, DustAlmConst=DustAlmConst, SyncAlmConst=SyncAlmConst, GMCA_mask=GMCA_mask
;
; INPUTS:
;   Data -- 2D IDL float array: Multichannel data on the Sphere (Healpix nested format). Data[*,i] is the ith channel.
;   NbrSource -- int: Number of sources. The number of source must be smaller or equal to the number of channels
;   NbrIter -- int: Number of iterations
;
; OUTPUTS:
;   RecSource -- 2D IDL float array: Multichannel reconstructed sources. RecSource[*,i]  (i=0,NbrSource-1) is the ith sources
;   MixintMat -- 2D IDL float array: Estimated mixing matrix ( Data = MixintMat # RecSource ).
;                                    MixintMat is of size NumberChannels x NbrSource ( MixintMat[0:NumberChannels-1, 0:NbrSource-1] )
;             
;  INPUT KEYWORDS:
;   SpecConst -- 2D IDL float array  : Initial mixing matrix with CMB and SZ spectral contraints
;   NbrScale -- int : Number of scales (default is 4) in the wavelet transform.
;	VARMEAN	 -- 2D IDL float array : noise covariance matrix given on input, if not set, it is assumed it is identity matrix
;   A0: 1D IDL float array [0..NumberChannels-1] : if set, then the first column of the matrix is considered
;                               as known and does not need to be estimated. We have MixintMat[*, 0] = A0
;   ForceGalactic : if SpecConst is set, it fixes the vales of the Dust and Synchrotron spectra 
;   DustAlmConst : if set, it constraints the first 769 alm  coefficients of the 2nd estimated sources (being the estimated dust)
;   SyncAlmConst : if set, it constraints the first 769 alm  coefficients of the 3rd estimated sources (being the estimated synchrotron)
;   GMCA_mask : if set the input data are masked
;
; EXTERNAL CALLS:
;       mrs_owttrans  
;   	get_all_faces
;   	mrs_allbandcleaning_perband
;   	put_all_faces
;
; EXAMPLE:
;       Compute the GMCA on a 2D data set considering 4 sources, 100 iterations and 5 wavelet scales
;          mrs_gmca, Data, 4, Sources, mat, 100, NbrScale=5
;         
; HISTORY:
;	Written: Jerome Bobin, 2007
;	April, 2007 File creation
;--------------------------------------------------------------------------------------------------------


pro mrs_gmca, x, ns, s_pi, AA, nit_max, NbrScale=NbrScale, InitMixMat=InitMixMat, BandThres=BandThres, A0=A0, VARMEAN=VARMEAN, SpecConst=SpecConst, ForceGalactic=ForceGalactic, DustAlmConst=DustAlmConst, SyncAlmConst=SyncAlmConst, GMCA_mask=GMCA_mask

;### FAST "GMCA" - Version du 26/03/07 with tuned thresholding
;### ns : nb de sources
;### s : sources estimees
;### A : matrice de melange estimee
;### nit_max : nombre d'iterations

if not keyword_set(NbrScale) then NbrScale=4

nx = size(x)
n = nx[1]
npix = n
nside = npix2nside(npix)
nc = nx[2]
s = dblarr(npix,ns)

if keyword_set(VARMEAN) then begin 
	print,'Noise Variance'
	iRn = invert(VARMEAN,/double)
endif

if not keyword_set(VARMEAN) then iRn = identity(nc)

if not keyword_set(InitMixMat) then begin
	AA = randomn(seed,nc,ns)
end else begin
	AA = InitMixMat
end

if keyword_set(SpecConst) then begin 

	AA = SpecConst
	print,'With CMB and SZ constraint'
	
	if keyword_set(ForceGalactic) then print,'With Dust and Synchrotron constraints'

endif

if keyword_set(DustAlmConst) then begin
	
	print,'Dust Alm Constraint'
	
	restore,'AlmDust.xdr',/verb
	
	ForceLmax = (size(AlmDust))[1]
			 
endif
		 
if keyword_set(SyncAlmConst) then begin
	
	print,'Synchrotron Alm Constraint'
	
	restore,'AlmSync.xdr',/verb
	
	ForceLmax = (size(AlmSync))[1]
			 
endif

if keyword_set(GMCA_mask) then print,'Mask is used'

if keyword_set(A0) then begin
	AA[*,0] = A0  ;--- Dans le cas ou l'on ne force qu'une colonne
end

Sc_Factor = 0.*dblarr(1,ns)

for pp=0,ns-1 do begin

        Sc_Factor[pp] = sqrt(norm(transpose(AA(*,pp))#AA(*,pp)))

        AA[*,pp] = AA[*,pp]/Sc_Factor[pp]

end

;moy = 0.*dblarr(1,nc)
;for uu=0,nc-1 do begin
;        moy(uu) = mean(x[*,uu])       
;        x[*,uu] = x[*,uu] - moy(uu)
;endfor

;#### Seuil en k*sigma_ns avec k decroissant lineairement

ts_max = 15.0  ;---- Regle le seuil initial
ts_min = 3.0
ts = ts_max

dts = (ts_max - ts_min)/(nit_max - 1.0)

nit = 0.0

;#### Transformation des données en ondelettes




            	    if not keyword_set(GMCA_mask) then mrs_owttrans, x[*,0],coefs1,  NbrScale=NbrScale
		    
		    if keyword_set(GMCA_mask) then mrs_owttrans, x[*,0]*GMCA_mask,coefs1,  NbrScale=NbrScale
				
		    nc1 = npix 
		    xcoeffs = dblarr(nc,nc1)
		    
		    put_all_faces,coefs1.coef,temp
		    
		    xcoeffs[0,*] = double(reform(temp,1,nc1));	float(reform(temp,1,nc1))
		    
		    ;xcoeffs[0,*] = float(reform(coefs1.coef,1,nc1))
		    
		    for tt=1,nc-1 do begin
		            
		            if not keyword_set(GMCA_mask) then mrs_owttrans, x[*,tt], coefs1,  NbrScale=NbrScale
			    
			    if keyword_set(GMCA_mask) then mrs_owttrans, x[*,tt]*GMCA_mask, coefs1,  NbrScale=NbrScale
		            
			    put_all_faces,coefs1.coef,temp
			    
			    xcoeffs[tt,*] = double(reform(temp,1,nc1));	float(reform(temp,1,nc1))
			    
		            ;xcoeffs[tt,*] = float(reform(coefs1.coef,1,nc1))
		            
		     endfor
		     
		     temp = 0

    print,'Starting FastGMCA'

    ;#### Corps du programme fast_gmca

    while (nit lt nit_max) do begin
               print,'niter ',nit  , ' ts - mad', ts
                nit = nit + 1.0

                           scoeffs = 0
                           ;W=fs_invert_2D_SVD(transpose(AA)#iRn#AA)
                           scoeffs = invert(transpose(AA)#iRn#AA,/DOUBLE)#transpose(AA)#iRn#xcoeffs;	W#transpose(AA)#iRn#xcoeffs;
                                                     
                           for pp =0,ns-1 do begin                     ;---- Inutile de seuiller le cmb                   

								if keyword_set(BandThres) then begin

									get_all_faces, transpose(scoeffs[pp,*]), temp
                            
									mrs_allbandcleaning_perband, temp, temp, KK = ts, NbrScale=NbrScale
                                     
									put_all_faces,temp,temp_i		
                                     
									scoeffs[pp,*] = transpose(temp_i)

									temp = 0
									
									temp_i = 0
									
								end else begin
								
									temp = scoeffs[pp,*]
									
									thrd = ts*mad( temp )
									
									ind = where( temp LT thrd, count )
									
									if count NE 0 then temp[ind] = 0.0
									
									scoeffs[pp,*] = temp
								
								end
        
                           endfor
        
		 ;##### Constraints on the Alm : (could be done in the wt domain)
		 
		 if keyword_set(DustAlmConst) then begin
		
			if norm(scoeffs[2,*]) gt 1e-6 then begin
		 
		 		get_all_faces,transpose(scoeffs[2,*]),temp
			
				coefs1.coef = temp
			
				mrs_owtrec,coefs1,DustMap
			
				mrs_almtrans,DustMap,DustTrans,/tab
			
				DustTrans.alm[0:ForceLmax-1,0:ForceLmax-1,*] = AlmDust
			
				mrs_almrec,DustTrans,DustMap
			
				DustTrans = 0
			
				mrs_owttrans, DustMap, coefs1,  NbrScale=NbrScale
			
				DustMap = 0
		            
				put_all_faces,coefs1.coef,temp
			    
				scoeffs[2,*] = double(reform(temp,1,nc1));	float(reform(temp,1,nc1))
			
				temp = 0
		 	
			endif
		 
		 endif
		 
		 if keyword_set(SyncAlmConst) then begin
		 
		 	if norm(scoeffs[3,*]) then begin
		 
		 		get_all_faces,transpose(scoeffs[3,*]),temp
			
				coefs1.coef = temp
			
				mrs_owtrec,coefs1,SyncMap
			
				mrs_almtrans,SyncMap,SyncTrans,/tab
			
				SyncTrans.alm[0:ForceLmax-1,0:ForceLmax-1,*] = AlmSync
			
				mrs_almrec,SyncTrans,SyncMap
			
				SyncTrans = 0
			
				mrs_owttrans, SyncMap, coefs1,  NbrScale=NbrScale
			
				SyncMap = 0
		            
				put_all_faces,coefs1.coef,temp
			    
				scoeffs[3,*] = double(reform(temp,1,nc1));	float(reform(temp,1,nc1))
			
				temp = 0
		 
		 	endif
		 
		 endif
	
                         
                 ;##### Estimation de AA :
                          
                 n_v = dblarr(1,ns)
                 for ll=0,ns-1 do n_v(ll) = norm(scoeffs[ll,*])
                
                if keyword_set(A0) then begin
					n_v[0] = 0.0 ;--- Dans le cas ou l'on ne force qu'une colonne
				end  
				
				if keyword_set(SpecConst) then begin
				
					n_v[1] = 0.0   ;---- Pour eviter de changer des valeurs de SZ_CL  

					if keyword_set(ForceGalactic) then n_v[2:3] = 0.	
				
				endif
				
                indix = where(n_v gt 1e-6)
                  
                if (max(n_v) gt 1e-6) then begin
                				;sz_id=size(indix)
                				;if sz_id[1] eq 1 then begin
                				;	W=scoeffs[indix,*]#transpose(scoeffs[indix,*])
                				;	W=1/double(W)
                				;end else begin
                             	;	W=fs_invert_2D_SVD(scoeffs[indix,*]#transpose(scoeffs[indix,*]))
                             	;end                   
                             tA = xcoeffs#transpose(scoeffs[indix,*])#invert(scoeffs[indix,*]#transpose(scoeffs[indix,*]),/DOUBLE)
                             ;xcoeffs#transpose(scoeffs[indix,*])#W
                                  
                            AA(*,indix) = tA
                                  
                endif          
                          
                for pp=0,ns-1 do begin
                          
                           na = sqrt(norm(transpose(AA(*,pp))#AA(*,pp)))

                            AA[*,pp] = AA[*,pp]/na
                                                                            
                endfor
                          
                  
                ts = ts - dts

    endwhile

    print,'Endind FastGMCA'
	;W=fs_invert_2D_SVD(transpose(AA)#AA)
    s_pi = transpose(invert(transpose(AA)#AA,/DOUBLE)#transpose(AA)#transpose(x));	transpose(W#transpose(AA)#transpose(x));

end


;--- cleaning in dwt at k*mad  / from the bands / Mad on all bands

pro mrs_allbandcleaning_perband,cc,cc_out,KK = kk,NbrScale=NbrScale

	KS = KK

	TransFace_in = cc

	N=double((size(TransFace_in))[1])

	TransFace_out = TransFace_in

	Nd = 0
	Nf = N/(2.^(NbrScale-1.))

	u = TransFace_in[Nd:Nf-1,Nd:Nf-1.,*]

	thrd = ks*mad(u)

	u = u*(abs(u) ge thrd)

	TransFace_out[Nd:Nf-1,Nd:Nf-1.,*] = u

	for pp=0,NbrScale - 2 do begin ;--- for each scale

        	Nd = N/(2.^(pp+1))
        	Nf = N/(2.^pp)-1.
		
		    udi = TransFace_in[Nd:Nf,Nd:Nf,*]
		    udr = TransFace_in[Nd:Nf,0:Nd-1.,*]
		    uga = TransFace_in[0:Nd-1.,Nd:Nf,*]
		
		    ut = [[[udi]],[[udr]],[[uga]]]
		
		    thrd = ks*mad(ut)

        	;--- The diagonal term
        
        	u = TransFace_in[Nd:Nf,Nd:Nf,*]
        	
        	;thrd = ks*mad(u)

		    u = u*(abs(u) ge thrd)
        
        	TransFace_out[Nd:Nf,Nd:Nf,*] = u
        
        	;--- The first off-diagonal term
        
        	u = TransFace_in[Nd:Nf,0:Nd-1.,*]
        	
        	;thrd = ks*mad(u)

	       	u = u*(abs(u) ge thrd)
        
        	TransFace_out[Nd:Nf,0:Nd-1.,*] = u
        
        	;--- The second off-diagonal term
        
        	u = TransFace_in[0:Nd-1.,Nd:Nf,*]
        	
        	;thrd = ks*mad(u)

    		u = u*(abs(u) ge thrd)

        	TransFace_out[0:Nd-1.,Nd:Nf,*] = u
				
        endfor
	
	cc_out = TransFace_out
	



end
