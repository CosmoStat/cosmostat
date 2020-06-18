


PRO invpwt1d, trans, recons
;Written August 2005, Yassir moudden & Ludovic Pourpard 
 
 typetrans = trans.typetrans
 n_band = trans.n_band
 to = trans.to
 from = trans.from
 nbscale = n_band
 
;**********************************************************************************
 
;test de puissance de 2

npoints = to(0)
pow2 = 0
if ( ceil(alog(npoints)/alog(2)) eq floor(alog(npoints)/alog(2))) then pow2 = 1
;**********************************************************************************
 
 
recons = trans.coef( 0:to(nbscale-1)-1, nbscale-1 )
 
 
 CASE 1 OF
 	
;**********************************************************************************

(typetrans eq 1) : BEGIN
 
 
 	for numscale = 0 , nbscale-2 do begin
        

		taille = size( recons )
 		taille = taille[1]
	
		compute_filter_htilde, taille, 1, htilde
 		compute_filter_gtilde, taille, 1, gtilde
 	
		tabfft = fft(recons,-1,/double)*N_ELEMENTS(recons) *htilde+fft(trans.coef(*,nbscale-2-numscale),-1,/double)*N_ELEMENTS(trans.coef(*,nbscale-2-numscale))* gtilde
		recons = real_part(fft(tabfft,1,/double)/N_ELEMENTS(tabfft))
	
 	endfor
 
 
 END
 
;**********************************************************************************

 
(typetrans eq 2) and ( pow2 eq 1) : BEGIN
 
 	
	for numscale = 0 , nbscale-2 do begin
        

		taille = size(recons)
 		taille = taille[1]
	
		recons2 = dblarr(2*taille)		;up-sampling
 		recons2( 2 * findgen(taille) ) = recons
	
		compute_filter_htilde, 2*taille, 1, htilde
 		compute_filter_gtilde, 2*taille, 1, gtilde
 	 
		tabfft = 2 * fft( recons2, -1, /double)*N_ELEMENTS(recons2) * htilde + fft(trans.coef(0:to(nbscale-2-numscale)-1,nbscale-2-numscale),-1,/double)*N_ELEMENTS(trans.coef(0:to(nbscale-2-numscale)-1,nbscale-2-numscale))* gtilde
		recons = real_part(fft(tabfft,1,/double)/N_ELEMENTS(tabfft))
	
 	endfor
 END
 
;**********************************************************************************
 
 (typetrans eq 2) and (pow2 eq 0) : BEGIN
 
 
	for numscale=0 , nbscale-2 do begin
        

		taille = size(recons)
 		taille = taille[1]
		taille2 = to(nbscale-2-numscale)
		taille2=fix(taille2)
		
		
		fft_recons = fft(recons,-1,/double) * N_ELEMENTS(recons)
		fft_recons2 = dcomplexarr(taille2)		;up-sampling
 		
		
		fft_recons2(0:floor(ceil(taille2/2.)/2.)-1) = fft_recons(0:floor(ceil(taille2/2.)/2.)-1)
		fft_recons2(taille2-ceil(taille2/4.):taille2-1) = fft_recons(floor(ceil(taille2/2.)/2.):taille-1)
		
		compute_filter_htilde, taille2, 1, htilde
 		compute_filter_gtilde, taille2, 1, gtilde
 		
		
		tabfft = 2 * fft_recons2 * htilde + fft(trans.coef(0:taille2-1,nbscale-2-numscale),-1,/double) * N_ELEMENTS(trans.coef(0:taille2-1,nbscale-2-numscale))* gtilde
		recons = real_part(fft(tabfft,1,/double) / N_ELEMENTS(tabfft))
	
	endfor
END

ENDCASE


end
