PRO GET_SPEC,rnu,spec

ImG=rim('ima_cmb.fits')
IM_ISOSPEC,ImG,spec,x=rnu,/nolog

END

FUNCTION C,y,N,rnu,spec

nuy = y/2.0/( N*sqrt(2.0)/2.0 )

jmax=size(rnu)
jmax=jmax[1]-1

loc = where( rnu GE nuy , c )

if c NE 0 then j=loc[0] else j=jmax

return, N*N*0.45*spec[j]

END

PRO CMB2D_SIM,F,V,N=N,rnu=rnu,spec=spec

;Genere une image carre NxN isotrope par l'intermediaire de sa transforme de Fourier
;
;Entrees:	N, entier, dimension de l'image
;			spec, tableau de float, spectre source de cmb spec[k] <==> spec( rnu[k] ) 
;			rnu, tableau de float, liste des frequences reduites
;
;Sorties:	F, tableau de float NxN, image
;			V, tableau de float NxN, carte des variances dans le domaine de Fourier

if not keyword_set(N) then N=128.

N=double(N)
U=dcindgen(N,N)
V=findgen(N,N)

if (not keyword_set(spec)) or (not keyword_set(rnu)) then GET_SPEC,rnu,spec

for I=0.,N-1 do begin
	for J=0.,N-1 do begin
  		K=sqrt( (I-N/2)^2 + (J-N/2)^2 );definition de l'image isotrope via sa FFT
  		H=C(K,N,rnu,spec)
  		RE=randomn(seed)*sqrt(H)
  		IM=randomn(seed)*sqrt(H)
  		U[I,J]=complex(RE,IM)
  		V[I,J]=H;Diag_Sigma_I2 domaine Fourier
	endfor
endfor

F=real_part(dfti(U))

F=F-mean(F)

END