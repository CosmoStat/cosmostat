;###################################################################################

function mrs_bandpass_sigma, Cl, Filter,  lmax=lmax, npix=npix, NormVal= NormVal

if not keyword_set(NormVal) then  NormVal = 1.
vs = size(filter)
lm = vs[1]
vs = size(Cl)
lm1 = vs[1]
if not keyword_set(lmax) then lmax = min([lm, lm1]) -1
 
 CoefN = NormVal* NormVal
 IndL = lindgen(LMAX+1)
IndL = 2. * IndL + 1.
Sig = total( IndL * Cl[0:Lmax]  * Filter[0:Lmax]^2 / CoefN)
Sig = Sig / (4. * !DPI)
return, sqrt(Sig)
end

;###################################################################################

function cmb_qualitymap,map,noisemap,nside=nside,BS=BS,NbrScale=NbrScale,Cl=Cl,coef=coef, TakeMin = TakeMin,Crit=Crit,LocalRMS=LocalRMS,Proba = Proba, lmax=lmax

if not keyword_set(Cl) then print,'We should estimate the power spectrum'

if not keyword_set(Crit) then Crit = 1

print,Crit

mrs_wttrans,map,out,NbrScale=NbrScale, lmax=lmax
c = out.coef
mrs_wttrans,noisemap,out,NbrScale=NbrScale, lmax=lmax
cn = out.coef
T = out.TabPsi

BSl = BS

coef = 0.*c
coef[*,NbrScale-1] = 1.

LocalRMS = 0.*c
Proba = LocalRMS

for s=0,NbrScale-2 do begin

	temp = c[*,s]
	tempn = cn[*,s]
	
	lrms = LocalRMS[*,s]
	
	CMBVar = mrs_bandpass_sigma(Cl, T[*,s])^2.
	
	for facenum=0,11 do begin ;--scan all faces
	
		f = get_one_face(temp,facenum)
		fn = get_one_face(tempn,facenum)
		g = 0.*f
		
		frms = get_one_face(lrms,facenum)
		
		nb_bin = double(nside)/double(BSl)
				
		for bx = 0,nb_bin-1. do begin ;--- scan bx
			for by = 0,nb_bin-1. do begin
			
				case Crit of
				
				1: begin
					p = max([CMBVar,sigma(f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^2. - sigma(fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^2.])
					g[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = CMBVar/p
					end
					
				2:  begin
					p = sigma(f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^2.
					g[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = (CMBVar+ sigma(fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^2.)/p
					end
			
				3:  begin
					p = total((f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]-mean(f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]))^4.)
					g[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = 3.*(CMBVar^2.+ sigma(fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^4.)/p   ;---- I should be 1 if Gaussian and noise + CMB only
					end
					
				4:  begin
					p = f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1]     ;----- Based on the l1 norm
					g[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = total(abs(p))/total(p*p)
					end
			
				endcase
				
				frms[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1] = sqrt(max([0.,sigma(f[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^2. - sigma(fn[bx*BSl:(bx+1)*BSl-1, by*BSl:(by+1)*BSl-1])^2.]))
				
			endfor
				
		endfor ;--- patch loops
	
		put_one_face,temp,g,facenum
		put_one_face,lrms,frms,facenum

	endfor ;--- face loop

	lpro = lrms
	mpro = mad(lpro)
	lpro = (lpro - median(lpro))/mpro
	lpro = 0.5*(1. - erf(lpro))

	coef[*,s] = temp
	LocalRMS[*,s] = lrms
	Proba[*,s] = lpro
		
	BSl = 2.*BSl

endfor

qmap = 0.*reform(coef[*,0])

if not keyword_set(TakeMin) then begin

	for s = 0,NbrScale-1 do qmap = qmap + coef[*,s] 
	qmap  = qmap/NbrScale

endif else begin

	qmap =coef[*,0]
	for s = 1,NbrScale-1 do qmap = min([[qmap],[coef[*,s]]],dim=2)

endelse

enside = nside/BS
					
filter = getidealbeam(0.*dblarr(4000.),lmin=0.,lmax=2.*enside,/tozero)

mrs_convol, qmap, Filter, conv
qmap = conv				
; qmap = convol_maps(qmap,filter)

return, qmap

end
