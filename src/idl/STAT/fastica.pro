;-------------------------------------------------------------------------------
;+
; NAME:
;	fastica
;	
; PURPOSE:
;	This is a very simple implementation of the fastica algorithm.
;   This code is a simplified adaptation of the code in Matlab freely
;    available at www.cis.hut.fi/projects/ica/fastica/
;
;    It uses a cubic non linearity. Some tuning of the parameters may be necessary to adapt to a particular 
;    application.
;
; CALLING SEQUENCE:
;	fastica, ObservSig, NbSource, MixingMat,  DeMixingMat, Sources, nonlinearity = nonlinearity
;
; INPUTS:
;		ObservSig : mixing of input signal (ObservSig = A # input signals)
;               ObservSig(i,*) = ith vector
;       NbSource  : number of sources in input signal
;
; OPTIONAL INPUTS:
;	nonlinearity: string = 'pow3',  'tanh', 'gaus', 'skew'  or 'sym' . Default is 'pow3'.
;	 
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;  MixingMat : demixing matrix (ObservSig = MixingMat #  Sources)
;  DeMixingMat : demixing matrix  (Sources = DeMixingMat # ObservSig)
;  Sources : Reconstructed sources = DeMixingMat # ObservSig
;                Sources(i, *) = ith source, with i = 0. NbSource - 1
;
; MODIFICATION HISTORY:
;Yassir Moudden & Jerome Bobin
;-------------------------------------------------------------------------------   


;---------- Fonction de blanchiment ----------------


pro whitenv_data, ns, data, w_data, whitmat, dewhitmat


	eigen_values = EIGENQL( data#transpose(data), /double, EIGENVECTORS=eigen_vectors )

	ww = eigen_values(0:ns-1.)					;keeping only the ns highest eigenvalues and corresponding vectors
	vv = eigen_vectors(*, 0:ns-1.)

	whitmat = diag_matrix(1.0/sqrt(ww)) #  transpose(vv)

	dewhitmat = vv # diag_matrix(sqrt(ww))						

	w_data = whitmat#data

return
end


;--------- Centrer les donnees ----------------------

pro center_data, data, cdata

	sn = size(data)
	cdata = dblarr(sn(1),sn(2))
	for k=0,sn(1)-1 do cdata(k,*) = data(k,*)-mean(data(k,*))

return
end


;====================================================

pro fastica, x_0, ns, A_est, W, s_est, nonlinearity = nonlinearity

; +++++++ Options ++++++ 

myy =1.
epsi = 0.00001
nmax = 200.0
a2 = 1.0 
a1 = 100.

if not keyword_set(nonlinearity) then begin
	nl_used = 'pow3'
endif else begin
	nl_used = nonlinearity
endelse

mode = 'sym'

if (nl_used EQ 'pow3') then nl = 1
if (nl_used EQ 'tanh') then nl = 2
if (nl_used EQ 'gaus') then nl = 3
if (nl_used EQ 'skew') then nl = 4
if (mode EQ 'sym') then mde = 1 

n = (size(x_0))(2)			; nombre d'echantillons par observation
nc = (size(x_0))(1)			; nombre de canaux dans les observations initiales

A_est = dblarr(nc,ns)
B_est = dblarr(ns,ns)
B_old = dblarr(ns,ns)
B_old2 = dblarr(ns,ns)
s_est = dblarr(ns,n)

; +++++++ Centrer les donnes +++++++++++++++++++++++++++++
center_data, x_0, cx_0

; +++++++ Blanchiement des donnees et reduction de donnees +++++++++++++++++++++++
whitenv_data, ns, cx_0, wx_0, white_mat, dewhite_mat
		;at this point,  (size(wx_0))(1) is equal to the number of ICs sought. 

; +++++++ Initialisation de B  et d'autres parametres ++++++++++++++++++++++++++++
if (mde EQ 1) then begin
	
	; Either some initial guess (make sure whitening, orthogonality , etc  are done correctly)
	; B = 				

	; Take random orthonormal initial vectors.
	B_est = randomu(seed,ns,ns)
	eigenvalues = eigenql(transpose(B_est)#B_est, EIGENVECTORS = evecs,RESIDUAL = residual, /ASCENDING)
	B_est = B_est # ( evecs # diag_matrix(1.0/sqrt(eigenvalues)) # transpose(evecs) )

	;la_svd,B_est#transpose(B_est),sb,ub,vb							
	;B_est = vb#diag_matrix(1.0/sqrt(sb))#transpose(ub)#B_est
	
	;Take the identity
	;B_est = identity( ns)
	
	tol = 0.0
	tol2 = 0.0
endif

; +++++++ Algorithme +++++++++++++++++++++++++++++++++++++

W = 0.*white_mat
A = 0.*dewhite_mat

n_pas = 0.0
	
;while (n_pas LE nmax) && ( (1. - tol) GT epsi) do begin
while (n_pas LE nmax) do begin

		n_pas=n_pas+1
	
		if n_pas eq 10 then nl = nl+10
		if n_pas eq 10 then myy = 1

		B_old = B_est
		
		CASE nl OF

		1: begin
			B_est = (wx_0#(transpose(wx_0)#B_est)^3.0)/double(n) -3.0*B_est
			end	
			

		11 : begin
				Y = transpose(wx_0)#B_est
				y3 = Y^3
				bb = total(Y * y3, 1)
				;dd = diag_matrix(1.0/ (bb - 3. * n))
				dd = diag_matrix(1.0 / bb)
				B_est = B_est + myy * ( B_est # dd # (transpose(Y) # y3 - diag_matrix(bb)) )
			end

		2: begin 
			 uns = dblarr(ns,1)+1.
			 hypTan =tanh(a1*(transpose(wx_0)#B_est)) 
			 mult = uns # total(1.0 - hypTan*hypTan, 1)
			 B_est = ( wx_0 # hypTan - a1*(mult*B_est)  )/double(n)
			end
		
		12 : begin
				Y = transpose(wx_0)#B_est
				hypTan = tanh(a1 * Y)
				bb = total(Y * hypTan, 1)
				;dd = diag_matrix(1.0/ (bb - a1 * total(1 - hypTan^2, 1) ))
				dd = diag_matrix(1.0/ bb)				
				B_est = B_est + myy * ( B_est # dd #(transpose(Y) # hypTan - diag_matrix(bb))  )
			end	



		3:begin 
			uns = dblarr(ns,1)+1.   
			argu = transpose(wx_0)#B_est
			argu_car = argu^2.0
			expo = exp(-0.5*a2*argu_car)
			gosse = argu*expo
			dgosse = (1.0 - a2*argu_car)*expo
			B_est = (  wx_0#gosse - (uns#total(dgosse,1))*B_est )/ double(n)
		end
		
		13 : begin
				Y = transpose(wx_0)#B_est
				y2 = Y^2.0
				expo = exp(-0.5*a2*y2)
				gosse = Y*expo
				bb = total(Y * gosse,1 )
				;dd = diag_matrix(1.0/ (bb - total( (1. - a2 * y2 ) * expo, 1  )))
				dd = diag_matrix(1.0/ bb) 
				B_est = B_est + myy * ( B_est # dd # (transpose(Y) # gosse - diag_matrix(bb)) )
			end

		4: B_est = (wx_0 #(transpose(wx_0)#B_est)^2.0)/double(n)

		14 : begin
				Y = transpose(wx_0)#B_est
				Gskew = Y ^ 2
				bb = total(Y * Gskew, 1)
				dd = diag_matrix(1.0/ bb )
				B_est = B_est + myy * ( B_est # dd  # (transpose ( Y )  # Gskew - diag_matrix(bb) ) )
			end

		ENDCASE

		; orthogonalisation
		;la_svd,transpose(B_est)#B_est,sb,ub,vb
		;B_est = B_est # ( transpose(ub) # diag_matrix(1. / sqrt(sb)) # vb)
		eigenvalues = eigenql(transpose(B_est)#B_est, EIGENVECTORS = evecs,RESIDUAL = residual, /ASCENDING)
		B_est = B_est # ( evecs # diag_matrix(1.0/sqrt(eigenvalues)) # transpose(evecs) )
				
		tol = min(abs(diag_matrix(transpose(B_est)#B_old)))
		
		;myy = 0.9*myy
		Y = transpose(wx_0)#B_est
		; plot, y(*,0), y(*,1), psym = 3
		
	endwhile
	
	
	W = transpose( B_est ) # white_mat
	A_est = dewhite_mat#B_est

	s_est = W#x_0

return
end



