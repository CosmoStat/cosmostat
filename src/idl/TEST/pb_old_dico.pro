
function ima2block, Image, BlockSize=BlockSize

B = -1
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'B = ima2block(Image,  BlockSize=BlockSize)'
   GOTO, CLOSING
 ENDIF
 
vsize = size(Image)
if vsize[0] NE 2 then BEGIN
   print, "Error: bad first parameter ..."
   print, "       Image has to be a 2 dimensional array"
   GOTO, CLOSING
 ENDIF
 
Nx = vsize[1]
Ny = vsize[2]

if keyword_set(BlockSize) then SizeBlock = BlockSize else SizeBlock = 8

Nxb = Nx / SizeBlock
Nyb = Ny / SizeBlock
Nb = Nxb*Nyb
TabBlock = fltarr(SizeBlock, SizeBlock, Nb)
TabVect = fltarr(Nb, SizeBlock*SizeBlock)
;------------------------------------------------------------
; function body
;------------------------------------------------------------

TabDebX = intarr(Nxb*Nyb)
TabDebY = intarr(Nxb*Nyb)
TabEndX = intarr(Nxb*Nyb)
TabEndY = intarr(Nxb*Nyb)

ind = 0
for i=0,Nxb-1 do BEGIN
for j=0,Nyb-1 do BEGIN
   Debx = i*SizeBlock
   Deby = j*SizeBlock
   Endx = (i+1)*SizeBlock-1
   Endy = (j+1)*SizeBlock-1
   TabBlock[*,*, ind] =  Image[Debx:Endx, Deby:Endy]
   TabVect[ind, *] = reform(Image[Debx:Endx, Deby:Endy], SizeBlock*SizeBlock)
   TabDebX[ind] = Debx
   TabDebY[ind] = Deby
   TabEndX[ind] = Endx
   TabEndY[ind] = Endy
   ind = ind + 1
   END
   END
 
 B = {Nx:Ny,  Ny: Ny, SizeBlock : SizeBlock , $
        Nxb: Nxb, Nyb:Nyb, Nb:Nb, $
        TabDebX: TabDebX, TabDebY: TabDebY, TabEndX: TabEndX, TabEndY: TabEndY, $
        TabBlock: TabBlock, TabVect:TabVect}


;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  return, B
 
 END
 ;===================================================================
 
 function block2ima, Block, norm=norm
Ima = -1
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ',   'Ima = block2ima(Block)'
   GOTO, CLOSING
 ENDIF
 
Ima = fltarr(Block.nx, Block.ny)
Nxb = Block.Nxb
Nyb = Block.Nyb
ind = 0
for i=0,Nxb-1 do BEGIN
for j=0,Nyb-1 do BEGIN
   Debx  = Block.TabDebX[ind]
   Deby  = Block.TabDebY[ind]
   Endx = Block.TabEndX[ind]
   Endy = Block.TabEndY[ind]
   B =  reform(  Block.TabVect[ind, *], Block.SizeBlock,   Block.SizeBlock)
   if keyword_set(norm) then BEGIN
     MinB = min(B)
     MaxB = Max(B)
     B = B - MinB
     if MaxB EQ MinB then B(*) = 1. $
     else  B = B  / (MaxB - MinB)
    END    
   Block.TabBlock[*,*, ind] = B
   Ima[Debx:Endx, Deby:Endy] =  B
   ind = ind + 1
   END
   END

 CLOSING:
 
  return, Ima
 
 END


pro tvatom, D, Z
  vs = size(D)
  na = vs(1)
  sizep= vs(2)
  w = sqrt(sizep)
  nna = floor(sqrt(na))
  BI = fltarr(nna * w, nna *  w)
  BIB = ima2block(BI,BlockSize=BlockSize)
  BIB.tabvect = D(0:nna*nna-1, *)
  Z = block2ima(BIB)
  load, z
end

  ;===================================================================
 function  perform_thresholding, x, t, type
; perform_thresholding - perform hard or soft thresholding
;
;   y = perform_thresholding(x, t, type);
;
;   type is either 'hard' or 'soft' or 'semisoft'
;   t is the threshold
;
;   if type is 'strict' then it keeps the t largest entry in each 
;   column of x.
;
;   Copyright (c) 2006 Gabriel Peyre

y = x
if  type EQ  'hard' then begin
   ind = where( abs(x) LT T, c)
   if c GT 0 then y(ind) = 0
end else if  type EQ  'soft' then begin
    s = abs(x) - t
    s = (s + abs(s))/2
    y = sign(x)*s
end else if  type EQ  'strict' then begin  
       vs = size(x)
       for i=0, vs(1)-1 do begin
          col = x(i,*)
          ind = sort(abs(col))
          v = abs(col(reverse(ind)))
          v = v(round(t)-1)
	       ind = where( abs(col) LT v, c)
	        if c GT 0 then col(ind) = 0
	        y(i, *) = col
	     end
end else begin
         print, 'Error: Unkwnown thresholding type.'
         y = -1
end

return, y
end

  ;===================================================================

 pro perform_dictionary_learning, Y,  D, X, E, niter_learning=niter_learning, learning_method=learning_method, redund=redund, thresh_type=thresh_type, lambda=lambda, lambda_min= lambda_min, lambda_max=lambda_max, coef_update_iter=coef_update_iter, mu=mu, init_ran_D=init_ran_D

;  function [D,X,E] = perform_dictionary_learning(Y,options)
; perform_dictionary_learning - learn a dictionnary using K-SVD algorithm
;
;   D = perform_dictionary_learning(Y,options)
;
;   Y is a matrix of size (n,m) containing m column examplar
;       vector in R^n.
;
;   D is a dictionnary matrix of size (n,K) of K vectors in R^n
;       that should approximate well the vectors of Y 
;       with few components.
;       K is given by options.K
;
;   The algorihtm perform the folowing optimization jointly on D and X
;       min_{D,X}  |Y-D*X|^2 + lambda * |X|_1
;           subject to columns of D having unit norm.
;   where |X|_1=\sum |X(i,j)| and lambda is given by options.lambda
;
;   It does this optimisation using block coordinate descent.
;       * Optimization of X knowing D amount to sparse coding, solved using
;       iterative thresholding:
;           X <- Thresh_lambda( X + D'*( Y - D*X ) )
;       This step is repeated options.coef_update_iter times.
;       * Optimization of D knowing X amount to L2 best fit:
;           D <- Y*X^+   where   X^+ = X'*(X*X')^{-1}
;
;   This algorithm is very much inspired by the MOD algorithm
;   for dictionary design.
;
;   The number of iteration of the algorithm is given in
;       options.niter_learning.
;   The number of iteration for sparse coding during each step of learning
;       is given by options.coef_update_iter.
;
;   options.options.learning_method can be set to:
;       'ksvd': in this case, the dictionary update is computed using the
;           K-SVD algorithm explained in 
;           "K-SVD: An Algorithm for Designing Overcomplete Dictionaries
;                   for Sparse Representation"
;           Michal Aharon Michael Elad Alfred Bruckstein, 2006
;       'mod': in this case, the dictionary update is computed using the L2
;           best fit as proposed by
;           "Method of optimal directions for frame design",
;           K. Engan, S. O. Aase, and J. H. Hus¿y,
;           in Proc. ICASSP, Phoenix, AZ, Mar. 1999, pp. 2443Ð2446. 
;       'randomized': apply randomly one or the other algorithm.
;
;   MOD is usualy faster but KSVD gives a better optimization of the energy.
;
;   Copyright (c) 2007 Gabriel Peyre


; n = WxW,  W = taille du patch 
; m = nombre de patch
vs = size(Y)
m = vs(1)
n = vs(2)

if keyword_set(niter_learning) then niter = niter_learning $
else niter = 30 

; K = number of atoms (number of rows) 
if keyword_set(redund) then k = double(redund*n) $
else k = double(2*n)

; coef_update_iter = niter for updating coefficients
if not keyword_set(coef_update_iter) then coef_update_iter = 20  
 
if not keyword_set(lambda) then  lambda = mean( sqrt( total( Y^2, 2 ) ) ) / 50.

; thresh_type = 'soft', 'hard', 'strict_sparsity'
if not keyword_set(thresh_type) then  thresh_type = 'strict'
  
; initialisation soit random, soit avec les patches
if not keyword_set(init_dico) then  init_dico = 'input'
 
 ; Sparsity level of the  estimated coefficients
 if not keyword_set(strict_sparsity) then  strict_sparsity = 5
 ; Descent step decreasing schedule: 1 ==> cst along inner iterations
 if not keyword_set(mu_damping) then  mu_damping = 1
 
 if not keyword_set(lambda_min) then  lambda_min = lambda
 if not keyword_set(lambda_max) then  lambda_max = lambda_min*4

if keyword_set(init_ran_D) then D = randomn(seed, K, n) $
else  begin
   D= Y
   if K LE m then D = D[0:long(k)-1, *]  $
   else D = [D, randomn(seed, K-m, n)]
end

for i=0,K-1 do D[i,*] = (D[i,*] - mean(D[i,*])) / sqrt( total(D[i,*]^2))
D[0,*] = 1. / sqrt(n)
 

 
 E = dblarr(niter)
 X = dblarr(m,K)
  
 print, ' niter = ', niter, ', coef_update_iter = ', coef_update_iter
 if  keyword_set(mu) then  print, 'mu'
 if  keyword_set(mu) then  MuVal = mu

disp, win=0
disp, win=1
 for iter=1,niter do BEGIN
     ;======================================
     ; Update coeff using BPDN
     ;======================================
     
        if not keyword_set(mu) then begin
    		LA_SVD, D, S, U, V, /DOUBLE
    		MuVal = 1. / max(S)^2
		end
        lambda = lambda_max - float(i-1)/float(niter-1)*(lambda_max-lambda_min)
    print, MuVal
     
       for inner=1,coef_update_iter do BEGIN
             X = X + MuVal * mu_damping * transpose(D)##( Y - D##X )
             if  thresh_type EQ  'strict' then  t = strict_sparsity $
             else if  thresh_type EQ  'soft' then  t = lambda*mu_damping* MuVal $
             ; for hard thresholding, the scaling is different
             else if thresh_type EQ   'hard' then t = lambda*sqrt(mu_damping* MuVal);
          
             X = perform_thresholding( X, t, thresh_type)
             
         END
        
         if thresh_type EQ 'strict' then  E(iter-1) = sqrt(total( (Y-D##X)^2)) $
        else E(iter-1) = 1/2*  total( (Y-D##X)^2)  + lambda * total( abs(X ) )
           
     ;======================================
     ; Update dictionary using MOD algorithm
     ;======================================
     ; dictionary fit
     ; D = Y * X'*(X*X')^(-1)
     pinvX = invert(transpose(X)#X,/DOUBLE)#transpose(X)
      D = Y## pinvX
     ; Normalize
    DD = sqrt( total(D^2., 2))
    ind = where( DD LT 1e-9, c)
    if c GT 0 then DD[ind] = 1.
    for i=0, c-1 do  D[ind[i], *] = randomn(seed, n)
    for i=0,K-1 do D[i,*] = D[i,*]  / sqrt( total(D[i,*]^2))
        
     ; sort the basis vector according to energy
     e = total( X^2, 1)
     ind = reverse(sort(e))
     TmpD = D
     TmpX = X
     for j=0,K-1 do BEGIN
            TmpD(j,*) = D(ind(j), *)
            TmpX(*, j) = X(*, ind(j))
     END
      X = TmpX
      D = TmpD
      wset, 0
      tvatom, D
      wset, 1
      load, x
      
END ; FOR 


END

  ;===================================================================

pro dico
BlockSize=8
i = rim('simu_sky.fits')
DIR='/Users/starck/Main/MRS/idl/TEST/dictionary-learning/images/'
FN = DIR + 'barb.png'
i = double(read_png(FN))
vs = size(i)
g = mygauss(vs(1), vs(2),15.)
g = g / total(g)
ig = conv(i, g)
load, i-ig
i1= i - ig
B = ima2block(i,BlockSize=BlockSize)
Y = B.TABVECT
perform_dictionary_learning, Y,  D, X, E
R = Y-D##X
tvatom, R 


vs = size(D)
na = vs(1)
sizep= vs(2)
w = sqrt(sizep)
BI = fltarr(floor(sqrt(na)) *  w, floor(sqrt(na)) *  w)
BIB = ima2block(BI,BlockSize=BlockSize)
BIB.tabvect = D(0:120, *)
Z = block2ima(BIB)


Z = block2ima(B)
hs, b
end
 