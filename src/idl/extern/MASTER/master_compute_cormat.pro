function master_compute_cormat,mll,matp,matq

matmult=mll*0.

for i=0,n_elements(mll(*,0))-1 do matmult(i,*)=mll(i,*)
kmat=(matp#matmult)#matq

return,kmat
end
