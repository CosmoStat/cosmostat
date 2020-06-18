pro master_coupling_matrix,ellw,well,lmax,deltal, $
                           mll,ell,ellbins,matp,matq, $
                           pqonly=pqonly,bintabmax=bintabmax, $
                           bintabmin=bintabmin, csbin=csbin


if not keyword_set(pqonly) then begin
;;;; Then compute the Mll matrix
    mll=master_make_mll(ellw,well,lmax)
    ell=ellw(0:lmax)
endif

;;;; Then compute the P and Q matrices and the binned ell
master_make_pq, deltal, lmax, matp, matq, ellbins, bintabmax=bintabmax, bintabmin=bintabmin, csbin=csbin
end
