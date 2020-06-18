function master_pseudcl2cl,matkinv,matp,cellmap

print, '!!!!!!!!!! spectre non moyenne par bande en entree !!!!!!!!!!!!!'

specbincorr=matkinv#matp#cellmap

help, cellmap
help, matkinv

return,specbincorr
end
