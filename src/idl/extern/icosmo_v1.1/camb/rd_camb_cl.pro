pro rd_camb_cl, fid,cosmo, cl, file = filecl

;Written by Anais Rassat, April 2009
;Find out current directory & change to directory where CAMB stuff is
cambpath = fid.calc.cambpath
cd, cambpath, current=old_dir   ;old_dir is the directory in which we are before changing to the CAMB directory
readcol,filecl, l, cll2
cl = cll2/l/(l+1)*2.d*!dpi
;Return to current directory
cd, old_dir

;Now make cl structure
ind = where(l ge fid.calc.l_ran[0] and l le fid.calc.l_ran[1])
l = l[ind]
cl = cl[ind]

cl = {l:l, cl:cl}
end
