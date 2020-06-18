pro mrsp_tqu2teb,tqu,teb,lmax=lmax
teb = tqu
mrsp_almtrans,tqu,alm,lmax=lmax

alm_e = alm
alm_e.alm(*,*,0) = alm_e.alm(*,*,1)
mrs_almrec,alm_e,e

alm_b = alm
alm_b.alm(*,*,0) = alm_b.alm(*,*,2)
mrs_almrec,alm_b,b

teb(*,1) =e
teb(*,2) =b
end

;============= 

function tqu2teb, D,lmax=lmax
mrsp_tqu2teb,D, teb,lmax=lmax
return, teb
end

;==============================

pro mrsp_teb2tqu,teb,tqu,lmax=lmax
mrs_almtrans,teb(*,0),alm_t,lmax=lmax
mrs_almtrans,teb(*,1),alm_e,lmax=lmax
mrs_almtrans,teb(*,2),alm_b,lmax=lmax
nalm = (size(alm_t.Alm))[1]
alm = dblarr(nalm,2,3)
alm(*,*,0) = alm_T.alm
alm(*,*,1) = alm_e.alm
alm(*,*,2) = alm_b.alm
; help,alm
out = {PixelType:alm_e.PixelType, tab: alm_e.tab, complex_alm: alm_e.complex_alm, nside : alm_e.nside, npix:alm_e.npix, ALM : Alm, norm: alm_e.norm, NormVal: alm_e.NormVal,lmin:alm_e.lmin,lmax:alm_e.lmax, TabNbrM: alm_e.TabNbrM, index:alm_e.index }
mrsp_almrec,out,tqu
end
;============= 

function teb2tqu, D,lmax=lmax
mrsp_teb2tqu, D, tqu,lmax=lmax
return, tqu
end

;============= 

function mrsp_gface, Map, Face
S1 = H2F(map[*,0])
S2 = H2F(map[*,1])
S3 = H2F(map[*,2])
sz = size(S1)
n1= sz[1]
n2= sz[2]
temp = S1[*,*,face]
Mxy = fltarr(n1,n2, 3)
Mxy[*,*,0] = S1[*,*,face]
Mxy[*,*,1] = S2[*,*,face]
Mxy[*,*,2] = S3[*,*,face]
return, Mxy
end

pro mrs_pface, Map, Face, FaceNumber
S1 = H2F(map[*,0])
S2 = H2F(map[*,1])
S3 = H2F(map[*,2])
sz = size(S1)
n1= sz[1]
n2= sz[2]
S1[*,*, FaceNumber] = Face[*,*,0]
S2[*,*, FaceNumber] = Face[*,*,1]
S3[*,*, FaceNumber] = Face[*,*,2]
map[*,0] = F2H(S1)
map[*,1] = F2H(S2)
map[*,2] = F2H(S3)
end


function mrsp_get_quface, Map, Face
S2 = H2F(map[*,1])
S3 = H2F(map[*,2])
sz = size(S2)
n1= sz[1]
n2= sz[2]
temp = S2[*,*,face]
Mxy = fltarr(n1,n2, 2)
Mxy[*,*,0] = S2[*,*,face]
Mxy[*,*,1] = S3[*,*,face]
return, Mxy
end

pro mrsp_put_quface, Map, Face, FaceNumber
S2 = H2F(map[*,1])
S3 = H2F(map[*,2])
sz = size(S2)
n1= sz[1]
n2= sz[2]
S2[*,*, FaceNumber] = Face[*,*,0]
S3[*,*, FaceNumber] = Face[*,*,1]
map[*,1] = F2H(S2)
map[*,2] = F2H(S3)
end



