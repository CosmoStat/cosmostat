

pro testpol


restore, /verb, 'map_cmb_nside256.xdr'
mollview, /on, /nest, pol=3, c, title='CMB'

restore, /verb, 'map_synchrotron_70Ghz_nside256.xdr'
mollview, /on, /nest, pol=3, is, title='Synchrotron 70 Ghz'

restore, /verb, 'map_thermal_dust_857Ghz_nside256.xdr'
ima = i1
mollview, /on, /nest, pol=3, i1, /log, vector_scale=1.e-20,  title='Dust 857 Ghz'

t =  mrsp_gface(ima, 4)
; t[*]= 0
i2 = ima
mrs_pface, i2, t, 4

ImaQU=  mrsp_get_quface(ima, 4)
; t[*]= 0
i2 = ima
mrsp_put_quface, i2, t, 4

NbrScale=5
Pol = t
mr_modphase_dwt_trans,  ImaQU,Trans, NbrScale= NbrScale



;================ FIG BACKPROJECT 2D OWT - manifold 
t = trans
x=4
y=4
z=0
t= trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 10
t.angCOEFF[x,y,z] =complex(0,1)
mr_modphase_dwt_rec,t,Rec
; tvxy, rec , skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_x3_y4.png',  /nobar, colt=20


x = 12
y = 11
z=0
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 10
t.angCOEFF[x,y,z] = complex(10,10)
mr_modphase_dwt_rec,t,Rec
; tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=1., png='fig_backproj_x12_y11.png', /nobar, colt=20

x = 25
y = 7
z=0
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 100
t.angCOEFF[x,y,z] = complex(30,1)
mr_modphase_dwt_rec,t,Rec
; tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=30., png='fig_backproj_x25_y7.png', /nobar, colt=20


x = 97
y = 32
z=0
a = 90.
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(cos(a),sin(a))
mr_modphase_dwt_rec,t,Rec
; tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=1., png='fig_backproj_x97_y32.png', /nobar, colt=20


; pas terrible
x = 32
y = 97
z=0
a = 90.
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(cos(a),sin(a))
mr_modphase_dwt_rec,t,Rec
; tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=1., png='fig_backproj_x32_y97.png', /nobar, colt=20



a=130.
x = 16
y = 48
z=0
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(cos(a),sin(a))
mr_modphase_dwt_rec,t,Rec
; tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000, png='fig_backproj_x16_y48.png', /nobar, colt=20,  vector_scale=0.8

a=0
x = 48
y = 16
z=0
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.angCOEFF[x,y,z] = complex(cos(a),sin(a))
t.MODCOEFF[x,y,z] = 3
mr_modphase_dwt_rec,t,Rec
;tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000, png='fig_backproj_x48_y16.png', /nobar, colt=20,  vector_scale=0.8

a=45
x = 48
y = 48
z=0
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.angCOEFF[x,y,z] = complex(cos(a),sin(a))
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = 20
mr_modphase_dwt_rec,t,Rec
;tvxy, rec, skip=2, length_factor=1
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000, png='fig_backproj_x48_y48.png', /nobar, colt=20,  vector_scale=0.8


;================ FIG BACKPROJECT 2D OWT - QU

t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x=7
y=8
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(-10, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale-1, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale-1, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=1., png='fig_backproj_qu_x7_y8.png', colt=-1


t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x=12
y=12
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(-10, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_qu_x12_y12.png', colt=-1

t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 97
y = 32
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(-10, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_qu_x97_y32.png', colt=-1


t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 16
y = 48
z=0
t.MODCOEFF[x,y,z] = 3
; t.angCOEFF[x,y,z] = complex(-10, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_q_x16_y48.png', colt=-1

t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 48
y = 48
z=0
t.MODCOEFF[x,y,z] = 3
; t.angCOEFF[x,y,z] = complex(-10, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_q_x48_y48.png', colt=-1

t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 48
y = 16
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(0, 0)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_q_x48_y16.png', colt=-1


t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 16
y = 48
z=0
t.angCOEFF[x,y,z] = complex(-10, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_u_x16_y48.png', colt=-1

t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 48
y = 48
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(40, 20)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_u_x48_y48.png', colt=-1

t = trans
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(0,0)
x = 48
y = 16
z=0
t.angCOEFF[x,y,z] = complex(-40, 0)
R1 = t.MODCOEFF
R1 =  bwt01_lift(R1, NbrScale, /inverse)
R2 = float(t.angCOEFF)
R2 =  bwt01_lift(R2, NbrScale, /inverse)
rec = fltarr(256,256,2)
Rec[*,*, 0] = R1
Rec[*,*, 1] = R2
tvxy, rec, skip=3, length=3
i2[*]=0
mrsp_put_quface, i2, Rec, 4
gnomview, i2, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.5, png='fig_backproj_u_x48_y16.png', colt=-1



 
;================ FIG BACKPROJECT 2D OWT - QU


;================== FIG DUST + Noise ==============
n = randomn(seed, 256L^2*12, 3) * sigma(ima) / 50
in = ima + n
mollview, /on, /nest, pol=3, in, /log, vector_scale=1.e-20,  title='Dust 857 Ghz + noise'

;================== FIG BACKPROJ E/B Wavelet ==============


mrsp_wttrans, Ima, WT, NbrScale=NbrScale, lmax=lmax, /DifInSH, MeyerWave=MeyerWave
W= wt
w.coef[*]=0
j = 3
w.coef[311296L, j, 2] = 100.
w.coef[311296L, j, 0] = 100.
mrsp_wtrec, w, rb
mollview, /on, /nest, pol=3, rb, vector_scale=0.001,  title='B-Wavelet Coefficient  Backprojection', png='fig_mol_backproj_ebwt_bj3.png'


gnomview, rb, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.001, title='B-Wavelet Coefficient  Backprojection', png='fig_backproj_ebwt_bj3.png', reso=20

w.coef[*]=0
j = 3
w.coef[311296L, j, 1] = 100.
w.coef[311296L, j, 0] = 100.
mrsp_wtrec, w, re
mollview, /on, /nest, pol=3, re, vector_scale=0.001,  title='E-Wavelet Coefficient  Backprojection (j=4)', png='fig_mol_backproj_ebwt_ej3.png'
gnomview, re, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.002, png='fig_backproj_ebwt_ej3.png', title='E-Wavelet Coefficient Backprojection (j=4)', reso=20


w.coef[*]=0
j = 2
w.coef[311296L, j, 2] = 100.
w.coef[311296L, j, 0] = 100.
mrsp_wtrec, w, rb
mollview, /on, /nest, pol=3, rb, title='B-Wavelet Coefficient  Backprojection (j=3)', png='fig_mol_backproj_ebwt_bj2.png',vector_scale=0.005
w.coef[*]=0
j = 2
w.coef[311296L, j, 1] = 100.
w.coef[311296L, j, 0] = 100.
mrsp_wtrec, w, re
mollview, /on, /nest, pol=3, re,  title='E-Wavelet Coefficient  Backprojection (j=3)', png='fig_mol_backproj_ebwt_ej2.png', vector_scale=0.005

 
w.coef[*]=0
j = 1
w.coef[311296L, j, 2] = 100.
w.coef[311296L, j, 0] = 100.
mrsp_wtrec, w, rb
mollview, /on, /nest, pol=3, rb, title='B-Wavelet Coefficient  Backprojection (j=2)', png='fig_mol_backproj_ebwt_bj1.png',vector_scale=0.01
w.coef[*]=0
j = 1
w.coef[311296L, j, 1] = 100.
w.coef[311296L, j, 0] = 100.
mrsp_wtrec, w, re
mollview, /on, /nest, pol=3, re,  title='E-Wavelet Coefficient  Backprojection (j=2)', png='fig_mol_backproj_ebwt_ej1.png', vector_scale=0.05
 





W= wt
w.coef[*]=0
j = 4
w.coef[311296L, j, 2] = 100.
mrsp_wtrec, w, rb
mollview, /on, /nest, pol=3, rb, vector_scale=0.1,  title='B-Coarse Coefficient  Backprojection', png='fig_mol_backproj_ebwt_bj4.png'

gnomview, rb, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.2, title='B-Coarse Coefficient  Backprojection', png='fig_backproj_ebwt_bj4.png'

w.coef[*]=0
j = 4
w.coef[311296L, j, 1] = 100.
mrsp_wtrec, w, re
mollview, /on, /nest, pol=3, re, vector_scale=0.1,  title='E-Coarse Coefficient  Backprojection', png='fig_mol_backproj_ebwt_ej4.png'
gnomview, re, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.2, png='fig_backproj_ebwt_ej4.png', title='E-Coarse Coefficient Backprojection'

NbrScale=7
mrsp_wttrans, Ima, WT, NbrScale=NbrScale, lmax=lmax, /DifInSH, MeyerWave=MeyerWave
W= wt
w.coef[*]=0
j = 5
w.coef[311296L, j, 2] = 100.
mrsp_wtrec, w, rb
mollview, /on, /nest, pol=3, rb, vector_scale=0.01,  title='B-Wavelet Coefficient  Backprojection', png='fig_mol_backproj_ebwt_bj5.png'

gnomview, rb, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.01, title='B-Wavelet Coefficient  Backprojection', png='fig_backproj_ebwt_bj5.png'

w.coef[*]=0
j = 5
w.coef[311296L, j, 1] = 100.
mrsp_wtrec, w, re
mollview, /on, /nest, pol=3, re, vector_scale=0.01,  title='E-Wavelet Coefficient  Backprojection', png='fig_mol_backproj_ebwt_ej5.png'
gnomview, re, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.01, png='fig_backproj_ebwt_ej5.png', title='E-Wavelet Coefficient Backprojection'


;================== FIG BACKPROJ Wavelet ==============

NbrScale=7
irec = ima
irec[*] = 0
mrs_wttrans, Ima[*,1], WT, NbrScale=NbrScale, lmax=lmax, /DifInSH, MeyerWave=MeyerWave
wt.coef[*]=0
j=5

wt.coef[311296L, j] = 100.
mrs_wtrec, wt, rec
irec[*,1] = rec
irec[*,0] = rec

mollview, /on, /nest, pol=3, irec, vector_scale=0.01 


wt.coef[*]=0
j=2
wt.coef[311296L, j] = 100.
mrs_wtrec, wt, rec
irec[*]=0
irec[*,1] = rec
irec[*,0] = rec
mollview, /on, /nest, pol=3, irec, vector_scale=0.01,  png='fig_mol_backproj_quwt_qj2.png', title='Q-Wavelet Coefficient Backprojection'


gnomview, irec, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.1,  png='fig_backproj_quwt_qj2.png', title='Q-Wavelet Coefficient Backprojection'

irec1 = irec

irec[*] = 0
irec[*,2] = rec
irec[*,0] = rec
mollview, /on, /nest, pol=3, irec, vector_scale=0.01,  png='fig_mol_backproj_quwt_uj2.png', title='U-Wavelet Coefficient Backprojection'

gnomview, irec, /on, /nest, pola=3, pxsize=1000,  vector_scale=0.1,  png='fig_backproj_quwt_uj2.png', title='U-Wavelet Coefficient Backprojection'

irec2 = irec

 
;================== FIG BACKPROJ E/B Curvelet ==============

mr1ok=1
z = fltarr(256L^2*12)
z1 = z
mrs_curtrans,z1,tzs,  nbrscale=6, /overlap, /ssr
tzs.ridscale4.coef[64,64,4]=1
mrs_currec, tzs, rs
tvs, rs
gnomview, rs, /nest, /on,   vect=0.005, RESO_ARCMIN=10

 
rr=fltarr(256L^2*12,3)
rsp_trans, rr, teb, /ebdec
teb.dec2[*]=0
teb.dec3[*]=rs
mrsp_rec, teb, rec1
rec1[*,0] = rs
;rec1[*,0] = 0
gnomview, rec1, /nest, /on, pol=3, vect=0.005, RESO_ARCMIN=5, title='B-Curvelet Coefficient  Backprojection', png='fig_backproj_ebcur_bj3.png'

mollview, rec1, /nest, /on, pol=3, vect=0.001,  title='B-Curvelet Coefficient  Backprojection', png='fig_mol_backproj_ebcur_bj3.png'


rr=fltarr(256L^2*12,3)
rsp_trans, rr, teb, /ebdec
teb.dec2[*]=rs
teb.dec3[*]=0
mrsp_rec, teb, rec2
rec2[*,0] = rs
;rec2[*,0] = 0

gnomview, rec2, /nest, /on, pol=3, vect=0.005, RESO_ARCMIN=5, title='E-Curvelet Coefficient  Backprojection', png='fig_backproj_ebcur_ej3.png'

mollview, rec2, /nest, /on, pol=3, vect=0.003,  title='E-Curvelet Coefficient  Backprojection', png='fig_mol_backproj_ebcur_ej3.png'

;==============FIG DENOISING ==============

restore, /verb, 'map_thermal_dust_857Ghz_nside256.xdr'
ima = i1
ima = ima / sigma(ima) * 100.
v = ima
v[*,0] = (v[*,0]>0) < 200
mollview, /on, /nest, pol=3, v, vector_scale=4,  title='Dust 857 Ghz', png='fig_dust.png'
mollview, /on, /nest, pol=3, ima, vector_scale=4,  /log, title='Dust 857 Ghz', png='fig_logdust.png'

n = randomn(seed, 256L^2*12, 3) * sigma(ima) / 10
in = ima + n
v = in
v[*,0] = (v[*,0]>0) < 200
mollview, /on, /nest, pol=3, v, title='Dust 857 Ghz + noise',vector_scale=8, png='fig_noisydust.png'
v[*,0] = ima[*,0]
mollview, /on, /nest, pol=3, v, /log, title='Dust 857 Ghz + noise',vector_scale=8, png='fig_lognoisydust.png'


mrsp_trans, /uwt, in, trans, /ebdec, /overlap, nbrscale=6
curt = trans
mrsp_threshold, curt, NSigma=5, Mad=Mad 
mrsp_rec, curt, rec
v = rec
v[*,0] = (v[*,0]>0) < 200
mollview, /on, /nest, pol=3, v, title='E-B Wavelet filtering',vector_scale=4, png='fig_ebwtfilter_noisydust.png'
v[*,0] = ima[*,0]

mollview, /on, /nest, pol=3, v, /log, title='E-B Wavelet filtering',vector_scale=4, png='fig_logebwtfilter_noisydust.png'



mrsp_trans, /uwt, in, trans,  /overlap, nbrscale=6
curt = trans
mrsp_threshold, curt, NSigma=5, Mad=Mad 
mrsp_rec, curt, rec
v = rec
v[*,0] = (v[*,0]>0) < 200
mollview, /on, /nest, pol=3, v, title='Q-U Wavelet filtering',vector_scale=4, png='fig_quwtfilter_noisydust.png'
v[*,0] = ima[*,0]
mollview, /on, /nest, pol=3, v, /log, title='Q-U Wavelet filtering',vector_scale=4, png='fig_logquwtfilter_noisydust.png'



mrsp_trans, /cur, in, trans, /ebdec, /overlap, nbrscale=6
curt = trans
mrsp_threshold, curt, NSigma=5, Mad=Mad 
mrsp_rec, curt, rec
v = rec
v[*,0] = (v[*,0]>0) < 200
mollview, /on, /nest, pol=3, v, title='E-B Curvelet Denoising',vector_scale=4, png='fig_ebcurfilter_noisydust.png' 
v[*,0] = ima[*,0]
mollview, /on, /nest, pol=3, v, /log, title='E-B Curvelet Denoising',vector_scale=4, png='fig_logebcurfilter_noisydust.png' 

mrsp_trans, /cur, in, trans,  /overlap, nbrscale=6
curt = trans
mrsp_threshold, curt, NSigma=5, Mad=Mad 
mrsp_rec, curt, rec
v = rec
v[*,0] = (v[*,0]>0) < 200
mollview, /on, /nest, pol=3, v, title='Q-U Curvelet Denoising',vector_scale=4, png='fig_qucurfilter_noisydust.png' 
v[*,0] = ima[*,0]
mollview, /on, /nest, pol=3, v, /log, title='Q-U Curvelet Denoising',vector_scale=4, png='fig_logqucurfilter_noisydust.png' 

;============= FIG PHASE Decomposition =======================

restore, /verb, 'map_thermal_dust_857Ghz_nside256.xdr'
ima = i1
mollview, /on, /nest, pol=3, i1, /log, vector_scale=1.e-20,  title='Dust 857 Ghz'
Pol =  mrsp_gface(ima, 4) * 1.e20
ta = alog(Pol[*,*,0])
; pa = imaginary(alog(z/abs(z)))

disp,win=0
disp, win=1
disp, win=2

NbrScale=5
; mr_modphase_dwt_trans, Pol,Trans, NbrScale= NbrScale


P1 =  Pol[*,*,1:2]  + randomn(seed, 256, 256,2) * 2.
wset, 0
tvxy, ima=ta,  Pol[*,*,1:2], /pol, skip=5, len=0.2

wset, 1
tvxy, ima=ta,  P1, /pol, skip=5, len=0.2

mr_modphase_uwt_trans, P1,  Tu, NbrScale= NbrScale
tu1 = tu
tu1.angcoeff[*,*,0:2] = 0
; tu1.modcoeff[*,*,0:2] = 0

mr_modphase_uwt_rec,  Tu1,  Ra1 

wset,2
tvxy, ima=ta, Ra1, /pol, skip=5, len=0.2


wset, 0
tvxy, Pol[*,*,1:2], /pol, skip=5, len=0.2
wset, 1
tvxy,  P1, /pol, skip=5, len=0.2
wset,2
tvxy, Ra1, /pol, skip=5, len=0.2


i = readfits('imapol256_face5.fits')
r = rim('recpol.fits')
ta = alog(i[*,*,0])
wset, 0
clear
tvxy, i[*,*,1:2], /pol, skip=5, len=20; , ima=ta
wset,1
tvxy, r[*,*,*,3], /pol, skip=5, len=20 ;, ima=ta


;============= TRASH =======================


mrsp_tvf, ima, 6, /log, length=30, /deg
a = mrsp_gface(ima, 6)
mrp_tvf, a, /log, length=30, /deg

NbrScale=5
Pol = a[*,*,1:2]
mrp_owttrans, Pol,Trans, NbrScale= NbrScale
mrp_owtrec,Trans,Rec

; mrp_wttrans, Pol, Trans, NbrScale=NbrScale
; mrp_wtrec, Trans,Rec

ar = a
ar[*,*,1:2] = Rec
mrp_tvf, ar, /log, length=30, /deg

t = trans
x = 11
y = 12

x=5
y=5
z=0
t.MODCOEFF[*] = 0
t.angCOEFF[*] = complex(1,0)
t.MODCOEFF[x,y,z] = 0.1
t.angCOEFF[x,y,z] = complex(10,10)

x = 11
y = 12
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(10,10)

; t.MODCOEFF[*] = 0
; t.angCOEFF[*] = complex(1,0)
x = 71
y = 32
z=0
t.MODCOEFF[x,y,z] = 5
t.angCOEFF[x,y,z] = complex(40, -20)

x = 71
y = 81
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(-40, 20)



;t.MODCOEFF[*] = 0
; t.angCOEFF[*] = complex(1,0)
x = 27
y = 3
z=0
t.MODCOEFF[x,y,z] = 3
t.angCOEFF[x,y,z] = complex(-40, 2)

 
; mrp_wtrec, t,Rec
mrp_owtrec,t,Rec

info, rec

ar = a
ar[*,*,0] = 1
ar[*,*,1:2] = Rec
mrp_tvf, ar, /log, length=10, /deg


;=============
  

mrsp_wttrans, Ima, WT, NbrScale=NbrScale, lmax=lmax, /DifInSH, MeyerWave=MeyerWave
W= wt
w.coef[*]=0
j = 3
w.coef[100, j, 1] = 100.
mrsp_wtrec, w, r
tvs, r, pol=3

w.coef[*]=0
j = 3
w.coef[100, j, 2] = 100.
mrsp_wtrec, w, r1
tvs, r1, pol=3


mrs_wttrans, Ima[*,0], WT1, NbrScale=NbrScale, lmax=lmax, /DifInSH, MeyerWave=MeyerWave
W= wt1
w.coef[*]=0
d = r1[*,0]
h = h2f(d)
h[128,128,4]=100
d = f2h(h)

j = 2
w.coef[*, j] = d
mrs_wtrec, w, x
r2 =  r1*0.
r2[*,1] = x
hr2 = mrsp_gface(r2, 4)
tvxy, hr2[*,*,1:2], skip=3, length=3

end




pro sim_make_pol
i = rims('/dsm/cosmo01/planck/wg2/real_map/planck-wg2.planck.fr/Challenge-2/PSM-maps_v0/PSM-maps_857GHz/i_dust_857GHz_milliKantenna_nside2048.fits')

i1 = mrs_resize(i, nside=256)

 ip = fltarray(786432, 3)
  g = rims('POL/g_dust_2048_5arcmin.fits')
   gamma = rims('gamma_dust_2048_5arcmin.fits')
    g1  = mrs_resize(g, nside=256)
     gamma1  = mrs_resize(gamma, nside=256)
      q = 0.1 * g1 * i1 * cos(gamma1*2)
      u =  0.1 * g1 * i1 * sin(2.*gamma1)
       ip[*,1]= q
       ip[*,2]= u
       
       save, filename='poldust256.xdr', ip
       
       end
