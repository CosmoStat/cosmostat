pro test_planck
; Oct 08 - Written by Anais Rassat
; PURPOSE: Just  a little test routine to learn how to use the Planck Fisher matrix
;code

f = set_fiducial()
c = mk_cosmo(f)
sv = mk_survey(f, 'sv1')

print, 'Calculating Lensing Fisher Matrix'
lensf = mk_fisher_lens(f, sv)

print, 'Calculating BAO Fisher Matrix'
fbao = set_fiducial('test2')
svbao = mk_survey(fbao, 'sv1')
baof = mk_fisher_bao(fbao, svbao)


print, 'Calculating Planck Fisher Matrix'
;Following line needed if you have changed the cosmology or want a fresh
;Planck Fisher matrix Calculation
cambpath = '/Users/arassat/Work/icosmo_v1.1/camb/CAMB_jw_Planck/CAMB/'
planckf = mk_fisher_planck(f,c, cambpath=cambpath)
;Following line is needed if you have previously calculated the Planck Fisher
;matrix and just want to read the output file
;rd_planck_fisher, c,planckf, file=file


test=comb_fisher(lensf, planckf)
testbao=comb_fisher(baof, planckf)
margin_fisher, lensf, lensf2, [0,1,1,0,0,0,0,0]
margin_fisher, planckf, planckf2, [0,1,1,0,0,0,0,0]
margin_fisher, test, test2, [0,1,1,0,0,0,0,0]   
margin_fisher, testbao, testbao2, [0,1,1,0,0,0,0,0]   
 margin_fisher, baof, baof2, [0,1,1,0,0,0]    
plt_fisher, lensf2, test2
plt_fisher, baof2, test2
stop
end
