function set_fid_stat, nsiderot=nsiderot, verb=verb, lmaxgen=lmaxgen,lowquad=lowquad, quadoct=quadoct, planaroct=planaroct, aoe=aoe, mparity=mparity, TQuad=TQuad

; Written by Anais Rassat, May 2013. 
; PURPOSE: Set a fiducial structure with parameters for Statistical
; 20/03/2014: JLS : change the default value for nsiderot to 512, previous one was 128.
; 20/03/2014: JLS : add a TheoryQuad keyword to the value can be changed from the calling routine

; Isotropy calculations. Used as input to code Stat_Anomalies.pro
; 
; LowQuad, QuadOct, ParityOct, AOE and MParity are structures.
if not keyword_set(LowQuad) then LowQuad={Test:1}
if not keyword_set(QuadOct) then QuadOct={Test:1} ; basic minimum is calculating QuadOct
if not keyword_set(PlanarOct) then PlanarOct={Test:1}
if not keyword_set(AOE) then AOE={Test:1}
if not keyword_set(Mparity) then Mparity={Test:1}
if not keyword_set(TQuad) then TQuad = 1161.3421d0   ; WMAP9 tquad

; First define general parameters 
if not keyword_set(nsiderot) then nsiderot = 512L 
if not keyword_set(verb) then verb = 1
if not keyword_set(lmaxgen) then lmaxgen = 5L

; Low Quadrupole options
if not tag_check(LowQuad,'Test',val=Test) then TestQuad = 0 else TestQuad = LowQuad.Test
if not tag_check(LowQuad,'TheoryQuad',val=TheoryQuad) then TheoryQuad = TQuad  else TheoryQuad = TheoryQuad.Test ; in units of muK^2

;if not keyword_set(TestQuad) then TestQuad = 0
LowQuad = {Test:TestQuad, TheoryQuad:TheoryQuad,Note:"Test of Low value of quadrupole."}


; Quadrupole Octopole Alignment options
if not tag_check(QuadOct,'Test',val=Test) then TestQuadOct = 0 else TestQuadOct = QuadOct.Test
if not tag_check(QuadOct,'nsiderot',val=nsiderotQO) then nsiderotQO = nsiderot else nsiderQO = QuadOct.nsiderot
QuadOct = {Test:TestQuadOct, nsiderot:nsiderotQO, Note:"Test of Quadrupole/Octopole ALignment"}

; Planar Octopole options
if not tag_check(PlanarOct,'Test',val=Test) then TestPlanarOct = 0 else TestPlanarOct = PlanarOct.Test
if not tag_check(PlanarOct,'Nsimu',val=NsimuPO) then NsimuPO = 100L else NsimuPO = PlanarOct.Nsimu
if not tag_check(PlanarOct,'Noprob',val=noprobPO) then NoprobPO = 1 else NoprobPO = PlanarOct.Noprob
if not tag_check(PlanarOct,'nsiderot',val=nsiderotPO) then nsiderotPO = nsiderot else nsiderotPO = PlanarOct.nsiderot
PlanarOct = {Test:TestPlanarOct, nsiderot:nsiderotPO, nsimu:NsimuPO,noprob:noprobPO,Note:"Test of planarity of Octopole."}

; Axis of Evil options
if not tag_check(AOE,'Test',val=Test) then TestAOE = 0 else TestAOE = AOE.Test
if not tag_check(AOE,'Nsimu',val=NsimuAOE) then NsimuAOE = 100L else NsimuAOE = AOE.Nsimu
if not tag_check(AOE,'Noprob',val=noprobAOE) then NoprobAOE = 1 else NoprobAOE = AOE.Noprob
if not tag_check(AOE,'nsiderot',val=nsiderotAOE) then nsiderotAOE = nsiderot else nsiderotAOE = AOE.nsiderot
if not tag_check(AOE,'lmax',val=lmaxAOE) then lmaxAOE = lmaxgen else lmaxAOE = AOE.lmax
AOE = {Test:TestAOE,nsiderot:nsiderotAOE,nsimu:nsimuAOE,noprob:noprobAOE, lmax:lmaxAOE, Note:"Test of Axis of Evil."}

; Mirror Parity options
if not tag_check(Mparity,'Test',val=Test) then TestMP = 0 else TestMP = MParity.Test
if not tag_check(Mparity,'nside_degrade',val=nside_degradeMP) then nside_degradeMP = 8L else nside_degradeMP = MParity.nside_degrade
if not tag_check(Mparity,'nsiderot',val=nsiderotMP) then nsiderotMP = 64L else nsiderotMP = MParity.nsiderot    
if not tag_check(Mparity,'Noprob',val=noprobMP) then NoprobMP = 1 else NoprobMP = MParity.Noprob
if not tag_check(Mparity,'Nsimu',val=NsimuMP) then NsimuMP = 100L else NsimuMP = Mparity.Nsimu
if not tag_check(Mparity,'lmax',val=lmaxMP) then lmaxMP = lmaxgen else lmaxMP = Mparity.lmax
MParity = {Test:TestMP,lmax:lmaxMP,nside_degrade:nside_degradeMP,nsiderot:nsiderotMP,noprob:noprobMP,nsimu:nsimuMP}

fidstat = {nsiderot:nsiderot,lmaxgen:lmaxgen,verb:verb,LowQuad:LowQuad, QuadOct:QuadOct, PlanarOct:PlanarOct, AOE:AOE, MParity:MParity}
return, fidstat
end
