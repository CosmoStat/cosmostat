;Compile necessary functions (change PATH according to the directory were the routines were installed)
.r /dsm/cosmo02/planck/sureau/code_idl/mrs_sparse_pointsource_removal.pro
.r /dsm/cosmo02/planck/sureau/code_idl/WMAP9_DATA_routines.pro

;LOCATION OF DATA (change this path too)
WMAP9_ROOT_DIR=getenv("DSKMSA")+'/WMAP9'
WMAP9_MAP_DIR= WMAP9_ROOT_DIR+'/data'
WMAP9_BEAM_DIR= WMAP9_ROOT_DIR +'/beams'
WMAP9_CTLG_DIR= WMAP9_ROOT_DIR +'/ptsrc_catalog'
WMAP9_GALMASK_DIR= WMAP9_ROOT_DIR +"/WMAP9_params/mask"


;;FROM THIS STEP, nothing to change

;Generate catalog
merged_ctlg_WMAP=WMAP9_PTSRCDATA_get_ctlg(WMAP9_CTLG_DIR)
print,(size(merged_ctlg_WMAP,/dim))[0]

;Example: PtSub for channel K, 1 iteration of SPSR
PtSub=WMAP9_PTSRCDATA_launch_SPSR(0, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=1)

;Example2: PtSub for channel W, 150 iteration of SPSR
PtSub=WMAP9_PTSRCDATA_launch_SPSR(4, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=150, RES_SPSR= RES_W_150,/verbose)
save,filename="ChannelW_test_fullpipe_150its.xdr", RES_W_150
PtSub=WMAP9_PTSRCDATA_launch_SPSR(4, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=300,INIT_SPSR= RES_W_150, RES_SPSR= RES_W_450,/verbose)


;Example 3: All iterations for all channels in a single step ... It takes a long time (better to subdivise)
PtSubK=WMAP9_PTSRCDATA_launch_SPSR(0, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=13350,/verbose)
PtSubKa=WMAP9_PTSRCDATA_launch_SPSR(1, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=9750,/verbose)
PtSubQ=WMAP9_PTSRCDATA_launch_SPSR(2, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=9750,/verbose)
PtSubV=WMAP9_PTSRCDATA_launch_SPSR(3, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=9750,/verbose)
PtSubW=WMAP9_PTSRCDATA_launch_SPSR(4, merged_ctlg_WMAP, WMAP9_BEAM_DIR, WMAP9_MAP_DIR, WMAP9_GALMASK_DIR,Niter=9750,/verbose)
