static mxArray *get_Mint_lo(const SERVO_T* Mint_lo);
static mxArray *get_aper(const APER_T* aper);
static mxArray *get_atmcfg(const GENATM_T* atmcfg);
static mxArray *get_cachedm_propdata(const PROPDATA_T* cachedm_propdata);
static mxArray *get_dmint(const SERVO_T* dmint);
static mxArray *get_evl_propdata_atm(const PROPDATA_T* evl_propdata_atm);
static mxArray *get_evl_propdata_dm(const PROPDATA_T* evl_propdata_dm);
static mxArray *get_parms(const PARMS_T* parms);
static mxArray *get_parms_aper(const APER_CFG_T* aper);
static mxArray *get_parms_atm(const ATM_CFG_T* atm);
static mxArray *get_parms_atmr(const ATMR_CFG_T* atmr);
static mxArray *get_parms_cn2(const CN2EST_CFG_T* cn2);
static mxArray *get_parms_dbg(const DBG_CFG_T* dbg);
static mxArray *get_parms_dm(const DM_CFG_T* dm);
static mxArray *get_parms_evl(const EVL_CFG_T* evl);
static mxArray *get_parms_fit(const FIT_CFG_T* fit);
static mxArray *get_parms_gpu(const GPU_CFG_T* gpu);
static mxArray *get_parms_load(const LOAD_CFG_T* load);
static mxArray *get_parms_lsr(const LSR_CFG_T* lsr);
static mxArray *get_parms_misreg(const MISREG_CFG_T* misreg);
static mxArray *get_parms_moao(const MOAO_CFG_T* moao);
static mxArray *get_parms_plot(const PLOT_CFG_T* plot);
static mxArray *get_parms_powfs(const POWFS_CFG_T* powfs, int nvar);
static mxArray *get_parms_powfs_llt(const LLT_CFG_T* llt);
static mxArray *get_parms_recon(const RECON_CFG_T* recon);
static mxArray *get_parms_sim(const SIM_CFG_T* sim);
static mxArray *get_parms_tomo(const TOMO_CFG_T* tomo);
static mxArray *get_parms_wfs(const WFS_CFG_T* wfs, int nvar);
static mxArray *get_parms_wfsr(const WFS_CFG_T* wfsr, int nvar);
static mxArray *get_powfs(const POWFS_T* powfs, int nvar);
static mxArray *get_powfs_dtf(const DTF_T* dtf);
static mxArray *get_powfs_etfprep(const ETF_T* etfprep);
static mxArray *get_powfs_etfsim(const ETF_T* etfsim);
static mxArray *get_powfs_intstat(const INTSTAT_T* intstat);
static mxArray *get_powfs_llt(const LLT_T* llt);
static mxArray *get_recon(const RECON_T* recon);
static mxArray *get_recon_FL(const MUV_T* FL);
static mxArray *get_recon_FR(const MUV_T* FR);
static mxArray *get_recon_LL(const MUV_T* LL);
static mxArray *get_recon_LR(const MUV_T* LR);
static mxArray *get_recon_RL(const MUV_T* RL);
static mxArray *get_recon_RR(const MUV_T* RR);
static mxArray *get_recon_cn2est(const CN2EST_T* cn2est);
static mxArray *get_recon_cn2est_pair(const CN2PAIR_T* pair);
static mxArray *get_recon_fdpcg(const FDPCG_T* fdpcg);
static mxArray *get_recon_fractal(const FRACTAL_T* fractal);
static mxArray *get_recon_invpsd(const INVPSD_T* invpsd);
static mxArray *get_recon_moao(const MOAO_T* moao);
static mxArray *get_recon_ngsmod(const NGSMOD_T* ngsmod);
static mxArray *get_simu(const SIM_T* simu);
static mxArray *get_uptint(const SERVO_T* uptint);
static mxArray *get_wfs_intsdata(const WFSINTS_T* wfs_intsdata);
static mxArray *get_wfs_propdata_atm(const PROPDATA_T* wfs_propdata_atm);
static mxArray *get_wfs_propdata_dm(const PROPDATA_T* wfs_propdata_dm);
static mxArray *get_Mint_lo(const SERVO_T* Mint_lo){
	if(!Mint_lo) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(Mint_lo->ap);i=mxAddField(tmp, "ap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(Mint_lo->ep);i=mxAddField(tmp, "ep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(Mint_lo->merrhist);i=mxAddField(tmp, "merrhist");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(Mint_lo->merrlast);i=mxAddField(tmp, "merrlast");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(Mint_lo->mint);i=mxAddField(tmp, "mint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(Mint_lo->mlead);i=mxAddField(tmp, "mlead");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(Mint_lo->mpreint);i=mxAddField(tmp, "mpreint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(Mint_lo->al);i=mxAddField(tmp, "al");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(Mint_lo->dt);i=mxAddField(tmp, "dt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(Mint_lo->initialized);i=mxAddField(tmp, "initialized");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_aper(const APER_T* aper){
	if(!aper) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(aper->amp);i=mxAddField(tmp, "amp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->amp1);i=mxAddField(tmp, "amp1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->ampground);i=mxAddField(tmp, "ampground");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->imcc);i=mxAddField(tmp, "imcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->locs);i=mxAddField(tmp, "locs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->locs_dm);i=mxAddField(tmp, "locs_dm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->mcc);i=mxAddField(tmp, "mcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->mod);i=mxAddField(tmp, "mod");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->opdadd);i=mxAddField(tmp, "opdadd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(aper->opdfloc);i=mxAddField(tmp, "opdfloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(aper->fcp);i=mxAddField(tmp, "fcp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(aper->ipcc);i=mxAddField(tmp, "ipcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(aper->sumamp2);i=mxAddField(tmp, "sumamp2");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_atmcfg(const GENATM_T* atmcfg){
	if(!atmcfg) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(atmcfg->screen);i=mxAddField(tmp, "screen");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atmcfg->spect);i=mxAddField(tmp, "spect");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->ilayer);i=mxAddField(tmp, "ilayer");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->l0);i=mxAddField(tmp, "l0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->method);i=mxAddField(tmp, "method");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->ninit);i=mxAddField(tmp, "ninit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->nlayer);i=mxAddField(tmp, "nlayer");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->nthread);i=mxAddField(tmp, "nthread");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->nx);i=mxAddField(tmp, "nx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->ny);i=mxAddField(tmp, "ny");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->r0);i=mxAddField(tmp, "r0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmcfg->share);i=mxAddField(tmp, "share");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_cachedm_propdata(const PROPDATA_T* cachedm_propdata){
	if(!cachedm_propdata) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(cachedm_propdata->locin);i=mxAddField(tmp, "locin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cachedm_propdata->locout);i=mxAddField(tmp, "locout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cachedm_propdata->mapin);i=mxAddField(tmp, "mapin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cachedm_propdata->mapout);i=mxAddField(tmp, "mapout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cachedm_propdata->ptsout);i=mxAddField(tmp, "ptsout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->alpha);i=mxAddField(tmp, "alpha");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->cubic_iac);i=mxAddField(tmp, "cubic_iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->index);i=mxAddField(tmp, "index");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->nooptim);i=mxAddField(tmp, "nooptim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cachedm_propdata->wrap);i=mxAddField(tmp, "wrap");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_dmint(const SERVO_T* dmint){
	if(!dmint) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(dmint->ap);i=mxAddField(tmp, "ap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dmint->ep);i=mxAddField(tmp, "ep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dmint->merrhist);i=mxAddField(tmp, "merrhist");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dmint->merrlast);i=mxAddField(tmp, "merrlast");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dmint->mint);i=mxAddField(tmp, "mint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dmint->mlead);i=mxAddField(tmp, "mlead");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dmint->mpreint);i=mxAddField(tmp, "mpreint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dmint->al);i=mxAddField(tmp, "al");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dmint->dt);i=mxAddField(tmp, "dt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dmint->initialized);i=mxAddField(tmp, "initialized");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_evl_propdata_atm(const PROPDATA_T* evl_propdata_atm){
	if(!evl_propdata_atm) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(evl_propdata_atm->locin);i=mxAddField(tmp, "locin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_atm->locout);i=mxAddField(tmp, "locout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_atm->mapin);i=mxAddField(tmp, "mapin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_atm->mapout);i=mxAddField(tmp, "mapout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_atm->ptsout);i=mxAddField(tmp, "ptsout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->alpha);i=mxAddField(tmp, "alpha");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->cubic_iac);i=mxAddField(tmp, "cubic_iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->index);i=mxAddField(tmp, "index");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->nooptim);i=mxAddField(tmp, "nooptim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_atm->wrap);i=mxAddField(tmp, "wrap");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_evl_propdata_dm(const PROPDATA_T* evl_propdata_dm){
	if(!evl_propdata_dm) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(evl_propdata_dm->locin);i=mxAddField(tmp, "locin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_dm->locout);i=mxAddField(tmp, "locout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_dm->mapin);i=mxAddField(tmp, "mapin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_dm->mapout);i=mxAddField(tmp, "mapout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl_propdata_dm->ptsout);i=mxAddField(tmp, "ptsout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->alpha);i=mxAddField(tmp, "alpha");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->cubic_iac);i=mxAddField(tmp, "cubic_iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->index);i=mxAddField(tmp, "index");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->nooptim);i=mxAddField(tmp, "nooptim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl_propdata_dm->wrap);i=mxAddField(tmp, "wrap");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms(const PARMS_T* parms){
	if(!parms) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(parms->dirs);i=mxAddField(tmp, "dirs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(parms->fdlock);i=mxAddField(tmp, "fdlock");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(parms->hipowfs);i=mxAddField(tmp, "hipowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(parms->lopowfs);i=mxAddField(tmp, "lopowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_aper(&(parms->aper));i=mxAddField(tmp, "aper");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_atm(&(parms->atm));i=mxAddField(tmp, "atm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_atmr(&(parms->atmr));i=mxAddField(tmp, "atmr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_cn2(&(parms->cn2));i=mxAddField(tmp, "cn2");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_dbg(&(parms->dbg));i=mxAddField(tmp, "dbg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_dm(parms->dm);i=mxAddField(tmp, "dm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_evl(&(parms->evl));i=mxAddField(tmp, "evl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_fit(&(parms->fit));i=mxAddField(tmp, "fit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_gpu(&(parms->gpu));i=mxAddField(tmp, "gpu");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_load(&(parms->load));i=mxAddField(tmp, "load");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_lsr(&(parms->lsr));i=mxAddField(tmp, "lsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_misreg(&(parms->misreg));i=mxAddField(tmp, "misreg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_moao(parms->moao);i=mxAddField(tmp, "moao");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_plot(&(parms->plot));i=mxAddField(tmp, "plot");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_powfs(parms->powfs, parms->npowfs);i=mxAddField(tmp, "powfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_recon(&(parms->recon));i=mxAddField(tmp, "recon");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_sim(&(parms->sim));i=mxAddField(tmp, "sim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_tomo(&(parms->tomo));i=mxAddField(tmp, "tomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_wfs(parms->wfs, parms->nwfs);i=mxAddField(tmp, "wfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms_wfsr(parms->wfsr, parms->nwfsr);i=mxAddField(tmp, "wfsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->ndm);i=mxAddField(tmp, "ndm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nhipowfs);i=mxAddField(tmp, "nhipowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nlopowfs);i=mxAddField(tmp, "nlopowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nlowfs);i=mxAddField(tmp, "nlowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nmoao);i=mxAddField(tmp, "nmoao");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nphypowfs);i=mxAddField(tmp, "nphypowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->npowfs);i=mxAddField(tmp, "npowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nsurf);i=mxAddField(tmp, "nsurf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->ntipowfs);i=mxAddField(tmp, "ntipowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->ntrpowfs);i=mxAddField(tmp, "ntrpowfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->ntsurf);i=mxAddField(tmp, "ntsurf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nwfs);i=mxAddField(tmp, "nwfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(parms->nwfsr);i=mxAddField(tmp, "nwfsr");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_aper(const APER_CFG_T* aper){
	if(!aper) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(aper->d);i=mxAddField(tmp, "d");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(aper->din);i=mxAddField(tmp, "din");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(aper->fnampuser);i=mxAddField(tmp, "fnampuser");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(aper->rotdeg);i=mxAddField(tmp, "rotdeg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(aper->fnamp);i=mxAddField(tmp, "fnamp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(aper->pupmask);i=mxAddField(tmp, "pupmask");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_atm(const ATM_CFG_T* atm){
	if(!atm) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(atm->ht);i=mxAddField(tmp, "ht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->ipsr);i=mxAddField(tmp, "ipsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->overx);i=mxAddField(tmp, "overx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->overy);i=mxAddField(tmp, "overy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->size);i=mxAddField(tmp, "size");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->wddeg);i=mxAddField(tmp, "wddeg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->ws);i=mxAddField(tmp, "ws");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atm->wt);i=mxAddField(tmp, "wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->evolve);i=mxAddField(tmp, "evolve");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->fractal);i=mxAddField(tmp, "fractal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->frozenflow);i=mxAddField(tmp, "frozenflow");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->hmax);i=mxAddField(tmp, "hmax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->iground);i=mxAddField(tmp, "iground");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->l0);i=mxAddField(tmp, "l0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->ninit);i=mxAddField(tmp, "ninit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->nps);i=mxAddField(tmp, "nps");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->nx);i=mxAddField(tmp, "nx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->nxm);i=mxAddField(tmp, "nxm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->nxn);i=mxAddField(tmp, "nxn");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->ny);i=mxAddField(tmp, "ny");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->nym);i=mxAddField(tmp, "nym");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->nyn);i=mxAddField(tmp, "nyn");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->r0);i=mxAddField(tmp, "r0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->r0z);i=mxAddField(tmp, "r0z");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->share);i=mxAddField(tmp, "share");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atm->wdrand);i=mxAddField(tmp, "wdrand");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_atmr(const ATMR_CFG_T* atmr){
	if(!atmr) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(atmr->ht);i=mxAddField(tmp, "ht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atmr->indps);i=mxAddField(tmp, "indps");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atmr->os);i=mxAddField(tmp, "os");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(atmr->wt);i=mxAddField(tmp, "wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->hmax);i=mxAddField(tmp, "hmax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->l0);i=mxAddField(tmp, "l0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->nps);i=mxAddField(tmp, "nps");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->r0);i=mxAddField(tmp, "r0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(atmr->r0z);i=mxAddField(tmp, "r0z");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_cn2(const CN2EST_CFG_T* cn2){
	if(!cn2) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(cn2->pair);i=mxAddField(tmp, "pair");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->hmax);i=mxAddField(tmp, "hmax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->keepht);i=mxAddField(tmp, "keepht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->moveht);i=mxAddField(tmp, "moveht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->nhtomo);i=mxAddField(tmp, "nhtomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->psol);i=mxAddField(tmp, "psol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->reset);i=mxAddField(tmp, "reset");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->saat);i=mxAddField(tmp, "saat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->step);i=mxAddField(tmp, "step");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->tomo);i=mxAddField(tmp, "tomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2->verbose);i=mxAddField(tmp, "verbose");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_dbg(const DBG_CFG_T* dbg){
	if(!dbg) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(dbg->tomo_maxit);i=mxAddField(tmp, "tomo_maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->annular_W);i=mxAddField(tmp, "annular_W");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->atm);i=mxAddField(tmp, "atm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->cmpgpu);i=mxAddField(tmp, "cmpgpu");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->deltafocus);i=mxAddField(tmp, "deltafocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->dmfullfov);i=mxAddField(tmp, "dmfullfov");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->ecovxx);i=mxAddField(tmp, "ecovxx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->fit);i=mxAddField(tmp, "fit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->force);i=mxAddField(tmp, "force");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->mvstlimit);i=mxAddField(tmp, "mvstlimit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->na_interp);i=mxAddField(tmp, "na_interp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->na_smooth);i=mxAddField(tmp, "na_smooth");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->nocgwarm);i=mxAddField(tmp, "nocgwarm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->psol);i=mxAddField(tmp, "psol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->pupmask);i=mxAddField(tmp, "pupmask");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->test);i=mxAddField(tmp, "test");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->tomo);i=mxAddField(tmp, "tomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->tomo_hxw);i=mxAddField(tmp, "tomo_hxw");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->usegwr);i=mxAddField(tmp, "usegwr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->useopdr);i=mxAddField(tmp, "useopdr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->wamethod);i=mxAddField(tmp, "wamethod");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dbg->wfslinearity);i=mxAddField(tmp, "wfslinearity");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_dm(const DM_CFG_T* dm){
	if(!dm) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(dm->dxcache);i=mxAddField(tmp, "dxcache");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dm->iastrokescale);i=mxAddField(tmp, "iastrokescale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->ar);i=mxAddField(tmp, "ar");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->dy);i=mxAddField(tmp, "dy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->guard);i=mxAddField(tmp, "guard");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->hist);i=mxAddField(tmp, "hist");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->histbin);i=mxAddField(tmp, "histbin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->histn);i=mxAddField(tmp, "histn");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->ht);i=mxAddField(tmp, "ht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->iac);i=mxAddField(tmp, "iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->iastroke);i=mxAddField(tmp, "iastroke");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->isground);i=mxAddField(tmp, "isground");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->ncache);i=mxAddField(tmp, "ncache");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->offset);i=mxAddField(tmp, "offset");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->order);i=mxAddField(tmp, "order");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->stroke);i=mxAddField(tmp, "stroke");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dm->vmisreg);i=mxAddField(tmp, "vmisreg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(dm->actfloat);i=mxAddField(tmp, "actfloat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(dm->actstuck);i=mxAddField(tmp, "actstuck");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(dm->hyst);i=mxAddField(tmp, "hyst");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(dm->iastrokefn);i=mxAddField(tmp, "iastrokefn");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_evl(const EVL_CFG_T* evl){
	if(!evl) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(evl->hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->psf);i=mxAddField(tmp, "psf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->psfgridsize);i=mxAddField(tmp, "psfgridsize");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->psfngsr);i=mxAddField(tmp, "psfngsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->psfr);i=mxAddField(tmp, "psfr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->psfsize);i=mxAddField(tmp, "psfsize");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->pttr);i=mxAddField(tmp, "pttr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->scalegroup);i=mxAddField(tmp, "scalegroup");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->sock);i=mxAddField(tmp, "sock");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->thetax);i=mxAddField(tmp, "thetax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->thetay);i=mxAddField(tmp, "thetay");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->wt);i=mxAddField(tmp, "wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(evl->wvl);i=mxAddField(tmp, "wvl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->cov);i=mxAddField(tmp, "cov");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->indoa);i=mxAddField(tmp, "indoa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->moao);i=mxAddField(tmp, "moao");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->nevl);i=mxAddField(tmp, "nevl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->nmod);i=mxAddField(tmp, "nmod");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->npsf);i=mxAddField(tmp, "npsf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->nthread);i=mxAddField(tmp, "nthread");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->nwvl);i=mxAddField(tmp, "nwvl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->psfhist);i=mxAddField(tmp, "psfhist");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->psfisim);i=mxAddField(tmp, "psfisim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->psfmean);i=mxAddField(tmp, "psfmean");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->psfol);i=mxAddField(tmp, "psfol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->rmax);i=mxAddField(tmp, "rmax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(evl->tomo);i=mxAddField(tmp, "tomo");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_fit(const FIT_CFG_T* fit){
	if(!fit) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(fit->hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fit->thetax);i=mxAddField(tmp, "thetax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fit->thetay);i=mxAddField(tmp, "thetay");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fit->wt);i=mxAddField(tmp, "wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->actinterp);i=mxAddField(tmp, "actinterp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->actslave);i=mxAddField(tmp, "actslave");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->assemble);i=mxAddField(tmp, "assemble");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->cachedm);i=mxAddField(tmp, "cachedm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->cachex);i=mxAddField(tmp, "cachex");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->indoa);i=mxAddField(tmp, "indoa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->lrt_piston);i=mxAddField(tmp, "lrt_piston");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->lrt_tt);i=mxAddField(tmp, "lrt_tt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->nfit);i=mxAddField(tmp, "nfit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->pos);i=mxAddField(tmp, "pos");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->precond);i=mxAddField(tmp, "precond");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->square);i=mxAddField(tmp, "square");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->svdthres);i=mxAddField(tmp, "svdthres");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fit->tikcr);i=mxAddField(tmp, "tikcr");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_gpu(const GPU_CFG_T* gpu){
	if(!gpu) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(gpu->evl);i=mxAddField(tmp, "evl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(gpu->fit);i=mxAddField(tmp, "fit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(gpu->lsr);i=mxAddField(tmp, "lsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(gpu->moao);i=mxAddField(tmp, "moao");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(gpu->psf);i=mxAddField(tmp, "psf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(gpu->tomo);i=mxAddField(tmp, "tomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(gpu->wfs);i=mxAddField(tmp, "wfs");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_load(const LOAD_CFG_T* load){
	if(!load) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(load->GS0);i=mxAddField(tmp, "GS0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(load->W);i=mxAddField(tmp, "W");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(load->fit);i=mxAddField(tmp, "fit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(load->i0);i=mxAddField(tmp, "i0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(load->mvst);i=mxAddField(tmp, "mvst");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(load->tomo);i=mxAddField(tmp, "tomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->GA);i=mxAddField(tmp, "GA");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->GP);i=mxAddField(tmp, "GP");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->HA);i=mxAddField(tmp, "HA");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->HXF);i=mxAddField(tmp, "HXF");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->HXW);i=mxAddField(tmp, "HXW");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->aloc);i=mxAddField(tmp, "aloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->atm);i=mxAddField(tmp, "atm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->cxx);i=mxAddField(tmp, "cxx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->floc);i=mxAddField(tmp, "floc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->locs);i=mxAddField(tmp, "locs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->mvm);i=mxAddField(tmp, "mvm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->mvmf);i=mxAddField(tmp, "mvmf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->mvmi);i=mxAddField(tmp, "mvmi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->ploc);i=mxAddField(tmp, "ploc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(load->xloc);i=mxAddField(tmp, "xloc");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_lsr(const LSR_CFG_T* lsr){
	if(!lsr) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(lsr->actslave);i=mxAddField(tmp, "actslave");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(lsr->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(lsr->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(lsr->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(lsr->svdthres);i=mxAddField(tmp, "svdthres");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(lsr->tikcr);i=mxAddField(tmp, "tikcr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(lsr->fnreg);i=mxAddField(tmp, "fnreg");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_misreg(const MISREG_CFG_T* misreg){
	if(!misreg) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(misreg->pupil);i=mxAddField(tmp, "pupil");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_moao(const MOAO_CFG_T* moao){
	if(!moao) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(moao->actslave);i=mxAddField(tmp, "actslave");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->ar);i=mxAddField(tmp, "ar");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->gdm);i=mxAddField(tmp, "gdm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->iac);i=mxAddField(tmp, "iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->lrt_ptt);i=mxAddField(tmp, "lrt_ptt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->order);i=mxAddField(tmp, "order");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->stroke);i=mxAddField(tmp, "stroke");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->used);i=mxAddField(tmp, "used");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(moao->actfloat);i=mxAddField(tmp, "actfloat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(moao->actstuck);i=mxAddField(tmp, "actstuck");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_plot(const PLOT_CFG_T* plot){
	if(!plot) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(plot->all);i=mxAddField(tmp, "all");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(plot->atm);i=mxAddField(tmp, "atm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(plot->opdx);i=mxAddField(tmp, "opdx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(plot->run);i=mxAddField(tmp, "run");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(plot->setup);i=mxAddField(tmp, "setup");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_powfs(const POWFS_CFG_T* powfs, int nvar){
	if(!powfs) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(nvar,1,0,0);
	int i;
	mxArray *tmp2;
for(int i_n=0; i_n<nvar; i_n++){
	tmp2=any2mx(powfs[i_n].scalegroup);i=mxAddField(tmp, "scalegroup");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].wfs);i=mxAddField(tmp, "wfs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].wfsind);i=mxAddField(tmp, "wfsind");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].wvl);i=mxAddField(tmp, "wvl");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].wvlwts);i=mxAddField(tmp, "wvlwts");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=get_parms_powfs_llt(powfs[i_n].llt);i=mxAddField(tmp, "llt");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].bkgrnd);i=mxAddField(tmp, "bkgrnd");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].bkgrndc);i=mxAddField(tmp, "bkgrndc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].cogoff);i=mxAddField(tmp, "cogoff");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].cogthres);i=mxAddField(tmp, "cogthres");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dfrs);i=mxAddField(tmp, "dfrs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dither);i=mxAddField(tmp, "dither");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dither_amp);i=mxAddField(tmp, "dither_amp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dither_gpll);i=mxAddField(tmp, "dither_gpll");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dither_npll);i=mxAddField(tmp, "dither_npll");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dither_nskip);i=mxAddField(tmp, "dither_nskip");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dither_nstat);i=mxAddField(tmp, "dither_nstat");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dtrat);i=mxAddField(tmp, "dtrat");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].embfac);i=mxAddField(tmp, "embfac");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].fieldstop);i=mxAddField(tmp, "fieldstop");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].gtype_recon);i=mxAddField(tmp, "gtype_recon");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].gtype_sim);i=mxAddField(tmp, "gtype_sim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].i0scale);i=mxAddField(tmp, "i0scale");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].idtrat);i=mxAddField(tmp, "idtrat");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].lo);i=mxAddField(tmp, "lo");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].moao);i=mxAddField(tmp, "moao");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].mtchadp);i=mxAddField(tmp, "mtchadp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].mtchcpl);i=mxAddField(tmp, "mtchcpl");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].mtchcr);i=mxAddField(tmp, "mtchcr");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].mtchcra);i=mxAddField(tmp, "mtchcra");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].mtchscl);i=mxAddField(tmp, "mtchscl");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].mtchstc);i=mxAddField(tmp, "mtchstc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].ncomp);i=mxAddField(tmp, "ncomp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].ncpa_method);i=mxAddField(tmp, "ncpa_method");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].neaphy);i=mxAddField(tmp, "neaphy");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].neasim);i=mxAddField(tmp, "neasim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].neaspeckle);i=mxAddField(tmp, "neaspeckle");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].needGS0);i=mxAddField(tmp, "needGS0");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].noisy);i=mxAddField(tmp, "noisy");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].nwfs);i=mxAddField(tmp, "nwfs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].nwfsr);i=mxAddField(tmp, "nwfsr");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].nwvl);i=mxAddField(tmp, "nwvl");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].order);i=mxAddField(tmp, "order");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].phystep);i=mxAddField(tmp, "phystep");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].phytype);i=mxAddField(tmp, "phytype");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].phytypesim);i=mxAddField(tmp, "phytypesim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].phyusenea);i=mxAddField(tmp, "phyusenea");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pistatout);i=mxAddField(tmp, "pistatout");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pistatstart);i=mxAddField(tmp, "pistatstart");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pistatstc);i=mxAddField(tmp, "pistatstc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixblur);i=mxAddField(tmp, "pixblur");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixoffx);i=mxAddField(tmp, "pixoffx");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixoffy);i=mxAddField(tmp, "pixoffy");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixpsa);i=mxAddField(tmp, "pixpsa");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixtheta);i=mxAddField(tmp, "pixtheta");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].psfout);i=mxAddField(tmp, "psfout");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].psol);i=mxAddField(tmp, "psol");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].radgx);i=mxAddField(tmp, "radgx");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].radpix);i=mxAddField(tmp, "radpix");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].radpixtheta);i=mxAddField(tmp, "radpixtheta");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].radrot);i=mxAddField(tmp, "radrot");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].rne);i=mxAddField(tmp, "rne");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].saat);i=mxAddField(tmp, "saat");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].safill2d);i=mxAddField(tmp, "safill2d");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].safocuspv);i=mxAddField(tmp, "safocuspv");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].saspherical);i=mxAddField(tmp, "saspherical");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].siglev);i=mxAddField(tmp, "siglev");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].sigscale);i=mxAddField(tmp, "sigscale");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].skip);i=mxAddField(tmp, "skip");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].trs);i=mxAddField(tmp, "trs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].usephy);i=mxAddField(tmp, "usephy");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].bkgrndfn);i=mxAddField(tmp, "bkgrndfn");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].bkgrndfnc);i=mxAddField(tmp, "bkgrndfnc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].fnllt);i=mxAddField(tmp, "fnllt");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].neareconfile);i=mxAddField(tmp, "neareconfile");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].neasimfile);i=mxAddField(tmp, "neasimfile");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].piinfile);i=mxAddField(tmp, "piinfile");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].saloc);i=mxAddField(tmp, "saloc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=str2mx(powfs[i_n].sninfile);i=mxAddField(tmp, "sninfile");mxSetFieldByNumber(tmp, i_n, i, tmp2);}
	return tmp;
}

static mxArray *get_parms_powfs_llt(const LLT_CFG_T* llt){
	if(!llt) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(llt->i);i=mxAddField(tmp, "i");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->misreg);i=mxAddField(tmp, "misreg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->ox);i=mxAddField(tmp, "ox");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->oy);i=mxAddField(tmp, "oy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->colprep);i=mxAddField(tmp, "colprep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->colsim);i=mxAddField(tmp, "colsim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->colsimdtrat);i=mxAddField(tmp, "colsimdtrat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->d);i=mxAddField(tmp, "d");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->n);i=mxAddField(tmp, "n");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->ttrat);i=mxAddField(tmp, "ttrat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(llt->widthp);i=mxAddField(tmp, "widthp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(llt->fnamp);i=mxAddField(tmp, "fnamp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(llt->fnprof);i=mxAddField(tmp, "fnprof");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(llt->fnrange);i=mxAddField(tmp, "fnrange");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(llt->fnsurf);i=mxAddField(tmp, "fnsurf");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_recon(const RECON_CFG_T* recon){
	if(!recon) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(recon->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->glao);i=mxAddField(tmp, "glao");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->mvm);i=mxAddField(tmp, "mvm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->sock);i=mxAddField(tmp, "sock");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->split);i=mxAddField(tmp, "split");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->warm_restart);i=mxAddField(tmp, "warm_restart");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_sim(const SIM_CFG_T* sim){
	if(!sim) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(sim->apdm);i=mxAddField(tmp, "apdm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->aplo);i=mxAddField(tmp, "aplo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->apupt);i=mxAddField(tmp, "apupt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->epdm);i=mxAddField(tmp, "epdm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->eplo);i=mxAddField(tmp, "eplo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->epupt);i=mxAddField(tmp, "epupt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->ncpa_hs);i=mxAddField(tmp, "ncpa_hs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->ncpa_thetax);i=mxAddField(tmp, "ncpa_thetax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->ncpa_thetay);i=mxAddField(tmp, "ncpa_thetay");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->ncpa_wt);i=mxAddField(tmp, "ncpa_wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(sim->seeds);i=mxAddField(tmp, "seeds");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->ahstfocus);i=mxAddField(tmp, "ahstfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->aldm);i=mxAddField(tmp, "aldm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->allo);i=mxAddField(tmp, "allo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->alupt);i=mxAddField(tmp, "alupt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->cachedm);i=mxAddField(tmp, "cachedm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->closeloop);i=mxAddField(tmp, "closeloop");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dmclip);i=mxAddField(tmp, "dmclip");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dmclipia);i=mxAddField(tmp, "dmclipia");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dmproj);i=mxAddField(tmp, "dmproj");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dt);i=mxAddField(tmp, "dt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dthi);i=mxAddField(tmp, "dthi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dtlo);i=mxAddField(tmp, "dtlo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dtrat_hi);i=mxAddField(tmp, "dtrat_hi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dtrat_lo);i=mxAddField(tmp, "dtrat_lo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->dtrat_skip);i=mxAddField(tmp, "dtrat_skip");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->ecnn);i=mxAddField(tmp, "ecnn");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->end);i=mxAddField(tmp, "end");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->evlol);i=mxAddField(tmp, "evlol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->fcfocus);i=mxAddField(tmp, "fcfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->fcttm);i=mxAddField(tmp, "fcttm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->fov);i=mxAddField(tmp, "fov");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->fuseint);i=mxAddField(tmp, "fuseint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->idealevl);i=mxAddField(tmp, "idealevl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->idealfit);i=mxAddField(tmp, "idealfit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->idealtomo);i=mxAddField(tmp, "idealtomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->idealwfs);i=mxAddField(tmp, "idealwfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->lpfocus);i=mxAddField(tmp, "lpfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->lpttm);i=mxAddField(tmp, "lpttm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->mffocus);i=mxAddField(tmp, "mffocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->mvmngpu);i=mxAddField(tmp, "mvmngpu");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->mvmport);i=mxAddField(tmp, "mvmport");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->mvmsize);i=mxAddField(tmp, "mvmsize");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->ncpa_calib);i=mxAddField(tmp, "ncpa_calib");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->ncpa_ndir);i=mxAddField(tmp, "ncpa_ndir");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->ncpa_ttr);i=mxAddField(tmp, "ncpa_ttr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->noatm);i=mxAddField(tmp, "noatm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->nseed);i=mxAddField(tmp, "nseed");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->pause);i=mxAddField(tmp, "pause");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->psfr);i=mxAddField(tmp, "psfr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->servotype_hi);i=mxAddField(tmp, "servotype_hi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->servotype_lo);i=mxAddField(tmp, "servotype_lo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->skysim);i=mxAddField(tmp, "skysim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->start);i=mxAddField(tmp, "start");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->uptideal);i=mxAddField(tmp, "uptideal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->wfsalias);i=mxAddField(tmp, "wfsalias");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->wsseq);i=mxAddField(tmp, "wsseq");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->za);i=mxAddField(tmp, "za");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->zoomdtrat);i=mxAddField(tmp, "zoomdtrat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->zoomgain);i=mxAddField(tmp, "zoomgain");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(sim->zoomshare);i=mxAddField(tmp, "zoomshare");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(sim->dmadd);i=mxAddField(tmp, "dmadd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(sim->mvmhost);i=mxAddField(tmp, "mvmhost");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=str2mx(sim->wspsd);i=mxAddField(tmp, "wspsd");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_tomo(const TOMO_CFG_T* tomo){
	if(!tomo) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(tomo->ahst_idealngs);i=mxAddField(tmp, "ahst_idealngs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->ahst_ttr);i=mxAddField(tmp, "ahst_ttr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->ahst_wt);i=mxAddField(tmp, "ahst_wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->assemble);i=mxAddField(tmp, "assemble");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->cgthres);i=mxAddField(tmp, "cgthres");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->cone);i=mxAddField(tmp, "cone");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->cxx);i=mxAddField(tmp, "cxx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->guard);i=mxAddField(tmp, "guard");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->iac);i=mxAddField(tmp, "iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->minwt);i=mxAddField(tmp, "minwt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->ninit);i=mxAddField(tmp, "ninit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->nxbase);i=mxAddField(tmp, "nxbase");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->piston_cr);i=mxAddField(tmp, "piston_cr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->pos);i=mxAddField(tmp, "pos");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->precond);i=mxAddField(tmp, "precond");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->predict);i=mxAddField(tmp, "predict");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->psol);i=mxAddField(tmp, "psol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->splitlrt);i=mxAddField(tmp, "splitlrt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->square);i=mxAddField(tmp, "square");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->svdthres);i=mxAddField(tmp, "svdthres");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(tomo->tikcr);i=mxAddField(tmp, "tikcr");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_parms_wfs(const WFS_CFG_T* wfs, int nvar){
	if(!wfs) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(nvar,1,0,0);
	int i;
	mxArray *tmp2;
for(int i_n=0; i_n<nvar; i_n++){
	tmp2=any2mx(wfs[i_n].wvlwts);i=mxAddField(tmp, "wvlwts");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].fitwt);i=mxAddField(tmp, "fitwt");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].powfs);i=mxAddField(tmp, "powfs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].siglev);i=mxAddField(tmp, "siglev");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].siglevsim);i=mxAddField(tmp, "siglevsim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].sock);i=mxAddField(tmp, "sock");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].thetax);i=mxAddField(tmp, "thetax");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs[i_n].thetay);i=mxAddField(tmp, "thetay");mxSetFieldByNumber(tmp, i_n, i, tmp2);}
	return tmp;
}

static mxArray *get_parms_wfsr(const WFS_CFG_T* wfsr, int nvar){
	if(!wfsr) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(nvar,1,0,0);
	int i;
	mxArray *tmp2;
for(int i_n=0; i_n<nvar; i_n++){
	tmp2=any2mx(wfsr[i_n].wvlwts);i=mxAddField(tmp, "wvlwts");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].fitwt);i=mxAddField(tmp, "fitwt");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].powfs);i=mxAddField(tmp, "powfs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].siglev);i=mxAddField(tmp, "siglev");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].siglevsim);i=mxAddField(tmp, "siglevsim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].sock);i=mxAddField(tmp, "sock");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].thetax);i=mxAddField(tmp, "thetax");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfsr[i_n].thetay);i=mxAddField(tmp, "thetay");mxSetFieldByNumber(tmp, i_n, i, tmp2);}
	return tmp;
}

static mxArray *get_powfs(const POWFS_T* powfs, int nvar){
	if(!powfs) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(nvar,1,0,0);
	int i;
	mxArray *tmp2;
for(int i_n=0; i_n<nvar; i_n++){
	tmp2=any2mx(powfs[i_n].GS0);i=mxAddField(tmp, "GS0");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].amp);i=mxAddField(tmp, "amp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].amp_tel);i=mxAddField(tmp, "amp_tel");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].bkgrnd);i=mxAddField(tmp, "bkgrnd");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].bkgrndc);i=mxAddField(tmp, "bkgrndc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].dtheta);i=mxAddField(tmp, "dtheta");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].focus);i=mxAddField(tmp, "focus");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].gamp);i=mxAddField(tmp, "gamp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].gloc);i=mxAddField(tmp, "gloc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].gradoff);i=mxAddField(tmp, "gradoff");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].gradphyoff);i=mxAddField(tmp, "gradphyoff");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].loc);i=mxAddField(tmp, "loc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].loc_dm);i=mxAddField(tmp, "loc_dm");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].loc_tel);i=mxAddField(tmp, "loc_tel");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].neasim);i=mxAddField(tmp, "neasim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].opdadd);i=mxAddField(tmp, "opdadd");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].opdbias);i=mxAddField(tmp, "opdbias");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].pts);i=mxAddField(tmp, "pts");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].realamp);i=mxAddField(tmp, "realamp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].realsaa);i=mxAddField(tmp, "realsaa");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].saa);i=mxAddField(tmp, "saa");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].saa_tel);i=mxAddField(tmp, "saa_tel");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].saimcc);i=mxAddField(tmp, "saimcc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].saloc);i=mxAddField(tmp, "saloc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].sodium);i=mxAddField(tmp, "sodium");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].sprint);i=mxAddField(tmp, "sprint");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].srot);i=mxAddField(tmp, "srot");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].srsa);i=mxAddField(tmp, "srsa");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].srsamax);i=mxAddField(tmp, "srsamax");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].sumamp);i=mxAddField(tmp, "sumamp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=any2mx(powfs[i_n].sumamp2);i=mxAddField(tmp, "sumamp2");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=get_powfs_dtf(powfs[i_n].dtf);i=mxAddField(tmp, "dtf");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=get_powfs_etfprep(powfs[i_n].etfprep);i=mxAddField(tmp, "etfprep");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=get_powfs_etfsim(powfs[i_n].etfsim);i=mxAddField(tmp, "etfsim");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=get_powfs_intstat(powfs[i_n].intstat);i=mxAddField(tmp, "intstat");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=get_powfs_llt(powfs[i_n].llt);i=mxAddField(tmp, "llt");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].areascale);i=mxAddField(tmp, "areascale");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].namp);i=mxAddField(tmp, "namp");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].ncompx);i=mxAddField(tmp, "ncompx");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].ncompy);i=mxAddField(tmp, "ncompy");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].npts);i=mxAddField(tmp, "npts");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].nsaimcc);i=mxAddField(tmp, "nsaimcc");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].nthread);i=mxAddField(tmp, "nthread");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].nwfs);i=mxAddField(tmp, "nwfs");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixpsax);i=mxAddField(tmp, "pixpsax");mxSetFieldByNumber(tmp, i_n, i, tmp2);
	tmp2=mxCreateDoubleScalar(powfs[i_n].pixpsay);i=mxAddField(tmp, "pixpsay");mxSetFieldByNumber(tmp, i_n, i, tmp2);}
	return tmp;
}

static mxArray *get_powfs_dtf(const DTF_T* dtf){
	if(!dtf) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(dtf->Ux);i=mxAddField(tmp, "Ux");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dtf->Uy);i=mxAddField(tmp, "Uy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dtf->nominal);i=mxAddField(tmp, "nominal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(dtf->si);i=mxAddField(tmp, "si");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(dtf->fused);i=mxAddField(tmp, "fused");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_powfs_etfprep(const ETF_T* etfprep){
	if(!etfprep) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(etfprep->p1);i=mxAddField(tmp, "p1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(etfprep->p2);i=mxAddField(tmp, "p2");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_powfs_etfsim(const ETF_T* etfsim){
	if(!etfsim) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(etfsim->p1);i=mxAddField(tmp, "p1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(etfsim->p2);i=mxAddField(tmp, "p2");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_powfs_intstat(const INTSTAT_T* intstat){
	if(!intstat) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(intstat->cogcoeff);i=mxAddField(tmp, "cogcoeff");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->fotf);i=mxAddField(tmp, "fotf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->gx);i=mxAddField(tmp, "gx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->gy);i=mxAddField(tmp, "gy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->i0);i=mxAddField(tmp, "i0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->i0sum);i=mxAddField(tmp, "i0sum");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->lotf);i=mxAddField(tmp, "lotf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->mtche);i=mxAddField(tmp, "mtche");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->otf);i=mxAddField(tmp, "otf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->saneaixy);i=mxAddField(tmp, "saneaixy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->saneaxy);i=mxAddField(tmp, "saneaxy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->saneaxyl);i=mxAddField(tmp, "saneaxyl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(intstat->sepsf);i=mxAddField(tmp, "sepsf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(intstat->nmtche);i=mxAddField(tmp, "nmtche");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(intstat->notf);i=mxAddField(tmp, "notf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(intstat->nsepsf);i=mxAddField(tmp, "nsepsf");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_powfs_llt(const LLT_T* llt){
	if(!llt) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(llt->amp);i=mxAddField(tmp, "amp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->imcc);i=mxAddField(tmp, "imcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->loc);i=mxAddField(tmp, "loc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->mcc);i=mxAddField(tmp, "mcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->ncpa);i=mxAddField(tmp, "ncpa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(llt->pts);i=mxAddField(tmp, "pts");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon(const RECON_T* recon){
	if(!recon) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(recon->DF);i=mxAddField(tmp, "DF");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GA);i=mxAddField(tmp, "GA");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GAhi);i=mxAddField(tmp, "GAhi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GAlo);i=mxAddField(tmp, "GAlo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GFall);i=mxAddField(tmp, "GFall");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GFlgs);i=mxAddField(tmp, "GFlgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GFngs);i=mxAddField(tmp, "GFngs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GP);i=mxAddField(tmp, "GP");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GP2);i=mxAddField(tmp, "GP2");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GWR);i=mxAddField(tmp, "GWR");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GX);i=mxAddField(tmp, "GX");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GXL);i=mxAddField(tmp, "GXL");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GXfocus);i=mxAddField(tmp, "GXfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GXlo);i=mxAddField(tmp, "GXlo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->GXtomo);i=mxAddField(tmp, "GXtomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->HA);i=mxAddField(tmp, "HA");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->HA_ncpa);i=mxAddField(tmp, "HA_ncpa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->HXF);i=mxAddField(tmp, "HXF");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->HXW);i=mxAddField(tmp, "HXW");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->HXWtomo);i=mxAddField(tmp, "HXWtomo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->L2);i=mxAddField(tmp, "L2");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->L2save);i=mxAddField(tmp, "L2save");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->MVA);i=mxAddField(tmp, "MVA");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->MVFM);i=mxAddField(tmp, "MVFM");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->MVGM);i=mxAddField(tmp, "MVGM");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->MVM);i=mxAddField(tmp, "MVM");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->MVModes);i=mxAddField(tmp, "MVModes");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->MVRngs);i=mxAddField(tmp, "MVRngs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->PDF);i=mxAddField(tmp, "PDF");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->PTT);i=mxAddField(tmp, "PTT");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->PTTF);i=mxAddField(tmp, "PTTF");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFdfa);i=mxAddField(tmp, "RFdfa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFdfx);i=mxAddField(tmp, "RFdfx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFlgsa);i=mxAddField(tmp, "RFlgsa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFlgsg);i=mxAddField(tmp, "RFlgsg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFlgsx);i=mxAddField(tmp, "RFlgsx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFngsa);i=mxAddField(tmp, "RFngsa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFngsg);i=mxAddField(tmp, "RFngsg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->RFngsx);i=mxAddField(tmp, "RFngsx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->TT);i=mxAddField(tmp, "TT");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->TTF);i=mxAddField(tmp, "TTF");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->W0);i=mxAddField(tmp, "W0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->W1);i=mxAddField(tmp, "W1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->ZZT);i=mxAddField(tmp, "ZZT");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->acmap);i=mxAddField(tmp, "acmap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->actcpl);i=mxAddField(tmp, "actcpl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->actfloat);i=mxAddField(tmp, "actfloat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->actinterp);i=mxAddField(tmp, "actinterp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->actslave);i=mxAddField(tmp, "actslave");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->actstuck);i=mxAddField(tmp, "actstuck");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->aimcc);i=mxAddField(tmp, "aimcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->aloc);i=mxAddField(tmp, "aloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->amap);i=mxAddField(tmp, "amap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->anloc);i=mxAddField(tmp, "anloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->anx);i=mxAddField(tmp, "anx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->any);i=mxAddField(tmp, "any");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->dm_ncpa);i=mxAddField(tmp, "dm_ncpa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->ecnn);i=mxAddField(tmp, "ecnn");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->fitNW);i=mxAddField(tmp, "fitNW");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->fitwt);i=mxAddField(tmp, "fitwt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->floc);i=mxAddField(tmp, "floc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->fmap);i=mxAddField(tmp, "fmap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->ht);i=mxAddField(tmp, "ht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->neam);i=mxAddField(tmp, "neam");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->ngrad);i=mxAddField(tmp, "ngrad");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->os);i=mxAddField(tmp, "os");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->ploc);i=mxAddField(tmp, "ploc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->ploc_tel);i=mxAddField(tmp, "ploc_tel");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->pmap);i=mxAddField(tmp, "pmap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->sanea);i=mxAddField(tmp, "sanea");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->saneai);i=mxAddField(tmp, "saneai");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->saneal);i=mxAddField(tmp, "saneal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->wt);i=mxAddField(tmp, "wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xcmap);i=mxAddField(tmp, "xcmap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xloc);i=mxAddField(tmp, "xloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xmap);i=mxAddField(tmp, "xmap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xmcc);i=mxAddField(tmp, "xmcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xnloc);i=mxAddField(tmp, "xnloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xnx);i=mxAddField(tmp, "xnx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(recon->xny);i=mxAddField(tmp, "xny");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_FL(&(recon->FL));i=mxAddField(tmp, "FL");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_FR(&(recon->FR));i=mxAddField(tmp, "FR");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_LL(&(recon->LL));i=mxAddField(tmp, "LL");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_LR(&(recon->LR));i=mxAddField(tmp, "LR");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_RL(&(recon->RL));i=mxAddField(tmp, "RL");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_RR(&(recon->RR));i=mxAddField(tmp, "RR");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_cn2est(recon->cn2est);i=mxAddField(tmp, "cn2est");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_fdpcg(recon->fdpcg);i=mxAddField(tmp, "fdpcg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_fractal(recon->fractal);i=mxAddField(tmp, "fractal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_invpsd(recon->invpsd);i=mxAddField(tmp, "invpsd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_moao(recon->moao);i=mxAddField(tmp, "moao");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_ngsmod(recon->ngsmod);i=mxAddField(tmp, "ngsmod");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->cxx);i=mxAddField(tmp, "cxx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->fitscl);i=mxAddField(tmp, "fitscl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->has_dfr);i=mxAddField(tmp, "has_dfr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->has_ttr);i=mxAddField(tmp, "has_ttr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->l0);i=mxAddField(tmp, "l0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->lowfs_gtilt);i=mxAddField(tmp, "lowfs_gtilt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->ndm);i=mxAddField(tmp, "ndm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->neamhi);i=mxAddField(tmp, "neamhi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->npsr);i=mxAddField(tmp, "npsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->nthread);i=mxAddField(tmp, "nthread");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->r0);i=mxAddField(tmp, "r0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(recon->sigmanlo);i=mxAddField(tmp, "sigmanlo");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_FL(const MUV_T* FL){
	if(!FL) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(FL->M);i=mxAddField(tmp, "M");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->MI);i=mxAddField(tmp, "MI");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->MIB);i=mxAddField(tmp, "MIB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->U);i=mxAddField(tmp, "U");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->Up);i=mxAddField(tmp, "Up");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->UpB);i=mxAddField(tmp, "UpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->V);i=mxAddField(tmp, "V");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->Vp);i=mxAddField(tmp, "Vp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FL->VpB);i=mxAddField(tmp, "VpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FL->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FL->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FL->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FL->nb);i=mxAddField(tmp, "nb");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FL->warm);i=mxAddField(tmp, "warm");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_FR(const MUV_T* FR){
	if(!FR) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(FR->M);i=mxAddField(tmp, "M");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->MI);i=mxAddField(tmp, "MI");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->MIB);i=mxAddField(tmp, "MIB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->U);i=mxAddField(tmp, "U");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->Up);i=mxAddField(tmp, "Up");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->UpB);i=mxAddField(tmp, "UpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->V);i=mxAddField(tmp, "V");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->Vp);i=mxAddField(tmp, "Vp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(FR->VpB);i=mxAddField(tmp, "VpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FR->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FR->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FR->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FR->nb);i=mxAddField(tmp, "nb");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(FR->warm);i=mxAddField(tmp, "warm");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_LL(const MUV_T* LL){
	if(!LL) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(LL->M);i=mxAddField(tmp, "M");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->MI);i=mxAddField(tmp, "MI");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->MIB);i=mxAddField(tmp, "MIB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->U);i=mxAddField(tmp, "U");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->Up);i=mxAddField(tmp, "Up");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->UpB);i=mxAddField(tmp, "UpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->V);i=mxAddField(tmp, "V");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->Vp);i=mxAddField(tmp, "Vp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LL->VpB);i=mxAddField(tmp, "VpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LL->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LL->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LL->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LL->nb);i=mxAddField(tmp, "nb");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LL->warm);i=mxAddField(tmp, "warm");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_LR(const MUV_T* LR){
	if(!LR) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(LR->M);i=mxAddField(tmp, "M");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->MI);i=mxAddField(tmp, "MI");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->MIB);i=mxAddField(tmp, "MIB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->U);i=mxAddField(tmp, "U");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->Up);i=mxAddField(tmp, "Up");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->UpB);i=mxAddField(tmp, "UpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->V);i=mxAddField(tmp, "V");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->Vp);i=mxAddField(tmp, "Vp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(LR->VpB);i=mxAddField(tmp, "VpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LR->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LR->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LR->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LR->nb);i=mxAddField(tmp, "nb");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(LR->warm);i=mxAddField(tmp, "warm");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_RL(const MUV_T* RL){
	if(!RL) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(RL->M);i=mxAddField(tmp, "M");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->MI);i=mxAddField(tmp, "MI");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->MIB);i=mxAddField(tmp, "MIB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->U);i=mxAddField(tmp, "U");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->Up);i=mxAddField(tmp, "Up");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->UpB);i=mxAddField(tmp, "UpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->V);i=mxAddField(tmp, "V");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->Vp);i=mxAddField(tmp, "Vp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RL->VpB);i=mxAddField(tmp, "VpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RL->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RL->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RL->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RL->nb);i=mxAddField(tmp, "nb");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RL->warm);i=mxAddField(tmp, "warm");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_RR(const MUV_T* RR){
	if(!RR) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(RR->M);i=mxAddField(tmp, "M");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->MI);i=mxAddField(tmp, "MI");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->MIB);i=mxAddField(tmp, "MIB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->U);i=mxAddField(tmp, "U");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->Up);i=mxAddField(tmp, "Up");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->UpB);i=mxAddField(tmp, "UpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->V);i=mxAddField(tmp, "V");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->Vp);i=mxAddField(tmp, "Vp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(RR->VpB);i=mxAddField(tmp, "VpB");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RR->alg);i=mxAddField(tmp, "alg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RR->bgs);i=mxAddField(tmp, "bgs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RR->maxit);i=mxAddField(tmp, "maxit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RR->nb);i=mxAddField(tmp, "nb");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(RR->warm);i=mxAddField(tmp, "warm");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_cn2est(const CN2EST_T* cn2est){
	if(!cn2est) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(cn2est->Pnk);i=mxAddField(tmp, "Pnk");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->cov1);i=mxAddField(tmp, "cov1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->cov2);i=mxAddField(tmp, "cov2");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->covc);i=mxAddField(tmp, "covc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->curi);i=mxAddField(tmp, "curi");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->dmht);i=mxAddField(tmp, "dmht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->dx);i=mxAddField(tmp, "dx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->embed);i=mxAddField(tmp, "embed");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->gxs);i=mxAddField(tmp, "gxs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->gys);i=mxAddField(tmp, "gys");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->ht);i=mxAddField(tmp, "ht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->htrecon);i=mxAddField(tmp, "htrecon");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->iPnk);i=mxAddField(tmp, "iPnk");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->mask);i=mxAddField(tmp, "mask");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->os);i=mxAddField(tmp, "os");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->overlap);i=mxAddField(tmp, "overlap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->r0);i=mxAddField(tmp, "r0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->wt);i=mxAddField(tmp, "wt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->wtconvert);i=mxAddField(tmp, "wtconvert");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(cn2est->wtrecon);i=mxAddField(tmp, "wtrecon");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon_cn2est_pair(cn2est->pair);i=mxAddField(tmp, "pair");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->count);i=mxAddField(tmp, "count");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->hmax);i=mxAddField(tmp, "hmax");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->l0);i=mxAddField(tmp, "l0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->nembed);i=mxAddField(tmp, "nembed");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->nsa);i=mxAddField(tmp, "nsa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->nwfs);i=mxAddField(tmp, "nwfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->nwfspair);i=mxAddField(tmp, "nwfspair");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->ovs);i=mxAddField(tmp, "ovs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(cn2est->r0m);i=mxAddField(tmp, "r0m");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_cn2est_pair(const CN2PAIR_T* pair){
	if(!pair) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=mxCreateDoubleScalar(pair->beta);i=mxAddField(tmp, "beta");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->dtheta);i=mxAddField(tmp, "dtheta");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->hsm);i=mxAddField(tmp, "hsm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->iht0);i=mxAddField(tmp, "iht0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->iht1);i=mxAddField(tmp, "iht1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->nsep);i=mxAddField(tmp, "nsep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->wfs0);i=mxAddField(tmp, "wfs0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->wfs1);i=mxAddField(tmp, "wfs1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->xstep);i=mxAddField(tmp, "xstep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(pair->ystep);i=mxAddField(tmp, "ystep");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_fdpcg(const FDPCG_T* fdpcg){
	if(!fdpcg) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(fdpcg->Mbinv);i=mxAddField(tmp, "Mbinv");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fdpcg->Minv);i=mxAddField(tmp, "Minv");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fdpcg->perm);i=mxAddField(tmp, "perm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fdpcg->permhf);i=mxAddField(tmp, "permhf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fdpcg->xloc);i=mxAddField(tmp, "xloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fdpcg->bs);i=mxAddField(tmp, "bs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fdpcg->nbx);i=mxAddField(tmp, "nbx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fdpcg->nby);i=mxAddField(tmp, "nby");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fdpcg->nxtot);i=mxAddField(tmp, "nxtot");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fdpcg->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fdpcg->square);i=mxAddField(tmp, "square");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_fractal(const FRACTAL_T* fractal){
	if(!fractal) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(fractal->xloc);i=mxAddField(tmp, "xloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(fractal->xopd);i=mxAddField(tmp, "xopd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fractal->l0);i=mxAddField(tmp, "l0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fractal->ninit);i=mxAddField(tmp, "ninit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fractal->r0);i=mxAddField(tmp, "r0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(fractal->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_invpsd(const INVPSD_T* invpsd){
	if(!invpsd) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(invpsd->fftxopd);i=mxAddField(tmp, "fftxopd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(invpsd->invpsd);i=mxAddField(tmp, "invpsd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(invpsd->xloc);i=mxAddField(tmp, "xloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(invpsd->square);i=mxAddField(tmp, "square");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_moao(const MOAO_T* moao){
	if(!moao) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(moao->HA);i=mxAddField(tmp, "HA");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->NW);i=mxAddField(tmp, "NW");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->W0);i=mxAddField(tmp, "W0");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->W1);i=mxAddField(tmp, "W1");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->actcpl);i=mxAddField(tmp, "actcpl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->actfloat);i=mxAddField(tmp, "actfloat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->actslave);i=mxAddField(tmp, "actslave");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->actstuck);i=mxAddField(tmp, "actstuck");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->aimcc);i=mxAddField(tmp, "aimcc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->aloc);i=mxAddField(tmp, "aloc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(moao->amap);i=mxAddField(tmp, "amap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(moao->used);i=mxAddField(tmp, "used");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_recon_ngsmod(const NGSMOD_T* ngsmod){
	if(!ngsmod) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(ngsmod->GM);i=mxAddField(tmp, "GM");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->IMCC);i=mxAddField(tmp, "IMCC");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->IMCC_TT);i=mxAddField(tmp, "IMCC_TT");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->MCC);i=mxAddField(tmp, "MCC");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->MCCP);i=mxAddField(tmp, "MCCP");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->Modes);i=mxAddField(tmp, "Modes");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->Pngs);i=mxAddField(tmp, "Pngs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->Ptt);i=mxAddField(tmp, "Ptt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->Rngs);i=mxAddField(tmp, "Rngs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(ngsmod->Wa);i=mxAddField(tmp, "Wa");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(ngsmod->ahstfocus);i=mxAddField(tmp, "ahstfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(ngsmod->aper_fcp);i=mxAddField(tmp, "aper_fcp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(ngsmod->hs);i=mxAddField(tmp, "hs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(ngsmod->ht);i=mxAddField(tmp, "ht");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(ngsmod->nmod);i=mxAddField(tmp, "nmod");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(ngsmod->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_simu(const SIM_T* simu){
	if(!simu) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(simu->cachedm);i=mxAddField(tmp, "cachedm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->cgres);i=mxAddField(tmp, "cgres");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->cle);i=mxAddField(tmp, "cle");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->cleNGSm);i=mxAddField(tmp, "cleNGSm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->cleNGSmp);i=mxAddField(tmp, "cleNGSmp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->clem);i=mxAddField(tmp, "clem");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->clemp);i=mxAddField(tmp, "clemp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->clep);i=mxAddField(tmp, "clep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->clmp);i=mxAddField(tmp, "clmp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->corrNGSm);i=mxAddField(tmp, "corrNGSm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->deltafocus);i=mxAddField(tmp, "deltafocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dm_evl);i=mxAddField(tmp, "dm_evl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dm_wfs);i=mxAddField(tmp, "dm_wfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmadd);i=mxAddField(tmp, "dmadd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmcmd);i=mxAddField(tmp, "dmcmd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmcmdlast);i=mxAddField(tmp, "dmcmdlast");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmfit);i=mxAddField(tmp, "dmfit");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmhist);i=mxAddField(tmp, "dmhist");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmproj);i=mxAddField(tmp, "dmproj");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmprojsq);i=mxAddField(tmp, "dmprojsq");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmpsol);i=mxAddField(tmp, "dmpsol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmreal);i=mxAddField(tmp, "dmreal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->dmrealsq);i=mxAddField(tmp, "dmrealsq");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->ecov);i=mxAddField(tmp, "ecov");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopd);i=mxAddField(tmp, "evlopd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopdcov);i=mxAddField(tmp, "evlopdcov");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopdcov_ngsr);i=mxAddField(tmp, "evlopdcov_ngsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopdcovol);i=mxAddField(tmp, "evlopdcovol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopdmean);i=mxAddField(tmp, "evlopdmean");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopdmean_ngsr);i=mxAddField(tmp, "evlopdmean_ngsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlopdmeanol);i=mxAddField(tmp, "evlopdmeanol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlpsfmean);i=mxAddField(tmp, "evlpsfmean");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlpsfmean_ngsr);i=mxAddField(tmp, "evlpsfmean_ngsr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->evlpsfolmean);i=mxAddField(tmp, "evlpsfolmean");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->gcov);i=mxAddField(tmp, "gcov");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->gngsmvst);i=mxAddField(tmp, "gngsmvst");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->gradacc);i=mxAddField(tmp, "gradacc");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->gradcl);i=mxAddField(tmp, "gradcl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->gradlastcl);i=mxAddField(tmp, "gradlastcl");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->gradlastol);i=mxAddField(tmp, "gradlastol");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->ints);i=mxAddField(tmp, "ints");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->lgsfocus);i=mxAddField(tmp, "lgsfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->lgsfocuslpf);i=mxAddField(tmp, "lgsfocuslpf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->ngsfocuslpf);i=mxAddField(tmp, "ngsfocuslpf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->ole);i=mxAddField(tmp, "ole");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->oleNGSm);i=mxAddField(tmp, "oleNGSm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->oleNGSmp);i=mxAddField(tmp, "oleNGSmp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->olep);i=mxAddField(tmp, "olep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->olmp);i=mxAddField(tmp, "olmp");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->opdevlground);i=mxAddField(tmp, "opdevlground");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->opdr);i=mxAddField(tmp, "opdr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->opdrhat);i=mxAddField(tmp, "opdrhat");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->opdrhatlast);i=mxAddField(tmp, "opdrhatlast");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->opdx);i=mxAddField(tmp, "opdx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->pistatout);i=mxAddField(tmp, "pistatout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->res);i=mxAddField(tmp, "res");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->surf);i=mxAddField(tmp, "surf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->telws);i=mxAddField(tmp, "telws");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->tsurf);i=mxAddField(tmp, "tsurf");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->ttmreal);i=mxAddField(tmp, "ttmreal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->uptcmds);i=mxAddField(tmp, "uptcmds");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->upterrs);i=mxAddField(tmp, "upterrs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->uptreal);i=mxAddField(tmp, "uptreal");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->wfspsfout);i=mxAddField(tmp, "wfspsfout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->winddir);i=mxAddField(tmp, "winddir");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->zoomavg);i=mxAddField(tmp, "zoomavg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->zoomerr);i=mxAddField(tmp, "zoomerr");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->zoomint);i=mxAddField(tmp, "zoomint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(simu->zoompos);i=mxAddField(tmp, "zoompos");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_Mint_lo(simu->Mint_lo);i=mxAddField(tmp, "Mint_lo");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_aper(simu->aper);i=mxAddField(tmp, "aper");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_atmcfg(simu->atmcfg);i=mxAddField(tmp, "atmcfg");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_cachedm_propdata(simu->cachedm_propdata);i=mxAddField(tmp, "cachedm_propdata");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_dmint(simu->dmint);i=mxAddField(tmp, "dmint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_evl_propdata_atm(simu->evl_propdata_atm);i=mxAddField(tmp, "evl_propdata_atm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_evl_propdata_dm(simu->evl_propdata_dm);i=mxAddField(tmp, "evl_propdata_dm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_parms(simu->parms);i=mxAddField(tmp, "parms");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_powfs(simu->powfs, simu->parms->npowfs);i=mxAddField(tmp, "powfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_recon(simu->recon);i=mxAddField(tmp, "recon");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_uptint(simu->uptint);i=mxAddField(tmp, "uptint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_wfs_intsdata(simu->wfs_intsdata);i=mxAddField(tmp, "wfs_intsdata");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_wfs_propdata_atm(simu->wfs_propdata_atm);i=mxAddField(tmp, "wfs_propdata_atm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=get_wfs_propdata_dm(simu->wfs_propdata_dm);i=mxAddField(tmp, "wfs_propdata_dm");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->cachedm_n);i=mxAddField(tmp, "cachedm_n");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->dt);i=mxAddField(tmp, "dt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->iseed);i=mxAddField(tmp, "iseed");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->isim);i=mxAddField(tmp, "isim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->last_report_time);i=mxAddField(tmp, "last_report_time");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->ngsfocus);i=mxAddField(tmp, "ngsfocus");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->perfevl_iground);i=mxAddField(tmp, "perfevl_iground");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->reconisim);i=mxAddField(tmp, "reconisim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->seed);i=mxAddField(tmp, "seed");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->telwsx);i=mxAddField(tmp, "telwsx");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->telwsy);i=mxAddField(tmp, "telwsy");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->tk_cache);i=mxAddField(tmp, "tk_cache");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->tk_eval);i=mxAddField(tmp, "tk_eval");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->tk_recon);i=mxAddField(tmp, "tk_recon");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->tk_wfs);i=mxAddField(tmp, "tk_wfs");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(simu->wfsints_isa);i=mxAddField(tmp, "wfsints_isa");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_uptint(const SERVO_T* uptint){
	if(!uptint) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(uptint->ap);i=mxAddField(tmp, "ap");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(uptint->ep);i=mxAddField(tmp, "ep");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(uptint->merrhist);i=mxAddField(tmp, "merrhist");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(uptint->merrlast);i=mxAddField(tmp, "merrlast");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(uptint->mint);i=mxAddField(tmp, "mint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(uptint->mlead);i=mxAddField(tmp, "mlead");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(uptint->mpreint);i=mxAddField(tmp, "mpreint");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(uptint->al);i=mxAddField(tmp, "al");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(uptint->dt);i=mxAddField(tmp, "dt");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(uptint->initialized);i=mxAddField(tmp, "initialized");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_wfs_intsdata(const WFSINTS_T* wfs_intsdata){
	if(!wfs_intsdata) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(wfs_intsdata->gradref);i=mxAddField(tmp, "gradref");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_intsdata->ints);i=mxAddField(tmp, "ints");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_intsdata->lltopd);i=mxAddField(tmp, "lltopd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_intsdata->opd);i=mxAddField(tmp, "opd");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_intsdata->pistatout);i=mxAddField(tmp, "pistatout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_intsdata->psfout);i=mxAddField(tmp, "psfout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_intsdata->iwfs);i=mxAddField(tmp, "iwfs");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_wfs_propdata_atm(const PROPDATA_T* wfs_propdata_atm){
	if(!wfs_propdata_atm) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(wfs_propdata_atm->locin);i=mxAddField(tmp, "locin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_atm->locout);i=mxAddField(tmp, "locout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_atm->mapin);i=mxAddField(tmp, "mapin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_atm->mapout);i=mxAddField(tmp, "mapout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_atm->ptsout);i=mxAddField(tmp, "ptsout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->alpha);i=mxAddField(tmp, "alpha");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->cubic_iac);i=mxAddField(tmp, "cubic_iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->index);i=mxAddField(tmp, "index");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->nooptim);i=mxAddField(tmp, "nooptim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_atm->wrap);i=mxAddField(tmp, "wrap");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static mxArray *get_wfs_propdata_dm(const PROPDATA_T* wfs_propdata_dm){
	if(!wfs_propdata_dm) return mxCreateDoubleMatrix(0,0,mxREAL);
	mxArray* tmp=mxCreateStructMatrix(1,1,0,0);
	int i;
	mxArray *tmp2;
	tmp2=any2mx(wfs_propdata_dm->locin);i=mxAddField(tmp, "locin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_dm->locout);i=mxAddField(tmp, "locout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_dm->mapin);i=mxAddField(tmp, "mapin");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_dm->mapout);i=mxAddField(tmp, "mapout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=any2mx(wfs_propdata_dm->ptsout);i=mxAddField(tmp, "ptsout");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->alpha);i=mxAddField(tmp, "alpha");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->cubic);i=mxAddField(tmp, "cubic");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->cubic_iac);i=mxAddField(tmp, "cubic_iac");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->index);i=mxAddField(tmp, "index");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->nooptim);i=mxAddField(tmp, "nooptim");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->scale);i=mxAddField(tmp, "scale");mxSetFieldByNumber(tmp, 0, i, tmp2);
	tmp2=mxCreateDoubleScalar(wfs_propdata_dm->wrap);i=mxAddField(tmp, "wrap");mxSetFieldByNumber(tmp, 0, i, tmp2);	return tmp;
}

static void get_data_help(){
	info2("powfs=maos('get','powfs')\n");
	info2("uptint=maos('get','uptint')\n");
	info2("parms=maos('get','parms')\n");
	info2("simu=maos('get','simu')\n");
	info2("atmcfg=maos('get','atmcfg')\n");
	info2("aper=maos('get','aper')\n");
	info2("dmint=maos('get','dmint')\n");
	info2("recon=maos('get','recon')\n");
}
static mxArray *get_data(SIM_T *simu, char *key){
	if(!strcmp(key, "powfs")) return get_powfs(simu->powfs, simu->parms->npowfs);
	else if(!strcmp(key, "uptint")) return get_uptint(simu->uptint);
	else if(!strcmp(key, "parms")) return get_parms(simu->parms);
	else if(!strcmp(key, "simu")) return get_simu(simu);
	else if(!strcmp(key, "atmcfg")) return get_atmcfg(simu->atmcfg);
	else if(!strcmp(key, "aper")) return get_aper(simu->aper);
	else if(!strcmp(key, "dmint")) return get_dmint(simu->dmint);
	else if(!strcmp(key, "recon")) return get_recon(simu->recon);
	else {get_data_help();return mxCreateDoubleMatrix(0,0,mxREAL);}
}
