mxArray *get_simu(SIM_T *simu){
if(!simu) return mxCreateDoubleMatrix(0,0,mxREAL);
int i;
int i_n;
mxArray *tmp=0, *tmp2=0;
int npowfs=simu->parms->npowfs;
mxArray *out=mxCreateStructMatrix(1,1,0,0);
//unknown1 type simu->dmprojsq:map_t*
//unknown1 type simu->init_rand:rand
//unknown1 type simu->dmint:SERVO_T
tmp=dcell2mx(simu->clemp);
i=mxAddField(out, "clemp");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->evlpsfolmean);
i=mxAddField(out, "evlpsfolmean");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->atmwd_rand:rand
//unknown1 type simu->evl_propdata_dm:PROPDATA_T
tmp=mxCreateDoubleScalar(simu->dt);
i=mxAddField(out, "dt");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->ttmreal);
i=mxAddField(out, "ttmreal");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateStructMatrix(1,1,0,0);
tmp2=mxCreateDoubleScalar(simu->aper->fcp);
i=mxAddField(tmp, "fcp");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->aper->ipcc);
i=mxAddField(tmp, "ipcc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->aper->imcc);
i=mxAddField(tmp, "imcc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->aper->opdadd);
i=mxAddField(tmp, "opdadd");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->aper->locs_dm:loc_t*
tmp2=d2mx(simu->aper->mcc);
i=mxAddField(tmp, "mcc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->aper->ampmask://map
tmp2=dcell2mx(simu->aper->opdfloc);
i=mxAddField(tmp, "opdfloc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=loc2mx(simu->aper->locs);
i=mxAddField(tmp, "locs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx((const dmat*)simu->aper->ampground);
i=mxAddField(tmp, "ampground");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->aper->amp);
i=mxAddField(tmp, "amp");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->aper->embed:locfft
tmp2=d2mx(simu->aper->amp1);
i=mxAddField(tmp, "amp1");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->aper->sumamp2);
i=mxAddField(tmp, "sumamp2");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->aper->mod);
i=mxAddField(tmp, "mod");
mxSetFieldByNumber(tmp, 0, i, tmp2);
i=mxAddField(out, "aper");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->perf_evl_pre:thread
tmp=d2mx(simu->cle);
i=mxAddField(out, "cle");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->wfspsfout:ccell*
tmp=dcell2mx(simu->lgsfocus);
i=mxAddField(out, "lgsfocus");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->gcov);
i=mxAddField(out, "gcov");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dmcmdlast);
i=mxAddField(out, "dmcmdlast");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dmcmd);
i=mxAddField(out, "dmcmd");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->cachedm_prop:thread_t*
tmp=mxCreateDoubleScalar(simu->wfsints_isa);
i=mxAddField(out, "wfsints_isa");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dmreal);
i=mxAddField(out, "dmreal");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->evlopd);
i=mxAddField(out, "evlopd");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->cachedm_n);
i=mxAddField(out, "cachedm_n");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->perfevl_iground);
i=mxAddField(out, "perfevl_iground");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->zoompos);
i=mxAddField(out, "zoompos");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->wfs_grad_pre:thread
//unknown1 type simu->wfs_ints:thread_t*
//unknown1 type simu->wfs_propdata_atm:PROPDATA_T
tmp=dcell2mx(simu->cgres);
i=mxAddField(out, "cgres");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->surf:map_t*
tmp=dcell2mx(simu->uptcmds);
i=mxAddField(out, "uptcmds");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->oleNGSmp);
i=mxAddField(out, "oleNGSmp");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->uptint:SERVO_T
//unknown1 type simu->cachedm_propdata:PROPDATA_T
//unknown1 type simu->atm2:map_t*
//unknown1 type simu->tsurf:rectmap_t*
tmp=dcell2mx(simu->clmp);
i=mxAddField(out, "clmp");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->lgsfocuslpf);
i=mxAddField(out, "lgsfocuslpf");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->clem);
i=mxAddField(out, "clem");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=ccell2mx(simu->opdrhat);
i=mxAddField(out, "opdrhat");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->cachedm:map_t**
//unknown1 type simu->wfs_grad_post:thread
tmp=dcell2mx(simu->dmadd);
i=mxAddField(out, "dmadd");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->res);
i=mxAddField(out, "res");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->opdevlground);
i=mxAddField(out, "opdevlground");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->telws_rand:rand
tmp=dcell2mx(simu->evlopdcov);
i=mxAddField(out, "evlopdcov");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->evlpsfmean_ngsr);
i=mxAddField(out, "evlpsfmean_ngsr");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->evl_propdata_atm:PROPDATA_T
tmp=mxCreateStructMatrix(1,1,0,0);
tmp2=dcell2mx(simu->wfs_intsdata->pistatout);
i=mxAddField(tmp, "pistatout");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->wfs_intsdata->gradref);
i=mxAddField(tmp, "gradref");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->wfs_intsdata->iwfs);
i=mxAddField(tmp, "iwfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=ccell2mx(simu->wfs_intsdata->psfout);
i=mxAddField(tmp, "psfout");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->wfs_intsdata->opd);
i=mxAddField(tmp, "opd");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->wfs_intsdata->lltopd);
i=mxAddField(tmp, "lltopd");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->wfs_intsdata->ints);
i=mxAddField(tmp, "ints");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->wfs_intsdata->parms:PARMS_T
//unknown1 type simu->wfs_intsdata->powfs:POWFS_T
i=mxAddField(out, "wfs_intsdata");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->evlopdmean);
i=mxAddField(out, "evlopdmean");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->winddir);
i=mxAddField(out, "winddir");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->wfs_prop_atm:thread_t*
tmp=dcell2mx(simu->gngsmvst);
i=mxAddField(out, "gngsmvst");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->pistatout:dcell*
tmp=mxCreateDoubleScalar(simu->reconisim);
i=mxAddField(out, "reconisim");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dmfit);
i=mxAddField(out, "dmfit");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->evlopdcovol);
i=mxAddField(out, "evlopdcovol");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateStructMatrix(1,1,0,0);
tmp2=d2mx(simu->recon->dx);
i=mxAddField(tmp, "dx");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->simu:SIM_T
//unknown1 type simu->recon->ngsmod:NGSMOD_T
//unknown1 type simu->recon->aper:APER_T
//unknown2 type simu->recon->RR:MUV_T
tmp2=mxCreateDoubleScalar(simu->recon->nthread);
i=mxAddField(tmp, "nthread");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->cxx);
i=mxAddField(tmp, "cxx");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->recon->RL:MUV_T
//unknown2 type simu->recon->LR:MUV_T
//unknown1 type simu->recon->GP2:spcell
//unknown1 type simu->recon->GP:spcell
//unknown1 type simu->recon->GX:spcell
tmp2=dcell2mx(simu->recon->RFdfa);
i=mxAddField(tmp, "RFdfa");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->GA:spcell
//unknown1 type simu->recon->fractal:FRACTAL_T
tmp2=mxCreateDoubleScalar(simu->recon->has_ttr);
i=mxAddField(tmp, "has_ttr");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->GXtomo:spcell
tmp2=dcell2mx(simu->recon->GXL);
i=mxAddField(tmp, "GXL");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->ndm);
i=mxAddField(tmp, "ndm");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->PTTF);
i=mxAddField(tmp, "PTTF");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->cn2est:CN2EST_T
tmp2=d2mx(simu->recon->W1);
i=mxAddField(tmp, "W1");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->W0:dsp
//unknown1 type simu->recon->xmap:map_t*
tmp2=dcell2mx(simu->recon->aimcc);
i=mxAddField(tmp, "aimcc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->actcpl);
i=mxAddField(tmp, "actcpl");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->neamhi);
i=mxAddField(tmp, "neamhi");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->HA:spcell
tmp2=dcell2mx(simu->recon->GFall);
i=mxAddField(tmp, "GFall");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->npsr);
i=mxAddField(tmp, "npsr");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->saneai:spcell
//unknown1 type simu->recon->saneal:spcell
tmp2=mxCreateDoubleScalar(simu->recon->lowfs_gtilt);
i=mxAddField(tmp, "lowfs_gtilt");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->parms:PARMS_T
tmp2=dcell2mx(simu->recon->ecnn);
i=mxAddField(tmp, "ecnn");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->MVGM);
i=mxAddField(tmp, "MVGM");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->PTT);
i=mxAddField(tmp, "PTT");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->TTF);
i=mxAddField(tmp, "TTF");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx((const dmat*)simu->recon->fmap);
i=mxAddField(tmp, "fmap");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->l0);
i=mxAddField(tmp, "l0");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->acmap:map_t*
tmp2=dcell2mx(simu->recon->RFngsa);
i=mxAddField(tmp, "RFngsa");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->anloc:long
tmp2=dcell2mx(simu->recon->RFngsg);
i=mxAddField(tmp, "RFngsg");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->fitscl);
i=mxAddField(tmp, "fitscl");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->RFngsx);
i=mxAddField(tmp, "RFngsx");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->recon->wt);
i=mxAddField(tmp, "wt");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->MVRngs);
i=mxAddField(tmp, "MVRngs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->GAlo:spcell
//unknown1 type simu->recon->aloc:loc_t*
tmp2=mxCreateDoubleScalar(simu->recon->r0);
i=mxAddField(tmp, "r0");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->xny:long
//unknown1 type simu->recon->xnx:long
//unknown1 type simu->recon->actfloat:icell
//unknown1 type simu->recon->invpsd:INVPSD_T
//unknown1 type simu->recon->L2:spcell
//unknown1 type simu->recon->xcmap:map_t*
tmp2=d2mx(simu->recon->os);
i=mxAddField(tmp, "os");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->GFlgs);
i=mxAddField(tmp, "GFlgs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->recon->FR:MUV_T
//unknown1 type simu->recon->L2save:spcell
//unknown1 type simu->recon->actstuck:icell
//unknown1 type simu->recon->GWR:spcell
//unknown1 type simu->recon->xnloc:long
//unknown2 type simu->recon->FL:MUV_T
//unknown1 type simu->recon->moao:MOAO_T
//unknown1 type simu->recon->HXWtomo:spcell
//unknown1 type simu->recon->ngrad:long
//unknown1 type simu->recon->ZZT:spcell
tmp2=dcell2mx(simu->recon->MVModes);
i=mxAddField(tmp, "MVModes");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=loc2mx(simu->recon->ploc);
i=mxAddField(tmp, "ploc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->fdpcg:FDPCG_T
//unknown1 type simu->recon->sanea:spcell
tmp2=d2mx((const dmat*)simu->recon->pmap);
i=mxAddField(tmp, "pmap");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->recon->ht);
i=mxAddField(tmp, "ht");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->MVA);
i=mxAddField(tmp, "MVA");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->dm_ncpa);
i=mxAddField(tmp, "dm_ncpa");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->sigmanlo);
i=mxAddField(tmp, "sigmanlo");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->PDF);
i=mxAddField(tmp, "PDF");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->ploc_tel:loc_t*
tmp2=d2mx(simu->recon->MVM);
i=mxAddField(tmp, "MVM");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->recon->neam);
i=mxAddField(tmp, "neam");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->MVFM);
i=mxAddField(tmp, "MVFM");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->DF);
i=mxAddField(tmp, "DF");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->xloc:loc_t*
tmp2=dcell2mx(simu->recon->GFngs);
i=mxAddField(tmp, "GFngs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->recon->has_dfr);
i=mxAddField(tmp, "has_dfr");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->anx:long
//unknown1 type simu->recon->any:long
//unknown1 type simu->recon->powfs:POWFS_T
//unknown2 type simu->recon->LL:MUV_T
tmp2=dcell2mx(simu->recon->TT);
i=mxAddField(tmp, "TT");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->fitNW);
i=mxAddField(tmp, "fitNW");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->GXfocus:spcell
//unknown1 type simu->recon->HXF:spcell
//unknown1 type simu->recon->GXlo:spcell
//unknown1 type simu->recon->HA_ncpa:spcell
//unknown1 type simu->recon->GAhi:spcell
tmp2=dcell2mx(simu->recon->RFdfx);
i=mxAddField(tmp, "RFdfx");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->HXW:spcell
//unknown1 type simu->recon->amap:map_t*
tmp2=dcell2mx(simu->recon->RFlgsa);
i=mxAddField(tmp, "RFlgsa");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->RFlgsx);
i=mxAddField(tmp, "RFlgsx");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->xmcc);
i=mxAddField(tmp, "xmcc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=loc2mx(simu->recon->floc);
i=mxAddField(tmp, "floc");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=d2mx(simu->recon->fitwt);
i=mxAddField(tmp, "fitwt");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=dcell2mx(simu->recon->RFlgsg);
i=mxAddField(tmp, "RFlgsg");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->recon->actslave:spcell
//unknown1 type simu->recon->actinterp:spcell
i=mxAddField(out, "recon");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->dmpsol:dcell*
//unknown1 type simu->atm:map_t*
tmp=dcell2mx(simu->dmhist);
i=mxAddField(out, "dmhist");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->gradlastcl);
i=mxAddField(out, "gradlastcl");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->ngsfocus);
i=mxAddField(out, "ngsfocus");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->tk_eval);
i=mxAddField(out, "tk_eval");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->ints:dcell*
tmp=mxCreateDoubleScalar(simu->iseed);
i=mxAddField(out, "iseed");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dm_wfs);
i=mxAddField(out, "dm_wfs");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->misc_rand:rand
tmp=mxCreateStructMatrix(1,1,0,0);
//unknown2 type simu->parms->load:LOAD_CFG_T
tmp2=d2mx(simu->parms->dirs);
i=mxAddField(tmp, "dirs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->parms->ntrpowfs);
i=mxAddField(tmp, "ntrpowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->parms->nwfs);
i=mxAddField(tmp, "nwfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->parms->fdlock:imat
//unknown2 type simu->parms->dbg:DBG_CFG_T
//unknown1 type simu->parms->hipowfs:imat
//unknown2 type simu->parms->lsr:LSR_CFG_T
//unknown1 type simu->parms->wfsr:WFS_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->nsurf);
i=mxAddField(tmp, "nsurf");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->parms->nlopowfs);
i=mxAddField(tmp, "nlopowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->parms->powfs:POWFS_CFG_T
//unknown1 type simu->parms->moao:MOAO_CFG_T
//unknown2 type simu->parms->plot:PLOT_CFG_T
//unknown2 type simu->parms->aper:APER_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->nphypowfs);
i=mxAddField(tmp, "nphypowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->parms->fit:FIT_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->nlowfs);
i=mxAddField(tmp, "nlowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->parms->misreg:MISREG_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->ndm);
i=mxAddField(tmp, "ndm");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->parms->gpu:GPU_CFG_T
//unknown2 type simu->parms->save:SAVE_CFG_T
//unknown2 type simu->parms->sim:SIM_CFG_T
//unknown1 type simu->parms->surf:char*
//unknown1 type simu->parms->dm:DM_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->ntsurf);
i=mxAddField(tmp, "ntsurf");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->parms->nhipowfs);
i=mxAddField(tmp, "nhipowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
tmp2=mxCreateDoubleScalar(simu->parms->nwfsr);
i=mxAddField(tmp, "nwfsr");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->parms->atmr:ATMR_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->ntipowfs);
i=mxAddField(tmp, "ntipowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown1 type simu->parms->lopowfs:imat
//unknown2 type simu->parms->cn2:CN2EST_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->npowfs);
i=mxAddField(tmp, "npowfs");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->parms->evl:EVL_CFG_T
tmp2=mxCreateDoubleScalar(simu->parms->nmoao);
i=mxAddField(tmp, "nmoao");
mxSetFieldByNumber(tmp, 0, i, tmp2);
//unknown2 type simu->parms->recon:RECON_CFG_T
//unknown2 type simu->parms->atm:ATM_CFG_T
//unknown1 type simu->parms->tsurf:char*
//unknown2 type simu->parms->tomo:TOMO_CFG_T
//unknown1 type simu->parms->wfs:WFS_CFG_T
i=mxAddField(out, "parms");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->ole);
i=mxAddField(out, "ole");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->hyst:HYST_T*
tmp=d2mx(simu->cleNGSm);
i=mxAddField(out, "cleNGSm");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->cleNGSmp);
i=mxAddField(out, "cleNGSmp");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->evlopdmeanol);
i=mxAddField(out, "evlopdmeanol");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->clep);
i=mxAddField(out, "clep");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->zoomavg);
i=mxAddField(out, "zoomavg");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=ccell2mx(simu->opdrhatlast);
i=mxAddField(out, "opdrhatlast");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->perf_evl_post:thread
tmp=dcell2mx(simu->ecov);
i=mxAddField(out, "ecov");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->oleNGSm);
i=mxAddField(out, "oleNGSm");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->olmp);
i=mxAddField(out, "olmp");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->olep);
i=mxAddField(out, "olep");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->dmrealsq:map_t*
tmp=dcell2mx(simu->evlpsfmean);
i=mxAddField(out, "evlpsfmean");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->isim);
i=mxAddField(out, "isim");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->dither:DITHER_T*
tmp=mxCreateDoubleScalar(simu->telwsy);
i=mxAddField(out, "telwsy");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->telwsx);
i=mxAddField(out, "telwsx");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dm_evl);
i=mxAddField(out, "dm_evl");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->gradcl);
i=mxAddField(out, "gradcl");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->zoomerr);
i=mxAddField(out, "zoomerr");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->telws);
i=mxAddField(out, "telws");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->seed);
i=mxAddField(out, "seed");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->wfs_rand:rand
tmp=dcell2mx(simu->ngsfocuslpf);
i=mxAddField(out, "ngsfocuslpf");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->evl_prop_dm:thread_t*
tmp=mxCreateDoubleScalar(simu->tk_recon);
i=mxAddField(out, "tk_recon");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->tk_cache);
i=mxAddField(out, "tk_cache");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateStructMatrix(2,1,0,0);
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].loc_tel:loc_t*
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].nwfs);
i=mxAddField(tmp, "nwfs");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].loc_dm:loc_t*
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].realamp);
i=mxAddField(tmp, "realamp");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].saa);
i=mxAddField(tmp, "saa");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].saa_tel);
i=mxAddField(tmp, "saa_tel");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].focus);
i=mxAddField(tmp, "focus");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].llt:LLT_T
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].nsaimcc);
i=mxAddField(tmp, "nsaimcc");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].sprint);
i=mxAddField(tmp, "sprint");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].pts:pts
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].saimcc:dcell*
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].sumamp);
i=mxAddField(tmp, "sumamp");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=loc2mx(simu->powfs[i_n].loc);
i=mxAddField(tmp, "loc");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].opdadd);
i=mxAddField(tmp, "opdadd");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].sodium);
i=mxAddField(tmp, "sodium");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].nthread);
i=mxAddField(tmp, "nthread");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].etfprep:ETF_T
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].realsaa);
i=mxAddField(tmp, "realsaa");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].srsa);
i=mxAddField(tmp, "srsa");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].fieldstop:locfft
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].ncompy);
i=mxAddField(tmp, "ncompy");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].ncompx);
i=mxAddField(tmp, "ncompx");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].gradoff);
i=mxAddField(tmp, "gradoff");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].srot);
i=mxAddField(tmp, "srot");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].amp_tel);
i=mxAddField(tmp, "amp_tel");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].bkgrnd);
i=mxAddField(tmp, "bkgrnd");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].sumamp2);
i=mxAddField(tmp, "sumamp2");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].bkgrndc);
i=mxAddField(tmp, "bkgrndc");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].neasim);
i=mxAddField(tmp, "neasim");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].intstat:INTSTAT_T
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].gamp);
i=mxAddField(tmp, "gamp");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=loc2mx(simu->powfs[i_n].gloc);
i=mxAddField(tmp, "gloc");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].srsamax);
i=mxAddField(tmp, "srsamax");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].dtheta);
i=mxAddField(tmp, "dtheta");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=d2mx(simu->powfs[i_n].amp);
i=mxAddField(tmp, "amp");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=loc2mx(simu->powfs[i_n].saloc);
i=mxAddField(tmp, "saloc");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].GS0:spcell
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].gradphyoff);
i=mxAddField(tmp, "gradphyoff");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].etfsim:ETF_T
}
for(i_n=0; i_n<npowfs; i_n++){
//unknown1 type simu->powfs[i_n].dtf:DTF_T
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=dcell2mx(simu->powfs[i_n].opdbias);
i=mxAddField(tmp, "opdbias");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].namp);
i=mxAddField(tmp, "namp");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].pixpsax);
i=mxAddField(tmp, "pixpsax");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].pixpsay);
i=mxAddField(tmp, "pixpsay");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].areascale);
i=mxAddField(tmp, "areascale");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
for(i_n=0; i_n<npowfs; i_n++){
tmp2=mxCreateDoubleScalar(simu->powfs[i_n].npts);
i=mxAddField(tmp, "npts");
mxSetFieldByNumber(tmp, i_n, i, tmp2);
}
i=mxAddField(out, "powfs");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->dmproj);
i=mxAddField(out, "dmproj");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=mxCreateDoubleScalar(simu->last_report_time);
i=mxAddField(out, "last_report_time");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->genscreen:GENSCREEN_T
tmp=mxCreateDoubleScalar(simu->tk_wfs);
i=mxAddField(out, "tk_wfs");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->Mint_lo:SERVO_T
tmp=dcell2mx(simu->evlopdcov_ngsr);
i=mxAddField(out, "evlopdcov_ngsr");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->evlopdmean_ngsr);
i=mxAddField(out, "evlopdmean_ngsr");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->atm_rand:rand
//unknown1 type simu->status:STATUS_T
tmp=dcell2mx(simu->opdx);
i=mxAddField(out, "opdx");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->wfs_propdata_dm:PROPDATA_T
tmp=dcell2mx(simu->gradacc);
i=mxAddField(out, "gradacc");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->opdr);
i=mxAddField(out, "opdr");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->wfs_prop_dm:thread_t*
tmp=dcell2mx(simu->upterrs);
i=mxAddField(out, "upterrs");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->gradlastol);
i=mxAddField(out, "gradlastol");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=dcell2mx(simu->uptreal);
i=mxAddField(out, "uptreal");
mxSetFieldByNumber(out, 0, i, tmp);
//unknown1 type simu->evl_prop_atm:thread_t*
tmp=dcell2mx(simu->deltafocus);
i=mxAddField(out, "deltafocus");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->zoomint);
i=mxAddField(out, "zoomint");
mxSetFieldByNumber(out, 0, i, tmp);
tmp=d2mx(simu->corrNGSm);
i=mxAddField(out, "corrNGSm");
mxSetFieldByNumber(out, 0, i, tmp);
return out;
}

