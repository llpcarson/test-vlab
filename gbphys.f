!> @file gbphys.f This file contains the gbphys subroutine.

!> @defgroup gbphys GFS Physics Implementation Layer
!> @brief Layer that invokes individual GFS physics routines
!> @{
             !
!at tune step========================================================= !
!  description:                                                         !
!
!                                                                      !
!  usage:                                                               !
!                                                                       !
!    call gbphys                                                        !
!       inputs:                                                         !
!         ( im,ix,levs,lsoil,lsm,ntrac,ncld,ntoz,ntcw,ntke,             !
!           ntiw,ntlnc,ntinc,                                           !
!           nmtvr,nrcm,ko3,lonr,latr,jcap,                              !
!           num_p3d,num_p2d,npdf3d,ncnvcld3d,                           !
!           kdt,lat,me,pl_coeff,nlons,ncw,flgmin,crtrh,cdmbgwd,         !
!           ccwf,dlqf,ctei_rm,clstp,cgwf,prslrd0,ral_ts,dtp,dtf,fhour,  !
!           solhr,slag,sdec,cdec,sinlat,coslat,pgr,ugrs,vgrs,           !
!           tgrs,qgrs,vvel,prsi,prsl,prslk,prsik,phii,phil,             !
!           rann,prdoz,poz,dpshc,hprime,xlon,xlat,                      !
!           h2o_phys,levh2o,h2opl,h2o_pres,h2o_coeff,                   !
!           isot,ivegsrc,                                               !
!           slope,shdmin,shdmax,snoalb,tg3,slmsk,vfrac,                 !
!           vtype,stype,uustar,oro,oro_uf,coszen,sfcdsw,sfcnsw,         !
!           sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                    !
!           sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                    !
!           slimskin_cpl,ulwsfcin_cpl,                                  !
!           dusfcin_cpl,dvsfcin_cpl,dtsfcin_cpl,dqsfcin_cpl,            !
!           sfcdlw,tsflw,sfcemis,sfalb,swh,swhc,hlw,hlwc,hlwd,lsidea,   !
!           ras,pre_rad,ldiag3d,lgocart,lssav,lssav_cpl                 !
!           xkzm_m,xkzm_h,xkzm_s,psautco,prautco,evpco,wminco,          !
!           pdfcld,shcnvcw,sup,redrag,hybedmf,dspheat,                  !
!           flipv,old_monin,cnvgwd,shal_cnv,                            !
!           imfshalcnv,imfdeepcnv,cal_pre,aero_in,                      !
!           mom4ice,mstrat,trans_trac,nstf_name,moist_adj,               !
!           thermodyn_id, sfcpress_id, gen_coord_hybrid,levr,adjtrc,nnp,!
!           cscnv,nctp,do_shoc,shocaftcnv,ntot3d,ntot2d,                !
!       input/outputs:                                                  !
!           hice,fice,tisfc,tsea,tprcp,cv,cvb,cvt,                      !
!           srflag,snwdph,weasd,sncovr,zorl,canopy,                     !
!           ffmm,ffhh,f10m,srunoff,evbsa,evcwa,snohfa,                  !
!           transa,sbsnoa,snowca,soilm,tmpmin,tmpmax,                   !
!           dusfc,dvsfc,dtsfc,dqsfc,totprcp,gflux,                      !
!           dlwsfc,ulwsfc,suntim,runoff,ep,cldwrk,                      !
!           dugwd,dvgwd,psmean,cnvprcp,spfhmin,spfhmax,rain,rainc,      !
!           dt3dt,dq3dt,du3dt,dv3dt,dqdt_v,cnvqc_v,acv,acvb,acvt,       !
!           slc,smc,stc,upd_mf,dwn_mf,det_mf,phy_f3d,phy_f2d,           !
!           dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl,                 !
!           dlwsfc_cpl,dswsfc_cpl,dnirbm_cpl,dnirdf_cpl,                !
!           dvisbm_cpl,dvisdf_cpl,rain_cpl,  nlwsfc_cpl,nswsfc_cpl,     !
!           nnirbm_cpl,nnirdf_cpl,nvisbm_cpl,nvisdf_cpl,snow_cpl,       !
!           xt,xs,xu,xv,xz,zm,xtts,xzts,d_conv,ifd,dt_cool,Qrain,       !
!           phy_fctd,                                                   !
!       outputs:                                                        !
!           gt0,gq0,gu0,gv0,t2m,q2m,u10m,v10m,                          !
!           zlvl,psurf,hpbl,pwat,t1,q1,u1,v1,                           !
!           chh,cmm,dlwsfci,ulwsfci,dswsfci,uswsfci,dusfci,dvsfci,      !
!           dtsfci,dqsfci,gfluxi,epi,smcwlt2,smcref2,wet1,              !
!           dusfci_cpl,dvsfci_cpl,dtsfci_cpl,dqsfci_cpl,                !
!           dlwsfci_cpl,dswsfci_cpl,                                    !
!           dnirbmi_cpl,dnirdfi_cpl,dvisbmi_cpl,dvisdfi_cpl,            !
!           nlwsfci_cpl,nswsfci_cpl,                                    !
!           nnirbmi_cpl,nnirdfi_cpl,nvisbmi_cpl,nvisdfi_cpl,            !
!           t2mi_cpl,q2mi_cpl,                                          !
!           u10mi_cpl,v10mi_cpl,tseai_cpl,psurfi_cpl,                   !
!           tref, z_c, c_0, c_d, w_0, w_d, rqtk,                        !
!           hlwd,lsidea                         )                       !
!                                                                       !
!  subprograms called:                                                  !
!                                                                       !
!     get_prs,  dcyc2t2_pre_rad (testing),    dcyc2t3,  sfc_diff,       !
!     sfc_ocean,sfc_drv,  sfc_land, sfc_sice, sfc_diag, moninp1,        !
!     moninp,   moninq1,  moninq,   gwdps,    ozphys,   get_phi,        !
!     sascnv,   sascnvn,  rascnv,   cs_convr, gwdc,     shalcvt3,shalcv,!
!     shalcnv,  cnvc90,   lrgscl,   gsmdrive, gscond,   precpd,         !
!     progt2.                                                           !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!           19xx  - ncep mrf/gfs                                        !
!           2002  - s. moorthi  modify and restructure and add Ferrier  !
!                               microphysics as an option               !
!           200x  - h. juang    modify (what?)                                  !
!      nov  2004  - x. wu       modify sea-ice model                    !
!      may  2005  - s. moorthi  modify and restructure                  !
!           2005  - s. lu       modify to include noah lsm              !
!      oct  2006  - h. wei      modify lsm options to include both      !
!                               noah and osu lsms.                      !
!           2006  - s. moorthi  added a. johansson's convective gravity !
!                               wave parameterization code              !
!           2007  - s. moorthi  added j. han's modified pbl/sas options !
!      dec  2007  - xu li       modified the operational version for    !
!                               nst model                               !
!           2008  - s. moorthi  applied xu li's nst model to new gfs    !
!      mar  2008  - y.-t. hou   added sunshine duration var (suntim) as !
!                     an input/output argument.                         !
!           2008  - jun wang    added spfhmax/spfhmin as input/output.  !
!      apr  2008  - y.-t. hou   added lw sfc emissivity var (sfcemis),  !
!                     define the lw sfc dn/up fluxes in two forms: atmos!
!                     and ground. also changed sw sfc net flux direction!
!                     (positive) from ground -> atmos to the direction  !
!                     of atmos -> ground. recode the program and add    !
!                     program documentation block.
!           2008/ - s. moorthi and y.t. hou upgraded the code to more   !
!           2009    modern form and changed all the inputs to MKS units.!
!      feb  2009  - s. moorthi  upgraded to add Hochun's gocart changes !
!      jul  2009  - s. moorthi  added rqtk for sela's semi-lagrangian   !
!      aug  2009  - s. moorthi  added j. han and h. pan updated shallow !
!                               convection package                      !
!      sep  2009  - s. moorthi  updated for the mcica (rrtm3) radiation !
!      dec  2010  - sarah lu    lgocart added to input arg;             !
!                               compute dqdt_v if inline gocart is on   !
!      feb  2011  - sarah lu    add the option to update surface diag   !
!                               fields (t2m,q2m,u10m,v10m) at the end   !
!      Jun  2011  - s. moorthi and Xu Li - updated the nst model        !
!                               !
!      sep  2011  - sarah lu    correct dqdt_v calculations             !
!      apr  2012  - henry juang add idea                                !
!      sep  2012  - s. moorthi  merge with operational version          !
!      Mar  2013  - Jun Wang    set idea heating rate to tmp tendency   !
!      May  2013  - Jun Wang    tmp updated after idea phys             !
!      Jun  2013  - s. moorthi  corrected a bug in 3d diagnostics for T !
!      Aug  2013  - s. moorthi updating J. Whitekar's changes related   !
!                              to stochastic physics perturnbation      !
!      Oct  2013  - Xingren Wu  add dusfci/dvsfci                       !
!      Mar  2014  - Xingren Wu  add "_cpl" for coupling                 !
!      Mar  2014  - Xingren Wu  add "nir/vis beam and nir/vis diff"     !
!      Apr  2014  - Xingren Wu  add "NET LW/SW including nir/vis"       !
!      Jan  2014  - Jun Wang    merge Moorthi's gwdc change and H.Juang !
!                               and F. Yang's energy conversion from GWD!
!      jan  2014  - y-t hou     revised sw sfc spectral component fluxes!
!                     for coupled mdl, added estimation of ocean albedo !
!                     without ice contamination.                        !
!      Jun  2014  - Xingren Wu  update net SW fluxes over the ocean     !
!                               (no ice contamination)                  !
!      Jul  2014  - Xingren Wu  add Sea/Land/Ice Mask - slmsk_cpl       !
!      Jul  2014  - s. moorthi  merge with standalone gfs and cleanup   !
!      Aug  2014  - s. moorthi  add tracer fixer                        !
!      Sep  2014  - Sarah Lu    disable the option to compute tracer    !
!                               scavenging in GFS phys (set fscav=0.)   !
!      Dec  2014  - Jun Wang    add cnvqc_v for gocart                  !

!  ====================  defination of variables  ====================  !
!      ---  2014  - D. Dazlich  Added Chikira-Sugiyama (CS) convection  !
!                               as an option in opr GFS.                !
!      Apr  2015    S. Moorthi  Added CS scheme to NEMS/GSM             !
!      Jun  2015    S. Moorthi  Added SHOC  to NEMS/GSM                 !
!      Aug  2015  - Xu  Li      change nst_fcst to be nstf_name         !
!                               and introduce depth mean SST            !
!      Sep  2015  - Xingren Wu  remove oro_cpl & slmsk_cpl              !
!      Sep  2015  - Xingren Wu  add sfc_cice                            !
!      Sep  2015  - Xingren Wu  connect CICE output to sfc_cice         !
!      Jan  2016  - P. Tripp    NUOPC/GSM merge                         !
!      Mar  2016  - J. Han  -   add ncnvcld3d integer                   !
!                               for convective cloudiness enhancement   !
!      Mar  2016  - J. Han  -   change newsas & sashal to imfdeepcnv    !
!                            &  imfshalcnv, respectively                !
!      Mar  2016    F. Yang     add pgr to rayleigh damping call        !
!      Mar  2016    S. Moorthi  add ral_ts                              !
!      Mar  2016    Ann Cheng   add morrison 2m microphysics (gsfc)     !
!      May  2016    S. Moorthi  cleanup 2m microphysics implementation  !
!      Jun  2016    X. Li       change all nst_fld as inout             !
!      jul  2016    S. Moorthi  fix some bugs in shoc/2m microphysics   !
!
!  ====================    end of description    =====================
!  ====================  definition of variables  ====================  !

!> @details This subroutine is the original driver for the GFS atmospheric physics.
!! It is responsible for calculating and applying tendencies of the atmospheric state
!! variables due to the atmospheric physics and due to the surface layer scheme. In addition,
!! this routine applies radiative heating rates that were calculated during the
!! antecedent call to the radiation scheme. Code within this subroutine is executed on the
!! physics sub-timestep. The sub-timestep loop is executed in the subroutine gloopb.
!!
!! Parameter descriptions include intent, name, description, and size
!! @param [in]     ix(1), im(1)   - integer, horiz dimension and num of used pts
!! @param [in]     levs(1)     - integer, vertical layer dimension
!! @param [in]     lsoil(1)    - integer, number of soil layers
!! @param [in]     lsm(1)      - integer, flag for land surface model to use
!!                            =0  for osu lsm; =1  for noah lsm
!! @param [in]     ntrac(1)    - integer, number of tracers
!! @param [in]     ncld(1)     - integer, number of cloud species
!! @param [in]     ntoz(1)     - integer, ozone location in the tracer array
!! @param [in]     ntcw(1)     - integer, cloud condensate location in the tracer
!!                                     array
!! @param [in]     ntke(1)     - integer, tke location in the tracer array
!! @param [in]     nmtvr(1)    - integer, number of topographic variables such as
!!                                     variance etc used in the GWD parameterization
!! @param [in]     nrcm(1)     - integer, second dimension for the random number
!!                                     array rann
!! @param [in]     ko3(1)      - integer, number of layers for ozone data
!! @param [in]     lonr(1),latr(1) - integer, number of lon/lat points
!! @param [in]     jcap(1)     - integer, number of spectral wave trancation
!!                                     used only by sascnv shalcnv
!! @param [in]     num_p3d(1)  - integer, number of 3D arrays needed for
!!                                      microphysics
!! @param [in]     num_p2d(1)  - integer, number of 2D arrays needed for
!!                                     microphysics
!! @param [in]     npdf3d(1)   - integer, number of 3d arrays associated with pdf
!!                                     based clouds/microphysics
!! @param [in]     ncnvcld3d(1) - integer, number of 3d arrays associated with
!!                                     convective cloudiness enhancement
!! @param [in]     kdt(1)       - integer, number of the current time step
!! @param [in]     lat(1)       - integer, latitude index - used for debug prints
!! @param [in]     me(1)        - integer, pe number - used for debug prints
!! @param [in]     pl_coeff(1) - integer, number coefficients in ozone forcing
!! @param [in]     nlons(im)    - integer, number of total grid points in a latitude
!!                                     circle through a point
!! @param [in]     ncw(2)      - integer, range of droplet number concentrations for
!!                                     Ferrier microphysics
!! @param [in]     flgmin(2)   - real, range of  minimum large ice fraction for
!!                                     Ferrier microphys
!! @param [in]     crtrh(3)    - real, critical relative humidity at the surface, PBL
!!                                  top and at the top of the atmosphere
!! @param [in]     cdmbgwd(2)  - real, multiplication factors for cdmb and gwd
!! @param [in]     ccwf(2)     - real, multiplication factor for critical cloud
!!                                  workfunction for RAS
!! @param [in]     dlqf(2)     - real, factor for cloud condensate detrainment from
!!                                  cloud edges (RAS)
!! @param [in]     ctei_rm(2)  - real, critical cloud top entrainment instability
!!                                  criteria (used if mstrat=.true.)
!! @param [in]     clstp(1)    - real, index used by cnvc90 (for convective clouds)
!!                                  legacy stuff - does not affect forecast
!! @param [in]     cgwf(2)     - real, multiplication factor for convective GWD
!! @param [in]     prslrd0(1)  - real, pressure level (Pa) from which Rayleigh Damping
!!                                  is applied
!! @param [in]     ral_ts(1)   - real  time scale for Rayleigh damping in days
!! @param [in]     dtp(1)      - real, physics time step in seconds
!! @param [in]     dtf(1)      - real, dynamics time step in seconds
!! @param [in]     fhour(1)    - real, forecast hour
!! @param [in]     solhr(1)    - real, fcst hour at the end of prev time step
!! @param [in]     slag(1)     - real, equation of time ( radian )
!! @param [in]     sdec(1),cdec(1) - real, sin and cos of the solar declination angle
!! @param [in]     sinlat(im)   - real, sin of latitude
!! @param [in]     coslat(im)   - real, cos of latitude
!! @param [in]     pgr(im)      - real, surface pressure (Pa)
!! @param [in]     ugrs(ix,levs),vgrs(ix,levs) - real, u/v component of layer wind
!! @param [in]     tgrs(ix,levs)     - real, layer mean temperature ( k )
!! @param [in]     qgrs(ix,levs,ntrac)     - real, layer mean tracer concentration
!! @param [in]     vvel(ix,levs)     - real, layer mean vertical velocity (Pa/s)
!! @param [in]     prsi(ix,levs+1)     - real, pressure at layer interfaces
!! @param [in]     prsl(ix,levs)     - real, mean layer presure
!! @param [in]     prsik(ix,levs+1)    - real, Exner function at layer interface
!! @param [in]     prslk(ix,levs)    - real, Exner function at layer
!! @param [in]     phii(ix,levs+1)     - real, interface geopotential (m^2/s^2)
!! @param [in]     phil(ix,levs)     - real, layer geopotential     (m^2/s^2)
!! @param [in]     rann(ix,nrcm)     - real, random number array (0-1)
!! @param [in]     prdout(ix,ko3,pl_coeff)   - real, ozone forcing data
!! @param [in]     poz(ko3)      - real, ozone forcing data level pressure (ln(Pa))
!! @param [in]     dpshc(im)    - real, maximum pressure depth for shallow convection
!! @param [in]     hprime(ix,nmtvr)   - real, orographic std dev
!! @param [in]     xlon(im),xlat(im) - real, longitude and latitude ( radian )
!! @param [in]     slope(im)    - real, sfc slope type for lsm
!! @param [in]     shdmin(im1)   - real, min fractional coverage of green veg
!! @param [in]     shdmax(im1)   - real, max fractnl cover of green veg (not used)
!! @param [in]     snoalb(im1)   - real, max snow albedo over land (for deep snow)
!! @param [in]     tg3(im)      - real, deep soil temperature
!! @param [in]     slmsk(im    - real, sea/land/ice mask (=0/1/2)
!! @param [in]     vfrac(im)    - real, vegetation fraction
!! @param [in]     vtype(im)    - real, vegetation type
!! @param [in]     stype(im)    - real, soil type
!! @param [in]     uustar(im)   - real, boundary layer parameter
!! @param [in]     oro(im)      - real, orography
!! @param [in]     oro_uf(im)   - real, unfiltered orography
!! @param [in]     coszen(im)   - real, avg cosz over daytime sw radiation interval
!! @param [in]     sfcdsw(im)   - real, total sky sfc downward sw flux ( w/m**2 )
!! @param [in]     sfcnsw(im)   - real, total sky sfc netsw flx into ground(w/m**2)
!! @param [in]     sfcdlw(im)   - real, total sky sfc downward lw flux ( w/m**2 )
!! @param [in]     tsflw(im)    - real, sfc air (layer 1) temp over lw interval (k)
!! @param [in]     sfcemis(im)  - real, sfc lw emissivity ( fraction )
!! @param [in]     sfalb(im)    - real, mean sfc diffused sw albedo
!! @param [in]     sfcnirbmu(im) - real, sfc nir-beam sw upward flux (w/m2)
!! @param [in]     sfcnirdfu(im) - real, sfc nir-diff sw upward flux (w/m2)
!! @param [in]     sfcvisbmu(im) - real, sfc uv+vis-beam sw upward flux (w/m2)
!! @param [in]     sfcvisdfu(im) - real, sfc uv+vis-diff sw upward flux (w/m2)
!! @param [in]     sfcnirbmd(im) - real, sfc nir-beam sw downward flux (w/m2)
!! @param [in]     sfcnirdfd(im) - real, sfc nir-diff sw downward flux (w/m2)
!! @param [in]     sfcvisbmd(im) - real, sfc uv+vis-beam sw downward flux (w/m2)
!! @param [in]     sfcvisdfd(im) - real, sfc uv+vis-diff sw downward flux (w/m2)
!! @param [in]     slimskin_cpl(im) - real,
!! @param [in]     ulwsfcin_cpl(im) - real,
!! @param [in]     dusfcin_cpl(im)  - real,
!! @param [in]     dvsfcin_cpl(im)  - real,
!! @param [in]     dtsfcin_cpl(im)  - real,
!! @param [in]     dqsfcin_cpl(im)  - real,
!! @param [in]     swh(ix,levs)      - real, total sky sw heating rates ( k/s )
!! @param [in]     swhc(ix,levs)     - real, clear sky sw heating rates ( k/s )
!! @param [in]     hlw(ix,levs)      - real, total sky lw heating rates ( k/s )
!! @param [in]     hlwc(ix,levs)     - real, clear sky lw heating rates ( k/s )
!! @param [in]     hlwd(ix,levs)     - real, idea  sky lw heating rates ( k/s )
!! @param [in]     ras(1)      - logical, flag for ras convection scheme
!! @param [in]     cscnv(1)    - logical, flag for Chikira-Sugiyama convection
!! @param [in]     nctp(1)     - integer, number of cloud types in CS scheme
!! @param [in]     do_shoc(1)  - logical, flag for SHOC
!! @param [in]     shocaftcnv(1) - logical, flag for SHOC
!! @param [in]     ntot3d(1)   - integer, number of total 3d fields for phy_f3d
!! @param [in]     ntot2d(1)   - integer, number of total 2d fields for phy_f2d
!! @param [in]     pre_rad(1)  - logical, flag for testing purpose
!! @param [in]     ldiag3d(1)  - logical, flag for 3d diagnostic fields
!! @param [in]     lgocart(1)  - logical, flag for 3d diagnostic fields for gocart
!! @param [in]     lssav(1)    - logical, flag controls data store and output
!! @param [in]     lssav_cpl(1) - logical, flag for save data for A/O/I coupling
!! @param [in]     flipv(1)    - logical, flag for vertical direction flip (ras)
!! @param [in]     xkzm_m(1)   - real, background vertical diffusion for momentum
!! @param [in]     xkzm_h(1)   - real, background vertical diffusion for heat, q
!! @param [in]     xkzm_s(1)   - real, sigma threshold for background mom. diffusn
!! @param [in]     psautco(2)  - real, auto conversion coeff from ice to snow
!! @param [in]     prautco(2)  - real, auto conversion coeff from cloud to rain
!! @param [in]     evpco(1)    - real, coeff for evaporation of largescale rain
!! @param [in]     wminco(1)   - real, water and ice minimum threshold for Zhao
!! @param [in]     pdfcld(1)   - logical, flag for pdfcld
!! @param [in]     shcnvcw(1)  - logical, flag for shallow convective cloud
!! @param [in]     sup(1)      - real, supsaturation for ho. nucleation of ice
!! @param [in]     redrag(1)   - logical, flag for reduced drag coeff. over sea
!! @param [in]     hybedmf(1)  - logical, flag for hybrid edmf pbl scheme
!! @param [in]     dspheat(1)  - logical, flag for tke dissipative heating
!! @param [in]     old_monin(1) - logical, flag for diff monin schemes
!! @param [in]     cnvgwd(1)   - logical, flag for conv gravity wave drag
!! @param [in]     shal_cnv(1) - logical, flag for calling shallow convection
!! @param [in]     imfshalcnv(1) - integer, flag for mass-flux shallow conv scheme
!!<PRE>
!!     1: July 2010 version of mass-flux shallow conv scheme
!!         current operational version as of 2016
!!     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
!!     0: modified Tiedtke's eddy-diffusion shallow conv scheme
!!    -1: no shallow convection used
!!</PRE>
!! @param [in]     imfdeepcnv(1) - integer, flag for mass-flux deep conv scheme
!!<PRE>
!!     1: July 2010 version of SAS conv scheme
!!           current operational version as of 2016
!!     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
!!     0: old SAS Convection scheme before July 2010
!!</PRE>
!! @param [in]     cal_pre(1)  - logical, flag controls precip type algorithm
!! @param [in]     mom4ice(1)  - logical, flag controls mom4 sea-ice
!! @param [in]     mstrat(1)   - logical, flag for moorthi approach for stratus
!! @param [in]     trans_trac(1) - logical, flag for convective transport of tracers
!! @param [in]     nstf_name(1)   - integer array, NSST related flag parameters
!!<PRE>
!!     nstf_name(1) : 0 = NSSTM off
!!                    1 = NSSTM on but uncoupled
!!                    2 = NSSTM on and coupled
!!     nstf_name(2) : 1 = NSSTM spin up on
!!                    0 = NSSTM spin up off
!!     nstf_name(3) : 1 = NSST analysis on
!!                    0 = NSSTM analysis off
!!     nstf_name(4) : zsea1 in mm
!!     nstf_name(5) : zsea2 in mm
!!</PRE>
!! @param [in]     moist_adj(1) - logical, flag for moist convective adjustment
!! @param [in]     thermodyn_id(1) - integer, valid for GFS only for get_prs/phi
!! @param [in]     sfcpress_id(1)  - integer, valid for GFS only for get_prs/phi
!! @param [in]     gen_coord_hybrid(1) - logical for pk=ak+bk*ps+ck*theta (Henry)
!! @param [in]     levr(1)  - integer, the number of layers GFS Radiative heating calculted at 1
!! @param [in]     adjtrc(ntrac)   - real, dynamics adjustments to tracers
!! @param [in]     nnp(1)      - integer, physics substep number
!! @param [in,out]     hice(im)      - real, sea-ice thickness
!! @param [in,out]     fice(im)      - real, sea-ice concentration
!! @param [in,out]     tisfc(im)     - real, sea-ice temperature
!! @param [in,out]     tsea(im)      - real, ground surface temperature ( k )
!! @param [in,out]     tprcp(im)     - real, total precipitation
!!                     the following three variables do not affect the forecast
!! @param [in,out]     cv(im)         - real, convective clouds amountt
!! @param [in,out]     cvb(im)        - real, convective clouds base pressure (kPa)
!! @param [in,out]     cvt(im)        - real, convective clouds top  pressure (kPa)
!! @param [in,out]     srflag(im)    - real, snow/rain flag for precipitation
!! @param [in,out]     snwdph(im)    - real, actual snow depth (mm) over land/sea ice
!! @param [in,out]     weasd(im)     - real, water equiv of accumulated  snow depth (kg/m**2)
!!                                      over land and sea ice
!! @param [in,out]     sncovr(im)    - real, snow cover over land
!! @param [in,out]     zorl(im)      - real, surface roughness
!! @param [in,out]     canopy(im)    - real, canopy water
!! @param [in,out]     ffmm(im)      - real, fm parameter from PBL scheme
!! @param [in,out]     ffhh(im)      - real, fh parameter from PBL scheme
!! @param [in,out]     f10m(im)      - real, fm at 10m
!! @param [in,out]     srunoff(im)   - real, surface water runoff (from lsm)
!! @param [in,out]     evbsa(im)     - real, noah lsm diagnostics
!! @param [in,out]     evcwa(im)     - real, noah lsm diagnostics
!! @param [in,out]     snohfa(im)    - real, noah lsm diagnostics
!! @param [in,out]     transa(im)    - real, noah lsm diagnostics
!! @param [in,out]     sbsnoa(im)    - real, noah lsm diagnostics
!! @param [in,out]     snowca(im)    - real, noah lsm diagnostics
!! @param [in,out]     soilm(im)     - real, soil moisture
!! @param [in,out]     tmpmin(im)    - real, min temperature at 2m height (k)
!! @param [in,out]     tmpmax(im)    - real, max temperature at 2m height (k)
!! @param [in,out]     dusfc(im)     - real, u component of surface stress
!! @param [in,out]     dvsfc(im)     - real, v component of surface stress
!! @param [in,out]     dtsfc(im)     - real, sensible heat flux (w/m2)
!! @param [in,out]     dqsfc(im)     - real, latent heat flux (w/m2)
!! @param [in,out]     totprcp(im)   - real, accumulated total precipitation (kg/m2)
!! @param [in,out]     gflux(im)     - real, groud conductive heat flux
!! @param [in,out]     dlwsfc(im)    - real, time accumulated sfc dn lw flux ( w/m**2 )
!! @param [in,out]     ulwsfc(im)    - real, time accumulated sfc up lw flux ( w/m**2 )
!! @param [in,out]     suntim(im)    - real, sunshine duration time (s)
!! @param [in,out]     runoff(im)    - real, total water runoff
!! @param [in,out]     ep(im)        - real, potential evaporation
!! @param [in,out]     cldwrk(im)    - real, cloud workfunction (valid only with sas)
!! @param [in,out]     dugwd(im)     - real, vertically integrated u change by OGWD
!! @param [in,out]     dvgwd(im)     - real, vertically integrated v change by OGWD
!! @param [in,out]     psmean(im)    - real, surface pressure (kPa)
!! @param [in,out]     cnvprcp(im)   - real, accumulated convective precipitation (kg/m2)
!! @param [in,out]     spfhmin(im)   - real, minimum specific humidity
!! @param [in,out]     spfhmax(im)   - real, maximum specific humidity
!! @param [in,out]     rain(im)      - real, total rain at this time step
!! @param [in,out]     rainc(im)     - real, convective rain at this time step
!! @param [in,out]     dt3dt(ix,les,6)     - real, temperature change due to physics
!! @param [in,out]     dq3dt(ix,levs,5+pl_coeff)     - real, moisture change due to physics
!! @param [in,out]     du3dt(ix,levs,4)     - real, u momentum change due to physics
!! @param [in,out]     dv3dt(ix,levs,4)     - real, v momentum change due to physics
!! @param [in,out]     dqdt_v(ix,levs)    - real, total moisture tendency (kg/kg/s)
!! @param [in,out]     cnvqc_v(ix,levs)   - real, total convective conensate (kg/kg)
!! @param [in,out]     acv(im)       - real,  array containing accumulated convective clouds
!! @param [in,out]     acvb(im),acvt(im) - real,  arrays used by cnvc90
!! @param [in,out]     slc(ix,soil)       - real, liquid soil moisture
!! @param [in,out]     smc(ix,soil)       - real, total soil moisture
!! @param [in,out]     stc(ix,soil)       - real, soil temperature
!! @param [in,out]     upd_mf(ix,levs)    - real, convective updraft mass flux
!! @param [in,out]     dwn_mf(ix,levs)    - real, convective downdraft mass flux
!! @param [in,out]     det_mf(ix,levs)    - real, convective detrainment mass flux
!                   ------- not used below -----------
!> @param [in,out]     dkh(ix,levs)       - real, vertical diffusion coefficient (gocart), not used
!! @param [in,out]     rnp(ix,levs)       - real, n cloud precipitation rate     (gocart), not used
!                   ------- not used  above -----------
!> @param [in,out]     phy_f3d(ix,levs,num_p3d)   - real, 3d arrays saved for restart
!! @param [in,out]     phy_f2d(ix,num_p2d)   - real, 2d arrays save for restart
!! @param [in,out]     dusfc_cpl(im) - real, sfc u-momentum flux       for A/O/I coupling
!! @param [in,out]     dvsfc_cpl(im) - real, sfc v-momentum flux       for A/O/I coupling
!! @param [in,out]     dtsfc_cpl(im) - real, sfc sensible heat flux    for A/O/I coupling
!! @param [in,out]     dqsfc_cpl(im) - real, sfc latent heat flux      for A/O/I coupling
!! @param [in,out]     dlwsfc_cpl(im) - real, sfc dnwd lw flux (w/m**2) for A/O/I coupling
!! @param [in,out]     dswsfc_cpl(im) - real, sfc dnwd sw flux (w/m**2) for A/O/I coupling
!! @param [in,out]     dnirbm_cpl(im) - real, sfc nir beam dnwd sw rad flux (w/m**2)
!! @param [in,out]     dnirdf_cpl(im) - real, sfc nir diff dnwd sw rad flux (w/m**2)
!! @param [in,out]     dvisbm_cpl(im) - real, sfc uv+vis beam dnwd sw rad flux (w/m**2)
!! @param [in,out]     dvisdf_cpl(im) - real, sfc uv+vis diff dnwd sw rad flux (w/m**2)
!! @param [in,out]     nlwsfc_cpl(im) - real, net dnwd lw flux (w/m**2) for A/O/I coupling
!! @param [in,out]     nswsfc_cpl(im) - real, net dnwd sw flux (w/m**2) for A/O/I coupling
!! @param [in,out]     nnirbm_cpl(im) - real, net nir beam dnwd sw rad flux (w/m**2)
!! @param [in,out]     nnirdf_cpl(im) - real, net nir diff dnwd sw rad flux (w/m**2)
!! @param [in,out]     nvisbm_cpl(im) - real, net uv+vis beam dnwd sw rad flux (w/m**2)
!! @param [in,out]     nvisdf_cpl(im) - real, net uv+vis diff dnwd sw rad flux (w/m**2)
!! @param [in,out]     rain_cpl(im)  - real, total precipitation  rain for A/O/I coupling
!! @param [in,out]     snow_cpl(im)  - real, total precipitation  snow for A/O/I coupling
!! @param [in,out]     xt(im)        - real, heat content in DTL
!! @param [in,out]     xs(im)        - real, salinity  content in DTL
!! @param [in,out]     xu(im)        - real, u-current content in DTL
!! @param [in,out]     xv(im)        - real, v-current content in DTL
!! @param [in,out]     xz(im)        - real, DTL thickness
!! @param [in,out]     zm(im)        - real, MXL thickness
!! @param [in,out]     xtts(im)      - real, d(xt)/d(ts)
!! @param [in,out]     xzts(im)      - real, d(xz)/d(ts)
!! @param [in,out]     d_conv(im)    - real, thickness of Free Convection Layer (FCL)
!! @param [in,out]     ifd(im)       - real, index to start DTM run or not
!! @param [in,out]     dt_cool(im)   - real, Sub-layer cooling amount
!! @param [in,out]     Qrain(im)     - real, sensible heat flux due to rainfall (watts)
!! @param [in,out]     phy_fctd(ix,nctp) - real, cloud base mass flux for CScnv
!! @param [out]     gt0(ix,levs)       - real, updated temperature
!! @param [out]     gq0(ix,levs,ntrac)       - real, updated tracers
!! @param [out]     gu0(ix,levs)       - real, updated zonal wind
!! @param [out]     gv0(ix,levs)       - real, update meridional wind
!! @param [out]     t2m(im),q2m(im)   - real, 2 meter temperature and humidity
!! @param [out]     u10m(im),v10m(im) - real, 10 meater u/v wind speed
!! @param [out]     zlvl(im)      - real, layer 1 height (m)
!! @param [out]     psurf(im)     - real, surface pressure (Pa)
!! @param [out]     hpbl(im)      - real, pbl height (m)
!! @param [out]     pwat(im)      - real, precipitable water
!! @param [out]     t1(im)        - real, layer 1 temperature (K)
!! @param [out]     q1(im)        - real, layer 1 specific humidity (kg/kg)
!! @param [out]     u1(im)        - real, layer 1 zonal wind (m/s)
!! @param [out]     v1(im)        - real, layer 1 merdional wind (m/s)
!! @param [out]     chh(im)       - real, thermal exchange coefficient
!! @param [out]     cmm(im)       - real, momentum exchange coefficient
!! @param [out]     dlwsfci(im)   - real, instantaneous sfc dnwd lw flux ( w/m**2 )
!! @param [out]     ulwsfci(im)   - real, instantaneous sfc upwd lw flux ( w/m**2 )
!! @param [out]     dswsfci(im)   - real, instantaneous sfc dnwd sw flux ( w/m**2 )
!! @param [out]     uswsfci(im)   - real, instantaneous sfc upwd sw flux ( w/m**2 )
!! @param [out]     dusfci(im)    - real, instantaneous u component of surface stress
!! @param [out]     dvsfci(im)    - real, instantaneous v component of surface stress
!! @param [out]     dtsfci(im)    - real, instantaneous sfc sensible heat flux
!! @param [out]     dqsfci(im)    - real, instantaneous sfc latent heat flux
!! @param [out]     gfluxi(im)    - real, instantaneous sfc ground heat flux
!! @param [out]     epi(im)       - real, instantaneous sfc potential evaporation
!! @param [out]     smcwlt2(im)   - real, wilting point (volumetric)
!! @param [out]     smcref2(im)   - real, soil moisture threshold (volumetric)
!! @param [out]     dusfci_cpl(im)  - real, sfc u-momentum flux at time step AOI cpl
!! @param [out]     dvsfci_cpl(im)  - real, sfc v-momentum flux at time step AOI cpl
!! @param [out]     dtsfci_cpl(im)  - real, sfc sensib heat flux at time step AOI cpl
!! @param [out]     dqsfci_cpl(im)  - real, sfc latent heat flux at time step AOI cpl
!! @param [out]     dlwsfci_cpl(im) - real, sfc dnwd lw flux at time step AOI cpl
!! @param [out]     dswsfci_cpl(im) - real, sfc dnwd sw flux at time step AOI cpl
!! @param [out]     dnirbmi_cpl(im) - real, sfc nir beam dnwd sw flx rad at time step
!! @param [out]     dnirdfi_cpl(im) - real, sfc nir diff dnwd sw flx rad at time step
!! @param [out]     dvisbmi_cpl(im) - real, sfc uv+vis beam dnwd sw flx at time step
!! @param [out]     dvisdfi_cpl(im) - real, sfc uv+vis diff dnwd sw flx at time step
!! @param [out]     nlwsfci_cpl(im) - real, net sfc dnwd lw flux at time step AOI cpl
!! @param [out]     nswsfci_cpl(im) - real, net sfc dnwd sw flux at time step AOI cpl
!! @param [out]     nnirbmi_cpl(im) - real, net nir beam dnwd sw flx rad at time step
!! @param [out]     nnirdfi_cpl(im) - real, net nir diff dnwd sw flx rad at time step
!! @param [out]     nvisbmi_cpl(im) - real, net uv+vis beam dnwd sw flx at time step
!! @param [out]     nvisdfi_cpl(im) - real, net uv+vis diff dnwd sw flx at time step
!! @param [out]     ocalnirbm_cpl(im) - real, ocean alb nir beam (no ice) at time step
!! @param [out]     ocalnirdf_cpl(im) - real, ocean alb nir diff (no ice) at time step
!! @param [out]     ocalvisbm_cpl(im) - real, ocean alb vis beam (no ice) at time step
!! @param [out]     ocalvisdf_cpl(im) - real, ocean alb vis diff (no ice) at time step
!! @param [out]     t2mi_cpl(im)    - real, T2m at time step AOI cpl
!! @param [out]     q2mi_cpl(im)    - real, Q2m at time step AOI cpl
!! @param [out]     u10mi_cpl(im)   - real, U10m at time step AOI cpl
!! @param [out]     v10mi_cpl(im)   - real, V10m at time step AOI cpl
!! @param [out]     tseai_cpl(im)   - real, sfc temp at time step AOI cpl
!! @param [out]     psurfi_cpl(im)  - real, sfc pressure at time step AOI cpl
!! @param [out]     tref(im)        - real, Reference Temperature
!! @param [out]     z_c(im)         - real, Sub-layer cooling thickness
!! @param [out]     c_0(im)         - real, coefficient1 to calculate d(Tz)/d(Ts)
!! @param [out]     c_d(im)         - real, coefficient2 to calculate d(Tz)/d(Ts)
!! @param [out]     w_0(im)         - real, coefficient3 to calculate d(Tz)/d(Ts)
!! @param [out]     w_d(im)         - real, coefficient4 to calculate d(Tz)/d(Ts)
!! @param [out]     rqtk(im)        - real, mass change due to moisture variation
!! @param [out]     dtdtr(ix,levs)       - real, temperature change due to radiative heating
!!                                      per time step (K)
!!  \section general General Algorithm
!!  -# Prepare input variables for calling individual parameterizations.
!!  -# Using a two-iteration loop, calculate the state variable tendencies for the surface layer.
!!  -# Calculate the state variable tendencies due to the PBL (vertical diffusion) scheme.
!!  -# Calculate the state variable tendencies due to orographic gravity wave drag and Rayleigh damping.
!!  -# Apply tendencies to the state variables calculated so far:
!!   - for temperature: radiation, surface, PBL, oro. GWD, Rayleigh damping
!!   - for momentum: surface, PBL, oro. GWD, Rayleigh damping
!!   - for water vapor: surface, PBL
!!  -# If SHOC is active and is supposed to be called before convection, call it and update the state variables within.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to deep convection.
!!  -# Calculate the state variable tendencies due to convective gravity wave drag and apply them afterwards.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to shallow convection.
!!  -# If SHOC is active and is supposed to be called after convection, call it and update the state variables within.
!!  -# If necessary, call the moist convective adjustment subroutine and update the state temperature and moisture variable within.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to microphysics.
!!  -# Determine the precipitation type and update land surface properties if necessary.
!!  -# Fill the output variables from the local variables as necessary and deallocate allocatable arrays.
!!  \section detailed Detailed Algorithm
!!  ## Prepare input variables for calling individual parameterizations.
!!  Before calling any parameterizations, there is a section at the beginning of the subroutine for
!!  preparing input arguments to the various schemes based on general input to the driver and initializing
!!  variables used throughout the driver.
!!  - General initialization:
!!   - set a flag for running in debug mode and the horizontal index of the column to print
!!   - calculate the pressure at layer centers, the exner function at layer centers and interfaces,
!!     geopotential at layer centers and interfaces, and the layer-centered pressure difference
!!   - calculate the ratio of dynamics time step to physics time step for applying tendencies
!!   - initialize local tendency arrays to zero
!!  - Radiation:
!!   - adjust radiative fluxes and heating rates to the shorter physics time step (from the longer radiation time step),
!!    unless idealized physics is true (lsidea) where radiative heating rates are set to 0
!!   - compute diagnostics from the radiation scheme needed for other schemes (e.g., downward longwave flux absorbed by the surface)
!!   - accumulate the upward and downward longwave fluxes at the surface
!!  - Surface:
!!   - set NOAH and OSU scheme variables from gbphys input arguments and initialize local soil moisture variables
!!   - set local sea ice variables from gbphys arguments
!!   - set up A/O/I coupling variables from gbphys arguments
!!  - PBL:
!!   - set the number of tracers that are diffused vertically
!!  - SHOC:
!!   - determine the index of TKE (ntk) in the convectively transported tracer array (clw)
!!   - allocate precipitation mixing ratio cloud droplet number concentration arrays
!!  - Deep Convection:
!!   - determine which tracers in the tracer input array undergo convective transport (valid only for the RAS and Chikira-Sugiyama schemes) and allocate a local convective transported tracer array (clw)
!!   - apply an adjustment to the tracers from the dynamics
!!   - calculate horizontal grid-related parameters needed for some parameterizations
!!   - calculate the maxiumum cloud base updraft speed for the Chikira-Sugiyama scheme
!!   - allocate array for cloud water and cloud cover (for non-RAS and non-Chikira-Sugiyama deep convective schemes)
!!  - Shallow Convection:
!!   - when using the Tiedtke shallow convection scheme with the stratus modifications, find the lowest
!!     model level where a temperature inversion exists in the absence of CTEI
!!  - Microphysics:
!!   - for the Morrison (MGB) scheme, calculate 'FRLAND' if the grid point is over land
!!   - allocate arrays associated with the Morrison scheme
!!   - assign the local critical relative humidity variables from the gbphys arguments
!!  - Gravity Wave Drag:
!!   - calculate the deep convective cloud fraction at cloud top for the convective GWD scheme
!>  ## Using a two-iteration loop, calculate the state variable tendencies for the surface layer.
!!   - each iteration of the loop calls the following:
!!    - 'sfc_diff' to calculate surface exchange coefficients and near-surface wind
!!    - surface energy balances routines are called regardless of surface type; the surface type is checked within each to determine whether the routine is "active"
!!    - for the surface energy balance over the ocean, call 'sfc_nst' if NSST is on, otherwise, call 'sfc_ocean'
!!    - for the surface energy balance over the land, call 'sfc_drv' for the NOAH model and 'sfc_land' for the OSU model
!!    - for the surface energy balance over sea ice, call sfc_sice; if A/O/I coupling, call sfc_cice
!!   - the initial iteration has flag_guess = F unless wind < 2 m/s; flag_iter = T
!!   - after the initial iteration, flag_guess = F and flag_iter = F (unless wind < 2 m/s and over a land surface or an ocean surface with NSST on)
!!   - the following actions are performed after the iteration to calculate surface energy balance:
!!    - set surface output variables from their local values
!!    - call 'sfc_diag' to calculate state variable values at 2 and 10 m as appropriate from near-surface model levels and the surface exchange coefficients
!!    - if A/O/I coupling, set coupling variables from local variables and calculate the open water albedo
!!    - finally, accumulate surface-related diagnostics and calculate the max/min values of T and q at 2 m height.
!>  ## Calculate the state variable tendencies due to the PBL (vertical diffusion) scheme.
      subroutine gbphys                                                 &
!  ---  inputs:
     &    ( im,ix,levs,lsoil,lsm,ntrac,ncld,ntoz,ntcw,ntke,             &
     &      ntiw,ntlnc,ntinc,                                           &
     &      nmtvr,nrcm,ko3,lonr,latr,jcap,                              &
     &      num_p3d,num_p2d,npdf3d,ncnvcld3d,                           &
     &      kdt,lat,me,pl_coeff,nlons,ncw,flgmin,crtrh,cdmbgwd,         &
     &      ccwf,dlqf,ctei_rm,clstp,cgwf,prslrd0,ral_ts,dtp,dtf,fhour,  &
     &      solhr,slag,sdec,cdec,sinlat,coslat,pgr,ugrs,vgrs,           &
     &      tgrs,qgrs,vvel,prsi,prsl,prslk,prsik,phii,phil,             &
     &      rann,prdoz,poz,dpshc,fscav,fswtr,hprime,xlon,xlat,          &
     &      h2o_phys,levh2o,h2opl,h2o_pres,h2o_coeff,                   &
     &      isot, ivegsrc,                                              &
     &      slope,shdmin,shdmax,snoalb,tg3,slmsk,vfrac,                 &
     &      vtype,stype,uustar,oro,oro_uf,coszen,sfcdsw,sfcnsw,         &

     &      sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                    &
     &      sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                    &
     &      slimskin_cpl,ulwsfcin_cpl,                                  &
     &      dusfcin_cpl,dvsfcin_cpl,dtsfcin_cpl,dqsfcin_cpl,            &
     &      sfcdlw,tsflw,sfcemis,sfalb,swh,swhc,hlw,hlwc,hlwd,lsidea,   &
     &      ras,pre_rad,ldiag3d,lgocart,lssav,lssav_cpl,                &

     &      xkzm_m,xkzm_h,xkzm_s,psautco,prautco,evpco,wminco,          &
     &      pdfcld,shcnvcw,sup,redrag,hybedmf,dspheat,                  &
     &      flipv,old_monin,cnvgwd,shal_cnv,                            &
     &      imfshalcnv,imfdeepcnv,cal_pre,aero_in,                      &
     &      mom4ice,mstrat,trans_trac,nstf_name,moist_adj,              &
     &      thermodyn_id, sfcpress_id, gen_coord_hybrid,levr,adjtrc,nnp,&
     &      cscnv,nctp,do_shoc,shocaftcnv,ntot3d,ntot2d,                &
!  ---  input/outputs:
     &      hice,fice,tisfc,tsea,tprcp,cv,cvb,cvt,                      &
     &      srflag,snwdph,weasd,sncovr,zorl,canopy,                     &
     &      ffmm,ffhh,f10m,srunoff,evbsa,evcwa,snohfa,                  &
     &      transa,sbsnoa,snowca,soilm,tmpmin,tmpmax,                   &
     &      dusfc,dvsfc,dtsfc,dqsfc,totprcp,gflux,                      &
     &      dlwsfc,ulwsfc,suntim,runoff,ep,cldwrk,                      &
     &      dugwd,dvgwd,psmean,cnvprcp,spfhmin,spfhmax,rain,rainc,      &

     &      dt3dt,dq3dt,du3dt,dv3dt,dqdt_v,cnvqc_v,acv,acvb,acvt,       &
     &      slc,smc,stc,upd_mf,dwn_mf,det_mf,phy_f3d,phy_f2d,           &
!    &      slc,smc,stc,upd_mf,dwn_mf,det_mf,dkh,rnp,phy_f3d,phy_f2d,   &
     &      dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl,                 &
     &      dlwsfc_cpl,dswsfc_cpl,dnirbm_cpl,dnirdf_cpl,                &
     &      dvisbm_cpl,dvisdf_cpl,rain_cpl,  nlwsfc_cpl,nswsfc_cpl,     &
     &      nnirbm_cpl,nnirdf_cpl,nvisbm_cpl,nvisdf_cpl,snow_cpl,       &

     &      xt,xs,xu,xv,xz,zm,xtts,xzts,d_conv,ifd,dt_cool,Qrain,       &
     &      tref, z_c, c_0, c_d, w_0, w_d,                              &
     &      phy_fctd,                                                   &
!  ---  outputs:
     &      gt0,gq0,gu0,gv0,t2m,q2m,u10m,v10m,                          &
     &      zlvl,psurf,hpbl,pwat,t1,q1,u1,v1,                           &
     &      chh,cmm,dlwsfci,ulwsfci,dswsfci,uswsfci,dusfci,dvsfci,      &
     &      dtsfci,dqsfci,gfluxi,epi,smcwlt2,smcref2,wet1,sr,           &

     &      rqtk,                                                       &
!       Stochastic physics perturnbation
     &      dtdtr,                                                      &

     &      dusfci_cpl,  dvsfci_cpl, dtsfci_cpl, dqsfci_cpl,            &
     &      dlwsfci_cpl, dswsfci_cpl,                                   &
     &      dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl,         &
     &      nlwsfci_cpl, nswsfci_cpl,                                   &
     &      nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl,         &
     &      t2mi_cpl,    q2mi_cpl,    u10mi_cpl,   v10mi_cpl,           &
     &      tseai_cpl,   psurfi_cpl                                     &
     &      )
!
      use machine ,   only : kind_phys
      use physcons,   only : con_cp, con_fvirt, con_g, con_rd, con_rv,  &
     &                       con_hvap, con_hfus, con_rerth, con_pi
     &,                      rhc_max, dxmin, dxinv, pa2mb, rlapse
      use module_nst_water_prop, only: get_dtzm_2d
      use cs_conv, only : cs_convr

      implicit none
!
!  ---  some constant parameters:

      real(kind=kind_phys), parameter :: hocp    = con_hvap/con_cp
!    &,                                  fhourpr = 0.0
     &,                                  qmin    = 1.0e-10
     &,                                  p850    = 85000.0
     &,                                  epsq    = 1.e-20
     &,                                  hsub    = con_hvap+con_hfus
     &,                                  czmin   = 0.0001      ! cos(89.994)

!  ---  inputs:

!  note: lgocart is the logical flag for in-line gocart;

      integer, intent(in) :: ix,   im,   levs, lsoil,   lsm,     ntrac, &
     &                       ncld, ntiw, ntlnc, ntinc,                  &
     &                       ntoz, ntcw, nmtvr,   nrcm,    ko3,         &
     &                       lonr, latr, jcap, num_p3d, num_p2d, kdt,   &
     &                       me,   pl_coeff, lat, npdf3d, ncnvcld3d,    &
     &                       thermodyn_id, sfcpress_id, levr, nnp, nctp,&
     &                       ntke, ntot3d, ntot2d,  h2o_coeff, levh2o,  &
     &                       isot, ivegsrc


      integer, intent(in) :: nlons(im), ncw(2)
      integer, intent(in) :: nstf_name(5)
      integer, intent(in) :: imfshalcnv, imfdeepcnv

      logical, intent(in) :: ras,        pre_rad,   ldiag3d, flipv,     &
     &                       old_monin,  cnvgwd,    aero_in,            &
     &                       redrag,     hybedmf,   dspheat,            &
     &                       lssav,                 mom4ice, mstrat,    &
     &                       trans_trac, moist_adj, cal_pre, cscnv,     &
     &                       shal_cnv,   gen_coord_hybrid,   lgocart,   &
     &                       lsidea,     lssav_cpl, pdfcld, shcnvcw,    &
     &                       do_shoc,    shocaftcnv, h2o_phys

      real(kind=kind_phys) :: adjtrc(ntrac)

      real(kind=kind_phys), dimension(im),            intent(in) ::     &
     &      sinlat, coslat, pgr,    dpshc,  xlon,   xlat,               &
     &      slope,  shdmin, shdmax, snoalb, tg3,    slmsk,  vfrac,      &
     &      vtype,  stype,  uustar, oro,    coszen, sfcnsw, sfcdsw,     &
     &      sfcnirbmu,      sfcnirdfu,      sfcvisbmu,      sfcvisdfu,  &
     &      sfcnirbmd,      sfcnirdfd,      sfcvisbmd,      sfcvisdfd,  &
     &      slimskin_cpl,   dusfcin_cpl,    dvsfcin_cpl,                &
     &      dtsfcin_cpl,    dqsfcin_cpl,    ulwsfcin_cpl,               &
     &      sfcdlw, tsflw,  sfalb,  sfcemis, oro_uf

      real(kind=kind_phys), dimension(ix,levs),       intent(in) ::     &
     &   ugrs, vgrs, tgrs, vvel, prsl, prslk, phil, swh, swhc, hlw, hlwc

!idea add by hmhj
      real(kind=kind_phys), intent(in) ::  hlwd(ix,levs,6)

      real(kind=kind_phys), intent(inout) ::  qgrs(ix,levs,ntrac)

      real(kind=kind_phys), dimension(ix,levs+1),     intent(in) ::     &
     &      prsi, prsik, phii

      real(kind=kind_phys), intent(in) ::  hprime(ix,nmtvr),            &
     &      prdoz(ix,ko3,pl_coeff),        rann(ix,nrcm), poz(ko3),     &
     &      h2opl(ix,levh2o,h2o_coeff),    h2o_pres(levh2o)

      real(kind=kind_phys), intent(in) ::  dtp,   dtf, fhour, solhr,    &
     &      slag,    sdec,      cdec,      ctei_rm(2), clstp,           &
     &      ccwf(2), crtrh(3),  flgmin(2), dlqf(2),    cdmbgwd(2),      &
     &      xkzm_m,  xkzm_h,    xkzm_s,    psautco(2), prautco(2),      &
     &      evpco,   wminco(2), cgwf(2),   prslrd0,    sup, ral_ts

!  ---  input/output:
      real(kind=kind_phys), dimension(im),            intent(inout) ::  &
     &      hice,   fice,    tisfc,  tsea,   tprcp,  cv,     cvb,  cvt, &
     &      srflag, snwdph,  weasd,  sncovr, zorl,   canopy, ffmm, ffhh,&
     &      f10m,   srunoff, evbsa,  evcwa,  snohfa, transa, sbsnoa,    &
     &      snowca, soilm,   tmpmin, tmpmax, dusfc,  dvsfc,  dtsfc,     &
     &      dqsfc,  totprcp, gflux,  dlwsfc, ulwsfc, suntim, runoff, ep,&
     &      cldwrk, dugwd,   dvgwd,  psmean, cnvprcp,spfhmin,spfhmax,   &
     &      rain,   rainc,   acv,    acvb,   acvt
      real(kind=kind_phys), dimension(im), optional,  intent(inout) ::  &
! for A/O/I coupling
     &      dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl,                 &
     &      dlwsfc_cpl,dswsfc_cpl,rain_cpl,  snow_cpl,                  &
     &      dnirbm_cpl,dnirdf_cpl,dvisbm_cpl,dvisdf_cpl,                &
     &      nlwsfc_cpl,nswsfc_cpl,                                      &
     &      nnirbm_cpl,nnirdf_cpl,nvisbm_cpl,nvisdf_cpl,                &
! for nst
     &      xt, xs, xu, xv, xz, zm, xtts, xzts, d_conv, ifd, dt_cool,
     &      Qrain, tref, z_c, c_0, c_d, w_0, w_d

!
      real(kind=kind_phys), dimension(ix,lsoil),      intent(inout) ::  &
     &      smc, stc, slc

      real(kind=kind_phys), dimension(ix,levs),       intent(inout) ::  &
     &      upd_mf, dwn_mf, det_mf, dqdt_v, cnvqc_v
!    &      upd_mf, dwn_mf, det_mf, dkh, rnp

      real(kind=kind_phys),                           intent(inout) ::  &
     &      phy_f3d(ix,levs,ntot3d),  phy_f2d(ix,ntot2d),               &
     &      dt3dt(ix,levs,6), du3dt(ix,levs,4), dv3dt(ix,levs,4),       &
     &      dq3dt(ix,levs,5+pl_coeff)

      real(kind=kind_phys),                           intent(inout) ::  &
     &      phy_fctd(ix,nctp)                                           &
      real(kind=kind_phys), dimension(ntrac-ncld+2) ::  fscav, fswtr

!  ---  output:
      real(kind=kind_phys), dimension(im),            intent(out) ::    &
     &      t2m,     q2m,     u10m,    v10m,    zlvl,   psurf, hpbl,    &
     &      pwat,    t1,      q1,      u1,      v1,     chh,   cmm,     &
     &      dlwsfci, ulwsfci, dswsfci, uswsfci,                         &
     &      dusfci,  dvsfci,  dtsfci,  dqsfci,                          &
     &      gfluxi,  epi,     smcwlt2, smcref2, wet1, sr

      real(kind=kind_phys), dimension(im), optional,  intent(out) ::    &
     &      dusfci_cpl,dvsfci_cpl,dtsfci_cpl,dqsfci_cpl,                &
     &      dlwsfci_cpl,dswsfci_cpl,                                    &
     &      dnirbmi_cpl,dnirdfi_cpl,dvisbmi_cpl,dvisdfi_cpl,            &
     &      nlwsfci_cpl,nswsfci_cpl,                                    &
     &      nnirbmi_cpl,nnirdfi_cpl,nvisbmi_cpl,nvisdfi_cpl,            &
     &      t2mi_cpl,q2mi_cpl,                                          &
     &      u10mi_cpl,v10mi_cpl,tseai_cpl,psurfi_cpl,                   &
     &      rqtk

      real(kind=kind_phys), dimension(ix,levs),       intent(out) ::    &
     &      gt0, gu0, gv0

      real(kind=kind_phys), dimension(ix,levs,ntrac), intent(out) ::    &
     &      gq0

!  ---  local:
!     real(kind=kind_phys) ::  fscav(ntrac-ncld-1)
!     real(kind=kind_phys),allocatable  ::  fscav(:), fswtr(:)
      real(kind=kind_phys), dimension(im)          :: ccwfac, garea,    &
     &      dlength, xncw,   cumabs, qmax,   cice,    zice,   tice,     &
!    &      gflx,    rain,   rainc,  rainl,  rain1,   raincs, evapc,    &
     &      gflx,                            rain1,   raincs,           &
     &      snowmt,  cd,     cdq,    qss,    dusfcg,  dvsfcg, dusfc1,   &
     &      dvsfc1,  dtsfc1, dqsfc1, rb,     rhscnpy, drain,  cld1d,    &
     &      evap,    hflx,   stress, t850,   ep1d,    gamt,   gamq,     &
     &      sigmaf,                  oc,     theta,   gamma,  sigma,    &
     &      elvmax,  wind,   work1,  work2,  runof,   xmu,              &
!    &      elvmax,  wind,   work1,  work2,  runof,   xmu,    oro_land, &
     &      fm10,    fh2,    tsurf,  tx1,    tx2,     ctei_r, flgmin_l, &
     &      evbs,    evcw,   trans,  sbsno,  snowc,   frland,           &
     &      adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw,gabsbdlw,        &
     &      adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                 &
     &      adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd,                 &
     &      xcosz,  tseal,  snohf,  dlqfac,  work3, ctei_rml, cldf,     &
     &      domr,   domzr,  domip,  doms,   psautco_l, prautco_l

      real(kind=kind_phys), dimension(im)          :: ocalnirbm_cpl,    &
     &      ocalnirdf_cpl,ocalvisbm_cpl,ocalvisdf_cpl

!    &      dswsfc, radsl,                                              &
!    &      dlwsf1,  ulwsf1, xcosz,  tseal,  snohf,   dlqfac,           &
!    &      domr,    domzr,  domip,  doms

!     real(kind=kind_phys), dimension(ix,levs)     :: ud_mf, dd_mf,     &
!    &      dt_mf, del
      real(kind=kind_phys), dimension(ix,levs)     :: del, dtdtr
      real(kind=kind_phys), dimension(im,levs-1)   :: dkt

      real(kind=kind_phys), dimension(im,levs)     :: rhc, dtdt,        &
     &      dudt, dvdt, gwdcu, gwdcv, dtdtc, dmdt,                      &
!    &      diagn1, diagn2, cuhr, cumchr,                               &
     &      qr_col, fc_ice, rainp, ud_mf, dd_mf, dt_mf, prnum
!    &      qr_col, fc_ice, rainp, ud_mf, dd_mf, dt_mf, shoc_cld, prnum

      real(kind=kind_phys), dimension(im,lsoil)    :: smsoil, stsoil,   &
     &      ai, bi, cci, rhsmc, zsoil, slsoil

      real(kind=kind_phys)                :: zsea1,zsea2
      real(kind=kind_phys), dimension(im) :: dtzm

      real(kind=kind_phys) :: rhbbot, rhbtop, rhpbl, frain, f_rain,     &
     &      f_ice, qi, qw, qr, wc, tem, tem1, tem2,  sume,  sumr, sumq, &
     &      dqdt(im,levs,ntrac), oa4(im,4), clx(im,4), albbm, xcosz_loc
      real(kind=kind_phys), parameter :: albdf=0.06


!           in clw, the first two varaibles are cloud water and ice.
!           from third to ntrac are convective transportable tracers,
!           third being the ozone, when ntrac=3 (valid only with ras)

      real(kind=kind_phys), allocatable :: clw(:,:,:), qpl(:,:),qpi(:,:)
     &,                                    ncpl(:,:),  ncpi(:,:)

      integer, dimension(im) :: kbot, ktop, kcnv, soiltyp, vegtype,     &
     &          kpbl, slopetyp, kinver, lmh, levshc, islmsk

      integer :: i, nvdiff, kk, ic, k, n, ipr, lv, k1, iter, levshcm,   &
     &           tracers, trc_shft, tottracer, num2, num3               &
     &,          nshocm, nshoc, ntk, ntln, ntin

      logical, dimension(im) :: flag_iter, flag_guess, invrsn           &
     &,                         skip_macro

      real(kind=kind_phys), dimension(im)          :: dtsfc_cice,       &
     &    dqsfc_cice, dusfc_cice, dvsfc_cice, ulwsfc_cice, tisfc_cice,  &
     &    tsea_cice, hice_cice, fice_cice

      integer, dimension(im) :: islmsk_cice
      logical, dimension(im) :: flag_cice

      logical :: lprnt, revap

      real(kind=kind_phys), allocatable :: cnvc(:,:),cnvw(:,:)
      real(kind=kind_phys) eng0, eng1, dtshoc
!
! for CS-convection
!     real(kind=kind_phys), parameter :: wcbmax1=2.8, wcbmax2=1.4
      real(kind=kind_phys), parameter :: wcbmax1=2.5, wcbmax2=1.5
!     real(kind=kind_phys), parameter :: wcbmax1=1.4, wcbmax2=1.4
      real(kind=kind_phys)  wcbmax(im)
!
      real(kind=kind_phys) tf, tcr, tcrf
!     parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
      parameter (tf=258.16, tcr=273.16, tcrf=1.0/(tcr-tf))
!
!     for 2 M microphysics
      real(kind=kind_phys), allocatable, dimension(:,:) :: qlcn, qicn
     &,      w_upi,cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,CLCN,        !Acheng
     &       CNV_FICE,CNV_NDROP,CNV_NICE
      real(kind=kind_phys), allocatable, dimension(:) :: cn_prc,cn_snr
!
!
!===> ...  begin here
!
! The option to compute tracer scavenging in GSM is disabled
!     do i=1, ntrac-ncld-1
!       fscav(i) = 0.
!     enddo


!  --- ...  set up check print point (for debugging)
!
!*************************************************************************
      lprnt = .false.
!     lprnt = me == 0 .and. kdt < 10 .and. lat == 22
      ipr = 1

!     if (me == 0 .and. kdt < 5)
!    &write(0,*)' In gbphys:', im,ix,levs,lsoil,lsm,
!    &  ntrac,ncld,ntoz,ntcw,
!    &  nmtvr,nrcm,ko3,lonr,latr,jcap,num_p3d,num_p2d,npdf3d,
!    &      kdt,lat,me,pl_coeff,ncw,flgmin,crtrh,cdmbgwd
!    &,' ccwf=',ccwf,' dlqf=',dlqf,' ras=',ras,
!    & ' evpco=',evpco,' wminco=',wminco,' levr=',levr
!     do i = 1, im
!       lprnt = kdt >= 256 .and. abs(xlon(i)*57.29578-136.07) < 0.08    &
!    &                   .and. abs(xlat(i)*57.29578+1.1)  < 0.101
!       lprnt = kdt >= 0 .and. abs(xlon(i)*57.29578-13.125) < 0.501     &
!    &                   .and. abs(xlat(i)*57.29578-79.043)  < 0.501
!       lprnt = kdt >= 0 .and. abs(xlon(i)*57.29578-8.4375) < 0.501     &
!    &                   .and. abs(xlat(i)*57.29578+1.4)  < 0.501
!       lprnt = kdt >= 40 .and. abs(xlon(i)*57.29578-288.75) < 1.501    &
!    &                   .and. abs(xlat(i)*57.29578+12.38)  < 1.501
!       lprnt = kdt >= 0 .and. abs(xlon(i)*57.29578-135.0) < 0.201      &
!    &                   .and. abs(xlat(i)*57.29578-10.476)  < 0.201
!       lprnt = kdt >= 0 .and. abs(xlon(i)*57.29578-110.3) < 0.201      &
!    &                   .and. abs(xlat(i)*57.29578-2.0)   < 0.201
!       if (kdt == 10)
!    &  print *,' i=',i,' xlon=',xlon(i)*57.29578,                      &
!    &                  ' xlat=',xlat(i)*57.29578,' i=',i,' me=',me
!       if (lprnt) then
!         ipr = i
!         exit
!       endif
!     enddo

!     lprnt = .false.
!     if(lprnt) then
!       write(0,*)' im=',im,' ix=',ix,' levs=',levs,' lsoil=',lsoil,    &
!    &   ' ntrac=',ntrac,' ntoz=',ntoz,' ntcw=',ntcw,' me=',me,         &
!    &   ' xlat=',xlat(ipr),' kdt=',kdt,' slmsk=',slmsk(ipr),           &
!    & ' ntke=',ntke,' num_p3d=',num_p3d,' xlon=',xlon(ipr)
!    &   ' tsea=',tsea(ipr),' tref=',tref(ipr),' dt_cool=',dt_cool(ipr),&
!    &   ' dt_warm=',2.0*xt(ipr)/xz(ipr),' nrcm=',nrcm,' xlon=',
!    &    xlon(ipr),                                                    &
!    &   ' dt_warm=',dt_warm(ipr),' nrcm=',nrcm,' xlon=',xlon(ipr),     &
!    &   ' sfalb=',sfalb(ipr),' kdt=',kdt
!       write(0,*) ' pgr=',pgr(ipr),' kdt=',kdt,' ipr=',ipr
!       write(0,*)' ipr=',ipr,' phy_f2d=',phy_f2d(ipr,1:num_p2d)
!       write(0,*),' ugrs=',ugrs(ipr,:)
!       write(0,*)' vgrs=',vgrs(ipr,:)
!       write(0,*)' tgrs=',tgrs(ipr,:),' kdt=',kdt,' ipr=',ipr
!    &,  ' xlon=',xlon(ipr),' xlat=',xlat(ipr)
!       write(0,*)' qgrs=',qgrs(ipr,:,1)
!       write(0,*)' ozg=',qgrs(ipr,:,2)
!       write(0,*)' tke=',qgrs(ipr,:,4)
!       print *,' clwb=',qgrs(ipr,:,ntiw)
!    &,  ' xlon=',xlon(ipr),' xlat=',xlat(ipr)
!     endif
!
!*************************************************************************
!
      do i = 1, im
        if(nint(slmsk(i)) == 1) then
          frland(i) = 1.0
        else
          frland(i) = 0.
        endif
      enddo

      nvdiff = ntrac           ! vertical diffusion of all tracers!
!
!  --- ...                       figure out number of extra tracers
!
      tottracer = 0            ! no convective transport of tracers
      if (trans_trac) then
        if (ntcw > 0) then
          if (ntoz < ntcw) then
            trc_shft = ntcw + ncld - 1
          else
            trc_shft = ntoz
          endif
        elseif (ntoz > 0) then
          trc_shft = ntoz
        else
          trc_shft = 1
        endif

        tracers   = ntrac - trc_shft
        tottracer = tracers
        if (ntoz > 0) tottracer = tottracer + 1  ! ozone is added separately
      endif
      if (ntke > 0) ntk = ntke - trc_shft + 3

!     if (lprnt) write(0,*)' trans_trac=',trans_trac,' tottracer=',     &
!                write(0,*)' trans_trac=',trans_trac,' tottracer=',     &
!    &                   tottracer,' trc_shft=',trc_shft,' kdt=',kdt
!    &,                  ntrac-ncld+2,' clstp=',clstp,' kdt=',kdt
!    &,' ntk=',ntk,' lat=',lat

      skip_macro = .false.
      allocate ( clw(ix,levs,tottracer+2) )
      if (do_shoc) then
        allocate (qpl(im,levs),  qpi(im,levs)
     &,           ncpl(im,levs), ncpi(im,levs))
        do k=1,levs
          do i=1,im
            ncpl(i,k) = 0.0
            ncpi(i,k) = 0.0
          enddo
        enddo
      endif

      if (.not. ras .or. .not. cscnv) then
        allocate ( cnvc(ix,levs), cnvw(ix,levs))
      endif
!     allocate (fscav(tottracer+3), fswtr(tottracer+3))

      if (ncld == 2) then         ! For MGB double moment microphysics
        allocate (qlcn(im,levs),      qicn(im,levs),   w_upi(im,levs)
     &,           cf_upi(im,levs),    CNV_MFD(im,levs),CNV_PRC3(im,levs)
     &,           CNV_DQLDT(im,levs), clcn(im,levs),   cnv_fice(im,levs)
     &,           cnv_ndrop(im,levs), cnv_nice(im,levs))
        allocate(cn_prc(im), cn_snr(im))
      else
        allocate (qlcn(1,1),      qicn(1,1),   w_upi(1,1)
     &,           cf_upi(1,1),    CNV_MFD(1,1),CNV_PRC3(1,1)
     &,           CNV_DQLDT(1,1), clcn(1,1),   cnv_fice(1,1)
     &,           cnv_ndrop(1,1), cnv_nice(1,1))
      endif

! The option to compute tracer scavenging in GSM is disabled
      do i=1, tottracer+3
        fscav(i) = 0.
        fswtr(i) = 0.
      enddo

      if (nnp == 1) then
        do n=1,ntrac
          if (abs(1.0-adjtrc(n)) > 1.0e-7) then
            do k=1,levs
              do i=1,im
                qgrs(i,k,n) = qgrs(i,k,n) * adjtrc(n)
              enddo
            enddo
          endif
        enddo
      endif
!
      call get_prs(im,ix,levs,ntrac,tgrs,qgrs,                          &
     &             thermodyn_id, sfcpress_id,                           &
     &             gen_coord_hybrid,                                    &
     &             prsi,prsik,prsl,prslk,phii,phil,del)
!    &             prsi,prsik,prsl,prslk,phii,phil,del,lprnt)
!
!     if (lprnt) then
!       write(0,*)' prsi=',prsi(ipr,:)
!       write(0,*)' prsik=',prsik(ipr,:),' me=',me,' kdt=',kdt
!       write(0,*)' prslk=',prslk(ipr,:),' me=',me,' kdt=',kdt
!       write(0,*)' phil=',phil(ipr,:),' me=',me,' kdt=',kdt,' ipr=',ipr
!    &,' lat=',lat,' prsi=',prsi(ipr,1)
!       write(0,*)' prsl=',prsl(ipr,:),' kdt=',kdt
!       print *,' del=',del(ipr,:)
!     endif
!
      rhbbot = crtrh(1)
      rhpbl  = crtrh(2)
      rhbtop = crtrh(3)
!
!  --- ...  frain=factor for centered difference scheme correction of rain amount.

      frain = dtf / dtp

      do i = 1, im
        sigmaf(i)   = max( vfrac(i),0.01 )
!       sigmaf(i)   = max( vfrac(i),0.3 )
        if (lsm == 0) sigmaf(i)   =  0.5 + vfrac(i) * 0.5

        islmsk(i)   = nint(slmsk(i))

        if (islmsk(i) == 2) then
         if (isot == 1) then
          soiltyp(i)  = 16
         else
          soiltyp(i)  = 9
         endif
         if (ivegsrc == 1) then
          vegtype(i)  = 15
         elseif(ivegsrc == 2) then
          vegtype(i)  = 13
         endif
          slopetyp(i) = 9                      !! clu: qa(slopetyp)
        else
          soiltyp(i)  = int( stype(i)+0.5 )
          vegtype(i)  = int( vtype(i)+0.5 )
          slopetyp(i) = int( slope(i)+0.5 )    !! clu: slope -> slopetyp
        endif
!  --- ...  xw: transfer ice thickness & concentration from global to local variables
        zice(i) = hice(i)
        cice(i) = fice(i)
        tice(i) = tisfc(i)

        if (lssav_cpl) then
          islmsk_cice(i) = nint(slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)

          ulwsfc_cice(i) = ulwsfcin_cpl(i)
          dusfc_cice(i)  = dusfcin_cpl(i)
          dvsfc_cice(i)  = dvsfcin_cpl(i)
          dtsfc_cice(i)  = dtsfcin_cpl(i)
          dqsfc_cice(i)  = dqsfcin_cpl(i)
          tisfc_cice(i)  = tisfc(i)
          tsea_cice(i)   = tsea(i)
          fice_cice(i)   = fice(i)
          hice_cice(i)   = hice(i)
        endif

        work1(i)   = (log(coslat(i) / (nlons(i)*latr)) - dxmin) * dxinv
        work1(i)   = max(0.0, min(1.0,work1(i)))
        work2(i)   = 1.0 - work1(i)
        psurf(i)   = pgr(i)
        work3(i)   = prsik(i,1) / prslk(i,1)
        tem1       = con_rerth * (con_pi+con_pi)*coslat(i)/nlons(i)
        tem2       = con_rerth * con_pi / latr
        garea(i)   = tem1 * tem2
        dlength(i) = sqrt( tem1*tem1+tem2*tem2 )
        cldf(i)    = cgwf(1)*work1(i) + cgwf(2)*work2(i)
        wcbmax(i)  = wcbmax1*work1(i) + wcbmax2*work2(i)
      enddo

!     if (lprnt) write(0,*)' in gbphys work1&2=',work1(ipr),work2(ipr)
!    &,' dxmin=',dxmin,' dxinv=',dxinv,' dx=',
!    &           log(coslat(ipr) / (nlons(ipr)*latr))
!    &,' coslat=',coslat(ipr),' nlons=',nlons(ipr),' latr=',latr

!  --- ...  transfer soil moisture and temperature from global to local variables

      do k = 1, lsoil
        do i = 1, im
          smsoil(i,k) = smc(i,k)
          stsoil(i,k) = stc(i,k)
          slsoil(i,k) = slc(i,k)          !! clu: slc -> slsoil
        enddo
      enddo

      do k = 1, levs
        do i = 1, im
          dudt(i,k)  = 0.
          dvdt(i,k)  = 0.
          dtdt(i,k)  = 0.
          dtdtc(i,k) = 0.
        enddo
      enddo

      do n = 1, ntrac
        do k = 1, levs
          do i = 1, im
            dqdt(i,k,n) = 0.
          enddo
        enddo
      enddo

!  --- ...  initialize dtdt with heating rate from dcyc2

!     if (lprnt) then
!       do ipr=1,im
!         write(0,*)' before dcyc2: im=',im,' lsoil=',lsoil,' levs=',   &
!    &                                                             levs &
!    &,    ' sde=',sdec,' cdec=',cdec,' tsea=',tsea(ipr),' ipr=',ipr    &
!    &,    ' lat=',lat,' me=',me,' kdt=',kdt                            &
!    &,    ' sfcdlw=',sfcdlw(ipr),' sfcnsw=',sfcnsw(ipr)
!         print *,' hlw=',hlw(ipr,:),' me=',me,' lat=',lat,xlon(ipr)
!         print *,' swh=',swh(ipr,:),' me=',me,' lat=',lat,xlon(ipr)
!       enddo
!     endif

!  --- ...  adjust mean radiation fluxes and heating rates to fit for
!           faster model time steps.
!      sw:  using cos of zenith angle as scaling factor
!      lw:  using surface air skin temperature as scaling factor

      if (pre_rad) then

        call dcyc2t3_pre_rad                                            &
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsea,tgrs(1,1),tgrs(1,1),                      &
     &       sfcdsw,sfcnsw,sfcdlw,swh,hlw,                              &
     &       sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   &
     &       sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   &
     &       ix, im, levs,                                              &
!  ---  input/output:
     &       dtdt,                                                      &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz,         &
     &       adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   &
     &       adjnirbmd,adjnirdfd,adjvisbmd,adjvisdfd                    &
     &     )

      else

        call dcyc2t3                                                    &
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsea,tgrs(1,1),tsflw,sfcemis,                  &
     &       sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    &
     &       sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   &
     &       sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   &
     &       ix, im, levs,                                              &
!  ---  input/output:
     &       dtdt,dtdtc,                                                &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz,         &
     &       adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   &
     &       adjnirbmd,adjnirdfd,adjvisbmd,adjvisdfd                    &
     &     )

!
! save temp change due to radiation - need for sttp stochastic physics
!---------------------------------------------------------------------
        do k=1,levs
          do i=1,im
             dtdtr(i,k) = dtdtr(i,k) + dtdtc(i,k)*dtf
          enddo
        enddo

      endif
!
      if (lsidea) then                       !idea jw
        do k = 1, levs
          do i = 1, im
!           dtdt(i,k) = hlwd(i,k,2)
            dtdt(i,k) = 0.
          enddo
        enddo
      endif

!  ---  convert lw fluxes for land/ocean/sea-ice models
!  note: for sw: adjsfcdsw and adjsfcnsw are zenith angle adjusted downward/net fluxes.
!        for lw: adjsfcdlw is (sfc temp adjusted) downward fluxe with no emiss effect.
!                adjsfculw is (sfc temp adjusted) upward fluxe including emiss effect.
!        one needs to be aware that that the absorbed downward lw flux (used by land/ocean
!        models as downward flux) is not the same as adjsfcdlw but a value reduced by
!        the factor of emissivity.  however, the net effects are the same when seeing
!        it either above the surface interface or below.
!
!   - flux above the interface used by atmosphere model:
!        down: adjsfcdlw;    up: adjsfculw = sfcemis*sigma*T**4 + (1-sfcemis)*adjsfcdlw
!        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
!   - flux below the interface used by lnd/oc/ice models:
!        down: sfcemis*adjsfcdlw;  up: sfcemis*sigma*T**4
!        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)

!  --- ...  define the downward lw flux absorbed by ground

      do i = 1, im
        gabsbdlw(i) = sfcemis(i) * adjsfcdlw(i)
      enddo

!     if( lsidea ) then     ! idea : moved temp adjust to idea_phys
!       print *,' in gbphys: lsidea is true '
!       DTDT = 0.
!     endif

      if (lssav) then      !  --- ...  accumulate/save output variables

!  --- ...  sunshine duration time is defined as the length of time (in mdl output
!           interval) that solar radiation falling on a plane perpendicular to the
!           direction of the sun >= 120 w/m2

        do i = 1, im
          if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
            tem1 = adjsfcdsw(i) / xcosz(i)

            if ( tem1 >= 120.0 ) then
              suntim(i) = suntim(i) + dtf
            endif
          endif
        enddo

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output

        do i = 1, im
          dlwsfc(i) = dlwsfc(i) + adjsfcdlw(i)*dtf
          if (lssav_cpl) then
            if (flag_cice(i)) adjsfculw(i) = ulwsfc_cice(i)
          endif
          ulwsfc(i) = ulwsfc(i) + adjsfculw(i)*dtf
          psmean(i) = psmean(i) + pgr(i)*dtf        ! mean surface pressure
        enddo

        if (ldiag3d) then
          if( lsidea ) then
            do k = 1, levs
              do i = 1, im
                dt3dt(i,k,1) = dt3dt(i,k,1) + hlwd(i,k,1)*dtf
                dt3dt(i,k,2) = dt3dt(i,k,2) + hlwd(i,k,2)*dtf
                dt3dt(i,k,3) = dt3dt(i,k,3) + hlwd(i,k,3)*dtf
                dt3dt(i,k,4) = dt3dt(i,k,4) + hlwd(i,k,4)*dtf
                dt3dt(i,k,5) = dt3dt(i,k,5) + hlwd(i,k,5)*dtf
                dt3dt(i,k,6) = dt3dt(i,k,6) + hlwd(i,k,6)*dtf
              enddo
            enddo
          else
            do k = 1, levs
              do i = 1, im
                dt3dt(i,k,1) = dt3dt(i,k,1) + hlw(i,k)*dtf
                dt3dt(i,k,2) = dt3dt(i,k,2) + swh(i,k)*dtf*xmu(i)
              enddo
            enddo
          endif
        endif

      endif    ! end if_lssav_block

      do i = 1, im
        kcnv(i)   = 0
        kinver(i) = levs
        invrsn(i) = .false.
        tx1(i)    = 0.0
        tx2(i)    = 10.0
        ctei_r(i) = 10.0
      enddo

!    Only used for old shallow convection with mstrat=.true.

      if (((imfshalcnv == 0 .and. shal_cnv) .or. old_monin)             &
     &                                     .and. mstrat) then
        do i = 1, im
          ctei_rml(i) = ctei_rm(1)*work1(i) + ctei_rm(2)*work2(i)
        enddo
        do k = 1, levs/2
          do i = 1, im
            if (prsi(i,1)-prsi(i,k+1) < 0.35*prsi(i,1)                  &
     &          .and. (.not. invrsn(i))) then
              tem = (tgrs(i,k+1)-tgrs(i,k)) / (prsl(i,k)-prsl(i,k+1))

              if ((tem > 0.00010 .and. tx1(i) < 0.0) .or.
     &            (tem-abs(tx1(i)) > 0.0 .and. tx2(i) < 0.0)) then
                invrsn(i) = .true.

                if (qgrs(i,k,1) > qgrs(i,k+1,1)) then
                  tem1 = tgrs(i,k+1) + hocp*max(qgrs(i,k+1,1),qmin)
                  tem2 = tgrs(i,k)   + hocp*max(qgrs(i,k,1),qmin)

                  tem1 = tem1 / prslk(i,k+1) - tem2 / prslk(i,k)

!  --- ...  (cp/l)(deltathetae)/(deltatwater) > ctei_rm -> conditon for CTEI
                  ctei_r(i) = (1.0/hocp)*tem1/(qgrs(i,k+1,1)-qgrs(i,k,1)&
     &                      + qgrs(i,k+1,ntcw)-qgrs(i,k,ntcw))
                else
                  ctei_r(i) = 10
                endif

                if ( ctei_rml(i) > ctei_r(i) ) then
                  kinver(i) = k
                else
                  kinver(i) = levs
                endif
              endif

              tx2(i) = tx1(i)
              tx1(i) = tem
            endif
          enddo
        enddo
      endif

!  --- ...  check print

!     ipr = 1
!     if (lprnt) then
!       write(0,*)' before progtm: im=',im,' lsoil=',lsoil              &
!    &,         ' nvdiff=',nvdiff,' adjsfcnsw=',adjsfcnsw(ipr)          &
!    &,         ' adjsfcdlw=',adjsfcdlw(ipr),'adjsfculw=',adjsfculw(ipr)&
!    &,         ' sfcemis=',sfcemis(ipr),' tsea2=',tsea(ipr)            &
!    &,         ' ipr=',ipr,' me=',me,' lat=',lat,' xlon=',xlon(ipr)    &
!    &,         ' kdt=',kdt
!       write(0,*)' dtdth=',dtdt(ipr,:),' kdt=',kdt
!     endif


!  --- ...  lu: initialize flag_guess, flag_iter, tsurf

      do i = 1, im
        tsurf(i)      = tsea(i)
        flag_guess(i) = .false.
        flag_iter(i)  = .true.
        drain(i)      = 0.0
        ep1d(i)       = 0.0
        runof(i)      = 0.0
        hflx(i)       = 0.0
        evap(i)       = 0.0

        evbs(i)       = 0.0
        evcw(i)       = 0.0
        trans(i)      = 0.0
        sbsno(i)      = 0.0
        snowc(i)      = 0.0
        snohf(i)      = 0.0
        zlvl(i)       = phil(i,1) / con_g
        smcwlt2(i)    = 0.0
        smcref2(i)    = 0.0
      enddo

!  --- ...  lu: iter-loop over (sfc_diff,sfc_drv,sfc_ocean,sfc_sice)

      do iter = 1, 2

!  --- ...  surface exchange coefficients
!
!     if (lprnt) write(0,*)' tsea=',tsea(ipr),' tsurf=',tsurf(ipr),iter
        call sfc_diff(im,pgr,ugrs,vgrs,tgrs,qgrs,zlvl,                  &
     &         snwdph,tsea,zorl,cd,cdq,rb,                              &
     &                prsl(1,1),work3,islmsk,                           &
     &                stress,ffmm,ffhh,                                 &
     &                uustar,wind,phy_f2d(1,num_p2d),fm10,fh2,          &
     &                sigmaf,vegtype,shdmax,ivegsrc,                    &
     &                tsurf, flag_iter, redrag)

!       if (lprnt) write(0,*)' cdq=',cdq(ipr),' iter=',iter             &
!    &,   ' wind=',wind(ipr),'phy_f2d=',phy_f2d(ipr,num_p2d),' ugrs='   &
!    &,   ugrs(ipr,1),' vgrs=',vgrs(ipr,1)

!  --- ...  lu: update flag_guess

        do i = 1, im
          if (iter == 1 .and. wind(i) < 2.0) then
            flag_guess(i) = .true.
          endif
        enddo

        if ( nstf_name(1) > 0 ) then

          do i = 1, im
            if ( islmsk(i) == 0 ) then
              tem      = (oro(i)-oro_uf(i)) * rlapse
              tseal(i) = tsea(i)  + tem
              tsurf(i) = tsurf(i) + tem
            endif
          enddo

!       if (lprnt) write(0,*)' tseaz1=',tsea(ipr),' tref=',tref(ipr)
!    &,   ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*(xt(ipr)/xz(ipr)
!    &,   ' kdt=',kdt
!    &,   ' tgrs=',tgrs(ipr,1),' prsl=',prsl(ipr,1)
!    &,   ' work3=',work3(ipr),' kdt=',kdt

          call sfc_nst                                                  &
     &       ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,tref,cd,cdq,            &
     &         prsl(1,1),work3,islmsk,xlon,sinlat,stress,               &
     &         sfcemis,gabsbdlw,adjsfcnsw,tprcp,dtf,kdt,solhr,xcosz,    &
     &         phy_f2d(1,num_p2d),flag_iter,flag_guess,nstf_name,       &
     &         lprnt,ipr,                                               &
!  --- Input/output
     &         tseal,tsurf,xt,xs,xu,xv,xz,zm,xtts,xzts,dt_cool,         &
     &         z_c,c_0,c_d,w_0,w_d,d_conv,ifd,qrain,                    &
!  ---  outputs:
     &         qss, gflx, cmm, chh, evap, hflx, ep1d)

!         if (lprnt) print *,' tseaz2=',tseal(ipr),' tref=',tref(ipr),
!    &     ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*xt(ipr)/xz(ipr),
!    &     ' kdt=',kdt

          do i = 1, im
            if ( islmsk(i) == 0 ) then
              tsurf(i) = tsurf(i) - (oro(i)-oro_uf(i)) * rlapse
            endif
          enddo

!  --- ...  run nsst model  ... ---

          if ( nstf_name(1) > 1 ) then
            zsea1 = 0.001*real(nstf_name(4))
            zsea2 = 0.001*real(nstf_name(5))
            call get_dtzm_2d(xt,xz,dt_cool,z_c,slmsk,
     &                       zsea1,zsea2,im,1,dtzm)
            do i = 1, im
              if ( islmsk(i) == 0 ) then
              tsea(i) = max(271.2,tref(i) + dtzm(i))
     &                      -(oro(i)-oro_uf(i))*rlapse
              endif
            enddo
          endif

!         if (lprnt) print *,' tseaz2=',tsea(ipr),' tref=',tref(ipr),   &
!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

        else

!  --- ...  surface energy balance over ocean

          call sfc_ocean                                                &
!  ---  inputs:
     &     ( im,pgr,ugrs,vgrs,tgrs,qgrs,tsea,cd,cdq,                    &
     &       prsl(1,1),work3,islmsk,phy_f2d(1,num_p2d),flag_iter,       &
!  ---  outputs:
     &       qss,cmm,chh,gflx,evap,hflx,ep1d                            &
     &     )

        endif       ! if ( nstf_name(1) > 0 ) then

!       if (lprnt) write(0,*)' sfalb=',sfalb(ipr),' ipr=',ipr           &
!    &,   ' weasd=',weasd(ipr),' snwdph=',snwdph(ipr)                   &
!    &,   ' tprcp=',tprcp(ipr),' kdt=',kdt,' iter=',iter
!    &,' tseabefland=',tsea(ipr)

!  --- ...  surface energy balance over land
!
        if (lsm == 1) then                          ! noah lsm call

!     if (lprnt) write(0,*)' tsead=',tsea(ipr),' tsurf=',tsurf(ipr),iter
!    &,' pgr=',pgr(ipr),' sfcemis=',sfcemis(ipr)

          call sfc_drv                                                  &

!  ---  inputs:
     &     ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,soiltyp,vegtype,sigmaf,   &
     &       sfcemis,gabsbdlw,adjsfcdsw,adjsfcnsw,dtf,tg3,cd,cdq,       &
     &       prsl(1,1),work3,zlvl,islmsk,phy_f2d(1,num_p2d),slopetyp,   &
     &       shdmin,shdmax,snoalb,sfalb,flag_iter,flag_guess,           &
     &       isot,ivegsrc,                                              &
!  ---  in/outs:
     &       weasd,snwdph,tsea,tprcp,srflag,smsoil,stsoil,slsoil,       &
     &       canopy,trans,tsurf,zorl,                                   &
!  ---  outputs:
     &       sncovr,qss,gflx,drain,evap,hflx,ep1d,runof,                &
     &       cmm,chh,evbs,evcw,sbsno,snowc,soilm,snohf,                 &
     &       smcwlt2,smcref2,wet1                                       &
     &     )

!     if (lprnt) write(0,*)' tseae=',tsea(ipr),' tsurf=',tsurf(ipr),iter
!    &,' phy_f2d=',phy_f2d(ipr,num_p2d)

        else                                       ! osu lsm call

          call sfc_land                                                 &
!  ---  inputs:
     &     ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,smsoil,soiltyp,           &
     &       sigmaf,vegtype,sfcemis,adjsfcdlw,adjsfcnsw,dtf,            &
     &            tg3,cd,cdq,prsl(1,1),work3,islmsk,                    &
!    &       zorl,tg3,cd,cdq,prsl(1,1),work3,islmsk,                    &
     &       phy_f2d(1,num_p2d),flag_iter,flag_guess,                   &
!  ---  input/outputs:
     &       weasd,tsea,tprcp,srflag,stsoil,canopy,tsurf,               &
!  ---  outputs:
     &       qss,snowmt,gflx,zsoil,rhscnpy,rhsmc,                       &
     &       ai,bi,cci,drain,evap,hflx,ep1d,cmm,chh,                    &
     &       evbs,evcw,trans,sbsno,snowc,soilm,                         &
     &       snohf,smcwlt2,smcref2                                      &
     &     )

        endif

!       if (lprnt) write(0,*)' tseabeficemodel =',tsea(ipr),' me=',me   &
!    &,   ' kdt=',kdt

!  --- ...  surface energy balance over seaice

        if (lssav_cpl) then
          do i = 1, im
            if (flag_cice(i)) then
               islmsk (i) = islmsk_cice(i)
            endif
          enddo
        endif

        call sfc_sice                                                   &
!  ---  inputs:
     &     ( im,lsoil,pgr,ugrs,vgrs,tgrs,qgrs,dtf,                      &
     &       sfcemis,gabsbdlw,adjsfcnsw,adjsfcdsw,srflag,               &
     &       cd,cdq,prsl(1,1),work3,islmsk,phy_f2d(1,num_p2d),          &
     &       flag_iter,mom4ice,lsm, lprnt,ipr,                          &
!    &       flag_iter,mom4ice,lsm,                                     &
!  ---  input/outputs:
     &       zice,cice,tice,weasd,tsea,tprcp,stsoil,ep1d,               &
!  ---  outputs:
     &       snwdph,qss,snowmt,gflx,cmm,chh,evap,hflx                   &
     &     )

        if (lssav_cpl) then
          do i = 1, im
            if (flag_cice(i)) then
               islmsk(i) = nint(slmsk(i))
            endif
          enddo

          call sfc_cice                                                 &
!  ---     inputs:
     &     ( im,ugrs,vgrs,tgrs,qgrs,cd,cdq,prsl(1,1),work3,             &
     &       islmsk_cice,phy_f2d(1,num_p2d),flag_iter,                  &
     &       dqsfc_cice,dtsfc_cice,                                     &
!  ---     outputs:
     &       qss,cmm,chh,evap,hflx                                      &
     &     )
        endif

!  --- ...  lu: update flag_iter and flag_guess

        do i = 1, im
          flag_iter(i)  = .false.
          flag_guess(i) = .false.

          if(islmsk(i) == 1 .and. iter == 1) then
            if (wind(i) < 2.0) flag_iter(i) = .true.
          elseif (islmsk(i) == 0 .and. iter == 1                        &
     &                           .and. nstf_name(1) > 0) then
            if (wind(i) < 2.0) flag_iter(i) = .true.
          endif
        enddo

      enddo   ! end iter_loop

      do i = 1, im
        epi(i)     = ep1d(i)
        dlwsfci(i) = adjsfcdlw(i)
        ulwsfci(i) = adjsfculw(i)
        uswsfci(i) = adjsfcdsw(i) - adjsfcnsw(i)
        dswsfci(i) = adjsfcdsw(i)
        gfluxi(i)  = gflx(i)
        t1(i)      = tgrs(i,1)
        q1(i)      = qgrs(i,1,1)
        u1(i)      = ugrs(i,1)
        v1(i)      = vgrs(i,1)
      enddo

      if (lsm == 0) then                          ! osu lsm call
        do i = 1, im
         sncovr(i) = 0.0
         if (weasd(i) > 0.0) sncovr(i) = 1.0
        enddo
      endif

!  --- ...  update near surface fields

      call sfc_diag(im,pgr,ugrs,vgrs,tgrs,qgrs,                         &
     &              tsea,qss,f10m,u10m,v10m,t2m,q2m,work3,              &
     &              evap,ffmm,ffhh,fm10,fh2)

      do i = 1, im
        phy_f2d(i,num_p2d) = 0.0
      enddo

      if (lssav_cpl) then
        do i = 1, im
          dlwsfci_cpl(i)   = adjsfcdlw(i)
          dswsfci_cpl(i)   = adjsfcdsw(i)
          dlwsfc_cpl(i)    = dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
          dswsfc_cpl(i)    = dswsfc_cpl(i) + adjsfcdsw(i)*dtf
          dnirbmi_cpl(i)   = adjnirbmd(i)
          dnirdfi_cpl(i)   = adjnirdfd(i)
          dvisbmi_cpl(i)   = adjvisbmd(i)
          dvisdfi_cpl(i)   = adjvisdfd(i)
          dnirbm_cpl(i)    = dnirbm_cpl(i) + adjnirbmd(i)*dtf
          dnirdf_cpl(i)    = dnirdf_cpl(i) + adjnirdfd(i)*dtf
          dvisbm_cpl(i)    = dvisbm_cpl(i) + adjvisbmd(i)*dtf
          dvisdf_cpl(i)    = dvisdf_cpl(i) + adjvisdfd(i)*dtf
          nlwsfci_cpl(i)   = adjsfcdlw(i)  - adjsfculw(i)
          nlwsfc_cpl(i)    = nlwsfc_cpl(i) + nlwsfci_cpl(i)*dtf
          t2mi_cpl(i)      = t2m(i)
          q2mi_cpl(i)      = q2m(i)
          u10mi_cpl(i)     = u10m(i)
          v10mi_cpl(i)     = v10m(i)
          tseai_cpl(i)     = tsea(i)
          psurfi_cpl(i)    = pgr(i)
        enddo

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

        do i = 1, im
          if (islmsk(i) /= 1) then  ! not a land point
!  ---  compute open water albedo
            xcosz_loc = max( 0.0, min( 1.0, xcosz(i) ))
            ocalnirdf_cpl(i) = 0.06
            ocalnirbm_cpl(i) = max(albdf, 0.026/(xcosz_loc**1.7+0.065)  &
     &                       + 0.15 * (xcosz_loc-0.1) * (xcosz_loc-0.5) &
     &                       * (xcosz_loc-1.0))
            ocalvisdf_cpl(i) = 0.06
            ocalvisbm_cpl(i) = ocalnirbm_cpl(i)

            nnirbmi_cpl(i) = adjnirbmd(i)-adjnirbmd(i)*ocalnirbm_cpl(i)
            nnirdfi_cpl(i) = adjnirdfd(i)-adjnirdfd(i)*ocalnirdf_cpl(i)
            nvisbmi_cpl(i) = adjvisbmd(i)-adjvisbmd(i)*ocalvisbm_cpl(i)
            nvisdfi_cpl(i) = adjvisdfd(i)-adjvisdfd(i)*ocalvisdf_cpl(i)
          else
            nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
            nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
            nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
            nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
          endif
          nswsfci_cpl(i) = nnirbmi_cpl(i) + nnirdfi_cpl(i)              &
     &                   + nvisbmi_cpl(i) + nvisdfi_cpl(i)
          nswsfc_cpl(i)  = nswsfc_cpl(i)  + nswsfci_cpl(i)*dtf
          nnirbm_cpl(i)  = nnirbm_cpl(i)  + nnirbmi_cpl(i)*dtf
          nnirdf_cpl(i)  = nnirdf_cpl(i)  + nnirdfi_cpl(i)*dtf
          nvisbm_cpl(i)  = nvisbm_cpl(i)  + nvisbmi_cpl(i)*dtf
          nvisdf_cpl(i)  = nvisdf_cpl(i)  + nvisdfi_cpl(i)*dtf
        enddo
      endif

      if (lssav) then
        do i = 1, im
          gflux(i)   = gflux(i)  + gflx(i)  * dtf
          evbsa(i)   = evbsa(i)  + evbs(i)  * dtf
          evcwa(i)   = evcwa(i)  + evcw(i)  * dtf
          transa(i)  = transa(i) + trans(i) * dtf
          sbsnoa(i)  = sbsnoa(i) + sbsno(i) * dtf
          snowca(i)  = snowca(i) + snowc(i) * dtf
          snohfa(i)  = snohfa(i) + snohf(i) * dtf
          ep(i)      = ep(i)     + ep1d(i)  * dtf

          tmpmax(i)  = max(tmpmax(i),t2m(i))
          tmpmin(i)  = min(tmpmin(i),t2m(i))

          spfhmax(i) = max(spfhmax(i),q2m(i))
          spfhmin(i) = min(spfhmin(i),q2m(i))
        enddo
      endif

!!!!!!!!!!!!!!!!!Commented by Moorthi on July 18, 2012 !!!!!!!!!!!!!!!!!!!
!     do i = 1, im
!  --- ...  compute coefficient of evaporation in evapc
!
!       if (evapc(i) > 1.0e0) evapc(i) = 1.0e0
!  --- ...  over snow cover or ice or sea, coef of evap =1.0e0
!       if (weasd(i) > 0.0 .or. slmsk(i) /= 1.0) evapc(i) = 1.0e0
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  --- ...  Boundary Layer and Free atmospheic turbulence parameterization

!     if (lprnt) write(0,*)' tsea3=',tsea(ipr),' slmsk=',slmsk(ipr)     &
!    &, ' kdt=',kdt,' evap=',evap(ipr)
!     if (lprnt)  write(0,*)' dtdtb=',dtdt(ipr,1:10)

!     do i = 1, im
!       if (islmsk(i) == 0) then
!         oro_land(i) = 0.0
!       else
!         oro_land(i) = oro(i)
!       endif
!     enddo

!     write(0,*)' before monin clstp=',clstp,' kdt=',kdt,' lat=',lat

      if (do_shoc) then
        call moninshoc(ix,im,levs,ntrac,ntcw,dvdt,dudt,dtdt,dqdt,       &
     &                 ugrs,vgrs,tgrs,qgrs,phy_f3d(1,1,ntot3d-1),       &  ! tkh
     &                 prnum,ntke,                                      &
     &                 prsik(1,1),rb,zorl,u10m,v10m,ffmm,ffhh,          &
     &                 tsea,hflx,evap,stress,wind,kpbl,                 &
     &                 prsi,del,prsl,prslk,phii,phil,dtp,               &
     &                 dusfc1,dvsfc1,dtsfc1,dqsfc1,dkt,hpbl,            &
     &                 kinver, xkzm_m, xkzm_h, xkzm_s, lprnt, ipr,me)
      else

        if (hybedmf) then

          call moninedmf(ix,im,levs,nvdiff,ntcw,dvdt,dudt,dtdt,dqdt,    &
     &       ugrs,vgrs,tgrs,qgrs,swh,hlw,xmu,                           &
     &       prsik(1,1),rb,zorl,u10m,v10m,ffmm,ffhh,                    &
     &       tsea,qss,hflx,evap,stress,wind,kpbl,                       &
     &       prsi,del,prsl,prslk,phii,phil,dtp,dspheat,                 &
     &       dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,dkt,            &
     &       kinver, xkzm_m, xkzm_h, xkzm_s, lprnt, ipr)

        elseif (.not. old_monin) then

          call moninq(ix,im,levs,nvdiff,ntcw,dvdt,dudt,dtdt,dqdt,       &
     &       ugrs,vgrs,tgrs,qgrs,swh,hlw,xmu,                           &
     &       prsik(1,1),rb,ffmm,ffhh,                                   &
     &       tsea,qss,hflx,evap,stress,wind,kpbl,                       &
     &       prsi,del,prsl,prslk,phii,phil,dtp,dspheat,                 &
     &       dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,dkt,            &
     &       kinver, xkzm_m, xkzm_h, xkzm_s, lprnt, ipr)

        else

          if (mstrat) then
            call moninp1(ix,im,levs,nvdiff,dvdt,dudt,dtdt,dqdt,         &
     &       ugrs,vgrs,tgrs,qgrs,                                       &
     &       prsik(1,1),rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,    &
     &       kpbl,prsi,del,prsl,prslk,phii,phil,dtp,                    &
     &       dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,dkt,            &
     &       kinver, xkzm_m, xkzm_h)
!    &       kinver, oro_land, ctei_r, ctei_rm, xkzm_m, xkzm_h)
          else
            call moninp(ix,im,levs,nvdiff,dvdt,dudt,dtdt,dqdt,          &
     &       ugrs,vgrs,tgrs,qgrs,                                       &
     &       prsik(1,1),rb,ffmm,ffhh,tsea,qss,hflx,evap,stress,wind,    &
     &       kpbl,prsi,del,prsl,phii,phil,dtp,                          &
     &       dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl,gamt,gamq,dkt,            &
     &       xkzm_m,xkzm_h)
          endif

        endif   ! end if_hybedmf
      endif   ! end if_do_shoc

      if (lssav_cpl) then
        do i = 1, im
          if (flag_cice(i)) then
             cice(i)   = fice_cice(i)
             tsea(i)   = tsea_cice(i)
             dusfc1(i) = dusfc_cice(i)
             dvsfc1(i) = dvsfc_cice(i)
             dqsfc1(i) = dqsfc_cice(i)
             dtsfc1(i) = dtsfc_cice(i)
          endif
        enddo
      endif

!     if (lprnt) then
!       write(0,*) ' dusfc1=',dusfc1(ipr),' kdt=',kdt,' lat=',lat
!       write(0,*)' dtsfc1=',dtsfc1(ipr)
!       write(0,*)' dqsfc1=',dqsfc1(ipr)
!       write(0,*)' dtdt=',dtdt(ipr,1:10)
!       print *,' dudtm=',dudt(ipr,:)
!     endif

!  --- ...  coupling insertion

      if (lssav_cpl) then
        do i=1, im
          dusfc_cpl(i)  = dusfc_cpl(i) + dusfc1(i)*dtf
          dvsfc_cpl(i)  = dvsfc_cpl(i) + dvsfc1(i)*dtf
          dtsfc_cpl(i)  = dtsfc_cpl(i) + dtsfc1(i)*dtf
          dqsfc_cpl(i)  = dqsfc_cpl(i) + dqsfc1(i)*dtf
          dusfci_cpl(i) = dusfc1(i)
          dvsfci_cpl(i) = dvsfc1(i)
          dtsfci_cpl(i) = dtsfc1(i)
          dqsfci_cpl(i) = dqsfc1(i)
        enddo
      endif
!-------------------------------------------------------lssav if loop ----------
      if (lssav) then
        do i = 1, im
          dusfc(i)  = dusfc(i) + dusfc1(i)*dtf
          dvsfc(i)  = dvsfc(i) + dvsfc1(i)*dtf
          dtsfc(i)  = dtsfc(i) + dtsfc1(i)*dtf
          dqsfc(i)  = dqsfc(i) + dqsfc1(i)*dtf
          dusfci(i) = dusfc1(i)
          dvsfci(i) = dvsfc1(i)
          dtsfci(i) = dtsfc1(i)
          dqsfci(i) = dqsfc1(i)
        enddo
!       if (lprnt) then
!         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
!    &     dtf,' kdt=',kdt,' lat=',lat
!       endif

        if (ldiag3d) then

          if (lsidea) then
            do k = 1, levs
              do i = 1, im
                dt3dt(i,k,3) = dt3dt(i,k,3) + dtdt(i,k)*dtf
              enddo
            enddo
          else
            do k = 1, levs
              do i = 1, im
                tem          = dtdt(i,k) - (hlw(i,k)+swh(i,k)*xmu(i))
                dt3dt(i,k,3) = dt3dt(i,k,3) + tem*dtf
              enddo
            enddo
          endif
          do k = 1, levs
            do i = 1, im
              du3dt(i,k,1) = du3dt(i,k,1) + dudt(i,k) * dtf
              du3dt(i,k,2) = du3dt(i,k,2) - dudt(i,k) * dtf
              dv3dt(i,k,1) = dv3dt(i,k,1) + dvdt(i,k) * dtf
              dv3dt(i,k,2) = dv3dt(i,k,2) - dvdt(i,k) * dtf
            enddo
          enddo
! update dqdt_v to include moisture tendency due to vertical diffusion
!         if (lgocart) then
!           do k = 1, levs
!             do i = 1, im
!               dqdt_v(i,k)  = dqdt(i,k,1) * dtf
!             enddo
!           enddo
!         endif
          do k = 1, levs
            do i = 1, im
              tem  = dqdt(i,k,1) * dtf
              dq3dt(i,k,1) = dq3dt(i,k,1) + tem
!             dqdt_v(i,k)  = tem
            enddo
          enddo
          if (ntoz > 0) then
            do k = 1, levs
              do i = 1, im
                dq3dt(i,k,5) = dq3dt(i,k,5) + dqdt(i,k,ntoz) * dtf
              enddo
            enddo
          endif
        endif

      endif   ! end if_lssav
!-------------------------------------------------------lssav if loop ----------
!
!            Orographic gravity wave drag parameterization
!            ---------------------------------------------

      if (nmtvr == 14) then         ! current operational - as of 2014

        do i = 1, im
          oc(i) = hprime(i,2)
        enddo
        do k = 1, 4
          do i = 1, im
            oa4(i,k) = hprime(i,k+2)
            clx(i,k) = hprime(i,k+6)
          enddo
        enddo
        do i = 1, im
          theta(i)  = hprime(i,11)
          gamma(i)  = hprime(i,12)
          sigma(i)  = hprime(i,13)
          elvmax(i) = hprime(i,14)
        enddo

      elseif (nmtvr == 10) then

        do i = 1, im
          oc(i) = hprime(i,2)
        enddo
        do k = 1, 4
          do i = 1, im
            oa4(i,k) = hprime(i,k+2)
            clx(i,k) = hprime(i,k+6)
          enddo
        enddo

      elseif (nmtvr == 6) then

        do i = 1, im
          oc(i) = hprime(i,2)
        enddo
        do k = 1, 4
          do i = 1, im
            oa4(i,k) = hprime(i,k+2)
            clx(i,k) = 0.0
          enddo
        enddo

      else

        oc = 0 ; oa4 = 0 ; clx = 0 ; theta = 0 ; gamma = 0 ; sigma = 0
        elvmax = 0

      endif   ! end if_nmtvr

!     write(0,*)' before gwd clstp=',clstp,' kdt=',kdt,' lat=',lat
      call gwdps(im, ix, im,   levs,  dvdt, dudt, dtdt,                 &
     &           ugrs,   vgrs, tgrs,  qgrs,                             &
     &           kpbl,   prsi, del,   prsl, prslk,                      &
     &           phii,   phil, dtp,                                     &
     &           kdt,    hprime(1,1), oc, oa4, clx,                     &
     &           theta,sigma,gamma,elvmax,dusfcg, dvsfcg,               &
     &           con_g,con_cp,con_rd,con_rv, lonr, nmtvr, cdmbgwd,      &
     &           me, lprnt,ipr)

!     if (lprnt)  print *,' dudtg=',dudt(ipr,:)

      if (lssav) then
        do i = 1, im
          dugwd(i) = dugwd(i) + dusfcg(i)*dtf
          dvgwd(i) = dvgwd(i) + dvsfcg(i)*dtf
        enddo

!       if (lprnt) print *,' dugwd=',dugwd(ipr),' dusfcg=',dusfcg(ipr)
!       if (lprnt) print *,' dvgwd=',dvgwd(ipr),' dvsfcg=',dvsfcg(ipr)

        if (ldiag3d) then
          do k = 1, levs
            do i = 1, im
              du3dt(i,k,2) = du3dt(i,k,2) + dudt(i,k) * dtf
              dv3dt(i,k,2) = dv3dt(i,k,2) + dvdt(i,k) * dtf
              dt3dt(i,k,2) = dt3dt(i,k,2) + dtdt(i,k) * dtf
            enddo
          enddo
        endif
      endif

!    Rayleigh damping  near the model top
      if( .not. lsidea .and. ral_ts > 0.0) then
!        call rayleigh_damp_mesopause(im, ix, im, levs, dvdt, dudt, dtdt,
!     &                   ugrs, vgrs, dtp, con_cp, levr, prsl, prslrd0)
!      else
        call rayleigh_damp(im, ix, im, levs, dvdt, dudt, dtdt, ugrs,
     &                     vgrs, dtp, con_cp, levr, pgr, prsl,
     &                     prslrd0, ral_ts)
      endif

!     if (lprnt) then
!       write(0,*)' tgrs1=',tgrs(ipr,1:10)
!       write(0,*)' dtdt=',dtdt(ipr,1:10)
!     endif

      do  k = 1, levs
        do i = 1, im
          gt0(i,k) = tgrs(i,k) + dtdt(i,k) * dtp
          gu0(i,k) = ugrs(i,k) + dudt(i,k) * dtp
          gv0(i,k) = vgrs(i,k) + dvdt(i,k) * dtp
        enddo
      enddo

      do n = 1, ntrac
        do k = 1, levs
          do i = 1, im
            gq0(i,k,n) = qgrs(i,k,n) + dqdt(i,k,n) * dtp
          enddo
        enddo
      enddo
!     if(lprnt) write(0,*)' gq0i=',gq0(ipr,:,ntiw)

      if( lsidea ) then            ! idea convective adjustment
        call ideaca_up(prsi,gt0,ix,im,levs+1)
      endif

!  --- ...  check print

!     if (me == 0) then
!       sumq = 0.0
!       do k = 1, levs
!         do i = 1, im
!           sumq = sumq + (dqdt(i,k,1)+dqdt(i,k,ntcw)) * del(i,k)
!         enddo
!       enddo

!       sume = 0.0
!       do i = 1, im
!         sume = sume + dqsfc1(i)
!       enddo

!       sumq = sumq * 1000.0 / con_g
!       sume = sume / con_hvap
!       print *,' after mon: sumq=',sumq,' sume=',sume, ' kdt=',kdt
!     endif

!  --- ...  ozone physics

      if (ntoz > 0 .and. ntrac >= ntoz) then

        if (pl_coeff > 4) then

         call ozphys_2015(ix,im,levs,ko3,dtp,gq0(1,1,ntoz),gq0(1,1,ntoz)&
     &,                   gt0, poz, prsl, prdoz, pl_coeff, del, ldiag3d &
     &,                   dq3dt(1,1,6), me)
        else

         call ozphys(ix,im,levs,ko3,dtp,gq0(1,1,ntoz),gq0(1,1,ntoz)     &
     &,              gt0, poz, prsl, prdoz, pl_coeff, del, ldiag3d      &
     &,              dq3dt(1,1,6), me)

        endif
      endif


      if (h2o_phys) then

        call h2ophys(ix,im,levs,levh2o,dtp,gq0(1,1,1),gq0(1,1,1)
     &,              h2o_pres,prsl,h2opl,h2o_coeff,ldiag3d
     &,              dq3dt(1,1,1), me)
      endif

!  --- ...  to side-step the ozone physics

!      if (ntrac >= 2) then
!        do k = 1, levs
!          gq0(k,ntoz) = qgrs(k,ntoz)
!        enddo
!      endif

!     if (lprnt) then
!       write(0,*) ' levs=',levs,' jcap=',jcap,' dtp',dtp               &
!    &,  ' slmsk=',slmsk(ilon,ilat),' kdt=',kdt
!       print *,' rann=',rann,' ncld=',ncld,' iq=',iq,' lat=',lat
!       print *,' pgr=',pgr
!       print *,' del=',del(ipr,:)
!       print *,' prsl=',prsl(ipr,:)
!       print *,' prslk=',prslk(ipr,:)
!       print *,' rann=',rann(ipr,1)
!       write(0,*)' gt0=',gt0(ipr,:)                                    &
!    &,         ' kdt=',kdt,' xlon=',xlon(ipr),' xlat=',xlat(ipr)
!       print *,' dtdt=',dtdt(ipr,:)
!       print *,' gu0=',gu0(ipr,:)
!       print *,' gv0=',gv0(ipr,:)
!       write(0,*)' gq0=',(gq0(ipr,k,1),k=1,levs),' lat=',lat
!       write(0,*)' gq0i2=',(gq0(ipr,k,ntiw),k=1,levs),' lat=',lat
!       print *,' gq1=',(gq0(ipr,k,ntcw),k=1,levs)
!       print *,' vvel=',vvel
!     endif

      if (ldiag3d) then

        do k = 1, levs
          do i = 1, im
            dtdt(i,k)   = gt0(i,k)
!           dqdt(i,k,1) = gq0(i,k,1)
            dudt(i,k)   = gu0(i,k)
            dvdt(i,k)   = gv0(i,k)
          enddo
        enddo

      elseif (cnvgwd) then

        do k = 1, levs
          do i = 1, im
            dtdt(i,k)   = gt0(i,k)
          enddo
        enddo

      endif   ! end if_ldiag3d/cnvgwd

      if (ldiag3d .or. lgocart) then
        do k = 1, levs
          do i = 1, im
            dqdt(i,k,1)   = gq0(i,k,1)
          enddo
        enddo
      endif   ! end if_ldiag3d/lgocart

      call get_phi(im,ix,levs,ntrac,gt0,gq0,                            &
     &             thermodyn_id, sfcpress_id,                           &
     &             gen_coord_hybrid,                                    &
     &             prsi,prsik,prsl,prslk,phii,phil)

!     if (lprnt) then
!       print *,' phii2=',phii(ipr,:)
!       print *,' phil2=',phil(ipr,:)
!     endif

      do k = 1, levs
        do i = 1, im
          clw(i,k,1) = 0.0
          clw(i,k,2) = -999.9
        enddo
      enddo
      if (.not. ras .or. .not. cscnv) then
        do k = 1, levs
          do i = 1, im
            cnvc(i,k)  = 0.0
            cnvw(i,k)  = 0.0
          enddo
        enddo
      endif

!     write(0,*)' before cnv clstp=',clstp,' kdt=',kdt,' lat=',lat

!  --- ...  for convective tracer transport (while using ras)

      if (ras .or. cscnv) then
        if (tottracer > 0) then

          if (ntoz > 0) then
            do k=1,levs
              do i=1,im
                clw(i,k,3) = gq0(i,k,ntoz)
              enddo
            enddo

            if (tracers > 0) then
              do n = 1, tracers
                do k=1,levs
                  do i=1,im
                    clw(i,k,3+n) = gq0(i,k,n+trc_shft)
                  enddo
                enddo
              enddo
            endif
          else
            do n = 1, tracers
                do k=1,levs
                  do i=1,im
                    clw(i,k,2+n) = gq0(i,k,n+trc_shft)
                  enddo
                enddo
            enddo
          endif

        endif
      endif   ! end if_ras or cfscnv

      do i = 1, im
        ktop(i)  = 1
        kbot(i)  = levs
      enddo

!  --- ...  calling condensation/precipitation processes
!           --------------------------------------------

      if (ntcw > 0) then

        do k = 1, levs
          do i = 1, im
            tem      = rhbbot - (rhbbot-rhbtop) * (1.0-prslk(i,k))
            tem      = rhc_max * work1(i) + tem * work2(i)
            rhc(i,k) = max(0.0, min(1.0,tem))
          enddo
        enddo

        if (ncld == 2) then
          do k = 1, levs
            do i = 1, im
              clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
              clw(i,k,2) = gq0(i,k,ntcw)                    ! water
            enddo
          enddo
        else

          if (num_p3d == 3) then    !  brad ferrier's microphysics

!  --- ...  algorithm to separate different hydrometeor species

            do k = 1, levs
              do i = 1, im
                wc     = gq0(i,k,ntcw)
                qi     = 0.
                qr     = 0.
                qw     = 0.
                f_ice  = max(0.0, min(1.0, phy_f3d(i,k,1)))
                f_rain = max(0.0, min(1.0, phy_f3d(i,k,2)))

                qi = f_ice*wc
                qw = wc-qi
                if (qw > 0.0) then
                  qr = f_rain*qw
                  qw = qw-qr
                endif

!               if (f_ice >= 1.0) then
!                 qi = wc
!               elseif (f_ice <= 0.0) then
!                 qw = wc
!               else
!                 qi = f_ice*wc
!                 qw = wc-qi
!               endif

!               if (qw > 0.0 .and. f_rain > 0.0) then
!                 if (f_rain >= 1.0) then
!                   qr = qw
!                   qw = 0.0
!                 else
!                   qr = f_rain*qw
!                   qw = qw-qr
!                 endif
!               endif

                qr_col(i,k) = qr
!               clw(i,k)    = qi + qw
                clw(i,k,1)  = qi
                clw(i,k,2)  = qw

!  --- ...  array to track fraction of "cloud" in the form of ice

!               if (qi+qw > epsq) then
!                 fc_ice(i,k) = qi / (qi+qw)
!               else
!                 fc_ice(i,k) = 0.0
!               endif
              enddo
            enddo

          elseif (num_p3d == 4) then   ! zhao-carr microphysics

            do i = 1, im
              psautco_l(i) = psautco(1)*work1(i) + psautco(2)*work2(i)
              prautco_l(i) = prautco(1)*work1(i) + prautco(2)*work2(i)
            enddo
            do k = 1, levs
              do i = 1, im
                clw(i,k,1) = gq0(i,k,ntcw)
              enddo
            enddo

          endif  ! end if_num_p3d
        endif    ! end if (ncld == 2)

      else    ! if_ntcw

        do i = 1, im
          psautco_l(i) = psautco(1)*work1(i) + psautco(2)*work2(i)
          prautco_l(i) = prautco(1)*work1(i) + prautco(2)*work2(i)
        enddo
        do k = 1, levs
          do i = 1, im
            rhc(i,k) = 1.0
          enddo
        enddo

      endif   ! end if_ntcw
!
!        Call SHOC if do_shoc is true and shocaftcnv is false
!
      if (do_shoc .and. .not. shocaftcnv) then

        if (ncld == 2) then
          skip_macro = do_shoc
          do k = 1, levs
            do i = 1, im
              clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
              clw(i,k,2) = gq0(i,k,ntcw)                    ! water
              ncpl(i,k)  = gq0(i,k,ntlnc)
              ncpi(i,k)  = gq0(i,k,ntinc)
            enddo
          enddo
        elseif (num_p3d == 4) then
          do k=1,levs
            do i=1,im
              qpl(i,k)   = 0.0
              qpi(i,k)   = 0.0
              tem = gq0(i,k,ntcw)                                       &
     &            * max(0.0, MIN(1.0, (TCR-gt0(i,k))*TCRF))
              clw(i,k,1) = tem                              ! ice
              clw(i,k,2) = gq0(i,k,ntcw) - tem              ! water
            enddo
          enddo
        endif

!       dtshoc = 60.0
!       dtshoc = 120.0
!       dtshoc = dtp
!       nshocm = (dtp/dtshoc) + 0.001
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!      if (lprnt) write(1000+me,*)' before shoc tke=',clw(ipr,:,ntk),
!    &' kdt=',kdt,' lat=',lat,'xlon=',xlon(ipr),' xlat=',xlat(ipr)

!     phy_f3d(1,1,ntot3d-2) - shoc determined sgs clouds
!     phy_f3d(1,1,ntot3d-1) - shoc determined diffusion coefficients
!     phy_f3d(1,1,ntot3d  ) - shoc determined  w'theta'
!
!         call shoc(ix, im, 1, levs, levs+1, dtshoc, me, lat,           &
          call shoc(ix, im, 1, levs, levs+1, dtp, me, lat,              &
     &              prsl(1,1), phii(1,1), phil(1,1),                    &
     &              gu0(1,1),gv0(1,1), vvel(1,1), gt0(1,1), gq0(1,1,1), &
     &              clw(1,1,1), clw(1,1,2), qpi, qpl,rhc, sup,          &
     &              phy_f3d(1,1,ntot3d-2), clw(1,1,ntk), hflx, evap,    &
     &              prnum, phy_f3d(1,1,ntot3d-1), phy_f3d(1,1,ntot3d),  &
     &              lprnt, ipr, ncpl, ncpi)

          if (ntlnc > 0 .and. ntinc > 0 .and. ncld >=2) then
            do k=1,levs
              do i=1,im
                gq0(i,k,ntlnc) = ncpl(i,k)
                gq0(i,k,ntinc) = ncpi(i,k)
              enddo
            enddo
          endif
!       do k=1,levs
!         do i=1,im
!           sgs_cld(i,k) = sgs_cld(i,k) + shoc_cld(i,k)
!         enddo
!       enddo
!     if (lprnt) write(0,*)' gt03=',gt0(ipr,1:10)
!     if (lprnt) write(0,*)' tke=',clw(ipr,1:10,ntk)

!      if (lprnt) write(1000+me,*)' after shoc tke=',clw(1,:,ntk),
!    &' kdt=',kdt
!       enddo
!
!      do k=1,levs
!      write(1000+me,*)' maxcld=',maxval(sgs_cld(1:im,k)),
!      write(1000+me,*)' maxtkh=',maxval(phy_f3d(1:im,k,ntot3d-1)),
!    &' k=',k,' kdt=',kdt,' lat=',lat
!      enddo

!     write(0,*)' aft shoc gt0=',gt0(1,:),' lat=',lat
!     write(0,*)' aft shoc gq0=',gq0(1,:,1),' lat=',lat
!     write(0,*)' aft shoc gu0=',gu0(1,:),' lat=',lat
!
      endif   ! if(do_shoc)


!  --- ...  calling convective parameterization
!
      if (.not. ras .and. .not. cscnv) then

        if (imfdeepcnv == 1) then             ! no random cloud top
          call sascnvn(im,ix,levs,jcap,dtp,del,prsl,pgr,phil,           &
     &                clw,gq0,gt0,gu0,gv0,cld1d,                        &
     &                rain1,kbot,ktop,kcnv,islmsk,                      &
     &                vvel,ncld,ud_mf,dd_mf,dt_mf,cnvw,cnvc)
        elseif (imfdeepcnv == 2) then
          call mfdeepcnv(im,ix,levs,dtp,del,prsl,pgr,phil,              &
     &                clw,gq0,gt0,gu0,gv0,cld1d,                        &
     &                rain1,kbot,ktop,kcnv,islmsk,garea,                &
     &                vvel,ncld,ud_mf,dd_mf,dt_mf,cnvw,cnvc)
!         if (lprnt) print *,' rain1=',rain1(ipr)
        elseif (imfdeepcnv == 0) then         ! random cloud top
          call sascnv(im,ix,levs,jcap,dtp,del,prsl,pgr,phil,            &
     &                clw,gq0,gt0,gu0,gv0,cld1d,                        &
     &                rain1,kbot,ktop,kcnv,islmsk,                      &
     &                vvel,rann,ncld,ud_mf,dd_mf,dt_mf,cnvw,cnvc)
!         if (lprnt) print *,' rain1=',rain1(ipr),' rann=',rann(ipr,1)
        endif

      else        ! ras or cscnv

        if(cscnv) then    ! Chikira-Sugiyama  convection scheme (via CSU)

!         fscav(:) = 0.0
!         fswtr(:) = 0.0
          call cs_convr(                                                & !DD
     &                  ix      ,im      ,levs    , tottracer+3 ,       & !DD
     &                  gt0     ,gq0     ,rain1   , clw       ,         & !DD
     &                  phil    ,phii    ,                              & !DD
     &                  prsl    ,prsi    ,                              & !DD
     &                  dtp     ,dtf     ,                              & !DD
     &                  ud_mf   ,dd_mf   ,dt_mf   ,                     & !DD
     &                  gu0     ,gv0     ,fscav, fswtr,                 & !DD
!    &                  phy_fctd, me         )                            !DD & moorthi
     &                  phy_fctd, me, wcbmax )                            !DD & moorthi
!    &                  phy_fctd     )                                    !DD
          do i = 1,im                                                     !DD
             rain1(i) = rain1(i) * (dtp*0.001)                            !DD
          enddo

        else      ! ras version 2

          if (ccwf(1) >= 0.0 .or. ccwf(2) >= 0 ) then
            do i=1,im
              ccwfac(i) = ccwf(1)*work1(i) + ccwf(2)*work2(i)
              dlqfac(i) = dlqf(1)*work1(i) + dlqf(2)*work2(i)
              lmh(i)    = levs
            enddo
          else
            do i=1,im
              ccwfac(i) = -999.0
              dlqfac(i) = 0.0
              lmh(i)    = levs
            enddo
          endif
!         if  (lprnt) write(0,*) ' calling ras for kdt=',kdt,' me=',me    &
!    &,                       ' lprnt=',lprnt,' ccwfac=',ccwfac(ipr)

          revap = .true.
!         if (ncld ==2) revap = .false.
          call rascnv(im,    ix,    levs,   dtp, dtf, rann              &
     &,               gt0,    gq0,   gu0,    gv0, clw, tottracer, fscav &
     &,               prsi,   prsl,   prsik,  prslk, phil,  phii        &
     &,               kpbl,   cd,     rain1,  kbot,  ktop,  kcnv        &
     &,               phy_f2d(1,num_p2d), flipv, pa2mb                  &
     &,               me, garea, lmh, ccwfac, nrcm, rhc                 &
     &,               ud_mf, dd_mf, dt_mf, dlqfac, lprnt, ipr, kdt,revap&
     &,               QLCN, QICN, w_upi,cf_upi, CNV_MFD, CNV_PRC3       &
     &,               CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,ncld )
        endif

!     if(lprnt) write(0,*)' after ras rain1=',rain1(ipr)
!    &,' cnv_prc3sum=',sum(cnv_prc3(ipr,1:levs))
!     if (lprnt) write(0,*)' gt04=',gt0(ipr,1:10)

        cld1d = 0

        if (lgocart) then
          do k = 1, levs
            do i = 1, im
              upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * frain
              dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * frain
              det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * frain
              cnvqc_v(i,k) = cnvqc_v(i,k) + (clw(i,k,1)+clw(i,k,2)-     &
     &                                         gq0(i,k,ntcw)) * frain
            enddo
          enddo
        endif ! if (lgocart)

!  --- ...  update the tracers due to convective transport

        if (tottracer > 0) then
          if (ntoz > 0) then                         ! for ozone
            do k=1,levs
              do i=1,im
                gq0(i,k,ntoz) = clw(i,k,3)
              enddo
            enddo

            if (tracers > 0) then                    ! for other tracers
              do n = 1, tracers
                do k=1,levs
                  do i=1,im
                    gq0(i,k,n+trc_shft) = clw(i,k,3+n)
                  enddo
                enddo
              enddo
            endif
          else
            do n = 1, tracers
              do k=1,levs
                do i=1,im
                  gq0(i,k,n+trc_shft) = clw(i,k,2+n)
                enddo
              enddo
            enddo
          endif
        endif
      endif   ! end if_not_ras
!
      do i = 1, im
        rainc(i) = frain * rain1(i)
      enddo

!
      if (lssav) then
        do i = 1, im
          cldwrk(i)  = cldwrk(i)  + cld1d(i) * dtf
          cnvprcp(i) = cnvprcp(i) + rainc(i)
        enddo

        if (ldiag3d) then
          do k = 1, levs
            do i = 1, im
              dt3dt(i,k,4) = dt3dt(i,k,4) + (gt0(i,k)-dtdt(i,k)) * frain
              dq3dt(i,k,2) = dq3dt(i,k,2) + (gq0(i,k,1)-dqdt(i,k,1))    &
     &                                                           * frain
              du3dt(i,k,3) = du3dt(i,k,3) + (gu0(i,k)-dudt(i,k)) * frain
              dv3dt(i,k,3) = dv3dt(i,k,3) + (gv0(i,k)-dvdt(i,k)) * frain

              upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * (con_g*frain)
              dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * (con_g*frain)
              det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * (con_g*frain)

            enddo
          enddo
        endif ! if (ldiag3d)

      endif   ! end if_lssav
!
!       update dqdt_v to include moisture tendency due to deep convection
      if (lgocart) then
        do k = 1, levs
          do i = 1, im
            dqdt_v(i,k)  = (gq0(i,k,1)-dqdt(i,k,1))  * frain
            upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * frain
            dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * frain
            det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * frain
            cnvqc_v(i,k) = cnvqc_v(i,k) + (clw(i,k,1)+clw(i,k,2))
     &                                                      *frain
          enddo
        enddo
      endif ! if (lgocart)
!
      if(npdf3d == 3  .and. num_p3d == 4) then
        num2 = num_p3d + 2
        num3 = num2 + 1
        do k = 1, levs
          do i = 1, im
            phy_f3d(i,k,num2) = cnvw(i,k)
            phy_f3d(i,k,num3) = cnvc(i,k)
          enddo
        enddo
      else if(npdf3d == 0  .and. ncnvcld3d == 1) then
        num2 = num_p3d + 1
        do k = 1, levs
          do i = 1, im
            phy_f3d(i,k,num2) = cnvw(i,k)
          enddo
        enddo
      endif

!
!----------------Convective gravity wave drag parameterization starting --------

      if (cnvgwd) then         !        call convective gravity wave drag

!  --- ...  calculate maximum convective heating rate            qmax [k/s]
!           cuhr = temperature change due to deep convection

        do i = 1, im
!         qmax(i)   = 0.
          cumabs(i) = 0.0
          work3(i)  = 0.0
        enddo

        do k = 1, levs
          do i = 1, im
!           cuhr(i,k) = (gt0(i,k)-dtdt(i,k)) / dtf
!           cuhr(i,k) = (gt0(i,k)-dtdt(i,k)) / dtp    ! moorthi

!           cumchr(i,k) = 0.
!           gwdcu(i,k)  = 0.
!           gwdcv(i,k)  = 0.
!           diagn1(i,k) = 0.
!           diagn2(i,k) = 0.

            if (k >= kbot(i) .and. k <= ktop(i)) then
!             qmax(i)   = max(qmax(i),cuhr(i,k))
!             cumabs(i) = cuhr(i,k) + cumabs(i)
              cumabs(i) = cumabs(i) + (gt0(i,k)-dtdt(i,k)) * del(i,k)
              work3(i)  = work3(i)  + del(i,k)
            endif
          enddo
        enddo
        do i=1,im
          if (work3(i) > 0.0) cumabs(i) = cumabs(i) / (dtp*work3(i))
        enddo

!       do i = 1, im
!         do k = kbot(i), ktop(i)
!           do k1 = kbot(i), k
!             cumchr(i,k) = cuhr(i,k1) + cumchr(i,k)
!           enddo
!           cumchr(i,k) = cumchr(i,k) / cumabs(i)
!         enddo
!       enddo

!  --- ...  begin check print ******************************************

!       if (lprnt) then
!         if (kbot(ipr) <= ktop(ipr)) then
!           write(*,*) 'kbot <= ktop     for (lat,lon) = ',             &
!    &            xlon(ipr)*57.29578,xlat(ipr)*57.29578
!           write(*,*) 'kcnv kbot ktop qmax dlength  ',kcnv(ipr),       &
!    &            kbot(ipr),ktop(ipr),(86400.*qmax(ipr)),dlength(ipr)
!           write(*,9000) kdt
!9000       format(/,3x,'k',5x,'cuhr(k)',4x,'cumchr(k)',5x,             &
!    &            'at kdt = ',i4,/)

!           do k = ktop(ipr), kbot(ipr),-1
!             write(*,9010) k,(86400.*cuhr(ipr,k)),(100.*cumchr(ipr,k))
!9010         format(2x,i2,2x,f8.2,5x,f6.0)
!           enddo
!         endif

!         if (fhour >= fhourpr) then
!           print *,' before gwdc in gbphys start print'
!           write(*,*) 'fhour ix im levs = ',fhour,ix,im,levs
!           print *,'dtp  dtf  = ',dtp,dtf

!           write(*,9100)
!9100       format(//,14x,'pressure levels',//                          &
!    &             ' ilev',7x,'prsi',8x,'prsl',8x,'delp',/)

!           k = levs + 1
!           write(*,9110) k,(10.*prsi(ipr,k))
!9110       format(i4,2x,f10.3)

!           do k = levs, 1, -1
!             write(*,9120) k,(10.*prsl(ipr,k)),(10.*del(ipr,k))
!             write(*,9110) k,(10.*prsi(ipr,k))
!           enddo
!9120       format(i4,12x,2(2x,f10.3))

!           write(*,9130)
!9130       format(//,10x,'before gwdc in gbphys',//,' ilev',6x,        &
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',8x,'dudt',8x,'dvdt',/)

!           do k = levs, 1, -1
!             write(*,9140) k,ugrs(ipr,k),gu0(ipr,k),                   &
!    &                        vgrs(ipr,k),gv0(ipr,k),                   &
!    &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
!    &                        dudt(ipr,k),dvdt(ipr,k)
!           enddo
!9140       format(i4,9(2x,f10.3))

!           print *,' before gwdc in gbphys end print'
!         endif
!       endif   ! end if_lprnt

!  --- ...  end check print ********************************************

        call gwdc(im, ix, im, levs, lat, ugrs, vgrs, tgrs, qgrs,        &
     &            prsl, prsi, del, cumabs, ktop, kbot, kcnv,cldf,       &
     &            con_g, con_cp, con_rd, con_fvirt, dlength,            &
     &            lprnt, ipr, fhour, gwdcu, gwdcv, dusfcg, dvsfcg)


!       if (lprnt) then
!         if (fhour >= fhourpr) then
!           print *,' after gwdc in gbphys start print'

!           write(*,9131)
!9131       format(//,10x,'after gwdc in gbphys',//,' ilev',6x,         &
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

!           do k = levs, 1, -1
!             write(*,9141) k,ugrs(ipr,k),gu0(ipr,k),                   &
!    &                        vgrs(ipr,k),gv0(ipr,k),                   &
!    &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
!    &                        gwdcu(ipr,k),gwdcv(ipr,k)
!           enddo
!9141       format(i4,9(2x,f10.3))

!           print *,' after gwdc in gbphys end print'
!         endif
!       endif

!  --- ...  write out cloud top stress and wind tendencies

        if (lssav) then
          do i = 1, im
            dugwd(i) = dugwd(i) + dusfcg(i)*dtf
            dvgwd(i) = dvgwd(i) + dvsfcg(i)*dtf
          enddo

          if (ldiag3d) then
            do k = 1, levs
              do i = 1, im
                du3dt(i,k,4) = du3dt(i,k,4) + gwdcu(i,k)  * dtf
                dv3dt(i,k,4) = dv3dt(i,k,4) + gwdcv(i,k)  * dtf
!               du3dt(i,k,2) = du3dt(i,k,2) + diagn1(i,k) * dtf
!               dv3dt(i,k,2) = dv3dt(i,k,2) + diagn2(i,k) * dtf
              enddo
            enddo
          endif
        endif   ! end if_lssav

!  --- ...  update the wind components with  gwdc tendencies

        do k = 1, levs
          do i = 1, im
            eng0      = 0.5*(gu0(i,k)*gu0(i,k)+gv0(i,k)*gv0(i,k))
            gu0(i,k)  = gu0(i,k) + gwdcu(i,k) * dtp
            gv0(i,k)  = gv0(i,k) + gwdcv(i,k) * dtp
            eng1      = 0.5*(gu0(i,k)*gu0(i,k)+gv0(i,k)*gv0(i,k))
            gt0(i,k)  = gt0(i,k) + (eng0-eng1)/(dtp*con_cp)
          enddo
        enddo

!       if (lprnt) then
!         if (fhour >= fhourpr) then
!           print *,' after tendency gwdc in gbphys start print'

!           write(*,9132)
!9132       format(//,10x,'after tendency gwdc in gbphys',//,' ilev',6x,&
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

!           do k = levs, 1, -1
!             write(*,9142) k,ugrs(ipr,k),gu0(ipr,k),vgrs(ipr,k),       &
!    &              gv0(ipr,k),tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),      &
!    &              gwdcu(ipr,k),gwdcv(ipr,k)
!           enddo
!9142       format(i4,9(2x,f10.3))

!           print *,' after tendency gwdc in gbphys end print'
!         endif
!       endif

      endif   ! end if_cnvgwd (convective gravity wave drag)

!----------------Convective gravity wave drag parameterization over --------

      if (ldiag3d) then
        do k = 1, levs
          do i = 1, im
            dtdt(i,k)   = gt0(i,k)
!           dqdt(i,k,1) = gq0(i,k,1)
          enddo
        enddo
      endif
      if (ldiag3d .or. lgocart) then
        do k = 1, levs
          do i = 1, im
            dqdt(i,k,1) = gq0(i,k,1)
          enddo
        enddo
      endif

!     write(0,*)' before do_shoc shal clstp=',clstp,' kdt=',kdt,
!    &         ' lat=',lat

      if (.not. do_shoc) then

        if (shal_cnv) then               ! Shallow convection parameterizations
!                                        --------------------------------------
          if (imfshalcnv == 1) then      ! opr option now at 2014
                                         !-----------------------
            call shalcnv(im,ix,levs,jcap,dtp,del,prsl,pgr,phil,         &
     &                   clw,gq0,gt0,gu0,gv0,                           &
     &                   rain1,kbot,ktop,kcnv,islmsk,                   &
     &                   vvel,ncld,hpbl,hflx,evap,ud_mf,dt_mf,          &
     &                   cnvw,cnvc)

            if (shcnvcw .and. num_p3d == 4 .and. npdf3d == 3 ) then
              do k = 1, levs
                do i = 1, im
                  phy_f3d(i,k,num2) = cnvw(i,k)
                  phy_f3d(i,k,num3) = cnvc(i,k)
!???              phy_f3d(i,k,num2) = phy_f3d(i,k,num2) + cnvw(i,k)
!???              phy_f3d(i,k,num3) = phy_f3d(i,k,num3) + cnvc(i,k)
                enddo
              enddo
            else if(npdf3d == 0  .and. ncnvcld3d == 1) then
              num2 = num_p3d + 1
              do k = 1, levs
                do i = 1, im
                  phy_f3d(i,k,num2) = cnvw(i,k)
                enddo
              enddo
            endif
            do i = 1, im
              raincs(i) = frain    * rain1(i)
              rainc(i)  = rainc(i) + raincs(i)
            enddo
            if (lssav) then
              do i = 1, im
                cnvprcp(i) = cnvprcp(i) + raincs(i)
              enddo
            endif

          elseif (imfshalcnv == 2) then
            call mfshalcnv(im,ix,levs,dtp,del,prsl,pgr,phil,            &
     &                     clw,gq0,gt0,gu0,gv0,                         &
     &                     rain1,kbot,ktop,kcnv,islmsk,garea,           &
     &                     vvel,ncld,hpbl,ud_mf,dt_mf,cnvw,cnvc)

            if (shcnvcw .and. num_p3d == 4 .and. npdf3d == 3 ) then
              do k = 1, levs
                do i = 1, im
                  phy_f3d(i,k,num2) = cnvw(i,k)
                  phy_f3d(i,k,num3) = cnvc(i,k)
!???              phy_f3d(i,k,num2) = phy_f3d(i,k,num2) + cnvw(i,k)
!???              phy_f3d(i,k,num3) = phy_f3d(i,k,num3) + cnvc(i,k)
                enddo
              enddo
            else if(npdf3d == 0  .and. ncnvcld3d == 1) then
              num2 = num_p3d + 1
              do k = 1, levs
                do i = 1, im
                  phy_f3d(i,k,num2) = cnvw(i,k)
                enddo
              enddo
            endif
            do i = 1, im
              raincs(i) = frain    * rain1(i)
              rainc(i)  = rainc(i) + raincs(i)
            enddo
            if (lssav) then
              do i = 1, im
                cnvprcp(i) = cnvprcp(i) + raincs(i)
              enddo
            endif

          elseif (imfshalcnv == 0) then    ! modified Tiedtke Shallow convecton
                                         !-----------------------------------
            do i = 1, im
              levshc(i) = 0
            enddo
            do k = 2, levs
              do i = 1, im
                if (prsi(i,1)-prsi(i,k) <= dpshc(i)) levshc(i) = k
              enddo
            enddo
            levshcm = 1
            do i = 1, im
              levshcm = max(levshcm, levshc(i))
            enddo

!           if (lprnt) print *,' levshcm=',levshcm,' gt0sh=',gt0(ipr,:)
!    &,    ' lat=',lat

            if (mstrat) then             !  As in CFSv2
              call shalcv(im,ix,levshcm,dtp,del,prsi,prsl,prslk,kcnv,   &
     &                    gq0,gt0,levshc,phil,kinver,ctei_r,ctei_rml    &
     &,                                                    lprnt,ipr)
            else
              call shalcvt3(im,ix,levshcm,dtp,del,prsi,prsl,prslk,      &
     &                      kcnv,gq0,gt0)
            endif
!           if (lprnt) print *,' levshcm=',levshcm,' gt0sha=',gt0(ipr,:)

          endif   ! end if_imfshalcnv
        endif     ! end if_shal_cnv

        if (lssav) then
!          update dqdt_v to include moisture tendency due to shallow convection
          if (lgocart) then
            do k = 1, levs
              do i = 1, im
                tem         = (gq0(i,k,1)-dqdt(i,k,1)) * frain
                dqdt_v(i,k) = dqdt_v(i,k)  + tem
              enddo
            enddo
          endif
          if (ldiag3d) then
            do k = 1, levs
              do i = 1, im
                dt3dt(i,k,5) = dt3dt(i,k,5) + (gt0(i,k)-dtdt(i,k))
     &                                                          * frain
                dq3dt(i,k,3) = dq3dt(i,k,3) + (gq0(i,k,1)-dqdt(i,k,1))  &
     &                                                          * frain
                dtdt(i,k)   = gt0(i,k)
                dqdt(i,k,1) = gq0(i,k,1)
              enddo
            enddo
          endif
        endif   ! end if_lssav
!
        do k = 1, levs
          do i = 1, im
            if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
          enddo
        enddo

!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' befshgt0=',gt0(ipr,:)
!         write(0,*) ' befshgq0=',gq0(ipr,:,1)
!       endif

      elseif (shocaftcnv) then ! if do_shoc is true and shocaftcnv is true call shoc
        if (ncld == 2) then
          skip_macro = do_shoc
          do k = 1, levs
            do i = 1, im
!             clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
!             clw(i,k,2) = gq0(i,k,ntcw)                    ! water
              ncpl(i,k)  = gq0(i,k,ntlnc)
              ncpi(i,k)  = gq0(i,k,ntinc)
            enddo
          enddo

!       else
!         if (clw(1,1,2) < -999.0) then ! if clw is not partitioned to ice and water
!           do k=1,levs
!             do i=1,im
!               tem = gq0(i,k,ntcw)                                     &
!    &              * max(0.0, MIN(1.0, (TCR-gt0(i,k))*TCRF))
!               clw(i,k,1) = tem                              ! ice
!               clw(i,k,2) = gq0(i,k,ntcw) - tem              ! water
!             enddo
!           enddo
!         endif     ! Anning ncld ==2
        endif
        do k=1,levs
          do i=1,im
            qpl(i,k)   = 0.0
            qpi(i,k)   = 0.0
          enddo
        enddo
!       dtshoc = 60.0
!       nshocm = (dtp/dtshoc) + 0.001
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!       call shoc(im, 1, levs, levs+1, dtp, me, lat,        &
!!       call shoc(im, 1, levs, levs+1, dtshoc, me, lat, &
!    &                       prsl(1:im,:), phii (1:im,:),  phil(1:im,:),&
!    &          gu0(1:im,:),gv0(1:im,:), vvel(1:im,:), gt0(1:im,:),     &
!    &                                                   gq0(1:im,:,1), &
!    &          clw(1:im,:,1), clw(1:im,:,2), qpi, qpl,  sgs_cld(1:im,:)&
!    &,         gq0(1:im,:,ntke),                                       &
!    &          phy_f3d(1:im,:,ntot3d-1), phy_f3d(1:im,:,ntot3d),       &
!    &          lprnt, ipr,                                             &
!    &          con_cp, con_g, con_hvap, con_hfus, con_hvap+con_hfus,   &
!    &          con_rv, con_rd, con_pi, con_fvirt)

!       call shoc(ix, im, 1, levs, levs+1, dtshoc, me, lat,             &
        call shoc(ix, im, 1, levs, levs+1, dtp, me, lat,                &
     &            prsl(1,1), phii(1,1), phil(1,1),                      &
     &            gu0(1,1),gv0(1,1), vvel(1,1), gt0(1,1), gq0(1,1,1),   &
     &            clw(1,1,1), clw(1,1,2), qpi, qpl,rhc, sup,            &
     &            phy_f3d(1,1,ntot3d-2),  gq0(1,1,ntke),hflx,evap,      &
     &            prnum, phy_f3d(1,1,ntot3d-1), phy_f3d(1,1,ntot3d),    &
     &            lprnt, ipr, ncpl, ncpi)

        if (ntlnc > 0 .and. ntinc > 0 .and. ncld >=2) then
          do k=1,levs
            do i=1,im
              gq0(i,k,ntlnc) = ncpl(i,k)
              gq0(i,k,ntinc) = ncpi(i,k)
            enddo
          enddo
        endif

!
!      do k=1,levs
!      write(1000+me,*)' maxtkh=',maxval(phy_f3d(1:im,k,ntot3d-1)),
!    &' k=',k,' kdt=',kdt,' lat=',lat
!      enddo

!     write(0,*)' aft shoc gt0=',gt0(1,:),' lat=',lat
!     write(0,*)' aft shoc gq0=',gq0(1,:,1),' lat=',lat
!     write(0,*)' aft shoc gu0=',gu0(1,:),' lat=',lat
!
      endif   ! if( .not. do_shoc)
!
!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' aftshgt0=',gt0(ipr,:)
!         write(0,*) ' aftshgq0=',gq0(ipr,:,1)
!       endif

      if (ntcw > 0) then

!  for microphysics

        if (ncld == 2) then           ! morrison microphysics
          do k = 1, levs
            do i = 1, im
              gq0(i,k,ntiw) = clw(i,k,1)                     ! ice
              gq0(i,k,ntcw) = clw(i,k,2)                     ! water
            enddo
          enddo

        elseif (num_p3d == 3) then    ! call brad ferrier's microphysics

          do k = 1, levs
            do i = 1, im
!             qi = clw(i,k)*fc_ice(i,k)
!             qw = clw(i,k) - qi
              qi = clw(i,k,1)
              qw = clw(i,k,2)

!  --- ...  algorithm to combine different hydrometeor species

!             gq0(i,k,ntcw) = max(epsq, qi+qw+qr_col(i,k))
              gq0(i,k,ntcw) = qi + qw + qr_col(i,k)

              if (qi <= epsq) then
                phy_f3d(i,k,1) = 0.
              else
                phy_f3d(i,k,1) = qi/gq0(i,k,ntcw)
              endif

              if (qr_col(i,k) <= epsq) then
                phy_f3d(i,k,2) = 0.
              else
                phy_f3d(i,k,2) = qr_col(i,k) / (qw+qr_col(i,k))
              endif

            enddo
          enddo

        elseif (num_p3d == 4) then    ! if_num_p3d

          do k = 1, levs
            do i = 1, im
              gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
            enddo
          enddo

        endif   ! end if_num_p3d

      else    ! if_ntcw

        do k = 1, levs
          do i = 1, im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo

      endif   ! end if_ntcw

!  Legacy routine which determines convectve clouds - should be removed at some point

      call cnvc90(clstp, im,   ix,   rainc, kbot, ktop, levs, prsi,     &
     &            acv,   acvb, acvt, cv,    cvb,  cvt)


      if (moist_adj) then       ! moist convective adjustment
!                                 ---------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       if (me == 0) then
!         sumq = 0.0
!         DO K=1,LEVS
!           do i=1,im
!             sumq = sumq - (gq0(i,k,1)+gq0(i,k,ntcw)) * del(i,k)
!           enddo
!         enddo
!       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       To call moist convective adjustment
!
!       if (lprnt) then
!         print *,' prsl=',prsl(ipr,:)
!         print *,' del=',del(ipr,:)
!         print *,' gt0b=',gt0(ipr,:)
!         print *,' gq0b=',gq0(ipr,:,1)
!       endif

        call mstcnv(im,ix,levs,dtp,gt0,gq0,prsl,del,prslk,rain1
     &,                          gq0(1,1,ntcw), rhc, lprnt,ipr)

!       if (lprnt) then
!         print *,' rain1=',rain1(ipr),' rainc=',rainc(ipr)
!         print *,' gt0a=',gt0(ipr,:)
!         print *,' gq0a=',gq0(ipr,:,1)
!       endif
        do i=1,im
          rainc(i) = rainc(i) + frain * rain1(i)
        enddo
        if(lssav) then
          do i=1,im
            cnvprcp(i) = cnvprcp(i) + rain1(i) * frain
          enddo

! update dqdt_v to include moisture tendency due to surface processes
! dqdt_v : instaneous moisture tendency (kg/kg/sec)
!          if (lgocart) then
!            do k=1,levs
!              do i=1,im
!                tem = (gq0(i,k,1)-dqdt(i,k,1)) * frain
!                dqdt_v(i,k) = dqdt_v(i,k) + tem
!                dqdt_v(i,k) = dqdt_v(i,k) / dtf
!              enddo
!            enddo
!          endif
          if (ldiag3d) then
            do k=1,levs
              do i=1,im
                dt3dt(i,k,4) = dt3dt(i,k,4) + (gt0(i,k)-dtdt(i,k))
     &                                      * frain
                dq3dt(i,k,2) = dq3dt(i,k,2) + (gq0(i,k,1)-dqdt(i,k,1))
     &                                      * frain
              enddo
            enddo
          endif
         endif
!
        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dtdt(i,k)   = gt0(i,k)
              dqdt(i,k,1) = gq0(i,k,1)
            enddo
          enddo
        endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       if (me == 0) then
!         DO K=1,LEVS
!           do i=1,im
!             sumq = sumq + (gq0(i,k,1)+gq0(i,k,ntcw)) * del(i,k)
!           enddo
!         enddo
!         sumr = 0.0
!         do i=1,im
!           sumr = sumr + rain1(i)
!         enddo
!         sumq = sumq * 1000.0 / grav
!         sumr = sumr *1000
!         print *,' after MCN: sumq=',sumq,' sumr=',sumr, ' kdt=',kdt
!       endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      endif               !       moist convective adjustment over

! dqdt_v : instaneous moisture tendency (kg/kg/sec)
      if (lgocart) then
        do k=1,levs
          do i=1,im
            dqdt_v(i,k) = dqdt_v(i,k) / dtf
          enddo
        enddo
      endif
!
!     grid-scale condensation/precipitations and microphysics parameterization
!     ------------------------------------------------------------------------

      if (ncld == 0) then           ! no cloud microphysics

        call lrgscl(ix,im,levs,dtp,gt0,gq0,prsl,del,prslk,rain1,clw)

      elseif (ncld == 1) then       ! microphysics with single cloud condensate

        if (num_p3d == 3) then      ! brad ferrier's microphysics
                                    ! ---------------------------
          do i = 1, im
            xncw(i)     = ncw(1)   * work1(i) + ncw(2)    * work2(i)
            flgmin_l(i) = flgmin(1)* work1(i) + flgmin(2) * work2(i)
          enddo

          if (kdt == 1 .and. abs(xlon(1)) < 0.0001) then
            write(0,*)' xncw=',xncw(1),' rhc=',rhc(1,1),' work1='       &
     &,          work1(1),' work2=',work2(1),' flgmin=',flgmin_l(1)     &
     &,         ' lon=',xlon(1) * 57.29578,' lat=',lat,' me=',me
          endif

          call gsmdrive(im, ix, levs, dtp, con_g, con_hvap, hsub, con_cp&
     &,                 me, lprnt, ipr                                  &
     &,                 prsl, del, rhc, xncw, flgmin_l                  &
     &,                 gt0, gq0(1,1,1), gq0(1,1,ntcw)                  &
     &,                 phy_f3d(1,1,1),  phy_f3d(1,1,2)                 &
     &,                 phy_f3d(1,1,3), rain1, sr)

        elseif (num_p3d == 4) then  ! call zhao/carr/sundqvist microphysics

          if (npdf3d /= 3) then               ! without pdf clouds

!           if (lprnt) then
!             write(0,*)' prsl=',prsl(ipr,:)
!             write(0,*) ' del=',del(ipr,:)
!             write(0,*) ' beflsgt0=',gt0(ipr,:),' kdt=',kdt
!             write(0,*) ' beflsgq0=',gq0(ipr,:,1),' kdt=',kdt
!           endif
                                              ! ------------------
            if (do_shoc) then
              call precpd_shoc(im, ix, levs, dtp, del, prsl,            &
     &                    gq0(1,1,1), gq0(1,1,ntcw), gt0, rain1, sr,    &
     &                    rainp, rhc, psautco_l, prautco_l, evpco,      &
     &                    wminco, phy_f3d(1,1,ntot3d-2), lprnt, ipr)
!    &                    wminco, sgs_cld(1:im,1:levs), lprnt, ipr)
!    &                    wminco, shoc_cld, lprnt, ipr)
            else
              call gscond(im, ix, levs, dtp, dtf, prsl, pgr,            &
     &                    gq0(1,1,1), gq0(1,1,ntcw), gt0,               &
     &                    phy_f3d(1,1,1), phy_f3d(1,1,2), phy_f2d(1,1), &
     &                    phy_f3d(1,1,3), phy_f3d(1,1,4), phy_f2d(1,2), &
     &                    rhc,lprnt, ipr)

              call precpd(im, ix, levs, dtp, del, prsl,                 &
     &                    gq0(1,1,1), gq0(1,1,ntcw), gt0, rain1, sr,    &
     &                    rainp, rhc, psautco_l, prautco_l, evpco,      &
     &                    wminco, lprnt, ipr)
            endif
!           if (lprnt) then
!             write(0,*)' prsl=',prsl(ipr,:)
!             write(0,*) ' del=',del(ipr,:)
!             write(0,*) ' aftlsgt0=',gt0(ipr,:),' kdt=',kdt
!             write(0,*) ' aftlsgq0=',gq0(ipr,:,1),' kdt=',kdt
!             write(0,*)' aft precpd rain1=',rain1(1:3),' lat=',lat
!           endif
          else                                ! with pdf clouds
                                              ! ---------------
            call gscondp(im, ix, levs, dtp, dtf, prsl, pgr,             &
     &                  gq0(1,1,1), gq0(1,1,ntcw), gt0,                 &
     &                  phy_f3d(1,1,1), phy_f3d(1,1,2), phy_f2d(1,1),   &
     &                  phy_f3d(1,1,3), phy_f3d(1,1,4), phy_f2d(1,2),   &
     &                  rhc,phy_f3d(1,1,num_p3d+1),sup,lprnt,           &
     &                  ipr,kdt)

            call precpdp(im, ix, levs, dtp, del, prsl, pgr,             &
     &                  gq0(1,1,1), gq0(1,1,ntcw), gt0, rain1,sr,       &
     &                  rainp, rhc, phy_f3d(1,1,num_p3d+1),             &
     &                  psautco_l, prautco_l, evpco, wminco,            &
     &                  lprnt, ipr)

          endif   ! end of grid-scale precip/microphysics options
        endif     ! end if_num_p3d

!     if (lprnt) write(0,*) ' rain1=',rain1(ipr),' rainc=',rainc(ipr)
!    &,' lat=',lat

      elseif (ncld == 2) then       ! MGB double-moment microphysics
!       Acheng used clw here for other code to run smoothly and minimum change
!       to make the code work. However, the nc and clw should be treated
!       in other procceses too.  August 28/2015; Hope that can be done next
!       year. I believe this will make the physical interaction more reasonable
!       Anning 12/5/2015 changed ntcw hold liquid only
        if (do_shoc) then
          do k=1,levs
            do i=1,im
              clw(i,k,1) = gq0(i,k,ntiw)             ! ice
              clw(i,k,2) = gq0(i,k,ntcw)             ! water
              phy_f3d(i,k,1) = phy_f3d(i,k,ntot3d-2) ! clouds from shoc
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              clw(i,k,1) = gq0(i,k,ntiw)             ! ice
              clw(i,k,2) = gq0(i,k,ntcw)             ! water
              phy_f3d(i,k,1) = min(1.0, phy_f3d(i,k,1)+cnvc(i,k))
                                                     ! clouds from t-dt and cnvc
            enddo
          enddo
        endif
!       notice clw ix instead of im
!       call m_micro_driver(im,ix,levs,flipv,del,dtp,prsl,prsi,
!    &    prslk,prsik,pgr,vvel,clw(1,1,2), QLCN, clw(1,1,1),QICN,
!     if (lprnt) write(0,*)' cnv_mfdbef=',cnv_mfd(ipr,:),' flipv=',flipv
!     if (lprnt) write(0,*)' clw1bef=',clw(ipr,:,1),' kdt=',kdt
        call m_micro_driver(im,   ix,   levs,  flipv,    dtp,
     &                      prsl, prsi, prslk, prsik,
     &                      vvel, clw(1,1,2), QLCN, clw(1,1,1),QICN,
     &                      hlw,  swh, w_upi, cf_upi,
     &                      FRLAND, HPBL, CNV_MFD, CNV_PRC3,
     &                      CNV_DQLDT, CLCN, gu0, gv0,
     &                      dusfc,  dvsfc, dusfc1, dvsfc1,
     &                      dusfc1, dvsfc1, CNV_FICE,
     &                      CNV_NDROP, CNV_NICE,  gq0(1,1,1),
     &                      gq0(1,1,ntcw), gq0(1,1,ntiw), gt0,
     &                      rain1, sr, gq0(1,1,ntlnc),
     &                      gq0(1,1,ntinc), phy_f3d(1,1,1), kbot,
     &                      aero_in, skip_macro, cn_prc, cn_snr,
     &                      lprnt, ipr, kdt)

!     if (lprnt) write(0,*) ' rain1=',rain1(ipr),' rainc=',rainc(ipr)
!    &,' cn_prc=',cn_prc(ipr),' cn_snr=',cn_snr(ipr)
!     if(lprnt) write(0,*) ' aftlsgq0=',gq0(ipr,:,1),' kdt=',kdt
!     if (lprnt) write(0,*)' clw1aft=',gq0(ipr,:,ntiw),' kdt=',kdt

      endif       ! end if_ncld

      do i = 1, im
        rain(i)  = rainc(i) + frain * rain1(i)
      enddo

      if (cal_pre) then       ! hchuang: add dominant precipitation type algorithm

        i = min(3,num_p3d)
        call calpreciptype(kdt,nrcm,im,ix,levs,levs+1,rann,
     &                     xlat,xlon,gt0,gq0,prsl,prsi,rain,
     &                     phii,num_p3d,tsea,sr,phy_f3d(1,1,i),             ! input
     &                     domr,domzr,domip,doms)                           ! output

!
!        if (lprnt) print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS '
!     &,DOMR(ipr),DOMZR(ipr),DOMIP(ipr),DOMS(ipr)
!        do i=1,im
!         if (abs(xlon(i)*57.29578-114.0) .lt. 0.2  .and.
!     &    abs(xlat(i)*57.29578-40.0) .lt. 0.2)
!     &    print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS ',
!     &    DOMR(i),DOMZR(i),DOMIP(i),DOMS(i)
!       end do
!       HCHUANG: use new precipitation type to decide snow flag for LSM snow accumulation

        do i=1,im
          if(doms(i) >0.0 .or. domip(i)>0.0)then
            srflag(i) = 1.
          else
            srflag(i) = 0.
          end if
        enddo
      endif

      if (lssav) then
        do i = 1, im
          totprcp(i) = totprcp(i) + rain(i)
!         cnvprcp(i) = cnvprcp(i) + rainc(i)
        enddo

        if (ldiag3d) then
          do k = 1, levs
            do i = 1, im
              dt3dt(i,k,6) = dt3dt(i,k,6) + (gt0(i,k)-dtdt(i,k)) * frain
              dq3dt(i,k,4) = dq3dt(i,k,4) + (gq0(i,k,1)-dqdt(i,k,1))    &
     &                                                           * frain
            enddo
          enddo
        endif
      endif

!  --- ...  estimate t850 for rain-snow decision

      do i = 1, im
        t850(i) = gt0(i,1)
      enddo

      do k = 1, levs-1
        do i = 1, im
          if (prsl(i,k) > p850 .and. prsl(i,k+1) <= p850) then
            t850(i) = gt0(i,k) - (prsl(i,k)-p850)                       &
     &              / (prsl(i,k)-prsl(i,k+1)) * (gt0(i,k)-gt0(i,k+1))
          endif
        enddo
      enddo

!  --- ...  lu: snow-rain detection is performed in land/sice module

      if (cal_pre) then ! hchuang: new precip type algorithm defines srflag
        do i = 1, im
          tprcp(i) = max(0.0, rain(i))  ! clu: rain -> tprcp
        enddo
      else
        do i = 1, im
          tprcp(i)  = max(0.0, rain(i) )! clu: rain -> tprcp
          srflag(i) = 0.                ! clu: default srflag as 'rain' (i.e. 0)

          if (t850(i) <= 273.16) then
            srflag(i) = 1.              ! clu: set srflag to 'snow' (i.e. 1)
          endif
        enddo
      endif

!  --- ...  coupling insertion

      if (lssav_cpl) then
        do i = 1, im
          if (t850(i) > 273.16) then
             rain_cpl(i) = rain_cpl(i) + rain(i)
          else
             snow_cpl(i) = snow_cpl(i) + rain(i)
          endif
        enddo
      endif

!  --- ...  end coupling insertion

      if (lsm == 0) then   ! for OSU land model

!  --- ...  wei: when calling osu lsm update soil moisture and canopy water
!                after precipitation computaion
        do i = 1, im
          if (t850(i) <= 273.16 .and. islmsk(i) /= 0) then
            weasd(i)  = weasd(i) + 1.e3*rain(i)
            tprcp(i)  = 0.
          endif
        enddo
        call progt2(im,lsoil,rhscnpy,rhsmc,ai,bi,cci,smsoil,            &
     &              islmsk,canopy,tprcp,runof,snowmt,                   &
     &              zsoil,soiltyp,sigmaf,dtf,me)

!  --- ...  wei: let soil liquid water equal to soil total water

        do k = 1, lsoil  ! let soil liquid water equal to soil total water
          do i = 1, im
            if (islmsk(i) == 1) then
              slsoil(i,k) = smsoil(i,k)
             endif
          enddo
        enddo

      endif   ! end if_lsm
!!!
!!! update surface diagnosis fields at the end of phys package
!!! this change allows gocart to use filtered wind fields
!!!
      if ( lgocart ) then
        call sfc_diag(im,pgr,gu0,gv0,gt0,gq0,                           &
     &                tsea,qss,f10m,u10m,v10m,t2m,q2m,work3,            &
     &                evap,ffmm,ffhh,fm10,fh2)

        if (lssav) then
          do i = 1, im
            tmpmax(i)  = max(tmpmax(i),t2m(i))
            tmpmin(i)  = min(tmpmin(i),t2m(i))

            spfhmax(i) = max(spfhmax(i),q2m(i))
            spfhmin(i) = min(spfhmin(i),q2m(i))
          enddo
        endif
      endif

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters

      if (lssav) then
        tem = dtf * 0.001
        do i = 1, im
          runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * tem
          srunoff(i) = srunoff(i) + runof(i) * tem
        enddo
      endif

!  --- ...  xw: return updated ice thickness & concentration to global array

      do i = 1, im
        if (islmsk(i) == 2) then
          hice(i)  = zice(i)
          fice(i)  = cice(i)
          tisfc(i) = tice(i)
        else
          hice(i)  = 0.0
          fice(i)  = 0.0
          tisfc(i) = tsea(i)
        endif
      enddo

!  --- ...  return updated smsoil and stsoil to global arrays

      do k = 1, lsoil
        do i = 1, im
          smc(i,k) = smsoil(i,k)
          stc(i,k) = stsoil(i,k)
          slc(i,k) = slsoil(i,k)
        enddo
      enddo

!  --- ...  calculate column precipitable water "pwat"

      do i = 1, im
        pwat(i)  = 0.
!       rqtk(i)  = 0.
!       work2(i) = 1.0 / pgr(i)
      enddo

      do k = 1, levs
        do i = 1, im
          work1(i) = 0.0
        enddo

        if (ncld > 0) then
          do ic = ntcw, ntcw+ncld-1
            do i = 1, im
              work1(i) = work1(i) + gq0(i,k,ic)
            enddo
          enddo
        endif

        do i = 1, im
          pwat(i) = pwat(i) + del(i,k)*(gq0(i,k,1)+work1(i))
          rqtk(i) = rqtk(i) + del(i,k)*(gq0(i,k,1)-qgrs(i,k,1))
        enddo
      enddo

      do i = 1, im
        pwat(i) = pwat(i) * (1.0/con_g)
!       rqtk(i) = rqtk(i) * work2(i)
      enddo
!
!     if (lat == 45) write(1000+me,*)' pwat=',pwat(1),' kdt=',kdt
!       if (lprnt) then
!         write(0,*) ' endgt0=',gt0(ipr,:),' kdt=',kdt
!         write(0,*) ' endgq0=',gq0(ipr,:,1),' kdt=',kdt
!       endif
      deallocate (clw)
      if (do_shoc) then
        deallocate (qpl, qpi, ncpl, ncpi)
      endif
      if (.not. ras .or. .not. cscnv) then
        deallocate (cnvc, cnvw)
      endif
!     deallocate (fscav, fswtr)
!
!     if (lprnt) write(0,*)' end of gbphys at kdt=',kdt,
!    &' rain=',rain(ipr),' rainc=',rainc(ipr)
!     if (lprnt) call mpi_quit(7)
!     if (kdt >2 ) call mpi_quit(7)
!     if (ncld == 2) then         ! For MGB double moment microphysics
        deallocate (qlcn,      qicn,    w_upi
     &,             cf_upi,    CNV_MFD, CNV_PRC3
     &,             CNV_DQLDT, clcn,    cnv_fice
     &,             cnv_ndrop, cnv_nice)
!     endif

      return
!...................................
      end subroutine gbphys
!-----------------------------------
!> @}
