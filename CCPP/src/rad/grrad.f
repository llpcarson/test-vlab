!> \file grrad.f This file is the radiation driver module. It prepares
!! the atmospheric profiles and invokes the main radiation calculation.

!> \defgroup rad RRTMG Shortwave/Longwave Radiation Scheme
!> @{
!!  \brief The GFS radiation scheme 
!!  \details Radiative processes are among the most complex and
!! computationally intensive parts of all model physics. As an 
!! essential component of modeling the atmosphere, radiation directly 
!! and indirectly connects all physics processes with model dynamics,
!! and it regulates the overall earth-atmosphere energy exchanges and 
!! transformations.
!!
!! The radiation package in GFS physics has standardized component 
!! modules (Table 1). The radiation driver module (\ref 
!! module_radiation_driver) is the interface with the Interoperable 
!! Physics Driver (IPD) for NGGPS, and it has three subroutines called
!! by IPD (Figure 1):
!! - radinit() is called in subroutine nuopc_phys_init to set up
!! radiation related fixed parameters.
!! - radupdate() is called in subroutine nuopc_rad_update to update
!! values between timesteps.
!! - grrad() is called in subroutine nuopc_rad_run, and it is the 
!! driver of radiation calculation.
!! \image html ipd_rad.png "Figure 1: Schematic illustration of the
!! communication between the GFS radiation package and IPD " width=10cm
!!
!! The schematic radiation module structure is shown in Table 1. 
!! \image html schematic_Rad_mod.png "Table 1: Schematic Radiation
!! Module Structure" width=10cm
!!
!> GFS radiation package is intended to provide a fast and accurate
!! method of determining the total radiative flux at any given 
!! location. These calculations provide both the total radiative flux
!! at the ground surface, which is needed to establish the surface 
!! energy budget, and the vertical radiative flux divergence, which is
!! used to calculate the radiative heating and cooling rates of a given
!! atmospheric layer. The magnitude of the terms in the surface energy
!! budget can set the stage for moist deep convection and are crucial
!! to the formation of low-level clouds. In addition, the vertical
!! radiative flux divergence can produce substantial cooling,
!! particularly at the tops of clouds, which can have strong dynamical
!! effects on cloud evolution.
!!
!! It uses a correlated-k distribution method and a transmittance lookup
!! table that is linearly scaled by optical depth to achieve high
!! accuracy and efficiency. The algorithm contains 140 unevenly
!! distributed quadrature points (reduced from the original set of 256)
!! to integrate the cumulative probability distribution functions of
!! absorption over 16 spectral bands. It employs the 
!! Clough-Kneizys-Davies (CKD_2.4) continuum model (Clough et al. 1992
!! \cite clough_et_al_1992) to compute absorption by water vapor at the
!! continuum band. Longwave cloud radiative properties external to the
!! RRTM depend on cloud liquid/ice water path and the effective radius
!! of ice particles and water droplets (Hu and Stamnes 1993 \cite 
!! hu_and_stamnes_1993; Ebert and Curry 1992 \cite ebert_and_curry_1992).
!!
!! Changes to Radiation Parameterization since 2007:
!! \n The longwave (LW) and the shortwave (SW) radiation
!! parameterizations in NCEP's operational GFS are both modified and 
!! optimized versions of the Rapid Radiative Transfer Model for GCMs 
!! (RRTMG_LW v2.3 and RRTMG_SW v2.3, respectively) developed at AER 
!! (Iacono et al. 2008 \cite iacono_et_al_2008,Mlawer et al. 1997 
!! \cite mlawer_et_al_1997, Iacono et al., 2000 \cite iacono_et_al_2000,
!! Clough et al., 2005 \cite clough_et_al_2005). The LW algorithm 
!! contains 140 unevenly distributed g-points (quadrature points) in 16
!! broad spectral bands, while the SW algorithm includes 112 g-points
!! in 14 bands. In addition to the major atmospheric absorbing gases of
!! ozone, water vapor, and carbon dioxide, the algorithm also includes
!! various minor absorbing species such as methane, nitrous oxide, 
!! oxygen, and in the longwave up to four types of halocarbons (CFCs).
!! To represent statistically the unresolved subgrid cloud variability
!! when dealing multi layered clouds, a Monte-Carlo Independent Column 
!! Approximation (McICA) method is used in the RRTMG radiative transfer.
!! A maximum-random cloud overlap method is used in both LW and SW 
!! radiation calculations. Cloud condensate path and effective radius
!! for water and ice are used for the calculation of cloud-radiative 
!! properties. Hu and Stamnes's method (1993) \cite hu_and_stamnes_1993
!! is used to treat water clouds in both LW and SW parameterizations. 
!! For ice clouds. Fu's parameterizations (1996,1998) \cite fu_1996 
!! \cite fu_et_al_1998 are used in the SW and LW, respectively.
!!
!! In the operational GFS, a climatological tropospheric aerosol with
!! a 5-degree horizontal resolution is used in both LW and SW 
!! radiations. A generalized spectral mapping formulation was developed
!! to compute radiative properties of various aerosol components for 
!! each of the radiation spectral bands. A separate stratospheric
!! volcanic aerosol parameterization was added that is capable of
!! handling volcanic events. In SW, a new table of incoming solar
!! constants is derived covering time period of 1850-2019 (Vandendool, 
!! personal communivation). An eleven-year solar cycle approximation
!! will be used for time out of the window period in long term climate
!! simulations. The SW albedo parameterization uses surface vegetation
!! type based seasonal climatology similar to that described in the
!! NCEP OFFICE Note 441 (Hou et al. 2002 \cite hou_et_al_2002) but with
!! a modification in the treatment of solar zenith angle dependency over
!! snow-free land surface (Yang et al. 2008 \cite yang_et_al_2008).
!! Similarly, vegetation type based non-black-body surface emissivity
!! is used for the LW radiation. Concentrations of atmospheric
!! greenhouse gases are either obtained from global network 
!! measurements, such as carbon dioxide (CO2), or taking the
!! climatological constants, the actual CO2 value for the forecast time
!! is an estimation based on the most recent five-year observations. In
!! the lower atmosphere (<3km) a monthly mean CO2 distribution in 15 
!! degree horizontal resolution is used, while a global mean monthly
!! value is used in the upper atmosphere.
!!
!> \defgroup module_radiation_driver module_radiation_driver
!> @{
!! \brief The GFS radiation driver module 
!! \details module_radiation_driver prepares the atmospheric profile, 
!! invokes the main radiation calculations, and computes radiative 
!! fluxes and heating rates for some arbitrary number of vertical
!! columns. This module also regulates the logistic running flow of
!! the computations, such as data initialization and update accordance
!! with forecast timing progress, the sequential order of subroutine
!! calls, and sorting results for final output.
!! There are three externally accessible subroutines: 
!! - radinit(): the initialization subroutine of radiation calculations
!! It is invoked by the model's initialization process and is independent
!! with forecat time progress.
!! - radupdate(): calls many update subroutines to check and update 
!! radiation required but time varying data sets and module variables.
!! It is placed inside a model's time advancing loop.
!! - grrad(): the driver of radiation calculation subroutines. It sets 
!! up profile variables for radiation input, including clouds, surface 
!! albedos, atmospheric aerosols, ozone, etc. It is located inside the
!! timing loop, and control the sequence of the radiative process 
!! calculations.
!! \version NCEP-Radiation_driver    v5.2  Jan 2013
! ==========================================================  !!!!!
!             'module_radiation_driver' descriptions           !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!   this is the radiation driver module.  it prepares atmospheric      !
!   profiles and invokes main radiation calculations.                  !
!                                                                      !
!   in module 'module_radiation_driver' there are twe externally       !
!   callable subroutine:                                               !
!                                                                      !
!      'radinit'    -- initialization routine                          !
!         input:                                                       !
!           ( si, NLAY, me )                                           !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'radupdate'  -- update time sensitive data used by radiations   !
!         input:                                                       !
!           ( idate,jdate,deltsw,deltim,lsswr, me )                    !
!         output:                                                      !
!           ( slag,sdec,cdec,solcon )                                  !
!                                                                      !
!      'grrad'      -- setup and invoke main radiation calls           !
!         input:                                                       !
!          ( prsi,prsl,prslk,tgrs,qgrs,tracer,vvl,slmsk,               !
!            xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,            !
!            alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,           !
!            sinlat,coslat,solhr,jdate,solcon,                         !
!            cv,cvt,cvb,fcice,frain,rrime,flgmin,                      !
!            icsdsw,icsdlw, ntcw,ncld,ntoz, NTRAC,NFXR,                !
!            dtlw,dtsw, lsswr,lslwr,lssav,                             !
!            IX, IM, LM, me, lprnt, ipt, kdt,deltaq,sup,cnvw,cnvc,     !
!         output:                                                      !
!            htrsw,topfsw,sfcfsw,dswcmp,uswcmp,sfalb,coszen,coszdg,    !
!            htrlw,topflw,sfcflw,tsflw,semis,cldcov,                   !
!         input/output:                                                !
!            fluxr                                                     !
!         optional output:                                             !
!            htrlw0,htrsw0,htrswb,htrlwb                               !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module physparam'                  in 'physparam.f'           !
!       'module funcphys'                   in 'funcphys.f'            !
!       'module physcons'                   in 'physcons.f'            !
!                                                                      !
!       'module module_radiation_gases'     in 'radiation_gases.f'     !
!       'module module_radiation_aerosols'  in 'radiation_aerosols.f'  !
!       'module module_radiation_surface'   in 'radiation_surface.f'   !
!       'module module_radiation_clouds'    in 'radiation_clouds.f'    !
!                                                                      !
!       'module module_radsw_cntr_para'     in 'radsw_xxxx_param.f'    !
!       'module module_radsw_parameters'    in 'radsw_xxxx_param.f'    !
!       'module module_radsw_main'          in 'radsw_xxxx_main.f'     !
!                                                                      !
!       'module module_radlw_cntr_para'     in 'radlw_xxxx_param.f'    !
!       'module module_radlw_parameters'    in 'radlw_xxxx_param.f'    !
!       'module module_radlw_main'          in 'radlw_xxxx_main.f'     !
!                                                                      !
!    where xxxx may vary according to different scheme selection       !
!                                                                      !
!                                                                      !
!   program history log:                                               !
!     mm-dd-yy    ncep         - created program grrad                 !
!     08-12-03    yu-tai hou   - re-written for modulized radiations   !
!     11-06-03    yu-tai hou   - modified                              !
!     01-18-05    s. moorthi   - NOAH/ICE model changes added          !
!     05-10-05    yu-tai hou   - modified module structure             !
!     12-xx-05    s. moorthi   - sfc lw flux adj by mean temperature   !
!     02-20-06    yu-tai hou   - add time variation for co2 data, and  !
!                                solar const. add sfc emiss change     !
!     03-21-06    s. Moorthi   - added surface temp over ice           !
!     07-28-06    yu-tai hou   - add stratospheric vocanic aerosols    !
!     03-14-07    yu-tai hou   - add generalized spectral band interp  !
!                                for aerosol optical prop. (sw and lw) !
!     04-10-07    yu-tai hou   - spectral band sw/lw heating rates     !
!     05-04-07    yu-tai hou   - make options for clim based and modis !
!                                based (h. wei and c. marshall) albedo !
!     09-05-08    yu-tai hou   - add the initial date and time 'idate' !
!                    and control param 'ICTM' to the passing param list!
!                    to handel different time/date requirements for    !
!                    external data (co2, aeros, solcon, ...)           !
!     10-10-08    yu-tai hou   - add the ICTM=-2 option for combining  !
!                    initial condition data with seasonal cycle from   !
!                    climatology.                                      !
!     03-12-09    yu-tai hou   - use two time stamps to keep tracking  !
!                    dates for init cond and fcst time. remove volcanic!
!                    aerosols data in climate hindcast (ICTM=-2).      !
!     03-16-09    yu-tai hou   - included sub-column clouds approx.    !
!                    control flags isubcsw/isubclw in initialization   !
!                    subroutine. passed auxiliary cloud control arrays !
!                    icsdsw/icsdlw (if isubcsw/isubclw =2, it will be  !
!                    the user provided permutation seeds) to the sw/lw !
!                    radiation calculation programs. also moved cloud  !
!                    overlapping control flags iovrsw/iovrlw from main !
!                    radiation routines to the initialization routines.!
!     04-02-09    yu-tai hou   - modified surface control flag iems to !
!                    have additional function of if the surface-air    !
!                    interface have the same or different temperature  !
!                    for radiation calculations.                       !
!     04-03-09    yu-tai hou   - modified to add lw surface emissivity !
!                    as output variable. changed the sign of sfcnsw to !
!                    be positive value denote solar flux goes into the !
!                    ground (this is needed to reduce sign confusion   !
!                    in other part of model)                           !
!     09-09-09    fanglin yang (thru s.moorthi) added QME5 QME6 to E-20!
!     01-09-10    sarah lu     - added gocart option, revised grrad for!
!                    gocart coupling. calling argument modifed: ldiag3 !
!                    removed; cldcov/fluxr sequence changed; cldcov is !
!                    changed from accumulative to instant field and    !
!                    from input/output to output field                 !
!     01-24-10    sarah lu     - added aod to fluxr, added prslk and   !
!                    oz to setaer input argument (for gocart coupling),!
!                    added tau_gocart to setaer output argument (for,  !
!                    aerosol diag by index of nv_aod)                  !
!     07-08-10    s.moorthi - updated the NEMS version for new physics !
!     07-28-10    yu-tai hou   - changed grrad interface to allow all  !
!                    components of sw/lw toa/sfc instantaneous values  !
!                    being passed to the calling program. moved the    !
!                    computaion of sfc net sw flux (sfcnsw) to the     !
!                    calling program. merged carlos' nmmb modification.!
!     07-30-10    s. moorthi - corrected some errors associated with   !
!                    unit changes                                      !
!     12-02-10    s. moorthi/y. hou - removed the use of aerosol flags !
!                    'iaersw' 'iaerlw' from radiations and replaced    !
!                    them by using the runtime variable laswflg and    !
!                    lalwflg defined in module radiation_aerosols.     !
!                    also replaced param nspc in grrad with the use of !
!                    max_num_gridcomp in module radiation_aerosols.    !
!     jun 2012    yu-tai hou   - added sea/land madk 'slmsk' to the    !
!                    argument list of subrotine setaer call for the    !
!                    newly modified horizontal bi-linear interpolation !
!                    in climatological aerosols schem. also moved the  !
!                    virtual temperature calculations in subroutines   !
!                    'radiation_clouds' and 'radiation_aerosols' to    !
!                    'grrad' to reduce repeat comps. renamed var oz as !
!                    tracer to reflect that it carries various prog    !
!                    tracer quantities.                                !
!                              - modified to add 4 compontents of sw   !
!                    surface downward fluxes to the output. (vis/nir;  !
!                    direct/diffused). re-arranged part of the fluxr   !
!                    variable fields and filled the unused slots for   !
!                    the new components.  added check print of select  !
!                    data (co2 value for now).                         !
!                              - changed the initialization subrution  !
!                    'radinit' into two parts: 'radinit' is called at  !
!                    the start of model run to set up radiation related!
!                    fixed parameters; and 'radupdate' is called in    !
!                    the time-loop to update time-varying data sets    !
!                    and module variables.                             !
!     sep 2012    h-m lin/y-t hou added option of extra top layer for  !
!                    models with low toa ceiling. the extra layer will !
!                    help ozone absorption at higher altitude.         !
!     nov 2012    yu-tai hou   - modified control parameters through   !
!                    module 'physparam'.                                !
!     jan 2013    yu-tai hou   - updated subr radupdate for including  !
!                    options of annual/monthly solar constant table.   !
!     mar 2013    h-m lin/y-t hou corrected a bug in extra top layer   !
!                    when using ferrier microphysics.                  !
!     may 2013    s. mooorthi - removed fpkapx                         !
!     jul 2013    r. sun - added pdf cld and convective cloud water and!
!                          cover for radiation                         !
!     aug 2013    s. moorthi  - port from gfs to nems                  !
!     13Feb2014   sarah lu - add aerodp to fluxr                       !
!     Apr 2014    Xingren Wu - add sfc SW downward fluxes nir/vis and  !
!                    sfcalb to export for A/O/I coupling               !
!     jun 2014    y-t hou    - revised code to include surface up and  !
!                    down spectral components sw fluxes as output.     !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radiation_driver     !
!........................................!
!
      use physparam
      use physcons,                 only : eps   => con_eps,            &
     &                                     epsm1 => con_epsm1,          &
     &                                     fvirt => con_fvirt           &
     &,                                    rocp  => con_rocp
      use funcphys,                 only : fpvs

      use module_radiation_astronomy,only: sol_init, sol_update, coszmn
      use module_radiation_gases,   only : NF_VGAS, getgases, getozn,   &
     &                                     gas_init, gas_update
      use module_radiation_aerosols,only : NF_AESW, NF_AELW, setaer,    &
     &                                     aer_init, aer_update,        &
     &                                     NSPC1                       
      use module_radiation_surface, only : NF_ALBD, sfc_init, setalb,   &
     &                                     setemis
      use module_radiation_clouds,  only : NF_CLDS, cld_init,           &
     &                                     progcld1, progcld2, progcld3,&
     &                                     diagcld1

      use module_radsw_parameters,  only : topfsw_type, sfcfsw_type,    &
     &                                     profsw_type,cmpfsw_type,NBDSW
      use module_radsw_main,        only : rswinit,  swrad

      use module_radlw_parameters,  only : topflw_type, sfcflw_type,    &
     &                                     proflw_type, NBDLW
      use module_radlw_main,        only : rlwinit,  lwrad
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGRAD='NCEP-Radiation_driver    v5.2  Jan 2013 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.1  Nov 2012 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.0  Aug 2012 '

!>\name Constant values

!> lower limit of saturation vapor pressure (=1.0e-10)
      real (kind=kind_phys) :: QMIN
!> lower limit of specific humidity (=1.0e-7)
      real (kind=kind_phys) :: QME5
!> lower limit of specific humidity (=1.0e-7)
      real (kind=kind_phys) :: QME6
!> EPSQ=1.0e-12
      real (kind=kind_phys) :: EPSQ
!     parameter (QMIN=1.0e-10, QME5=1.0e-5,  QME6=1.0e-6,  EPSQ=1.0e-12)
      parameter (QMIN=1.0e-10, QME5=1.0e-7,  QME6=1.0e-7,  EPSQ=1.0e-12)
!     parameter (QMIN=1.0e-10, QME5=1.0e-20, QME6=1.0e-20, EPSQ=1.0e-12)

!> lower limit of toa pressure value in mb 
      real, parameter :: prsmin = 1.0e-6

!> control flag for LW surface temperature at air/ground interface
!! (default=0, the value will be set in subroutine radinit)
      integer :: itsfc  =0

!> new data input control variables (set/reset in subroutines radinit/radupdate):
      integer :: month0=0,   iyear0=0,   monthd=0

!> control flag for the first time of reading climatological ozone data 
!! (set/reset in subroutines radinit/radupdate, it is used only if the
!! control parameter ioznflg=0)
      logical :: loz1st =.true.

!> optional extra top layer on top of low ceiling models
!!\n LTP=0: no extra top layer
      integer, parameter :: LTP = 0   ! no extra top layer
!     integer, parameter :: LTP = 1   ! add an extra top layer

!> control flag for extra top layer
      logical, parameter :: lextop = (LTP > 0)

!  ---  publicly accessible module programs:

      public radinit, radupdate, grrad


! =================
      contains
! =================

!> This subroutine initialize a model's radiation process through 
!! calling of specific initialization subprograms that directly
!! related to radiation calculations. This subroutine needs to be
!! invoked only once at the start stage of a model's run, and the
!! call is placed outside of both the time advancement loop and 
!! horizontal grid loop.
!> \param si       model vertical sigma interface
!> \param nlay     number of model vertical layers
!> \param me       print control flag
!> \section gen_radinit General Algorithm
!> @{
!-----------------------------------
      subroutine radinit( si, NLAY, me )
!...................................

!  ---  inputs:
!     &     ( si, NLAY, me )
!  ---  outputs:
!          ( none )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radinit     initialization of radiation calculations    !
!                                                                       !
! usage:        call radinit                                            !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   wcoss                                                   !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   si               : model vertical sigma interface                   !
!   NLAY             : number of model vertical layers                  !
!   me               : print control flag                               !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in module physparam)                     !
!   isolar   : solar constant cntrol flag                               !
!              = 0: use the old fixed solar constant in "physcon"       !
!              =10: use the new fixed solar constant in "physcon"       !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cycl apprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cycl apprx!
!   iaerflg  : 3-digit aerosol flag (abc for volc, lw, sw)              !
!              a:=0 use background stratospheric aerosol                !
!                =1 include stratospheric vocanic aeros                 !
!              b:=0 no topospheric aerosol in lw radiation              !
!                =1 compute tropspheric aero in 1 broad band for lw     !
!                =2 compute tropspheric aero in multi bands for lw      !
!              c:=0 no topospheric aerosol in sw radiation              !
!                =1 include tropspheric aerosols for sw                 !
!   ico2flg  : co2 data source control flag                             !
!              =0: use prescribed global mean co2 (old  oper)           !
!              =1: use observed co2 annual mean value only              !
!              =2: use obs co2 monthly data with 2-d variation          !
!   ictmflg  : =yyyy#, external data ic time/date control flag          !
!              =   -2: same as 0, but superimpose seasonal cycle        !
!                      from climatology data set.                       !
!              =   -1: use user provided external data for the          !
!                      forecast time, no extrapolation.                 !
!              =    0: use data at initial cond time, if not            !
!                      available, use latest, no extrapolation.         !
!              =    1: use data at the forecast time, if not            !
!                      available, use latest and extrapolation.         !
!              =yyyy0: use yyyy data for the forecast time,             !
!                      no further data extrapolation.                   !
!              =yyyy1: use yyyy data for the fcst. if needed, do        !
!                      extrapolation to match the fcst time.            !
!   ioznflg  : ozone data source control flag                           !
!              =0: use climatological ozone profile                     !
!              =1: use interactive ozone profile                        !
!   ialbflg  : albedo scheme control flag                               !
!              =0: climatology, based on surface veg types              !
!              =1: modis retrieval based surface albedo scheme          !
!   iemsflg  : emissivity scheme cntrl flag (ab 2-digit integer)        !
!              a:=0 set sfc air/ground t same for lw radiation          !
!                =1 set sfc air/ground t diff for lw radiation          !
!              b:=0 use fixed sfc emissivity=1.0 (black-body)           !
!                =1 use varying climtology sfc emiss (veg based)        !
!                =2 future development (not yet)                        !
!   icldflg  : cloud optical property scheme control flag               !
!              =0: use diagnostic cloud scheme                          !
!              =1: use prognostic cloud scheme (default)                !
!   icmphys  : cloud microphysics scheme control flag                   !
!              =1 zhao/carr/sundqvist microphysics scheme               !
!              =2 brad ferrier microphysics scheme                      !
!              =3 zhao/carr/sundqvist microphysics+pdf cloud & cnvc,cnvw!
!   iovrsw   : control flag for cloud overlap in sw radiation           !
!   iovrlw   : control flag for cloud overlap in lw radiation           !
!              =0: random overlapping clouds                            !
!              =1: max/ran overlapping clouds                           !
!   isubcsw  : sub-column cloud approx control flag in sw radiation     !
!   isubclw  : sub-column cloud approx control flag in lw radiation     !
!              =0: with out sub-column cloud approximation              !
!              =1: mcica sub-col approx. prescribed random seed         !
!              =2: mcica sub-col approx. provided random seed           !
!   lcrick   : control flag for eliminating CRICK                       !
!              =t: apply layer smoothing to eliminate CRICK             !
!              =f: do not apply layer smoothing                         !
!   lcnorm   : control flag for in-cld condensate                       !
!              =t: normalize cloud condensate                           !
!              =f: not normalize cloud condensate                       !
!   lnoprec  : precip effect in radiation flag (ferrier microphysics)   !
!              =t: snow/rain has no impact on radiation                 !
!              =f: snow/rain has impact on radiation                    !
!   ivflip   : vertical index direction control flag                    !
!              =0: index from toa to surface                            !
!              =1: index from surface to toa                            !
!                                                                       !
!  subroutines called: sol_init, aer_init, gas_init, cld_init,          !
!                      sfc_init, rlwinit, rswinit                       !
!                                                                       !
!  usage:       call radinit                                            !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, me

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none, to module variables)

!  ---  locals:

!
!===> ...  begin here
!
!> -# Set up control variables and external module variables in 
!!    module physparam
      itsfc  = iemsflg / 10             ! sfc air/ground temp control
      loz1st = (ioznflg == 0)           ! first-time clim ozone data read flag
      month0 = 0
      iyear0 = 0
      monthd = 0

      if (me == 0) then
!       print *,' NEW RADIATION PROGRAM STRUCTURES -- SEP 01 2004'
        print *,' NEW RADIATION PROGRAM STRUCTURES BECAME OPER. ',      &
     &          '  May 01 2007'
        print *, VTAGRAD                !print out version tag
        print *,' - Selected Control Flag settings: ICTMflg=',ictmflg,  &
     &    ' ISOLar =',isolar, ' ICO2flg=',ico2flg,' IAERflg=',iaerflg,  &
     &    ' IALBflg=',ialbflg,' IEMSflg=',iemsflg,' ICLDflg=',icldflg,  &
     &    ' ICMPHYS=',icmphys,' IOZNflg=',ioznflg
        print *,' IVFLIP=',ivflip,' IOVRSW=',iovrsw,' IOVRLW=',iovrlw,  &
     &    ' ISUBCSW=',isubcsw,' ISUBCLW=',isubclw
!       write(0,*)' IVFLIP=',ivflip,' IOVRSW=',iovrsw,' IOVRLW=',iovrlw,&
!    &    ' ISUBCSW=',isubcsw,' ISUBCLW=',isubclw
        print *,' LCRICK=',lcrick,' LCNORM=',lcnorm,' LNOPREC=',lnoprec
        print *,' LTP =',LTP,', add extra top layer =',lextop

        if ( ictmflg==0 .or. ictmflg==-2 ) then
          print *,'   Data usage is limited by initial condition!'
          print *,'   No volcanic aerosols'
        endif

        if ( isubclw == 0 ) then
          print *,' - ISUBCLW=',isubclw,' No McICA, use grid ',         &
     &            'averaged cloud in LW radiation'
        elseif ( isubclw == 1 ) then
          print *,' - ISUBCLW=',isubclw,' Use McICA with fixed ',       &
     &            'permutation seeds for LW random number generator'
        elseif ( isubclw == 2 ) then
          print *,' - ISUBCLW=',isubclw,' Use McICA with random ',      &
     &            'permutation seeds for LW random number generator'
        else
          print *,' - ERROR!!! ISUBCLW=',isubclw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw == 0 ) then
          print *,' - ISUBCSW=',isubcsw,' No McICA, use grid ',         &
     &            'averaged cloud in SW radiation'
        elseif ( isubcsw == 1 ) then
          print *,' - ISUBCSW=',isubcsw,' Use McICA with fixed ',       &
     &            'permutation seeds for SW random number generator'
        elseif ( isubcsw == 2 ) then
          print *,' - ISUBCSW=',isubcsw,' Use McICA with random ',      &
     &            'permutation seeds for SW random number generator'
        else
          print *,' - ERROR!!! ISUBCSW=',isubcsw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw /= isubclw ) then
          print *,' - *** Notice *** ISUBCSW /= ISUBCLW !!!',           &
     &            isubcsw, isubclw
        endif
      endif

!> -# Initialization
!! - astronomy initialization routine:
!! call module_radiation_astronomy::sol_init()
!! - aerosols initialization routine:
!! call module_radiation_aerosols::aer_init()
!! - CO2 and other gases intialization routine:
!! call module_radiation_gases::gas_init()
!! - surface intialization routine:
!! call module_radiation_surface::sfc_init()
!! - cloud initialization routine:
!! call module_radiation_clouds::cld_init()
!! - LW radiation initialization routine:
!! call module_radlw_main::rlwinit()
!! - SW radiation initialization routine:
!! call module_radsw_main::rswinit()
!     Initialization

      call sol_init ( me )          !  --- ...  astronomy initialization routine

      call aer_init ( NLAY, me )    !  --- ...  aerosols initialization routine

      call gas_init ( me )          !  --- ...  co2 and other gases initialization routine

      call sfc_init ( me )          !  --- ...  surface initialization routine

      call cld_init ( si, NLAY, me) !  --- ...  cloud initialization routine

      call rlwinit ( me )           !  --- ...  lw radiation initialization routine

      call rswinit ( me )           !  --- ...  sw radiation initialization routine
!
      return
!...................................
      end subroutine radinit
!-----------------------------------
!> @}

!> This subroutine checks and updates time sensitive data used by 
!! radiation computations. This subroutine needs to be placed inside
!! the time advancement loop but outside of the horizontal grid loop.
!! It is invoked at radiation calling frequncy but before any actual
!! radiative transfer computations.
!! \param idate          NCEP absolute date and time of intial condition 
!!                       (year,month,day,time-zone,hour,minute,second,
!!                        mil-second)
!! \param jdate          NCEP absolute date and time at forecast time
!!                       (year,month,day,time-zone,hour,minute,second,
!!                        mil-second)
!! \param deltsw         SW radiation calling time interval in seconds
!! \param deltim         model advancing time-step duration in seconds
!! \param lsswr          logical control flag for SW radiation calculations
!! \param me             print control flag
!! \param slag           equation of time in radians
!! \param sdec,cdec      sine and cosine of the solar declination angle
!! \param solcon         solar constant adjusted by sun-earth distance \f$(W/m^2)\f$
!> \section gen_radupdate General Algorithm
!> @{
!-----------------------------------
      subroutine radupdate( idate,jdate,deltsw,deltim,lsswr, me,    &
     &       slag,sdec,cdec,solcon)
!...................................

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radupdate   calls many update subroutines to check and  !
!   update radiation required but time varying data sets and module     !
!   variables.                                                          !
!                                                                       !
! usage:        call radupdate                                          !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   ibm sp                                                   !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   idate(8)       : ncep absolute date and time of initial condition   !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   jdate(8)       : ncep absolute date and time at fcst time           !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   deltsw         : sw radiation calling frequency in seconds          !
!   deltim         : model timestep in seconds                          !
!   lsswr          : logical flags for sw radiation calculations        !
!   me             : print control flag                                 !
!                                                                       !
!  outputs:                                                             !
!   slag           : equation of time in radians                        !
!   sdec, cdec     : sin and cos of the solar declination angle         !
!   solcon         : sun-earth distance adjusted solar constant (w/m2)  !
!                                                                       !
!  external module variables:                                           !
!   isolar   : solar constant cntrl  (in module physparam)               !
!              = 0: use the old fixed solar constant in "physcon"       !
!              =10: use the new fixed solar constant in "physcon"       !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cycl apprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cycl apprx!
!   ictmflg  : =yyyy#, external data ic time/date control flag          !
!              =   -2: same as 0, but superimpose seasonal cycle        !
!                      from climatology data set.                       !
!              =   -1: use user provided external data for the          !
!                      forecast time, no extrapolation.                 !
!              =    0: use data at initial cond time, if not            !
!                      available, use latest, no extrapolation.         !
!              =    1: use data at the forecast time, if not            !
!                      available, use latest and extrapolation.         !
!              =yyyy0: use yyyy data for the forecast time,             !
!                      no further data extrapolation.                   !
!              =yyyy1: use yyyy data for the fcst. if needed, do        !
!                      extrapolation to match the fcst time.            !
!                                                                       !
!  module variables:                                                    !
!   loz1st   : first-time clim ozone data read flag                     !
!                                                                       !
!  subroutines called: sol_update, aer_update, gas_update               !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: idate(:), jdate(:), me
      logical, intent(in) :: lsswr

      real (kind=kind_phys), intent(in) :: deltsw, deltim

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon

!  ---  locals:
      integer :: iyear, imon, iday, ihour
      integer :: kyear, kmon, kday, khour

      logical :: lmon_chg       ! month change flag
      logical :: lco2_chg       ! cntrl flag for updating co2 data
      logical :: lsol_chg       ! cntrl flag for updating solar constant
!
!===> ...  begin here
!
!> -# Set up time stamp at fcst time and that for green house gases
!! (currently co2 only)
!  --- ...  time stamp at fcst time

      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihour = jdate(5)

!  --- ...  set up time stamp used for green house gases (** currently co2 only)

      if ( ictmflg==0 .or. ictmflg==-2 ) then  ! get external data at initial condition time
        kyear = idate(1)
        kmon  = idate(2)
        kday  = idate(3)
        khour = idate(5)
      else                           ! get external data at fcst or specified time
        kyear = iyear
        kmon  = imon
        kday  = iday
        khour = ihour
      endif   ! end if_ictmflg_block

      if ( month0 /= imon ) then
        lmon_chg = .true.
        month0   = imon
      else
        lmon_chg = .false.
      endif

!> -# Call module_radiation_astronomy::sol_update(), yearly update, no
!! time interpolation.
      if (lsswr) then

        if ( isolar == 0 .or. isolar == 10 ) then
          lsol_chg = .false.
        elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
        else
          lsol_chg = ( isolar==4 .and. lmon_chg )
        endif
        iyear0 = iyear

        call sol_update                                                 &
!  ---  inputs:
     &     ( jdate,kyear,deltsw,deltim,lsol_chg, me,                    &
!  ---  outputs:
     &       slag,sdec,cdec,solcon                                      &
     &     )

      endif  ! end_if_lsswr_block

!> -# Call module_radiation_aerosols::aer_update(), monthly update, no
!! time interpolation
      if ( lmon_chg ) then
        call aer_update ( iyear, imon, me )
      endif

!> -# Call co2 and other gases update routine:
!! module_radiation_gases::gas_update()
      if ( monthd /= kmon ) then
        monthd = kmon
        lco2_chg = .true.
      else
        lco2_chg = .false.
      endif

      call gas_update ( kyear,kmon,kday,khour,loz1st,lco2_chg, me )

      if ( loz1st ) loz1st = .false.

!> -# Call surface update routine (currently not needed)
!     call sfc_update ( iyear, imon, me )

!> -# Call clouds update routine (currently not needed)
!     call cld_update ( iyear, imon, me )
!
      return
!...................................
      end subroutine radupdate
!-----------------------------------
!> @}

!> This subroutine is the driver of main radiation calculations. It
!! sets up column profiles, such as pressure, temperature, moisture,
!! gases, clouds, aerosols, etc., as well as surface radiative 
!! characteristics, such as surface albedo, and emissivity. The call
!! of this subroutine is placed inside both the time advancing loop
!! and the horizontal grid loop.
!! \param prsi       model level pressure in Pa
!! \param prsl       model layer mean pressure in Pa
!! \param prslk      exner function = \f$ (p/p0)^{rocp} \f$
!! \param tgrs       model layer mean temperature in K
!! \param qgrs       layer specific humidity in gm/gm
!! \param tracer     layer prognostic tracer amount mixing-ratio, 
!!                   including: ozone,cloud condensate,aerosols,etc
!! \param vvl        layer mean vertical velocity in pa/sec
!!                   (used only for the legacy diagnostic style of 
!!                    cloud scheme)
!! \param slmsk      sea/land mask array (sea:0,land:1,sea-ice:2)
!! \param xlon       grid longitude in radians,ok for both 0->2pi or 
!!                   -pi->+pi ranges
!! \param xlat       grid latitude in radians, default to pi/2->-pi/2
!!                   range, otherwise need to adjust in the called
!!                   subroutine
!! \param tsfc       surface temperature in K
!! \param snowd      snow depth water equivalent in mm (used when
!!                   control flag ialbflg=1)
!! \param sncovr     snow cover in fraction (used when contrl flag
!!                   ialbflg=1)
!! \param snoalb     maximum snow albedo in fraction (used when control
!!                   flag ialbflg=1)
!! \param zorl       surface roughness in cm
!! \param hprim      topographic standard deviation in m
!! \param alvsf      ialbflg=0: uv+visible albedo with strong cosz 
!!                              dependency (z=60)
!!\n                 ialbflg=1: uv+visible black sky albedo (z=60 degree)
!! \param alnsf      ialbflg=0: near IR albedo with strong cosz 
!!                              dependency (z=60)
!!\n                 ialbflg=1: near IR black sky albedo (z=60 degree)
!! \param alvwf      ialbflg=0: uv+visible albedo with weak cosz 
!!                              dependency (z=60)
!!\n                 ialbflg=1: uv+visible white sky albedo
!! \param alnwf      ialbflg=0: near IR albedo with weak cosz
!!                              dependency (z=60)
!!\n                 ialbflg=1: near IR white sky albedo
!! \param facsf      fractional coverage with strong cosz dependency
!! \param facwf      fractional coverage with weak cosz dependency
!! \param fice       fraction ice cover over open water grid
!! \param tisfc      surface temperature over ice cover in K
!! \param sinlat     sine of latitude for the model grid
!! \param coslat     cosine of latitude for the model grid
!! \param solhr      hour time after 00z at the current time-step
!! \param jdate      current forecast date and time (year, month, 
!!                   day,time-zone,hour, minute, second, mil-second)
!! \param solcon     solar constant (sun-earth distant adjusted) in \f$W/m^2\f$
!! \param cv         fraction of convective cloud cover 
!!                   (for diagnostic clouds only)
!! \param cvt,cvb    convective cloud top/bottom pressure in pa
!!                   (for diagnostic clouds only)
!! \param fcice      fraction of cloud ice content 
!!                   (for Ferrier microphysics scheme only)
!! \param frain      fraction of rain water 
!!                   (for Ferrier microphysics scheme only)
!! \param rrime      mass ratio of total to unrimed ice content 
!!                   (>= 1, for Ferrier microphysics scheme only)
!! \param flgmin     minimum large ice fraction 
!!                   (for Ferrier microphysics scheme only)
!! \param icsdsw,icsdlw    auxiliary cloud control arrays for radiations 
!!                   if isubcsw/isubclw (\ref physparam) are set to 2,
!!                   the arrays contains random seeds for the sub-column 
!!                   cloud overlap scheme, McICA, used in SW/LW radiations
!! \param ntcw       =0: no cloud condensate calculated; 
!!\n                 >0: tracer array location index for cloud condensate
!! \param ncld       only used when ntcw>0
!! \param ntoz       =0: use climatological ozone profile 
!!\n                 >0: use interactive ozone profile
!! \param NTRAC      number of tracers
!! \param NFXR       number of fields (second dimension) of I/O array fluxr
!! \param dtlw,dtsw  time durations for LW/SW radiation calls in second
!! \param lsswr,lslwr    logical control flags for SW/LW radiation calls
!! \param lssav      logical control flag for storing 3-d cloud field
!! \param IX,IM      horizontal dimension and number of used points
!! \param LM         vertical layer dimension
!! \param me         control flag for parallel process
!! \param lprnt      control flag for diagnostic printout
!! \param ipt        grid-point index for diagnostic printout (debugging)
!! \param kdt        time-step sequential number
!! \param deltaq     half width of pdf cloud uniform total water distribution
!!                   (for pdf cloud cover scheme)
!! \param sup        supersaturation in pdf cloud when t is very low
!!                   (for pdf cloud cover scheme)
!! \param cnvw       layer convective cloud water content
!!                   (for pdf cloud scheme)
!! \param cnvc       layer convective cloud cover
!!                   (for pdf cloud scheme)
!! \param htrsw      total sky SW heating rate in k/sec
!! \param topfsw     derived type, SW radiation fluxes at TOA, components:
!!                   (check module_radsw_parameters for definition)
!! \n          %upfxc   - total-sky upward SW flux at toa (\f$W/m^2\f$)
!! \n          %dnflx   - total-sky downward SW flux at toa (\f$W/m^2\f$)
!! \n          %upfx0   - clear-sky upward SW flux at toa (\f$W/m^2\f$)
!! \param sfcfsw     derived type, SW radiation fluxes at surface, components: 
!!                   (check module_radsw_parameters for definition)
!! \n          %upfxc   - total-sky upward SW flux at sfc (\f$W/m^2\f$)
!! \n          %dnfxc   - total-sky downward SW flux at sfc (\f$W/m^2\f$)
!! \n          %upfx0   - clear-sky upward SW flux at sfc (\f$W/m^2\f$)
!! \n          %dnfx0   - clear-sky downward SW flux at sfc (\f$W/m^2\f$)
!! \param dswcmp     downward surface SW spectral components:
!! \n          (:, 1)   - total-sky sfc downward nir direct flux
!! \n          (:, 2)   - total-sky sfc downward nir diffused flux
!! \n          (:, 3)   - total-sky sfc downward uv+vis direct flux
!! \n          (:, 4)   - total-sky sfc downward uv+vis diffused flux
!! \param uswcmp     upward surface SW spectral components:
!! \n          (:, 1)   - total-sky sfc upward nir direct flux
!! \n          (:, 2)   - total-sky sfc upward nir diffused flux
!! \n          (:, 3)   - total-sky sfc upward uv+vis direct flux
!! \n          (:, 4)   - total-sky sfc upward uv+vis diffused flux
!! \param sfalb      mean surface diffused albedo for SW radiation
!! \param coszen     mean cosine of solar zenith angle over radiation calling period
!! \param coszdg     daytime mean cosine of zenith angle over the radiation 
!!                   calling period
!! \param htrlw      total-sky LW heating rate in k/sec
!! \param topflw     derived type, LW radiation fluxes at TOA, component:
!!                   (check module_radlw_paramters for definition)
!! \n          %upfxc   - total-sky upward LW flux at toa (\f$W/m^2\f$)
!! \n          %upfx0   - clear-sky upward LW flux at toa (\f$W/m^2\f$)
!! \param sfcflw     derived type, LW radiation fluxes at surface, component:
!!                   (check module_radlw_paramters for definition)
!! \n          %upfxc   - total-sky upward LW flux at sfc (\f$W/m^2\f$)
!! \n          %upfx0   - clear-sky upward LW flux at sfc (\f$W/m^2\f$)
!! \n          %dnfxc   - total-sky downward LW flux at sfc (\f$W/m^2\f$)
!! \n          %dnfx0   - clear-sky downward LW flux at sfc (\f$W/m^2\f$)
!! \param tsflw      surface air temp during LW calculation call in K
!! \param semis      surface emissivity in fraction for LW radiation
!! \param cldcov     3-d cloud fraction
!! \param fluxr      array for saving time accumulated 2-d fields that are 
!!                   defined as:
!! \n          (:, 1)   - toa total-sky upward LW radiation flux
!! \n          (:, 2)   - toa total-sky upward SW radiation flux
!! \n          (:, 3)   - sfc total-sky upward SW radiation flux
!! \n          (:, 4)   - sfc total-sky downward SW radiation flux
!! \n          (:, 5)   - high domain cloud fraction
!! \n          (:, 6)   - mid  domain cloud fraction
!! \n          (:, 7)   - low  domain cloud fraction
!! \n          (:, 8)   - high domain mean cloud top pressure
!! \n          (:, 9)   - mid  domain mean cloud top pressure
!! \n          (:,10)   - low  domain mean cloud top pressure
!! \n          (:,11)   - high domain mean cloud base pressure
!! \n          (:,12)   - mid  domain mean cloud base pressure
!! \n          (:,13)   - low  domain mean cloud base pressure
!! \n          (:,14)   - high domain mean cloud top temperature
!! \n          (:,15)   - mid  domain mean cloud top temperature
!! \n          (:,16)   - low  domain mean cloud top temperature
!! \n          (:,17)   - total cloud fraction
!! \n          (:,18)   - boundary layer domain cloud fraction
!! \n          (:,19)   - sfc total-sky downward LW radiation flux
!! \n          (:,20)   - sfc total-sky upward LW radiation flux
!! \n          (:,21)   - sfc total-sky downward SW UV-B radiation flux
!! \n          (:,22)   - sfc clear-sky downward SW UV-B radiation flux
!! \n          (:,23)   - TOA incoming solar radiation flux
!! \n          (:,24)   - sfc UV+visible beam downward SW radiation flux
!! \n          (:,25)   - sfc UV+visible diffused downward SW radiation flux
!! \n          (:,26)   - sfc near-IR beam downward SW radiation flux
!! \n          (:,27)   - sfc near-IR diffused downward SW radiation flux
!! \n          (:,28)   - toa clear-sky upward LW radiation flux
!! \n          (:,29)   - toa clear-sky upward SW radiation flux
!! \n          (:,30)   - sfc clear-sky downward LW radiation flux
!! \n          (:,31)   - sfc clear-sky upward SW radiation flux
!! \n          (:,32)   - sfc clear-sky downward SW radiation flux
!! \n          (:,33)   - sfc clear-sky upward LW radiation flux
!! \n optional:
!! \n          (:,34)   - aerosol AOD at 550nm (all components)
!! \n          (:,35)   - aerosol AOD at 550nm for du component
!! \n          (:,36)   - aerosol AOD at 550nm for bc component
!! \n          (:,37)   - aerosol AOD at 550nm for oc component
!! \n          (:,38)   - aerosol AOD at 550nm for su component
!! \n          (:,39)   - aerosol AOD at 550nm for ss component
!! \param htrswb     spectral bands distributed total sky SW heating rate in k/sec
!! \param htrlwb     spectral bands distributed total sky LW heating rate in k/sec
!!
!> \section gen_grrad General Algorithm
!> @{
!-----------------------------------
      subroutine grrad                                                &                                            
     &     ( prsi,prsl,prslk,tgrs,qgrs,tracer,vvl,slmsk,              &      !  ---  inputs
     &       xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,           &  
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,          &  
     &       sinlat,coslat,solhr,jdate,solcon,                        &  
     &       cv,cvt,cvb,fcice,frain,rrime,flgmin,                     &  
     &       icsdsw,icsdlw, ntcw,ncld,ntoz, NTRAC,NFXR,               &  
     &       dtlw,dtsw, lsswr,lslwr,lssav, shoc_cld,lmfshal,lmfdeep2, &  
     &       IX, IM, LM, me, lprnt, ipt, kdt, deltaq,sup,cnvw,cnvc,   &  
     &       htrsw,topfsw,sfcfsw,dswcmp,uswcmp,sfalb,coszen,coszdg,   &      !  ---  outputs:
     &       htrlw,topflw,sfcflw,tsflw,semis,cldcov,                  &  
     &       fluxr                                                    &      !  ---  input/output:
     &,      htrlw0,htrsw0,htrswb,htrlwb                              &      !! ---  optional outputs:
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
!    this program is the driver of radiation calculation subroutines. * !
!    It sets up profile variables for radiation input, including      * !
!    clouds, surface albedos, atmospheric aerosols, ozone, etc.       * !
!                                                                     * !
!    usage:        call grrad                                         * !
!                                                                     * !
!    subprograms called:                                              * !
!                  setalb, setemis, setaer, getozn, getgases,         * !
!                  progcld1, progcld2, diagcds,                       * !
!                  swrad, lwrad, fpvs                                 * !
!                                                                     * !
!    attributes:                                                      * !
!      language:   fortran 90                                         * !
!      machine:    ibm-sp, sgi                                        * !
!                                                                     * !
!                                                                     * !
!  ====================  definition of variables  ====================  !
!                                                                       !
!    input variables:                                                   !
!      prsi  (IX,LM+1) : model level pressure in Pa                     !
!      prsl  (IX,LM)   : model layer mean pressure Pa                   !
!      prslk (IX,LM)   : exner function = (p/p0)**rocp                  !
!      tgrs  (IX,LM)   : model layer mean temperature in k              !
!      qgrs  (IX,LM)   : layer specific humidity in gm/gm               !
!      tracer(IX,LM,NTRAC):layer prognostic tracer amount/mixing-ratio  !
!                        incl: oz, cwc, aeros, etc.                     !
!      vvl   (IX,LM)   : layer mean vertical velocity in pa/sec         !
!      slmsk (IM)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
!      xlon  (IM)      : grid longitude in radians, ok for both 0->2pi  !
!                        or -pi -> +pi ranges                           !
!      xlat  (IM)      : grid latitude in radians, default to pi/2 ->   !
!                        -pi/2 range, otherwise adj in subr called      !
!      tsfc  (IM)      : surface temperature in k                       !
!      snowd (IM)      : snow depth water equivalent in mm              !
!      sncovr(IM)      : snow cover in fraction                         !
!      snoalb(IM)      : maximum snow albedo in fraction                !
!      zorl  (IM)      : surface roughness in cm                        !
!      hprim (IM)      : topographic standard deviation in m            !
!      alvsf (IM)      : mean vis albedo with strong cosz dependency    !
!      alnsf (IM)      : mean nir albedo with strong cosz dependency    !
!      alvwf (IM)      : mean vis albedo with weak cosz dependency      !
!      alnwf (IM)      : mean nir albedo with weak cosz dependency      !
!      facsf (IM)      : fractional coverage with strong cosz dependen  !
!      facwf (IM)      : fractional coverage with weak cosz dependency  !
!      fice  (IM)      : ice fraction over open water grid              !
!      tisfc (IM)      : surface temperature over ice fraction          !
!      sinlat(IM)      : sine of the grids' corresponding latitudes     !
!      coslat(IM)      : cosine of the grids' corresponding latitudes   !
!      solhr           : hour time after 00z at the t-stepe             !
!      jdate (8)       : current forecast date and time                 !
!                        (yr, mon, day, t-zone, hr, min, sec, mil-sec)  !
!      solcon          : solar constant (sun-earth distant adjusted)    !
!      cv    (IM)      : fraction of convective cloud                   !
!      cvt, cvb (IM)   : convective cloud top/bottom pressure in pa     !
!      fcice           : fraction of cloud ice  (in ferrier scheme)     !
!      frain           : fraction of rain water (in ferrier scheme)     !
!      rrime           : mass ratio of total to unrimed ice ( >= 1 )    !
!      flgmin          : minimim large ice fraction                     !
!      icsdsw/icsdlw   : auxiliary cloud control arrays passed to main  !
!           (IM)         radiations. if isubcsw/isubclw (input to init) !
!                        are set to 2, the arrays contains provided     !
!                        random seeds for sub-column clouds generators  !
!      lmfshal         : mass-flux shallow conv scheme flag             !
!      lmfdeep2        : scale-aware mass-flux deep conv scheme flag    !
!      ntcw            : =0 no cloud condensate calculated              !
!                        >0 array index location for cloud condensate   !
!      ncld            : only used when ntcw .gt. 0                     !
!      ntoz            : =0 climatological ozone profile                !
!                        >0 interactive ozone profile                   !
!      NTRAC           : dimension veriable for array oz                !
!      NFXR            : second dimension of input/output array fluxr   !
!      dtlw, dtsw      : time duration for lw/sw radiation call in sec  !
!      lsswr, lslwr    : logical flags for sw/lw radiation calls        !
!      lssav           : logical flag for store 3-d cloud field         !
!      IX,IM           : horizontal dimention and num of used points    !
!      LM              : vertical layer dimension                       !
!      me              : control flag for parallel process              !
!      lprnt           : control flag for diagnostic print out          !
!      ipt             : index for diagnostic printout point            !
!      kdt             : time-step number                               !
!      deltaq          : half width of uniform total water distribution !
!      sup             : supersaturation in pdf cloud when t is very low!
!      cnvw            : layer convective cloud water                   !
!      cnvc            : layer convective cloud cover                   !
!                                                                       !
!    output variables:                                                  !
!      htrsw (IX,LM)   : total sky sw heating rate in k/sec             !
!      topfsw(IM)      : sw radiation fluxes at toa, components:        !
!                      (check module_radsw_parameters for definition)   !
!       %upfxc           - total sky upward sw flux at toa (w/m**2)     !
!       %dnflx           - total sky downward sw flux at toa (w/m**2)   !
!       %upfx0           - clear sky upward sw flux at toa (w/m**2)     !
!      sfcfsw(IM)      : sw radiation fluxes at sfc, components:        !
!                      (check module_radsw_parameters for definition)   !
!       %upfxc           - total sky upward sw flux at sfc (w/m**2)     !
!       %dnfxc           - total sky downward sw flux at sfc (w/m**2)   !
!       %upfx0           - clear sky upward sw flux at sfc (w/m**2)     !
!       %dnfx0           - clear sky downward sw flux at sfc (w/m**2)   !
!      dswcmp(IX,4)    : dn sfc sw spectral components:                 !
!       ( :, 1)          -  total sky sfc downward nir direct flux      !
!       ( :, 2)          -  total sky sfc downward nir diffused flux    !
!       ( :, 3)          -  total sky sfc downward uv+vis direct flux   !
!       ( :, 4)          -  total sky sfc downward uv+vis diff flux     !
!      uswcmp(IX,4)    : up sfc sw spectral components:                 !
!       ( :, 1)          -  total sky sfc upward nir direct flux        !
!       ( :, 2)          -  total sky sfc upward nir diffused flux      !
!       ( :, 3)          -  total sky sfc upward uv+vis direct flux     !
!       ( :, 4)          -  total sky sfc upward uv+vis diff flux       !
!      sfalb (IM)      : mean surface diffused sw albedo                !
!      coszen(IM)      : mean cos of zenith angle over rad call period  !
!      coszdg(IM)      : daytime mean cosz over rad call period         !
!      htrlw (IX,LM)   : total sky lw heating rate in k/sec             !
!      topflw(IM)      : lw radiation fluxes at top, component:         !
!                        (check module_radlw_paramters for definition)  !
!       %upfxc           - total sky upward lw flux at toa (w/m**2)     !
!       %upfx0           - clear sky upward lw flux at toa (w/m**2)     !
!      sfcflw(IM)      : lw radiation fluxes at sfc, component:         !
!                        (check module_radlw_paramters for definition)  !
!       %upfxc           - total sky upward lw flux at sfc (w/m**2)     !
!       %upfx0           - clear sky upward lw flux at sfc (w/m**2)     !
!       %dnfxc           - total sky downward lw flux at sfc (w/m**2)   !
!       %dnfx0           - clear sky downward lw flux at sfc (w/m**2)   !
!      semis (IM)      : surface lw emissivity in fraction              !
!      cldcov(IX,LM)   : 3-d cloud fraction                             !
!      tsflw (IM)      : surface air temp during lw calculation in k    !
!                                                                       !
!    input and output variables:                                        !
!      fluxr (IX,NFXR) : to save time accumulated 2-d fields defined as:!
!                 1      - toa total sky upwd lw radiation flux         !
!                 2      - toa total sky upwd sw radiation flux         !
!                 3      - sfc total sky upwd sw radiation flux         !
!                 4      - sfc total sky dnwd sw radiation flux         !
!                 5      - high domain cloud fraction                   !
!                 6      - mid  domain cloud fraction                   !
!                 7      - low  domain cloud fraction                   !
!                 8      - high domain mean cloud top pressure          !
!                 9      - mid  domain mean cloud top pressure          !
!                10      - low  domain mean cloud top pressure          !
!                11      - high domain mean cloud base pressure         !
!                12      - mid  domain mean cloud base pressure         !
!                13      - low  domain mean cloud base pressure         !
!                14      - high domain mean cloud top temperature       !
!                15      - mid  domain mean cloud top temperature       !
!                16      - low  domain mean cloud top temperature       !
!                17      - total cloud fraction                         !
!                18      - boundary layer domain cloud fraction         !
!                19      - sfc total sky dnwd lw radiation flux         !
!                20      - sfc total sky upwd lw radiation flux         !
!                21      - sfc total sky dnwd sw uv-b radiation flux    !
!                22      - sfc clear sky dnwd sw uv-b radiation flux    !
!                23      - toa incoming solar radiation flux            !
!                24      - sfc vis beam dnwd sw radiation flux          !
!                25      - sfc vis diff dnwd sw radiation flux          !
!                26      - sfc nir beam dnwd sw radiation flux          !
!                27      - sfc nir diff dnwd sw radiation flux          !
!                28      - toa clear sky upwd lw radiation flux         !
!                29      - toa clear sky upwd sw radiation flux         !
!                30      - sfc clear sky dnwd lw radiation flux         !
!                31      - sfc clear sky upwd sw radiation flux         !
!                32      - sfc clear sky dnwd sw radiation flux         !
!                33      - sfc clear sky upwd lw radiation flux         !
!optional        34      - aeros opt depth at 550nm (all components)    !
!                35      - aeros opt depth at 550nm for du component    !
!                36      - aeros opt depth at 550nm for bc component    !
!                37      - aeros opt depth at 550nm for oc component    !
!                38      - aeros opt depth at 550nm for su component    !
!                39      - aeros opt depth at 550nm for ss component    !
!                                                                       !
!    optional output variables:                                         !
!      htrswb(IX,LM,NBDSW) : spectral band total sky sw heating rate    !
!      htrlwb(IX,LM,NBDLW) : spectral band total sky lw heating rate    !
!                                                                       !
!                                                                       !
!    definitions of internal variable arrays:                           !
!                                                                       !
!     1. fixed gases:         (defined in 'module_radiation_gases')     !
!          gasvmr(:,:,1)  -  co2 volume mixing ratio                    !
!          gasvmr(:,:,2)  -  n2o volume mixing ratio                    !
!          gasvmr(:,:,3)  -  ch4 volume mixing ratio                    !
!          gasvmr(:,:,4)  -  o2  volume mixing ratio                    !
!          gasvmr(:,:,5)  -  co  volume mixing ratio                    !
!          gasvmr(:,:,6)  -  cf11 volume mixing ratio                   !
!          gasvmr(:,:,7)  -  cf12 volume mixing ratio                   !
!          gasvmr(:,:,8)  -  cf22 volume mixing ratio                   !
!          gasvmr(:,:,9)  -  ccl4 volume mixing ratio                   !
!                                                                       !
!     2. cloud profiles:      (defined in 'module_radiation_clouds')    !
!                ---  for  prognostic cloud  ---                        !
!          clouds(:,:,1)  -  layer total cloud fraction                 !
!          clouds(:,:,2)  -  layer cloud liq water path                 !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud     !
!          clouds(:,:,4)  -  layer cloud ice water path                 !
!          clouds(:,:,5)  -  mean effective radius for ice cloud        !
!          clouds(:,:,6)  -  layer rain drop water path                 !
!          clouds(:,:,7)  -  mean effective radius for rain drop        !
!          clouds(:,:,8)  -  layer snow flake water path                !
!          clouds(:,:,9)  -  mean effective radius for snow flake       !
!                ---  for  diagnostic cloud  ---                        !
!          clouds(:,:,1)  -  layer total cloud fraction                 !
!          clouds(:,:,2)  -  layer cloud optical depth                  !
!          clouds(:,:,3)  -  layer cloud single scattering albedo       !
!          clouds(:,:,4)  -  layer cloud asymmetry factor               !
!                                                                       !
!     3. surface albedo:      (defined in 'module_radiation_surface')   !
!          sfcalb( :,1 )  -  near ir direct beam albedo                 !
!          sfcalb( :,2 )  -  near ir diffused albedo                    !
!          sfcalb( :,3 )  -  uv+vis direct beam albedo                  !
!          sfcalb( :,4 )  -  uv+vis diffused albedo                     !
!                                                                       !
!     4. sw aerosol profiles: (defined in 'module_radiation_aerosols')  !
!          faersw(:,:,:,1)-  sw aerosols optical depth                  !
!          faersw(:,:,:,2)-  sw aerosols single scattering albedo       !
!          faersw(:,:,:,3)-  sw aerosols asymmetry parameter            !
!                                                                       !
!     5. lw aerosol profiles: (defined in 'module_radiation_aerosols')  !
!          faerlw(:,:,:,1)-  lw aerosols optical depth                  !
!          faerlw(:,:,:,2)-  lw aerosols single scattering albedo       !
!          faerlw(:,:,:,3)-  lw aerosols asymmetry parameter            !
!                                                                       !
!     6. sw fluxes at toa:    (defined in 'module_radsw_main')          !
!        (topfsw_type -- derived data type for toa rad fluxes)          !
!          topfsw(:)%upfxc  -  total sky upward flux at toa             !
!          topfsw(:)%dnfxc  -  total sky downward flux at toa           !
!          topfsw(:)%upfx0  -  clear sky upward flux at toa             !
!                                                                       !
!     7. lw fluxes at toa:    (defined in 'module_radlw_main')          !
!        (topflw_type -- derived data type for toa rad fluxes)          !
!          topflw(:)%upfxc  -  total sky upward flux at toa             !
!          topflw(:)%upfx0  -  clear sky upward flux at toa             !
!                                                                       !
!     8. sw fluxes at sfc:    (defined in 'module_radsw_main')          !
!        (sfcfsw_type -- derived data type for sfc rad fluxes)          !
!          sfcfsw(:)%upfxc  -  total sky upward flux at sfc             !
!          sfcfsw(:)%dnfxc  -  total sky downward flux at sfc           !
!          sfcfsw(:)%upfx0  -  clear sky upward flux at sfc             !
!          sfcfsw(:)%dnfx0  -  clear sky downward flux at sfc           !
!                                                                       !
!     9. lw fluxes at sfc:    (defined in 'module_radlw_main')          !
!        (sfcflw_type -- derived data type for sfc rad fluxes)          !
!          sfcflw(:)%upfxc  -  total sky upward flux at sfc             !
!          sfcflw(:)%dnfxc  -  total sky downward flux at sfc           !
!          sfcflw(:)%dnfx0  -  clear sky downward flux at sfc           !
!                                                                       !
!! optional radiation outputs:                                          !
!!   10. sw flux profiles:    (defined in 'module_radsw_main')          !
!!       (profsw_type -- derived data type for rad vertical profiles)   !
!!         fswprf(:,:)%upfxc - total sky upward flux                    !
!!         fswprf(:,:)%dnfxc - total sky downward flux                  !
!!         fswprf(:,:)%upfx0 - clear sky upward flux                    !
!!         fswprf(:,:)%dnfx0 - clear sky downward flux                  !
!!                                                                      !
!!   11. lw flux profiles:    (defined in 'module_radlw_main')          !
!!       (proflw_type -- derived data type for rad vertical profiles)   !
!!         flwprf(:,:)%upfxc - total sky upward flux                    !
!!         flwprf(:,:)%dnfxc - total sky downward flux                  !
!!         flwprf(:,:)%upfx0 - clear sky upward flux                    !
!!         flwprf(:,:)%dnfx0 - clear sky downward flux                  !
!!                                                                      !
!!   12. sw sfc components:   (defined in 'module_radsw_main')          !
!!       (cmpfsw_type -- derived data type for component sfc fluxes)    !
!!         scmpsw(:)%uvbfc  -  total sky downward uv-b flux at sfc      !
!!         scmpsw(:)%uvbf0  -  clear sky downward uv-b flux at sfc      !
!!         scmpsw(:)%nirbm  -  total sky sfc downward nir direct flux   !
!!         scmpsw(:)%nirdf  -  total sky sfc downward nir diffused flux !
!!         scmpsw(:)%visbm  -  total sky sfc downward uv+vis direct flx !
!!         scmpsw(:)%visdf  -  total sky sfc downward uv+vis diff flux  !
!                                                                       !
!    external module variables:                                         !
!     ivflip           : control flag for in/out vertical indexing      !
!                        =0 index from toa to surface                   !
!                        =1 index from surface to toa                   !
!     icmphys          : cloud microphysics scheme control flag         !
!                        =1 zhao/carr/sundqvist microphysics scheme     !
!                        =2 brad ferrier microphysics scheme            !
!                        =3 zhao/carr/sundqvist microphysics +pdf cloud !
!                                                                       !
!    module variables:                                                  !
!     itsfc            : =0 use same sfc skin-air/ground temp           !
!                        =1 use diff sfc skin-air/ground temp (not yet) !
!                                                                       !
!  ======================  end of definitions  =======================  !
!
      implicit none

!  ---  inputs: (for rank>1 arrays, horizontal dimensioned by IX)
      integer,  intent(in) :: IX,IM, LM, NTRAC, NFXR, me,               
     &                        ntoz, ntcw, ncld, ipt, kdt
      integer,  intent(in) :: icsdsw(IM), icsdlw(IM), jdate(8)

      logical,  intent(in) :: lsswr, lslwr, lssav, lprnt, shoc_cld,    
     &                        lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(IX,LM+1), intent(in) ::  prsi

      real (kind=kind_phys), dimension(IX,LM),   intent(in) ::  prsl,   
     &       prslk, tgrs, qgrs, vvl, fcice, frain, rrime, deltaq, cnvw, 
     &       cnvc
      real (kind=kind_phys), dimension(IM), intent(in) :: flgmin
      real(kind=kind_phys), intent(in) :: sup

      real (kind=kind_phys), dimension(IM),      intent(in) ::  slmsk,  
     &       xlon, xlat, tsfc, snowd, zorl, hprim, alvsf, alnsf, alvwf, 
     &       alnwf, facsf, facwf, cv, cvt, cvb, fice, tisfc,            
     &       sncovr, snoalb, sinlat, coslat

      real (kind=kind_phys), intent(in) :: solcon, dtlw, dtsw, solhr,   
     &       tracer(IX,LM,NTRAC)

      real (kind=kind_phys), dimension(IX,LM),intent(inout):: cldcov

!  ---  outputs: (horizontal dimensioned by IX)
      real (kind=kind_phys), dimension(IX,LM),intent(out):: htrsw,htrlw

      real (kind=kind_phys), dimension(IX,4), intent(out) :: dswcmp,    
     &       uswcmp

      real (kind=kind_phys), dimension(IM),   intent(out):: tsflw,      
     &       sfalb, semis, coszen, coszdg

      type (topfsw_type), dimension(IM), intent(out) :: topfsw
      type (sfcfsw_type), dimension(IM), intent(out) :: sfcfsw

      type (topflw_type), dimension(IM), intent(out) :: topflw
      type (sfcflw_type), dimension(IM), intent(out) :: sfcflw

!  ---  variables are for both input and output:
      real (kind=kind_phys), intent(inout) :: fluxr(IX,NFXR)

!! ---  optional outputs:
      real (kind=kind_phys), dimension(IX,LM,NBDSW), optional,          
     &                       intent(out) :: htrswb
      real (kind=kind_phys), dimension(IX,LM,NBDLW), optional,         
     &                       intent(out) :: htrlwb
      real (kind=kind_phys), dimension(ix,lm), optional,                &
     &                       intent(out) :: htrlw0
      real (kind=kind_phys), dimension(ix,lm), optional,                &
     &                       intent(out) :: htrsw0

!  ---  local variables: (horizontal dimensioned by IM)
      real (kind=kind_phys), dimension(IM,LM+1+LTP):: plvl, tlvl

      real (kind=kind_phys), dimension(IM,LM+LTP)  :: plyr, tlyr, qlyr, &
     &       olyr, rhly, qstl, vvel, clw, prslk1, tem2da, tem2db, tvly

      real (kind=kind_phys), dimension(IM) :: tsfa, cvt1, cvb1, tem1d,  &
     &       sfcemis, tsfg, tskn

      real (kind=kind_phys), dimension(IM,LM+LTP,NF_CLDS) :: clouds
      real (kind=kind_phys), dimension(IM,LM+LTP,NF_VGAS) :: gasvmr
      real (kind=kind_phys), dimension(IM,       NF_ALBD) :: sfcalb
      real (kind=kind_phys), dimension(IM,       NSPC1)   :: aerodp     
      real (kind=kind_phys), dimension(IM,LM+LTP,NTRAC)   :: tracer1

      real (kind=kind_phys), dimension(IM,LM+LTP,NBDSW,NF_AESW)::faersw
      real (kind=kind_phys), dimension(IM,LM+LTP,NBDLW,NF_AELW)::faerlw

      real (kind=kind_phys), dimension(IM,LM+LTP) :: htswc
      real (kind=kind_phys), dimension(IM,LM+LTP) :: htlwc

      real (kind=kind_phys), dimension(IM,LM+LTP) :: gcice, grain, grime

!! ---  may be used for optional sw/lw outputs:
!!      take out "!!" as needed
      real (kind=kind_phys), dimension(IM,LM+LTP)   :: htsw0
!!    type (profsw_type),    dimension(IM,LM+1+LTP) :: fswprf
      type (cmpfsw_type),    dimension(IM)          :: scmpsw
      real (kind=kind_phys), dimension(IM,LM+LTP,NBDSW) :: htswb

      real (kind=kind_phys), dimension(IM,LM+LTP)   :: htlw0
!!    type (proflw_type),    dimension(IM,LM+1+LTP) :: flwprf
      real (kind=kind_phys), dimension(IM,LM+LTP,NBDLW) :: htlwb

      real (kind=kind_phys) :: raddt, es, qs, delt, tem0d, cldsa(IM,5)

      integer :: i, j, k, k1, lv, icec, itop, ibtc, nday, idxday(IM),   &
     &       mbota(IM,3), mtopa(IM,3), LP1, nb, LMK, LMP, kd, lla, llb, &
     &       lya, lyb, kt, kb

!  ---  for debug test use
!     real (kind=kind_phys) :: temlon, temlat, alon, alat
!     integer :: ipt
!     logical :: lprnt1

!
!===> ...  begin here
!
      LP1 = LM + 1               ! num of in/out levels

!  --- ...  set local /level/layer indexes corresponding to in/out variables

      LMK = LM + LTP             ! num of local layers
      LMP = LMK + 1              ! num of local levels

      if ( lextop ) then
        if ( ivflip == 1 ) then    ! vertical from sfc upward
          kd = 0                   ! index diff between in/out and local
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
          lla = LMK                ! local index at the 2nd level from top
          llb = LMP                ! local index at toa level
          lya = LM                 ! local index for the 2nd layer from top
          lyb = LP1                ! local index for the top layer
        else                       ! vertical from toa downward
          kd = 1                   ! index diff between in/out and local
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
          lla = 2                  ! local index at the 2nd level from top
          llb = 1                  ! local index at toa level
          lya = 2                  ! local index for the 2nd layer from top
          lyb = 1                  ! local index for the top layer
        endif                    ! end if_ivflip_block
      else
        kd = 0
        if ( ivflip == 1 ) then  ! vertical from sfc upward
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
        else                     ! vertical from toa downward
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
        endif                    ! end if_ivflip_block
      endif   ! end if_lextop_block

      raddt = min(dtsw, dtlw)

!  --- ...  for debug test
!     alon = 120.0
!     alat = 29.5
!     ipt = 0
!     do i = 1, IM
!       temlon = xlon(i) * 57.29578
!       if (temlon < 0.0) temlon = temlon + 360.0
!       temlat = xlat(i) * 57.29578
!       lprnt1 = abs(temlon-alon) < 1.1 .and. abs(temlat-alat) < 1.1
!       if ( lprnt1 ) then
!         ipt = i
!         exit
!       endif
!     enddo

!     print *,' in grrad : raddt=',raddt

      if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      else                              ! use diff sfc skin-air/ground temp
        do i = 1, IM
!!        tskn(i) = ta  (i)               ! not yet
!!        tsfg(i) = tg  (i)               ! not yet
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      endif

!> -# Prepare atmospheric profiles for radiation input.
!
!           convert pressure unit from pa to mb
      do k = 1, LM
        k1 = k + kd
        do i = 1, IM
!         plvl(i,k1)   = 10.0 * prsi(i,k)   ! cb (kpa) to mb (hpa)
!         plyr(i,k1)   = 10.0 * prsl(i,k)   ! cb (kpa) to mb (hpa)
          plvl(i,k1)   = 0.01 * prsi(i,k)   ! pa to mb (hpa)
          plyr(i,k1)   = 0.01 * prsl(i,k)   ! pa to mb (hpa)
          tlyr(i,k1)   = tgrs(i,k)
          prslk1(i,k1) = prslk(i,k)

!>  - Compute relative humidity.
!         es  = min( prsl(i,k), 0.001 * fpvs( tgrs(i,k) ) )   ! fpvs in pa
          es  = min( prsl(i,k),  fpvs( tgrs(i,k) ) )  ! fpvs and prsl in pa
          qs  = max( QMIN, eps * es / (prsl(i,k) + epsm1*es) )
          rhly(i,k1) = max( 0.0, min( 1.0, max(QMIN, qgrs(i,k))/qs ) )
          qstl(i,k1) = qs
        enddo
      enddo
      do j = 1, NTRAC
        do k = 1, LM
          k1 = k + kd
          do i = 1, IM
             tracer1(i,k1,j) = tracer(i,k,j)
          enddo
        enddo
      enddo

      do i = 1, IM
!       plvl(i,LP1+kd) = 10.0 * prsi(i,LP1)  ! cb (kpa) to mb (hpa
        plvl(i,LP1+kd) = 0.01 * prsi(i,LP1)  ! pa to mb (hpa)
      enddo

      if ( lextop ) then                 ! values for extra top layer
        do i = 1, IM
          plvl(i,llb) = prsmin
          if ( plvl(i,lla) <= prsmin ) plvl(i,lla) = 2.0*prsmin
          plyr(i,lyb)   = 0.5 * plvl(i,lla)
          tlyr(i,lyb)   = tlyr(i,lya)
!         prslk1(i,lyb) = (plyr(i,lyb)*0.001) ** rocp   ! plyr in hPa
          prslk1(i,lyb) = (plyr(i,lyb)*0.00001) ** rocp ! plyr in Pa

          rhly(i,lyb)   = rhly(i,lya)
          qstl(i,lyb)   = qstl(i,lya)
        enddo

        do j = 1, NTRAC
          do i = 1, IM
!  ---  note: may need to take care the top layer amount
             tracer1(i,lyb,j) = tracer1(i,lya,j)
          enddo
        enddo
      endif

!>  - Compute extra variables needed for Ferrier's microphysics

      if (icmphys == 2) then
        do k = 1, LM
          k1 = k + kd

          do i = 1, IM
            gcice(i,k1) = fcice(i,k)
            grain(i,k1) = frain(i,k)
            grime(i,k1) = rrime(i,k)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            gcice(i,lyb) = fcice(i,lya)
            grain(i,lyb) = frain(i,lya)
            grime(i,lyb) = rrime(i,lya)
          enddo
        endif
      endif   ! if_icmphys

!>  - Get layer ozone mass mixing ratio (if use ozone climatology data, 
!!    call getozn()).

      if (ntoz > 0) then            ! interactive ozone generation

        do k = 1, LMK
          do i = 1, IM
            olyr(i,k) = max( QMIN, tracer1(i,k,ntoz) )
          enddo
        enddo

      else                          ! climatological ozone

!     print *,' in grrad : calling getozn'
        call getozn                                                     &
!  ---  inputs:
     &     ( prslk1,xlat,                                               &
     &       IM, LMK,                                                   &
!  ---  outputs:
     &       olyr                                                       &
     &     )

      endif                            ! end_if_ntoz

!>  - Call coszmn(), to compute cosine of zenith angle.

      call coszmn                                                       &
!  ---  inputs:
     &     ( xlon,sinlat,coslat,solhr, IM, me,                          &
!  ---  outputs:
     &       coszen, coszdg                                             &
     &      )

!>  - Call getgases(), to set up non-prognostic gas volume mixing
!!    ratioes (gasvmr).
!  - gasvmr(:,:,1)  -  co2 volume mixing ratio
!  - gasvmr(:,:,2)  -  n2o volume mixing ratio
!  - gasvmr(:,:,3)  -  ch4 volume mixing ratio
!  - gasvmr(:,:,4)  -  o2  volume mixing ratio
!  - gasvmr(:,:,5)  -  co  volume mixing ratio
!  - gasvmr(:,:,6)  -  cf11 volume mixing ratio
!  - gasvmr(:,:,7)  -  cf12 volume mixing ratio
!  - gasvmr(:,:,8)  -  cf22 volume mixing ratio
!  - gasvmr(:,:,9)  -  ccl4 volume mixing ratio

!  --- ...  set up non-prognostic gas volume mixing ratioes

      call getgases                                                     &
!  ---  inputs:
     &    ( plvl, xlon, xlat,                                           &
     &      IM, LMK,                                                    &
!  ---  outputs:
     &      gasvmr                                                      &
     &     )

!>  - Get temperature at layer interface, and layer moisture.

      do k = 2, LMK
        do i = 1, IM
          tem2da(i,k) = log( plyr(i,k) )
          tem2db(i,k) = log( plvl(i,k) )
        enddo
      enddo

      if (ivflip == 0) then              ! input data from toa to sfc

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = 1.0
          tsfa  (i)   = tlyr(i,LMK)                  ! sfc layer air temp
          tlvl(i,1)   = tlyr(i,1)
          tlvl(i,LMP) = tskn(i)
        enddo

        do k = 1, LM
          k1 = k + kd

          do i = 1, IM
            qlyr(i,k1) = max( tem1d(i), qgrs(i,k) )
            tem1d(i)   = min( QME5, qlyr(i,k1) )
            tvly(i,k1) = tgrs(i,k) * (1.0 + fvirt*qlyr(i,k1)) ! virtual T (K)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 2, LMK
          do i = 1, IM
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

      else                               ! input data from sfc to toa

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,LMP) = tlyr(i,LMK)
        enddo

        do k = LM, 1, -1
          do i = 1, IM
            qlyr(i,k) = max( tem1d(i), qgrs(i,k) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = tgrs(i,k) * (1.0 + fvirt*qlyr(i,k)) ! virtual T (K)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 1, LMK-1
          do i = 1, IM
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

      endif                              ! end_if_ivflip

!>  - Check for daytime points for SW radiation.

      nday = 0
      do i = 1, IM
        if (coszen(i) >= 0.0001) then
          nday = nday + 1
          idxday(nday) = i
        endif
      enddo

!      write(0,*)' plvl=',plvl(ipt,1:65)
!      write(0,*)' plyr=',plyr(ipt,1:64)
!      write(0,*)' tlyr=',tlyr(ipt,1:64)
!      write(0,*)' tlvl=',tlvl(ipt,1:65)
!      write(0,*)' qlyr=',qlyr(ipt,1:10)*1000

!>  - Call module_radiation_aerosols::setaer(),to setup aerosols
!! property profile for radiation.

!check  print *,' in grrad : calling setaer '

      call setaer                                                       &
!  ---  inputs:
     &     ( plvl,plyr,prslk1,tvly,rhly,slmsk,tracer1,xlon,xlat,        &
     &       IM,LMK,LMP, lsswr,lslwr,                                   &
!  ---  outputs:
     &       faersw,faerlw,aerodp                                       &
     &     )

!>  - Obtain cloud information for radiation calculations
!!    (clouds,cldsa,mtopa,mbota)
!!\n   for  prognostic cloud:
!!    - For Zhao/Moorthi's prognostic cloud scheme,
!!      call module_radiation_clouds::progcld1()
!!    - For Ferrier's microphysics,
!!      call module_radiation_clouds::progcld2()
!!    - For Zhao/Moorthi's prognostic cloud+pdfcld,
!!      call module_radiation_clouds::progcld3()
!>  - If cloud condensate is not computed (ntcw=0), using the legacy 
!!   cloud scheme, compute cloud information based on Slingo's 
!!   diagnostic cloud scheme (call module_radiation_clouds::diagcld1())

!  --- ...  obtain cloud information for radiation calculations

      if (ntcw > 0) then                   ! prognostic cloud scheme

        do k = 1, LMK
          do i = 1, IM
            clw(i,k) = 0.0
          enddo

          do j = 1, ncld
            lv = ntcw + j - 1
            do i = 1, IM
              clw(i,k) = clw(i,k) + tracer1(i,k,lv)   ! cloud condensate amount
            enddo
          enddo
        enddo

        do k = 1, LMK
          do i = 1, IM
            if ( clw(i,k) < EPSQ ) clw(i,k) = 0.0
          enddo
        enddo
!
!  --- add suspended convective cloud water to grid-scale cloud water
!      only for cloud fraction & radiation computation
!      it is to enhance cloudiness due to suspended convec cloud water
!      for zhao/moorthi's (icmphys=1) & 
!          ferrier's (icmphys=2) microphysics schemes
!
        if (icmphys == 1 .or. icmphys == 2) then
          do k = 1, LMK
            do i = 1, IM
              clw(i,k) = clw(i,k) + cnvw(i,k)
            enddo
          enddo
        endif
!

        if (icmphys == 1) then           ! zhao/moorthi's prognostic cloud scheme

          call progcld1                                                 &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk, IM, LMK, LMP,                             &
     &       shoc_cld, lmfshal, lmfdeep2, cldcov(1:im,1:lm),            &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )

        elseif (icmphys == 2) then       ! ferrier's microphysics

!     print *,' in grrad : calling progcld2'
          call progcld2                                                 &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk, gcice,grain,grime,flgmin,                 &
     &       IM, LMK, LMP, lmfshal, lmfdeep2,                           &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )
!
        elseif(icmphys == 3) then      ! zhao/moorthi's prognostic cloud+pdfcld

          call progcld3                                                 &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &
     &       xlat,xlon,slmsk,                                           &
     &       im, lmk, lmp,                                              &
     &       deltaq, sup,kdt,me,                                        &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )

        endif                            ! end if_icmphys

      else                               ! diagnostic cloud scheme

        do i = 1, IM
!         cvt1(i) = 10.0 * cvt(i)
!         cvb1(i) = 10.0 * cvb(i)
          cvt1(i) = 0.01 * cvt(i)
          cvb1(i) = 0.01 * cvb(i)

        enddo

        do k = 1, LM
          k1 = k + kd

          do i = 1, IM
!           vvel(i,k1) = 10.0 * vvl(i,k)
            vvel(i,k1) = 0.01 * vvl(i,k)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            vvel(i,lyb) = vvel(i,lya)
          enddo
        endif

!  ---  compute diagnostic cloud related quantities

        call diagcld1                                                   &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,rhly,vvel,cv,cvt1,cvb1,                     &
     &       xlat,xlon,slmsk,                                           &
     &       IM, LMK, LMP,                                              &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )

      endif                                ! end_if_ntcw

!  --- ...  start radiation calculations 
!           remember to set heating rate unit to k/sec!
!> -# Start SW radiation calculations
      if (lsswr) then

!>  - Call module_radiation_surface::setalb() to setup surface albedo.
!!  for SW radiation.


        call setalb                                                     &
!  ---  inputs:
     &     ( slmsk,snowd,sncovr,snoalb,zorl,coszen,tsfg,tsfa,hprim,     &
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
     &       IM,                                                        &
!  ---  outputs:
     &       sfcalb                                                     &
     &     )


        do i = 1, IM
          sfalb(i) = max(0.01, 0.5 * (sfcalb(i,2) + sfcalb(i,4)))
        enddo

        if (nday > 0) then

!>  - Call module_radsw_main::swrad(), to compute SW heating rates and 
!!   fluxes.
!     print *,' in grrad : calling swrad'

          if ( present(htrswb) .and. present(htrsw0)) then

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      HSW0=htsw0,FLXPRF=fswprf                                   &
     &,      hsw0=htsw0,hswb=htswb,fdncmp=scmpsw                        &
     &     )

            do k = 1, LM
              k1 = k + kd

              do j = 1, NBDSW
                do i = 1, IM
                  htrswb(i,k,j) = htswb(i,k1,j)
                enddo
              enddo
            enddo

          else if ( present(htrswb) .and. .not. present(htrsw0) ) then

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      hsw0=htsw0,flxprf=fswprf                                   &
     &,      hswb=htswb,fdncmp=scmpsw                                   &
     &     )

            do k = 1, lm
              k1 = k + kd
              do j = 1, nbdsw
                do i = 1, im
                  htrswb(i,k,j) = htswb(i,k1,j)
                enddo
              enddo
            enddo

          else if ( present(htrsw0) .and. .not. present(htrswb) ) then

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      hsw0=htsw0,flxprf=fswprf                                   &
     &,      hsw0=htsw0,fdncmp=scmpsw                                   &
     &     )

          else

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      HSW0=htsw0,FLXPRF=fswprf,HSWB=htswb                        &
     &,      FDNCMP=scmpsw                                              &
     &     )

          endif

          do k = 1, LM
            k1 = k + kd

            do i = 1, IM
              htrsw(i,k) = htswc(i,k1)
            enddo
          enddo
          if (present(htrsw0)) then
             do k = 1, lm
               k1 = k + kd
               do i = 1, im
                 htrsw0(i,k) = htsw0(i,k1)
               enddo
             enddo
          endif

!  --- surface down and up spectral component fluxes
!>  - Save two spectral bands' surface downward and upward fluxes for 
!!    output.

          do i = 1, IM
            dswcmp(i,1) = scmpsw(i)%nirbm
            dswcmp(i,2) = scmpsw(i)%nirdf
            dswcmp(i,3) = scmpsw(i)%visbm
            dswcmp(i,4) = scmpsw(i)%visdf

            uswcmp(i,1) = scmpsw(i)%nirbm * sfcalb(i,1)
            uswcmp(i,2) = scmpsw(i)%nirdf * sfcalb(i,2)
            uswcmp(i,3) = scmpsw(i)%visbm * sfcalb(i,3)
            uswcmp(i,4) = scmpsw(i)%visdf * sfcalb(i,4)
          enddo

        else                   ! if_nday_block

          do k = 1, LM
            do i = 1, IM
              htrsw(i,k) = 0.0
            enddo
          enddo

          sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          topfsw = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

          do k = 1, 4
          do i = 1, IM
            dswcmp(i,k) = 0.0
            uswcmp(i,k) = 0.0
          enddo
          enddo

!! ---  optional:
!!        fswprf= profsw_type( 0.0, 0.0, 0.0, 0.0 )

          if ( present(htrswb) ) then
            do j = 1, NBDSW
              do k = 1, LM
                do i = 1, IM
                  htrswb(i,k,j) = 0.0
                enddo
              enddo
            enddo
          endif
          if ( present(htrsw0) ) then
              do k = 1, lm
                do i = 1, im
                  htrsw0(i,k) = 0.0
                enddo
              enddo
          endif

        endif                  ! end_if_nday

      endif                                ! end_if_lsswr

!      write(0,*)' htrsw=',htrsw(ipt,1:64)*86400
!> -# Start LW radiation calculations
      if (lslwr) then

!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.

        call setemis                                                    &
!  ---  inputs:
     &     ( xlon,xlat,slmsk,snowd,sncovr,zorl,tsfg,tsfa,hprim,         &
     &       IM,                                                        &
!  ---  outputs:
     &       sfcemis                                                    &
     &     )

!>  - Call module_radlw_main::lwrad(), to compute LW heating rates and 
!!    fluxes.
!     print *,' in grrad : calling lwrad'

        if ( present(htrlwb) .and. present(htrlw0) ) then

          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      HLW0=htlw0,FLXPRF=flwprf                                   &
     &,      hlw0=htlw0                                                 &
     &,      hlwb=htlwb                                                 &
     &     )

          do k = 1, LM
            k1 = k + kd

            do j = 1, NBDLW
              do i = 1, IM
                htrlwb(i,k,j) = htlwb(i,k1,j)
              enddo
            enddo
          enddo

        else if ( present(htrlwb) .and. .not. present(htrlw0) ) then

          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      hlw0=htlw0,flxprf=flwprf                                   &
     &,      hlwb=htlwb                                                 &
     &     )

          do k = 1, lm
            k1 = k + kd

            do j = 1, nbdlw
              do i = 1, im
                htrlwb(i,k,j) = htlwb(i,k1,j)
              enddo
            enddo
          enddo
        else if ( present(htrlw0) .and. .not. present(htrlwb) ) then

          !print *,'call lwrad saving clear sky component'
          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      hlw0=htlw0,flxprf=flwprf                                   &
     &,      hlw0=htlw0                                                 &
     &     )

        else

          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       IM, LMK, LMP, lprnt,                                       &
!  ---  outputs:
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      HLW0=htlw0,FLXPRF=flwprf,HLWB=htlwb                        &
     &     )

        endif

!> -# Save calculation results
        do i = 1, IM
          semis (i) = sfcemis(i)
!>  - Save surface air temp for diurnal adjustment at model t-steps
          tsflw (i) = tsfa(i)
        enddo

        do k = 1, LM
          k1 = k + kd

          do i = 1, IM
            htrlw(i,k) = htlwc(i,k1)
          enddo
        enddo

        if (present(htrlw0)) then
           do k = 1, lm
             k1 = k + kd
             do i = 1, im
               htrlw0(i,k) = htlw0(i,k1)
             enddo
           enddo
        endif

      endif                                ! end_if_lslwr

!>  - For time averaged output quantities (including total-sky and 
!!    clear-sky SW and LW fluxes at TOA and surface; conventional 
!!    3-domain cloud amount, cloud top and base pressure, and cloud top 
!!    temperature; aerosols AOD, etc.), store computed results in 
!!    corresponding slots of array fluxr with appropriate time weights.

!  --- ...  collect the fluxr data for wrtsfc

      if (lssav) then

         if ( lsswr ) then
          do i = 1, IM
             fluxr(i,34) = fluxr(i,34) + dtsw*aerodp(i,1)  ! total aod at 550nm
             fluxr(i,35) = fluxr(i,35) + dtsw*aerodp(i,2)  ! DU aod at 550nm
             fluxr(i,36) = fluxr(i,36) + dtsw*aerodp(i,3)  ! BC aod at 550nm
             fluxr(i,37) = fluxr(i,37) + dtsw*aerodp(i,4)  ! OC aod at 550nm
             fluxr(i,38) = fluxr(i,38) + dtsw*aerodp(i,5)  ! SU aod at 550nm
             fluxr(i,39) = fluxr(i,39) + dtsw*aerodp(i,6)  ! SS aod at 550nm
          enddo
        endif
!  ---  save lw toa and sfc fluxes

        if (lslwr) then
          do i = 1, IM
!  ---  lw total-sky fluxes
            fluxr(i,1 ) = fluxr(i,1 ) + dtlw * topflw(i)%upfxc   ! total sky top lw up
            fluxr(i,19) = fluxr(i,19) + dtlw * sfcflw(i)%dnfxc   ! total sky sfc lw dn
            fluxr(i,20) = fluxr(i,20) + dtlw * sfcflw(i)%upfxc   ! total sky sfc lw up
!  ---  lw clear-sky fluxes
            fluxr(i,28) = fluxr(i,28) + dtlw * topflw(i)%upfx0   ! clear sky top lw up
            fluxr(i,30) = fluxr(i,30) + dtlw * sfcflw(i)%dnfx0   ! clear sky sfc lw dn
            fluxr(i,33) = fluxr(i,33) + dtlw * sfcflw(i)%upfx0   ! clear sky sfc lw up
          enddo
        endif

!  ---  save sw toa and sfc fluxes with proper diurnal sw wgt. coszen=mean cosz over daylight
!       part of sw calling interval, while coszdg= mean cosz over entire interval

        if (lsswr) then
          do i = 1, IM
            if (coszen(i) > 0.) then
!  ---                                  sw total-sky fluxes
!                                       -------------------
              tem0d = dtsw * coszdg(i)  / coszen(i)
              fluxr(i,2 ) = fluxr(i,2)  + topfsw(i)%upfxc * tem0d  ! total sky top sw up
              fluxr(i,3 ) = fluxr(i,3)  + sfcfsw(i)%upfxc * tem0d  ! total sky sfc sw up
              fluxr(i,4 ) = fluxr(i,4)  + sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
!  ---                                  sw uv-b fluxes
!                                       --------------
              fluxr(i,21) = fluxr(i,21) + scmpsw(i)%uvbfc * tem0d  ! total sky uv-b sw dn
              fluxr(i,22) = fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d  ! clear sky uv-b sw dn
!  ---                                  sw toa incoming fluxes
!                                       ----------------------
              fluxr(i,23) = fluxr(i,23) + topfsw(i)%dnfxc * tem0d  ! top sw dn
!  ---                                  sw sfc flux components
!                                       ----------------------
              fluxr(i,24) = fluxr(i,24) + scmpsw(i)%visbm * tem0d  ! uv/vis beam sw dn
              fluxr(i,25) = fluxr(i,25) + scmpsw(i)%visdf * tem0d  ! uv/vis diff sw dn
              fluxr(i,26) = fluxr(i,26) + scmpsw(i)%nirbm * tem0d  ! nir beam sw dn
              fluxr(i,27) = fluxr(i,27) + scmpsw(i)%nirdf * tem0d  ! nir diff sw dn
!  ---                                  sw clear-sky fluxes
!                                       -------------------
              fluxr(i,29) = fluxr(i,29) + topfsw(i)%upfx0 * tem0d  ! clear sky top sw up
              fluxr(i,31) = fluxr(i,31) + sfcfsw(i)%upfx0 * tem0d  ! clear sky sfc sw up
              fluxr(i,32) = fluxr(i,32) + sfcfsw(i)%dnfx0 * tem0d  ! clear sky sfc sw dn
            endif
          enddo
        endif

!  ---  save total and boundary layer clouds

        if (lsswr .or. lslwr) then
          do i = 1, IM
            fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)
            fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)
          enddo

!  ---  save cld frac,toplyr,botlyr and top temp, note that the order
!       of h,m,l cloud is reversed for the fluxr output.
!  ---  save interface pressure (pa) of top/bot

          do j = 1, 3
            do i = 1, IM
              tem0d = raddt * cldsa(i,j)
              itop  = mtopa(i,j) - kd
              ibtc  = mbota(i,j) - kd
              fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
              fluxr(i,11-j) = fluxr(i,11-j) + tem0d * prsi(i,itop+kt)
              fluxr(i,14-j) = fluxr(i,14-j) + tem0d * prsi(i,ibtc+kb)
              fluxr(i,17-j) = fluxr(i,17-j) + tem0d * tgrs(i,itop)
            enddo
          enddo
        endif

        if (.not. shoc_cld) then
          do k = 1, LM
            k1 = k + kd
            do i = 1, IM
              cldcov(i,k) = clouds(i,k1,1)
            enddo
          enddo
        endif

!  ---  save optional vertically integrated aerosol optical depth at
!       wavelenth of 550nm aerodp(:,1), and other optional aod for
!       individual species aerodp(:,2:NSPC1)

!       if ( laswflg ) then
!         if ( NFXR > 33 ) then
!           do i = 1, IM
!             fluxr(i,34) = fluxr(i,34) + dtsw*aerodp(i,1)  ! total aod at 550nm (all species)
!           enddo

!           if ( lspcodp ) then
!             do j = 2, NSPC1
!               k = 33 + j

!               do i = 1, IM
!                 fluxr(i,k) = fluxr(i,k) + dtsw*aerodp(i,j) ! aod at 550nm for indiv species
!               enddo
!             enddo
!           endif     ! end_if_lspcodp
!         else
!           print *,'  !Error! Need to increase array fluxr size NFXR ',&
!    &              ' to be able to output aerosol optical depth'
!           stop
!         endif     ! end_if_nfxr
!       endif       ! end_if_laswflg

      endif                                ! end_if_lssav
!
      return
!...................................
      end subroutine grrad
!-----------------------------------


!
!> @}
!........................................!
      end module module_radiation_driver !
!========================================!
!> @}
!>@}
