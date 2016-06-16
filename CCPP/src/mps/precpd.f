!> \file precpd.f
!! This file contains the subroutine that calculates precipitation processes from suspended cloud water/ice

!> \ingroup MPscheme
!> \defgroup precip Precipitation (snow or rain) Production
!! @{
!> This subroutine computes the conversion from condensation to precipitation (snow or rain) or
!! evaporation of rain.
!! \brief The parameterization of precipitation is required in order to remove water from the atmosphere
!! and transport it to the ground. In the scheme discussed here, simplifications in the precipitation parameterization
!! are used due to computational limitations required by operational NWP models. First, consideration of particle
!! size and shape can be avoided by using the bulk
!! parameterization method introduced by Kessler (1969) \cite kessler_1969. Second, only two types of precipitation
!! , rain and snow, are considered in this scheme. Third, only the most important microphysical
!! processes associated with the formation of rain and snow are included (production of rain from
!! cloud water, production of snow from cloud ice, melting of snow to form rain below the
!! freezing level, and the evaporation of precipitation). Finally, the fourth simplification
!! is that precipitation is diagnostically calculated directly from the cloud mixing ratio.
!! \param[in] im        horizontal number of used pts
!! \param[in] ix        horizontal dimension
!! \param[in] km        vertical layer dimension
!! \param[in] dt        time step in seconds
!! \param[in] del       pressure layer thickness (bottom to top)
!! \param[in] prsl      pressure values for model layers (bottom to top)
!! \param[in,out] q         specific humidity (updated in the code)
!! \param[in,out] cwm       condensate mixing ratio (updated in the code)
!! \param[in,out] t         temperature (updated in the code)
!! \param[out] rn        precipitation over one time-step dt (m/dt)
!! \param[out] sr        "snow ratio", ratio of snow to total precipitation
!! \param[out] rainp    rainwater path
!! \param[in] u00k      the critical value of relative humidity for large-scale condensation
!! \param[in] psautco   auto conversion coeff from ice to snow
!! \n                   = 4.0E-4; defined in module_MP_GFS.F90
!! \param[in] prautco   auto conversion coeff from cloud to rain
!! \n                   = 1.0E-4; defined in module_MP_GFS.F90
!! \param[in] evpco     coeff for evaporation of largescale rain
!! \n                   = 2.0E-5; defined in module_MP_GFS.F90
!! \param[in] wminco    water and ice minimum threshold to conversion from condensate to precipitation
!! \n                   = \1.0E-5, 1.0E-5\; defined in module_MP_GFS.F90
!! \param[in] lprnt     logical print flag
!! \param[in] jpr       check print point for debugging
!> \section general General Algorithm
!> @{
       subroutine precpd (im,ix,km,dt,del,prsl,q,cwm,t,rn,sr
     &,                   rainp,u00k,psautco,prautco,evpco,wminco
     &,                   lprnt,jpr)
!
!
!     ******************************************************************
!     *                                                                *
!     *           subroutine for precipitation processes               *
!     *           from suspended cloud water/ice                       *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  originally created by  q. zhao                jan. 1995       *
!     *                         -------                                *
!     *  modified and rewritten by shrinivas moorthi   oct. 1998       *
!     *                            -----------------                   *
!     *  and                       hua-lu pan                          *
!     *                            ----------                          *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     *  zhao and carr (1997), monthly weather review (august)         *
!     *  sundqvist et al., (1989) monthly weather review. (august)     *
!     *  chuang 2013, modify sr to define frozen precipitation fraction*
!     ******************************************************************
!
!     in this code vertical indexing runs from surface to top of the
!     model
!
!     argument list:
!     --------------
!       im         : inner dimension over which calculation is made
!       ix         : maximum inner dimension
!       km         : number of vertical levels
!       dt         : time step in seconds
!       del(km)    : pressure layer thickness (bottom to top)
!       prsl(km)   : pressure values for model layers (bottom to top)
!       q(ix,km)   : specific humidity (updated in the code)
!       cwm(ix,km) : condensate mixing ratio (updated in the code)
!       t(ix,km)   : temperature       (updated in the code)
!       rn(im)     : precipitation over one time-step dt (m/dt)
!old      sr(im)     : index (=-1 snow, =0 rain/snow, =1 rain)
!new    sr(im)     : "snow ratio", ratio of snow to total precipitation
!       cll(ix,km) : cloud cover
!hchuang rn(im) unit in m per time step
!        precipitation rate conversion 1 mm/s = 1 kg/m2/s
!
      use machine , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g, hvap => con_hvap, hfus => con_hfus
     &,             ttp => con_ttp, cp => con_cp
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!     include 'constant.h'
!
      real (kind=kind_phys) g,      h1,    h1000
     &,                     d00
     &,                     elwv,   eliv,  row
     &,                     epsq,   eliw
     &,                     rcp,    rrow
       parameter (g=grav,         h1=1.e0,     h1000=1000.0
     &,           d00=0.e0
     &,           elwv=hvap,      eliv=hvap+hfus,   row=1.e3
     &,           epsq=2.e-12
     &,           eliw=eliv-elwv, rcp=h1/cp,   rrow=h1/row)
!
      real(kind=kind_phys), parameter :: cons_0=0.0,     cons_p01=0.01
     &,                                  cons_20=20.0
     &,                                  cons_m30=-30.0, cons_50=50.0
!
      integer im, ix, km, jpr
      real (kind=kind_phys) q(ix,km),   t(ix,km),    cwm(ix,km)
     &,                                 del(ix,km),  prsl(ix,km) !    &,                     cll(im,km), del(ix,km),  prsl(ix,km)
     &,                     rn(im),      sr(im)
     &,                     dt
     &,                     rainp(im,km), rnp(im), !hchuang code change [+1l] : add record to record information in vertical in addition to total column precrl
     &                      psautco(im), prautco(im), evpco, wminco(2)
!
!
      real (kind=kind_phys) err(im),      ers(im),     precrl(im)
     &,                     precsl(im),   precrl1(im), precsl1(im)
     &,                     rq(im),       condt(im)
     &,                     conde(im),    rconde(im),  tmt0(im)
     &,                     wmin(im,km),  wmink(im),   pres(im)
     &,                     wmini(im,km), ccr(im)
     &,                     tt(im),       qq(im),      ww(im)
     &,                     u00k(im,km)
     &,                     zaodt
       real (kind=kind_phys) cclim(km)
!
      integer iw(im,km), ipr(im), iwl(im),     iwl1(im)
!
       logical comput(im)
       logical lprnt
!
      real (kind=kind_phys) ke,   rdt,  us, climit, cws, csm1
     &,                     crs1, crs2, cr, aa2,     dtcp,   c00, cmr
     &,                     tem,  c1,   c2, wwn
!    &,                     tem,  c1,   c2, u00b,    u00t,   wwn
     &,                     precrk, precsk, pres1,   qk,     qw,  qi
     &,                     qint, fiw, wws, cwmk, expf
     &,                     psaut, psaci, amaxcm, tem1, tem2
     &,                     tmt0k, psm1, psm2, ppr
     &,                     rprs,  erk,   pps, sid, rid, amaxps
     &,                     praut, fi, qc, amaxrq, rqkll
      integer i, k, ihpr, n
!
!-----------------------preliminaries ---------------------------------
!
!     do k=1,km
!       do i=1,im
!         cll(i,k) = 0.0
!       enddo
!     enddo
!
      rdt     = h1 / dt
!     ke      = 2.0e-5  ! commented on 09/10/99  -- opr value
!     ke      = 2.0e-6
!     ke      = 1.0e-5
!!!   ke      = 5.0e-5
!!    ke      = 7.0e-5
      ke      = evpco
!     ke      = 7.0e-5
      us      = h1
      climit  = 1.0e-20
      cws     = 0.025
!
      zaodt   = 800.0 * rdt
!
      csm1    = 5.0000e-8   * zaodt
      crs1    = 5.00000e-6  * zaodt
      crs2    = 6.66600e-10 * zaodt
      cr      = 5.0e-4      * zaodt
      aa2     = 1.25e-3     * zaodt
!
      ke      = ke * sqrt(rdt)
!     ke      = ke * sqrt(zaodt)
!
      dtcp    = dt * rcp
!
!     c00 = 1.5e-1 * dt
!     c00 = 10.0e-1 * dt
!     c00 = 3.0e-1 * dt          !05/09/2000
!     c00 = 1.0e-4 * dt          !05/09/2000
!     c00 = prautco * dt         !05/09/2000
      cmr = 1.0 / 3.0e-4
!     cmr = 1.0 / 5.0e-4
!     c1  = 100.0
      c1  = 300.0
      c2  = 0.5
!
!
!--------calculate c0 and cmr using lc at previous step-----------------
!
      do k=1,km
        do i=1,im
          tem   = (prsl(i,k)*0.00001)
!         tem   = sqrt(tem)
          iw(i,k)    = 0.0
!         wmin(i,k)  = 1.0e-5 * tem
!         wmini(i,k) = 1.0e-5 * tem       ! testing for ras
!

          wmin(i,k)  = wminco(1) * tem
          wmini(i,k) = wminco(2) * tem


          rainp(i,k) = 0.0

        enddo
      enddo
      do i=1,im
!       c0(i)  = 1.5e-1
!       cmr(i) = 3.0e-4
!
        iwl1(i)    = 0
        precrl1(i) = d00
        precsl1(i) = d00
        comput(i)  = .false.
        rn(i)      = d00
        sr(i)      = d00
        ccr(i)     = d00
!
        rnp(i)     = d00
      enddo
!> -# Select columns where rain can be produced.
!------------select columns where rain can be produced--------------
      do k=1, km-1
        do i=1,im
          tem = min(wmin(i,k), wmini(i,k))
          if (cwm(i,k) > tem) comput(i) = .true.
        enddo
      enddo
      ihpr = 0
      do i=1,im
        if (comput(i)) then
           ihpr      = ihpr + 1
           ipr(ihpr) = i
        endif
      enddo
!***********************************************************************
!-----------------begining of precipitation calculation-----------------
!***********************************************************************
!     do k=km-1,2,-1
      do k=km,1,-1
        do n=1,ihpr
          precrl(n) = precrl1(n)
          precsl(n) = precsl1(n)
          err  (n)  = d00
          ers  (n)  = d00
          iwl  (n)  = 0
!
          i         = ipr(n)
          tt(n)     = t(i,k)
          qq(n)     = q(i,k)
          ww(n)     = cwm(i,k)
          wmink(n)  = wmin(i,k)
          pres(n)   = prsl(i,k)
!
          precrk = max(cons_0,    precrl1(n))
          precsk = max(cons_0,    precsl1(n))
          wwn    = max(ww(n), climit)
!         if (wwn .gt. wmink(n) .or. (precrk+precsk) .gt. d00) then
          if (wwn > climit .or. (precrk+precsk) > d00) then
            comput(n) = .true.
          else
            comput(n) = .false.
          endif
        enddo
!
!       es(1:ihpr) = fpvs(tt(1:ihpr))
        do n=1,ihpr
          if (comput(n)) then
            i = ipr(n)
            conde(n)  = (dt/g) * del(i,k)
            condt(n)  = conde(n) * rdt
            rconde(n) = h1 / conde(n)
            qk        = max(epsq,  qq(n))
            tmt0(n)   = tt(n) - 273.16
            wwn       = max(ww(n), climit)
!
!           pl = pres(n) * 0.01
!           call qsatd(tt(n), pl, qc)
!           rq(n) = max(qq(n), epsq) / max(qc, 1.0e-10)
!           rq(n) = max(1.0e-10, rq(n))           ! -- relative humidity---
!
!  the global qsat computation is done in pa
            pres1   = pres(n)
!           qw      = es(n)
            qw      = min(pres1, fpvs(tt(n)))
            qw      = eps * qw / (pres1 + epsm1 * qw)
            qw      = max(qw,epsq)
!
!           tmt15 = min(tmt0(n), cons_m15)
!           ai    = 0.008855
!           bi    = 1.0
!           if (tmt0(n) .lt. -20.0) then
!             ai = 0.007225
!             bi = 0.9674
!           endif
!           qi   = qw * (bi + ai*min(tmt0(n),cons_0))
!           qint = qw * (1.-0.00032*tmt15*(tmt15+15.))
!
            qi   = qw
            qint = qw
!           if (tmt0(n).le.-40.) qint = qi
!
!-------------------ice-water id number iw------------------------------
            if(tmt0(n) < -15.) then
               fi = qk - u00k(i,k)*qi
               if(fi > d00 .or. wwn > climit) then
                  iwl(n) = 1
               else
                  iwl(n) = 0
               endif
!           endif
            elseif (tmt0(n) >= 0.) then
               iwl(n) = 0
!
!           if(tmt0(n).lt.0.0.and.tmt0(n).ge.-15.0) then
            else
              iwl(n) = 0
              if(iwl1(n) == 1 .and. wwn > climit) iwl(n) = 1
            endif
!
!           if(tmt0(n).ge.0.) then
!              iwl(n) = 0
!           endif
!----------------the satuation specific humidity------------------------
            fiw   = float(iwl(n))
            qc    = (h1-fiw)*qint + fiw*qi
!----------------the relative humidity----------------------------------
            if(qc <= 1.0e-10) then
               rq(n) = d00
            else
               rq(n) = qk / qc
            endif
!----------------cloud cover ratio ccr----------------------------------
            if(rq(n) < u00k(i,k)) then
                   ccr(n) = d00
            elseif(rq(n) >= us) then
                   ccr(n) = us
            else
                 rqkll  = min(us,rq(n))
                 ccr(n) = h1-sqrt((us-rqkll)/(us-u00k(i,k)))
            endif
!
          endif
        enddo
!-------------------ice-water id number iwl------------------------------
!       do n=1,ihpr
!         if (comput(n) .and.  (ww(n) .gt. climit)) then
!           if (tmt0(n) .lt. -15.0
!    *         .or. (tmt0(n) .lt. 0.0 .and. iwl1(n) .eq. 1))
!    *                                      iwl(n) = 1
!             cll(ipr(n),k) = 1.0                           ! cloud cover!
!             cll(ipr(n),k) = min(1.0, ww(n)*cclim(k))      ! cloud cover!
!         endif
!       enddo
!
!> -# Precipitation production by auto conversion and accretion
!!    - The autoconversion of cloud ice to snow (\f$P_{saut}\f$) is simulated
!! using the equation from Lin et al. (1983) \cite lin_et_al_1983
!!\f[
!!   P_{saut}=a_{1}(m-m_{i0})
!!\f]
!! where \f$m_{i0}\f$ is the threshold of cloud ice mixing ratio for production of snow
!! from cloud ice and is set to a value of \f$1.0\times 10^{-5} (kg/kg)\f$. Since snow
!! production in this process is caused by the increase in size of cloud ice particles due
!! to depositional growth and aggregation of small ice particles, \f$P_{saut}\f$ is
!! a function of temperature as determined by coefficient \f$a_{1}\f$, given by
!! \f[
!!   a_{1}=10^{-3}exp\left[ 0.025\left(T-273.15\right)\right]
!! \f]
!!
!!    - The accretion of cloud ice by snow (\f$P_{saci}\f$) in the regions where  cloud ice
!! exists is simulated by
!!\f[
!!  P_{saci}=C_{s}mP_{s}
!!\f]
!! where \f$P_{s}\f$ is the precipitation rate of snow. The collection coefficient \f$C_{s}\f$
!! is a function of temperature since the open structures of ice crystals at relative warm
!! temperatures are more likely to stick, given a collision, than crystals of other shapes
!! (Rogers 1979 \cite rogers_1979). Above the freezing level, \f$C_{s}\f$ is expressed by
!!\f[
!!   C_{s}=c_{1}exp[c_{2}(T-273.15)]
!!\f]
!! where \f$c_{1}=1.25\times 10^{-3} m^{2}kg^{-1}s^{-1}\f$ and \f$c_{2}=0.025 K^{-1}\f$ are used.
!! \f$C_{s}\f$ is set to zero below the freezing level.
!!
!---   precipitation production --  auto conversion and accretion
!
        do n=1,ihpr
          if (comput(n) .and. ccr(n) > 0.0) then
            wws    = ww(n)
            cwmk   = max(cons_0, wws)
            i      = ipr(n)
!           amaxcm = max(cons_0, cwmk - wmink(n))
            if (iwl(n) == 1) then                 !  ice phase
               amaxcm = max(cons_0, cwmk - wmini(i,k))
               expf      = dt * exp(0.025*tmt0(n))
               psaut     = min(cwmk, psautco(i)*expf*amaxcm)
               ww(n)     = ww(n) - psaut
               cwmk      = max(cons_0, ww(n))
!              cwmk      = max(cons_0, ww(n)-wmini(i,k))
               psaci     = min(cwmk, aa2*expf*precsl1(n)*cwmk)

               ww(n)     = ww(n) - psaci

               precsl(n) = precsl(n) + (wws - ww(n)) * condt(n)
            else                                    !  liquid water
!
!>    - Following Sundqvist et al. (1989) \cite sundqvist_et_al_1989, the autoconversion of cloud water to rain
!! (\f$P_{raut}\f$) can be parameterized from the cloud water mixing ratio \f$m\f$ and cloud
!! coverage \f$b\f$, that is,
!!\f[
!!  P_{raut}=c_{0}m\left\{1-exp[-(\frac{m}{m_{r}b})^{2}]\right\}
!!\f]
!! where \f$c_{0}\f$ and \f$m_{r}\f$ are \f$1.0\times 10^{-4}s^{-1}\f$ and \f$3.0\times 10^{-4}\f$,
!! respectively.
!          for using sundqvist precip formulation of rain
!
               amaxcm    = max(cons_0, cwmk - wmink(n))
!!             amaxcm    = cwmk
               tem1      = precsl1(n) + precrl1(n)
               tem2      = min(max(cons_0, 268.0-tt(n)), cons_20)
               tem       = (1.0+c1*sqrt(tem1*rdt)) * (1+c2*sqrt(tem2))
!
               tem2      = amaxcm * cmr * tem / max(ccr(n),cons_p01)
               tem2      = min(cons_50, tem2*tem2)
!              praut     = c00  * tem * amaxcm * (1.0-exp(-tem2))
               praut     = (prautco(i)*dt) * tem * amaxcm
     &                                     * (1.0-exp(-tem2))
               praut     = min(praut, cwmk)
               ww(n)     = ww(n) - praut
!
!          below is for zhao's precip formulation (water)
!
!              amaxcm    = max(cons_0, cwmk - wmink(n))
!              praut     = min(cwmk, c00*amaxcm*amaxcm)
!              ww(n)     = ww(n) - praut
!
!              cwmk      = max(cons_0, ww(n))
!              tem1      = precsl1(n) + precrl1(n)
!              pracw     = min(cwmk, cr*dt*tem1*cwmk)
!              ww(n)     = ww(n) - pracw
!
               precrl(n) = precrl(n) + (wws - ww(n)) * condt(n)
!
!hchuang code change [+1l] : add record to record information in vertical
! turn rnp in unit of ww (cwm and q, kg/kg ???)
               rnp(n) = rnp(n) + (wws - ww(n))
            endif
          endif
        enddo
!> -# Evaporation of precipitation (\f$E_{rr}\f$ and \f$E_{rs}\f$)
!!\n Evaporation of precipitation is an important process that  moistens the layers below
!! cloud base. Through this process, some of the precipitating water is evaporated
!! back to the atmosphere and the precipitation efficiency is reduced.
!!    - Evaporation of rain is calculated using the equation (Sundqvist 1988 \cite sundqvist_1988):
!!\f[
!!   E_{rr}=k_{e}(f_{0}-f)(P_{r})^{\beta}
!!\f]
!! where \f$k_{e}\f$ and \f$\beta\f$ are empirical parameters.
!! In Zhao and Carr (1997) \cite zhao_and_carr_1997, \f$k_{e}=2\times 10^{-5}m^{1}kg^{-0.5}s^{-1}\f$ and \f$\beta=0.5\f$
!! gave reasonable values of rain evaporation.
!!
!!    - Evaporation of snow is calculated using the equation:
!!\f[
!!  E_{rs}=[C_{rs1}+C_{rs2}(T-273.15)](\frac{f_{0}-f}{f_{0}})P_{s}
!!\f]
!! where \f$C_{rs1}=5\times 10^{-6}m^{2}kg^{-1}s^{-1}\f$ and
!! \f$C_{rs2}=6.67\times 10^{-10}m^{2}kg^{-1}K^{-1}s^{-1}\f$. The evaporation of melting snow
!! below the freezing level is ignored in this scheme because of the difficulty in the latent
!! heat treatment since the surface of a melting snowflake is usually covered by a thin layer
!! of liquid water.
!
!-----evaporation of precipitation-------------------------
!**** err & ers positive--->evaporation-- negtive--->condensation
!
        do n=1,ihpr
          if (comput(n)) then
            i      = ipr(n)
            qk     = max(epsq,  qq(n))
            tmt0k  = max(cons_m30, tmt0(n))
            precrk = max(cons_0,    precrl(n))
            precsk = max(cons_0,    precsl(n))
            amaxrq = max(cons_0,    u00k(i,k)-rq(n)) * conde(n)
!----------------------------------------------------------------------
! increase the evaporation for strong/light prec
!----------------------------------------------------------------------
            ppr    = ke * amaxrq * sqrt(precrk)
!           ppr    = ke * amaxrq * sqrt(precrk*rdt)
            if (tmt0(n) .ge. 0.) then
              pps = 0.
            else
              pps = (crs1+crs2*tmt0k) * amaxrq * precsk / u00k(i,k)
            end if
!---------------correct if over-evapo./cond. occurs--------------------
            erk=precrk+precsk
            if(rq(n).ge.1.0e-10)  erk = amaxrq * qk * rdt / rq(n)
            if (ppr+pps .gt. abs(erk)) then
               rprs   = erk / (precrk+precsk)
               ppr    = precrk * rprs
               pps    = precsk * rprs
            endif
            ppr       = min(ppr, precrk)
            pps       = min(pps, precsk)
            err(n)    = ppr * rconde(n)
            ers(n)    = pps * rconde(n)
            precrl(n) = precrl(n) - ppr
!hchuang code change [+1l] : add record to record information in vertical
! use err for kg/kg/dt not the ppr (mm/dt=kg/m2/dt)
!
            rnp(n) = rnp(n) - err(n)
!
            precsl(n) = precsl(n) - pps
          endif
        enddo
!> -# Melting of snow (\f$P_{sm1}\f$ and \f$P_{sm2}\f$)
!!\n In this scheme, we allow snow melting to take place in certain temperature regions
!! below the freezing level in two ways. In both cases, the melted snow is assumed to become
!! raindrops.
!!    - One is the continuous melting of snow due to the increase in temperature as it falls
!! down through the freezing level. This process is parameterized as a function of temperature
!! and snow precipitation rate, that is,
!!\f[
!! P_{sm1}=C_{sm}(T-273.15)^{\alpha}P_{s}
!!\f]
!! In Zhao and Carr (1997) \cite zhao_and_carr_1997, parameter values of \f$C_{sm}=5\times 10^{-8}m^{2}kg^{-1}K^{-2}s^{-1}\f$
!! and \f$\alpha=2\f$ cause the falling snow to melt almost completely before it reaches the
!! \f$T=278.15 K\f$ level.
!!    - Another is the immediate melting of melting snow by collection of the cloud water
!! below the freezing level. In order to calculate the melting rate, the
!! collection rate of cloud water by melting snow is computed first. Similar to the collection of cloud water
!! by rain, the collection of cloud water by melting snow can be parameterized to be proportional
!! to the cloud water mixing ratio \f$m\f$ and the precipitation rate of snow \f$P_{s}\f$:
!!\f[
!!   P_{sacw}=C_{r}mP_{s}
!!\f]
!! where \f$C_{r}\f$ is the collection coefficient,\f$C_{r}=5.0\times 10^{-4}m^{2}kg^{-1}s^{-1}\f$
!!. The melting rate of snow then can be computed from
!!\f[
!!   P_{sm2}=C_{ws}P_{sacw}
!!\f]
!! where \f$C_{ws}=0.025\f$
!--------------------melting of the snow--------------------------------
        do n=1,ihpr
          if (comput(n)) then
            if (tmt0(n) .gt. 0.) then
               amaxps = max(cons_0,    precsl(n))
               psm1   = csm1 * tmt0(n) * tmt0(n) * amaxps
               psm2   = cws * cr * max(cons_0, ww(n)) * amaxps
               ppr    = (psm1 + psm2) * conde(n)
               if (ppr .gt. amaxps) then
                 ppr  = amaxps
                 psm1 = amaxps * rconde(n)
               endif
               precrl(n) = precrl(n) + ppr
!
!hchuang code change [+1l] : add record to record information in vertical
! turn ppr (mm/dt=kg/m2/dt) to kg/kg/dt -> ppr/air density (kg/m3)
               rnp(n) = rnp(n) + ppr * rconde(n)
!
               precsl(n) = precsl(n) - ppr
            else
               psm1 = d00
            endif
!
!---------------update t and q------------------------------------------
!>    - Update t and q.
            tt(n) = tt(n) - dtcp * (elwv*err(n)+eliv*ers(n)+eliw*psm1)
            qq(n) = qq(n) + dt * (err(n)+ers(n))
          endif
        enddo
!
        do n=1,ihpr
          iwl1(n)    = iwl(n)
          precrl1(n) = max(cons_0, precrl(n))
          precsl1(n) = max(cons_0, precsl(n))
          i          = ipr(n)
          t(i,k)     = tt(n)
          q(i,k)     = qq(n)
          cwm(i,k)   = ww(n)
          iw(i,k)    = iwl(n)
!hchuang code change [+1l] : add record to record information in vertical
! rnp = precrl1*rconde(n) unit in kg/kg/dt
!
          rainp(i,k) = rnp(n)
        enddo
!
!  move water from vapor to liquid should the liquid amount be negative
!
        do i = 1, im
          if (cwm(i,k) < 0.) then
            tem      = q(i,k) + cwm(i,k)
            if (tem >= 0.0) then
              q(i,k)   = tem
              t(i,k)   = t(i,k) - elwv * rcp * cwm(i,k)
              cwm(i,k) = 0.
            elseif (q(i,k) > 0.0) then
              cwm(i,k) = tem
              t(i,k)   = t(i,k) + elwv * rcp * q(i,k)
              q(i,k)   = 0.0
            endif
          endif
        enddo
!
      enddo                               ! k loop ends here!
!**********************************************************************
!-----------------------end of precipitation processes-----------------
!**********************************************************************
!
!> -# Compute precipitation at surface and determine fraction of frozen precipitation.
      do n=1,ihpr
        i = ipr(n)
        rn(i) = (precrl1(n)  + precsl1(n)) * rrow  ! precip at surface
!
!----sr=1 if sfc prec is rain ; ----sr=-1 if sfc prec is snow
!----sr=0 for both of them or no sfc prec
!
!        rid = 0.
!        sid = 0.
!        if (precrl1(n) .ge. 1.e-13) rid = 1.
!        if (precsl1(n) .ge. 1.e-13) sid = -1.
!        sr(i) = rid + sid  ! sr=1 --> rain, sr=-1 -->snow, sr=0 -->both
! chuang, june 2013: change sr to define fraction of frozen precipitation instead
! because wpc uses it in their winter experiment

        rid = precrl1(n) + precsl1(n)
        if (rid < 1.e-13) then
           sr(i) = 0.
        else
           sr(i) = precsl1(n)/rid
        endif
      enddo
!
      return
      end
!! @}
!> @}
