!!!!!  ==========================================================  !!!!!
!!!!!                    module physparam description              !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!     This module defines commonly used control variables/parameters   !
!     in physics related programs.                                     !
!                                                                      !
!     Section 1 contains control variables defined in the form of      !
!     parameter. They are pre-determined choices and not adjustable    !
!     during model's run-time.                                         !
!                                                                      !
!     Section 2 contains control variables defined as module variables.!
!     They are more flexible to be changed during run-time by either   !
!     through input namelist, or through model environment condition.  !
!     They are preassigned here as the default values.                 !
!                                                                      !
!!!!!  ==========================================================  !!!!!

!> This module defines commonly used control variables/parameters in physics related programs.
!! \n Section 1 contains control variables defined in the form of parameter. They are pre-determined
!! choices and not adjustable during model's run-time.
!! \n Section 2 contains control variables defined as module variables. They are more flexible to 
!! be changed during run-time by either through input namelist, or through model environment condition.
!! They are preassigned here as the default values.

!========================================!
      module physparam                   !
!........................................!
!
!     implicit   none

!  --- ...  define kind parameters here

!   ** if already exist, use the module containing kind definitions
      use machine

!   ** otherwise, define kind parameter here
!     implicit   none
!     integer, public, parameter :: kind_io4 = 4
!     integer, public, parameter :: kind_io8 = 8
!     integer, public, parameter :: kind_phys= selected_real_kind(13,60) ! the '60' maps to 64-bit real
!      .....

!     implicit   none
!
      public

!==================================================================================
!  Section - 1 -
!     control flags are pre-set as run-time non-adjuztable parameters.
!==================================================================================

! ............................................. !
!  -1.1- control flags for sw radiation         !
! ............................................. !
!> sw heating rate unit control flag
!! \n =1:k/day; =2:k/second.
      integer,parameter :: iswrate = 2  ! sw heating rate unit control flag
                                        ! =1:k/day; =2:k/second.
!> sw rare gases effect control flag (ch4,n2o,o2,...)
!! \n =0:no; =1:yes.
      integer,parameter :: iswrgas = 1  ! sw rare gases effect control flag (ch4,n2o,o2,...)
                                        ! =0:no; =1:yes.
!> sw optical property for liquid clouds
!! \n =0:input cld opt depth, ignoring iswcice setting
!! \n =1:input cwp,rew, use hu and stamnes (1993) method
!! \n =2:not defined yet.
      integer,save      :: iswcliq = 1  ! sw optical property for liquid clouds
                                        ! =0:input cld opt depth, ignoring iswcice setting
                                        ! =1:input cwp,rew, use hu and stamnes(1993) method
                                        ! =2:not defined yet
!> sw optical property for ice clouds (only iswcliq>0)
!! \n =0:not defined yet
!! \n =1:input cip,rei, use ebert and curry (1992) method
!! \n =2:input cip,rei, use streamer v3.0 (2001) method
!! \n =3:input cip,rei, use fu's method (1996)
      integer,save      :: iswcice = 3  ! sw optical property for ice clouds (only iswcliq>0)
                                        ! =0:not defined yet
                                        ! =1:input cip,rei, use ebert and curry (1992) method
                                        ! =2:input cip,rei, use streamer v3.0 (2001) method
                                        ! =3:input cip,rei, use fu's method (1996) method
!> sw control flag for 2-stream transfer scheme
!! \n =1:delta-eddington    (joseph et al., 1976)
!! \n =2:pifm               (zdunkowski et al., 1980)
!! \n =3:discrete ordinates (liou, 1973)
      integer,parameter :: iswmode = 2  ! sw control flag for 2-stream transfer scheme
                                        ! =1:delta-eddington    (joseph et al., 1976)
                                        ! =2:pifm               (zdunkowski et al., 1980)
                                        ! =3:discrete ordinates (liou, 1973)


! ............................................. !
!  -1.2- control flags for lw radiation         !
! ............................................. !
!> lw heating rate unit 
!! \n =1:k/day; =2:k/second.
      integer,parameter :: ilwrate = 2  ! lw heating rate unit (1:k/day; 2:k/second)
                                        ! =1:k/day; =2:k/second.
!> lw rare gases effect control flag (ch4,n2o,o2,cfcs,...)
!! \n =0:no; =1:yes.
      integer,parameter :: ilwrgas = 1  ! lw rare gases effect control flag (ch4,n2o,o2,cfcs...)
                                        ! =0:no; =1:yes.
!> lw optical property for liquid clouds
!! \n =0:input cld opt depth, ignoring ilwcice setting
!! \n =1:input cwp,rew,use hu and stamnes (1993) method
!! \n =2:not defined yet
      integer,save      :: ilwcliq = 1  ! lw optical property for liquid clouds
                                        ! =0:input cld opt depth, ignoring ilwcice setting
                                        ! =1:input cwp,rew, use hu and stamnes(1993) method
                                        ! =2:not defined yet
!> lw optical property for ice clouds (only ilwcliq>0)
!! \n =0:not defined yet
!! \n =1:input cip,rei, use ebert and curry (1992) method
!! \n =2:input cip,rei, use streamer (1996) method
!! \n =3:input cip,rei, use fu's method (1998)
      integer,save      :: ilwcice = 3  ! lw optical property for ice clouds (only ilwcliq>0)
                                        ! =0:not defined yet
                                        ! =1:input cip,rei, use ebert and curry (1992) method
                                        ! =2:input cip,rei, use streamer (1996) method
                                        ! =3:input cip,rei, use fu's method (1998) method

! ............................................. !
!  -1.3- control flag for lw aerosol property   !
!> control flag for lw aerosol property
!! \n =t: use 1 broad-band lw aeros properties
!! \n =f: use multi bands aeros properties
      logical,parameter :: lalw1bd =.false. ! =t: use 1 broad-band lw aeros properties
                                            ! =f: use multi bands aeros properites




!==================================================================================
!  Section - 2 -
!     values of control flags might be re-set in initialization subroutines
!       (may be adjusted at run time based on namelist input or run condition)
!==================================================================================

! ............................................. !
!  -2.1- for module radiation_astronomy         !
! ............................................. !
!> solar constant scheme control flag
      integer, save :: isolar  = 0      ! solar constant scheme control flag
!> external solar constant data table
      character, save :: solar_file*26  ! external solar constant data table
!     data solar_file   / 'solarconstantdata.txt     ' /
      data solar_file   / 'solarconstant_noaa_a0.txt ' /

! ............................................. !
!  -2.2- for module radiation_aerosols          !
! ............................................. !
!> aerosol model scheme control flag
      integer, save :: iaermdl = 0      ! aerosol model scheme control flag
!> aerosol effect control flag
      integer, save :: iaerflg = 0      ! aerosol effect control flag
!> lw aerosols effect control flag
      logical, save :: lalwflg = .true. ! lw aerosols effect control flag
!> sw aerosols effect control flag
      logical, save :: laswflg = .true. ! sw aerosols effect control flag
!> stratospheric volcanic effect flag
      logical, save :: lavoflg = .true. ! stratospheric volcanic effect flag
!> external aerosols data file
      character, save :: aeros_file*26  ! external aerosols data file
!     data aeros_file   / 'climaeropac_global.txt    ' /
      data aeros_file   / 'aerosol.dat               ' /

! ............................................. !
!  -2.3- for module radiation_gases             !
! ............................................. !
!> co2 data source control flag
      integer, save :: ico2flg = 0      ! co2 data source control flag
!> external data time/date control flag
      integer, save :: ictmflg = 0      ! external data time/date control flag
!> ozone dta source control flag
      integer, save :: ioznflg = 1      ! ozone data source control flag
!> external co2 2d monthly obsv data table
      character, save :: co2dat_file*26 ! external co2 2d monthly obsv data table
!> external co2 global annual mean data table
      character, save :: co2gbl_file*26 ! external co2 global annual mean data tb
!> external co2 user defined data table
      character, save :: co2usr_file*26 ! external co2 user defined data table
!> external co2 clim monthly cycle data table
      character, save :: co2cyc_file*26 ! external co2 clim monthly cycle data tb
      data co2dat_file   / 'co2historicaldata_2004.txt' /   !year is run-time selected
      data co2gbl_file   / 'co2historicaldata_glob.txt' /
      data co2usr_file   / 'co2userdata.txt           ' /
      data co2cyc_file   / 'co2monthlycyc.txt         ' /

! ............................................. !
!  -2.4- for module radiation_clouds            !
! ............................................. !
!> cloud optical property scheme control flag
      integer, save :: icldflg = 1      ! cloud optical property scheme control flag
!> cloud microphysics scheme control flag
      integer, save :: icmphys = 1      ! cloud microphysics scheme control flag
!> cloud overlapping control flag for sw
      integer, save :: iovrsw  = 1      ! cloud overlapping control flag for sw
!> cloud overlapping control flag for lw
      integer, save :: iovrlw  = 1      ! cloud overlapping control flag for lw
!> eliminating CRICK control flag
      logical, save :: lcrick  =.false. ! eliminating CRICK control flag
!> in-cld condensate control flag
      logical, save :: lcnorm  =.false. ! in-cld condensate control flag
!> precip effect on radiation flag (ferrier microphysics)
      logical, save :: lnoprec =.false. ! precip effect on radiation flag (ferrier microphysics)
!> shallow convection flag
      logical, save :: lsashal =.false. ! shallow convection flag

! ............................................. !
!  -2.5- for module radiation_surface           !
! ............................................. !
!> surface albedo scheme control flag
      integer, save :: ialbflg = 0      ! surface albedo scheme control flag
!> surface emissivity scheme control flag
      integer, save :: iemsflg = 0      ! surface emissivity scheme control flag
!> external sfc emissivity data table
      character, save :: semis_file*26  ! external sfc emissivity data table
      data semis_file   / 'sfc_emissivity_idx.txt    ' /

! ............................................. !
!  -2.6- general purpose                        !
! ............................................. !
!> vertical profile indexing flag
      integer, save :: ivflip  = 1      ! vertical profile indexing flag
!> sub-column cloud approx flag in sw radiation
      integer, save :: isubcsw = 0      ! sub-column cloud approx flag in sw radiation
!> sub-column cloud approx flag in lw radiation
      integer, save :: isubclw = 0      ! sub-column cloud approx flag in lw radiation
!> initial permutation seed for mcica radiation
      integer, save :: ipsd0   = 0      ! initial permutation seed for mcica radiation

!
!...................................!
      end module physparam          !
!===================================!
