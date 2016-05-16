!> \file physparam.f
!! This file contains module physparam.


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

!> \ingroup rad
!! \defgroup physparam physparam
!! @{
!> This module defines commonly used control variables/parameters in physics related programs.
!! -# contains control variables defined in the form of parameter. They are pre-determined 
!!    choices and not adjustable during model's run-time.
!! -# contains control variables defined as module variables. They are more flexible to be 
!!    changed during run-time by either through input namelist, or through model environment condition.
!!    They are preassigned here as the default values.
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
!> \name 1.1 control flags for sw radiation

!> sw heating rate unit control flag: =1:k/day; =2:k/second.
      integer,parameter :: iswrate = 2  ! sw heating rate unit control flag
                                        ! =1:k/day; =2:k/second.
!> sw rare gases effect control flag (ch4,n2o,o2,...): =0:no; =1:yes.
      integer,parameter :: iswrgas = 1  ! sw rare gases effect control flag (ch4,n2o,o2,...)
                                        ! =0:no; =1:yes.
!> sw optical property for liquid clouds
!!\n =0:input cld opt depth, ignoring iswcice setting
!!\n =1:input cwp,rew, use hu and stamnes(1993) method
!!\n =2:not defined yet
      integer,save      :: iswcliq = 1  ! sw optical property for liquid clouds
                                        ! =0:input cld opt depth, ignoring iswcice setting
                                        ! =1:input cwp,rew, use hu and stamnes(1993) method
                                        ! =2:not defined yet
!> sw optical property for ice clouds (only iswcliq>0)
!!\n =0:not defined yet
!!\n =1:input cip,rei, use ebert and curry (1992) method
!!\n =2:input cip,rei, use streamer v3.0 (2001) method
!!\n =3:input cip,rei, use fu's method (1996) method
      integer,save      :: iswcice = 3  ! sw optical property for ice clouds (only iswcliq>0)
                                        ! =0:not defined yet
                                        ! =1:input cip,rei, use ebert and curry (1992) method
                                        ! =2:input cip,rei, use streamer v3.0 (2001) method
                                        ! =3:input cip,rei, use fu's method (1996) method
!> sw control flag for 2-stream transfer scheme
!!\n =1:delta-eddington    (joseph et al., 1976)
!!\n =2:pifm               (zdunkowski et al., 1980)
!!\n =3:discrete ordinates (liou, 1973)
      integer,parameter :: iswmode = 2  ! sw control flag for 2-stream transfer scheme
                                        ! =1:delta-eddington    (joseph et al., 1976)
                                        ! =2:pifm               (zdunkowski et al., 1980)
                                        ! =3:discrete ordinates (liou, 1973)


! ............................................. !
!  -1.2- control flags for lw radiation         !
! ............................................. !
!> \name 1.2 control flags for lw radiation

!> lw heating rate unit (1:k/day; 2:k/second): =1:k/day; =2:k/second.
      integer,parameter :: ilwrate = 2  

!> lw rare gases effect control flag (ch4,n2o,o2,cfcs...): =0:no; =1:yes.
      integer,parameter :: ilwrgas = 1  

!> lw optical property for liquid clouds
!!\n =0:input cld opt depth, ignoring ilwcice setting
!!\n =1:input cwp,rew, use hu and stamnes(1993) method
!!\n =2:not defined yet
      integer,save      :: ilwcliq = 1 

!> lw optical property for ice clouds (only ilwcliq>0)
!!\n =0:not defined yet
!!\n =1:input cip,rei, use ebert and curry (1992) method
!!\n =2:input cip,rei, use streamer (1996) method
!!\n =3:input cip,rei, use fu's method (1998) method
      integer,save      :: ilwcice = 3  

! ............................................. !
!  -1.3- control flag for lw aerosol property   !

!>\name 1.3 control flag for lw aerosol property

!> =t: use 1 broad-band lw aeros properties
!!\n =f: use multi bands aeros properites
      logical,parameter :: lalw1bd =.false. 

!==================================================================================
!  Section - 2 -
!     values of control flags might be re-set in initialization subroutines
!       (may be adjusted at run time based on namelist input or run condition)
!==================================================================================

! ............................................. !
!  -2.1- for module radiation_astronomy         !
! ............................................. !
!> \name 2.1 for module radiation_astronomy 

!> solar constant scheme control flag
      integer, save :: isolar  = 0   

!> external solar constant data table,solarconstant_noaa_a0.txt
      character, save :: solar_file*26 
!     data solar_file   / 'solarconstantdata.txt     ' /
      data solar_file   / 'solarconstant_noaa_a0.txt ' /

! ............................................. !
!  -2.2- for module radiation_aerosols          !
! ............................................. !
!> \name 2.2 for module radiation_aerosols

!> aerosol model scheme control flag
      integer, save :: iaermdl = 0      
!> aerosol effect control flag
      integer, save :: iaerflg = 0     
!> lw aerosols effect control flag
      logical, save :: lalwflg = .true. 
!> sw aerosols effect control flag
      logical, save :: laswflg = .true. 
!> stratospheric volcanic effect flag
      logical, save :: lavoflg = .true. 
!> external aerosols data file: aerosol.dat
      character, save :: aeros_file*26  
!     data aeros_file   / 'climaeropac_global.txt    ' /
      data aeros_file   / 'aerosol.dat               ' /

! ............................................. !
!  -2.3- for module radiation_gases             !
! ............................................. !
!> \name 2.3 for module radiation_gases

!> co2 data source control flag
      integer, save :: ico2flg = 0    
!> external data time/date control flag
      integer, save :: ictmflg = 0     
!> ozone data source control flag
      integer, save :: ioznflg = 1      
!> external co2 2d monthly obsv data table: co2historicaldata_2004.txt
      character, save :: co2dat_file*26 
!> external co2 global annual mean data tb: co2historicaldata_glob.txt
      character, save :: co2gbl_file*26 
!> external co2 user defined data table: co2userdata.txt 
      character, save :: co2usr_file*26 
!> external co2 clim monthly cycle data tb: co2monthlycyc.txt
      character, save :: co2cyc_file*26
      data co2dat_file   / 'co2historicaldata_2004.txt' /   !year is run-time selected
      data co2gbl_file   / 'co2historicaldata_glob.txt' /
      data co2usr_file   / 'co2userdata.txt           ' /
      data co2cyc_file   / 'co2monthlycyc.txt         ' /

! ............................................. !
!  -2.4- for module radiation_clouds            !
! ............................................. !
!> \name 2.4 for module radiation_clouds

!> cloud optical property scheme control flag
      integer, save :: icldflg = 1      
!> cloud micorphysics scheme control flag
      integer, save :: icmphys = 1    
!> cloud overlapping control flag for SW
      integer, save :: iovrsw  = 1   
!> cloud overlapping control flag for LW
      integer, save :: iovrlw  = 1  
!> eliminating CRICK control flag
      logical, save :: lcrick  =.false. 
!> in-cld condensate control flag
      logical, save :: lcnorm  =.false.
!> precip effect on radiation flag (Ferrier microphysics)
      logical, save :: lnoprec =.false. 
!> shallow convetion flag
      logical, save :: lsashal =.false. 

! ............................................. !
!  -2.5- for module radiation_surface           !
! ............................................. !
!> \name 2.5 for module radiation_surface

!> surface albedo scheme control flag
      integer, save :: ialbflg = 0      
!> surface emissivity scheme control flag
      integer, save :: iemsflg = 0     

!> external sfc emissivity data table: sfc_emissivity_idx.txt
      character, save :: semis_file*26  
      data semis_file   / 'sfc_emissivity_idx.txt    ' /

! ............................................. !
!  -2.6- general purpose                        !
! ............................................. !
!> \name 2.6 general purpose 

!> vertical profile indexing flag
      integer, save :: ivflip  = 1      
!> sub-column cloud approx flag in SW radiation
      integer, save :: isubcsw = 0      
!> sub-column cloud approx flag in LW radiation
      integer, save :: isubclw = 0      
!> initial permutaion seed for mcica radiation
      integer, save :: ipsd0   = 0    

!
!...................................!
      end module physparam          !
!===================================!
!! @}
