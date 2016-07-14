!> \file physparam.f
!! This file contains module physparam.

!  ==========================================================  !!!!!
!                    module physparam description              !!!!!
!  ==========================================================  !!!!!
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
!> \defgroup physparam physparam
!! @{
!> This module defines commonly used control variables and parameters in physics related programs.
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
!> \name  -1.1- Control flags for SW radiation       
! ............................................. !

!> SW heating rate unit control flag: =1:k/day; =2:k/second.
      integer,parameter :: iswrate = 2  
                                       
!> SW rare gases effect control flag (ch4,n2o,o2,...): =0:no; =1:yes.
      integer,parameter :: iswrgas = 1
                                       
!> SW optical property for liquid clouds
!!\n =0:input cld opt depth, ignoring iswcice setting
!!\n =1:input cwp,rew, use Hu and Stamnes(1993) \cite hu_and_stamnes_1993 method
!!\n =2:not defined yet
      integer,save      :: iswcliq = 1  
                                       
!> SW optical property for ice clouds (only iswcliq>0)
!!\n =0:not defined yet
!!\n =1:input cip,rei, use Ebert and Curry (1992) \cite ebert_and_curry_1992 method
!!\n =2:input cip,rei, use Streamer v3.0 (2001) \cite key_2001 method
!!\n =3:input cip,rei, use Fu's method (1996) \cite fu_1996 method
      integer,save      :: iswcice = 3  
                                       
!> SW control flag for 2-stream transfer scheme
!!\n =1:delta-eddington    (Joseph et al. 1976 \cite joseph_et_al_1976)
!!\n =2:pifm               (Zdunkowski et al. 1980 \cite zdunkowski_et_al_1980)
!!\n =3:discrete ordinates (Liou, 1973 \cite liou_1973)
      integer,parameter :: iswmode = 2 

! ............................................. !
!> \name  -1.2- Control flags for LW radiation  
! ............................................. !

!> LW heating rate unit (1:k/day; 2:k/second): =1:k/day; =2:k/second.
      integer,parameter :: ilwrate = 2  

!> LW rare gases effect control flag (ch4,n2o,o2,cfcs...): =0:no; =1:yes.
      integer,parameter :: ilwrgas = 1  

!> LW optical property for liquid clouds
!!\n =0:input cld opt depth, ignoring ilwcice setting
!!\n =1:input cwp,rew, use Hu and Stamnes(1993) \cite hu_and_stamnes_1993 method
!!\n =2:not defined yet
      integer,save      :: ilwcliq = 1 

!> LW optical property for ice clouds (only ilwcliq>0)
!!\n =0:not defined yet
!!\n =1:input cip,rei, use Ebert and Curry (1992) \cite ebert_and_curry_1992 method
!!\n =2:input cip,rei, use Streamer (1996) method
!!\n =3:input cip,rei, use Fu's method (1998) \cite fu_et_al_1998 method
      integer,save      :: ilwcice = 3  

! ............................................. !
!>\name  -1.3- Control flag for LW aerosol property 

!> =t: use 1 broad-band LW aeros properties
!!\n =f: use multi bands aeros properites
      logical,parameter :: lalw1bd =.false. 

!==================================================================================
!  Section - 2 -
!     values of control flags might be re-set in initialization subroutines
!       (may be adjusted at run time based on namelist input or run condition)
!==================================================================================

! ............................................. !
!>\name  -2.1- For module radiation_astronomy 
! ............................................. !

!> solar constant scheme control flag
      integer, save :: isolar  = 0   

!> external solar constant data table,solarconstant_noaa_a0.txt
      character, save :: solar_file*26 
!     data solar_file   / 'solarconstantdata.txt     ' /
      data solar_file   / 'solarconstant_noaa_a0.txt ' /

! ............................................. !
!> \name  -2.2- For module radiation_aerosols   
! ............................................. !

!> aerosol model scheme control flag
      integer, save :: iaermdl = 0      
!> aerosol effect control flag
      integer, save :: iaerflg = 0     
!> LW aerosols effect control flag
      logical, save :: lalwflg = .true. 
!> SW aerosols effect control flag
      logical, save :: laswflg = .true. 
!> stratospheric volcanic effect flag
      logical, save :: lavoflg = .true. 
!> external aerosols data file: aerosol.dat
      character, save :: aeros_file*26  
!     data aeros_file   / 'climaeropac_global.txt    ' /
      data aeros_file   / 'aerosol.dat               ' /

! ............................................. !
!> \name  -2.3- For module radiation_gases 
! ............................................. !

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
!>\name  -2.4- For module radiation_clouds 
! ............................................. !

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
!>\name  -2.5- For module radiation_surface 
! ............................................. !

!> surface albedo scheme control flag
      integer, save :: ialbflg = 0      
!> surface emissivity scheme control flag
      integer, save :: iemsflg = 0     

!> external sfc emissivity data table: sfc_emissivity_idx.txt
      character, save :: semis_file*26  
      data semis_file   / 'sfc_emissivity_idx.txt    ' /

! ............................................. !
!> \name  -2.6- general purpose     
! ............................................. !

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
