This directory is to contain the source repositories. The doxygen configuration files expect the following structure:
/src
 /gsm
  /dev_branch
  /trunk

To check out the necessary code from the EMC SVN repo, execute the following commands while in the /src/gsm directory:
svn checkout --ignore-externals https://svnemc.ncep.noaa.gov/projects/gsm/branches/DTC/phys-doc dev_branch
svn checkout --ignore-externals https://svnemc.ncep.noaa.gov/projects/gsm/trunk trunk
