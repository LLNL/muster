# the name of the target operating system
set(CMAKE_SYSTEM_NAME BlueGeneP-dynamic)

# Set search paths to prefer local, admin-installed wrappers for the BG backend compilers
set(BGP_XL_COMPILER_SEARCH_PATHS /usr/local/bin /usr/bin)

# GNU C Compilers
find_program(CMAKE_C_COMPILER       bgxlc    ${BGP_XL_COMPILER_SEARCH_PATHS} /opt/ibmcmp/vac/bg/9.0/bin)
find_program(CMAKE_CXX_COMPILER     bgxlC    ${BGP_XL_COMPILER_SEARCH_PATHS} /opt/ibmcmp/vacpp/bg/9.0/bin)
find_program(CMAKE_Fortran_COMPILER bgxlf90  ${BGP_XL_COMPILER_SEARCH_PATHS} /opt/ibmcmp/xlf/bg/11.1/bin)

# Make sure MPI_COMPILER wrapper matches the gnu compilers.  
# Prefer local machine wrappers to driver wrappers here too.
find_program(MPI_COMPILER NAMES mpixlcxx mpixlc++ mpixlC mpixlc
  PATHS 
  /usr/local/bin
  /usr/bin
  /bgsys/drivers/ppcfloor/comm/bin
  /bgsys/drivers/ppcfloor/comm/default/bin)
