# the name of the target operating system
set(CMAKE_SYSTEM_NAME BlueGeneP)

# set the compiler
set(CMAKE_C_COMPILER   /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc )
set(CMAKE_CXX_COMPILER /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-g++ )

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /bgsys/drivers/ppcfloor/
    /bgsys/drivers/ppcfloor/gnu-linux/powerpc-bgp-linux/
    /bgsys/drivers/ppcfloor/runtime/
    /bgsys/drivers/ppcfloor/comm/
    /bgsys/drivers/ppcfloor/comm/default/
    /bgsys/drivers/ppcfloor/comm/sys/
)

# adjust the default behaviour of the FIND_XXX() commands:
# search headers and libraries in the target environment, search 
# programs in the host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
