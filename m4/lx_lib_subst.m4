##################################################################################################
# Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
# Produced at the Lawrence Livermore National Laboratory  
# Written by Todd Gamblin, tgamblin@llnl.gov.
# LLNL-CODE-433662
# All rights reserved.  
#
# This file is part of Muster. For details, see http://github.com/tgamblin/muster. 
# Please also read the LICENSE file for further information.
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the disclaimer below.
#  * Redistributions in binary form must reproduce the above copyright notice, this list of
#    conditions and the disclaimer (as noted below) in the documentation and/or other materials
#    provided with the distribution.
#  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##################################################################################################


#
#  LX_LIB_SUBST (name, symbol, NAME, [dir], [libs])
#  ------------------------------------------------------------------------
#  This tests for the presence of a library, given a libname and a symbol.
#  
#  If the library is found, it uses AC_SUBST to export the following:
#     <NAME>_LDFLAGS   standard linker args (-L/foo/lib -lname, etc.)
#     <NAME>_RPATH     rpath arguments for library locations (-R /foo/lib, etc.)
#  
#  Standard vars like LIBS, LDFLAGS, etc. are unmodified.  If the library
#  is not found, then have_<name> is set to "no".
#  
#  If dir is provided, -L<dir> will be added to the link line and an rpath argument
#  will be added to the command line.
#
#  You can optionally supply specific library names, if the name of the module
#  does not match the name of the library.  The library name defaults to -l<name>,
#  but can be overridden by supplying the libs parameter.  This is useful for 
#  packages with multiple libraries, e.g. "-lfoo -lbar -lbaz".
#
AC_DEFUN([LX_LIB_SUBST],
[
  # If dir param is there, add it to the link line.
  if test "x$4" != x; then
    $3_LDFLAGS="-L$4"
  else
    $3_LDFLAGS=""
  fi

  # If extra libraries are supplied, add those too.
  if test "x$5" != "x"; then
    $3_LDFLAGS="$$3_LDFLAGS $5"
  else
    $3_LDFLAGS="$$3_LDFLAGS -l$1"
  fi

  # search the link line for -Lwhatever and add rpath args for each one.
  $3_RPATH=""
  for elt in $$3_LDFLAGS; do
     if echo $elt | grep -q '^-L'; then
        rpath=`echo $elt | sed 's/^-L/-R /'`
        if [[ -z "$$3_RPATH" ]]; then
            $3_RPATH="$rpath"
        else
            $3_RPATH="$$3_RPATH $rpath"
        fi
     fi
  done

  # Do the actual check to determine if we can link against this library.
  OLD_LDFLAGS="$LDFLAGS"
  LDFLAGS="$$3_LDFLAGS $LDFLAGS"
  
  LX_SAFE_CHECK_LIB([$1], [$2],
                    [have_$1=yes],
                    [AC_MSG_NOTICE([Couldn't find lib$1.])
                     $3_LDFLAGS=""
                     have_$1=no])

  LDFLAGS="$OLD_LDFLAGS"

  AC_SUBST($3_LDFLAGS)
  AC_SUBST($3_RPATH)
])

