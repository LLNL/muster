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
#  LX_HEADER_SUBST (name, header, NAME, includes)
#   ------------------------------------------------------------------------
#  This tests for the presence of a header file, given its name, and a directory
#  to search.  If found, it uses AC_SUBST to export NAME_CPPFLAGS for the header.
#  Standard var CPPFLAGS is unmodified.  If the header is not found, then have_name
#  is set to "no".
# 
AC_DEFUN([LX_HEADER_SUBST],
[
  $3_CPPFLAGS="$4"
  OLD_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$$3_CPPFLAGS $CPPFLAGS"

  AC_CHECK_HEADER([$2],
                  [have_$1=yes],
                  [AC_MSG_NOTICE([Couldn't find $2.])
                   have_$1=no])
  AC_SUBST($3_CPPFLAGS)
  CPPFLAGS="$OLD_CPPFLAGS"
])

