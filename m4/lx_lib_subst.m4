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

