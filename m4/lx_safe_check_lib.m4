# 
#  LX_SAFE_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
# 
#  This is a safe wrapper for AC_CHECK_LIB, which will NOT append to 
#  $LIBS.  The default version will add the library to the link line
#  for every target, while this just makes sure that you have the lib.
#  
AC_DEFUN([LX_SAFE_CHECK_LIB],
[
    OLD_LIBS="$LIBS"
    AC_CHECK_LIB($1,$2,$3,$4,$5)
    LIBS="$OLD_LIBS"
])

