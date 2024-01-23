## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

N  checking package dependencies
   Imports includes 28 non-default packages.
   Importing from so many packages makes the package vulnerable to any of
   them becoming unavailable.  Move as many as possible to Suggests and
   use conditionally.
   


rhub::check_for_cran currently says that the build fails on ubuntu-gcc-release and fedora-clang-devel.
However, I checked the build logs, and all of them say that the package successfully built. Please let me know
which is the case here.

