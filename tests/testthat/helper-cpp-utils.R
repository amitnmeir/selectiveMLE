# Compile the C++ wrappers that expose sampler helper functions.
# The implementations mirror the internal routines in
# `src/lassoSamplerCpp.cpp` but are only built here for testing so that
# the low level logic can be exercised from R.
Rcpp::sourceCpp('cpp_utils.cpp')
