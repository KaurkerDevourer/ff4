#pragma once
#include <execinfo.h>
#include <iostream>
#include <cassert>

#define assertm(exp, msg) assert(((void)msg, exp))

void backtrace() {
  void* callstack[128];
  int i, frames = backtrace(callstack, 128);
  char** strs = backtrace_symbols(callstack, frames);
  for (i = 0; i < frames; ++i) {
    std::cout << strs[i] << std::endl;
  }
  free(strs);
}


#define ASSERT_EQUAL(A, B) \
    if (A != B) { \
        backtrace(); \
        std::cout << A << " != " << B << std::endl; \
    }

    