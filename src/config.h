#pragma once

#include <sofa/config.h>

// Used to disable warnings on external project headers (MSVC only since gcc and clang already do it)
#ifdef _MSC_VER
#define DISABLE_ALL_WARNINGS_BEGIN __pragma(warning(push, 0))
#define DISABLE_ALL_WARNINGS_END   __pragma(warning(pop))
#else
#define DISABLE_ALL_WARNINGS_BEGIN
#define DISABLE_ALL_WARNINGS_END
#endif