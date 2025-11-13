/*******************************************************************************
 * Supersite tracing helpers
 *
 * Light-weight utilities used across the supersite instrumentation to gate
 * verbose logging behind SHAPEIT5_TEST_TRACE (or the future dedicated env var).
 ******************************************************************************/
#pragma once

#include <cstdarg>
#include <cstdlib>
#include <cstdio>

inline bool supersite_trace_enabled() {
    static int trace_flag = -1;
    if (trace_flag < 0) {
        const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
        trace_flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return trace_flag == 1;
}

inline void supersite_trace_log(const char* fmt, ...) {
    if (!supersite_trace_enabled()) return;
    va_list args;
    va_start(args, fmt);
    std::vfprintf(stderr, fmt, args);
    va_end(args);
    std::fflush(stderr);
}
