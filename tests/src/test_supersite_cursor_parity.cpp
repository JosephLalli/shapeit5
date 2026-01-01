#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <unistd.h>

#include <libgen.h> // For dirname
#include <limits.h> // For PATH_MAX


#include "test_common.h"
int main() {
    TEST_INIT("test_supersite_cursor_parity");
    const char* home = getenv("HOME");
    std::string ld = (home ? std::string(home) + "/.linuxbrew/lib:/usr/local/lib" : std::string("/usr/local/lib"));
    const char* old_ld = getenv("LD_LIBRARY_PATH");
    if (old_ld) ld += std::string(":") + old_ld;

    // Get absolute path of current executable
    char self_path[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", self_path, PATH_MAX);
    if (count == -1) {
        std::perror("readlink");
        return 2;
    }
    self_path[count] = '\0';

    // Extract directory of current executable (tests/bin/)
    std::string self_dir = dirname(self_path);
    
    // Construct absolute path to test_supersite_expansion_parity
    std::string expansion_parity_path = self_dir + "/test_supersite_expansion_parity";

    std::string cmd = "LD_LIBRARY_PATH='" + ld + "' " + expansion_parity_path + " 2>&1";

    std::string output;
    const int bufsize = 4096;
    char buf[bufsize];
    
    // Try first command (tests/bin/ path)
    FILE* f = popen(cmd.c_str(), "r");
    if (!f) {
        std::perror("popen");
        return 2;
    }
    
    while (fgets(buf, bufsize, f)) {
        output += buf;
    }
    int rc = pclose(f);
    if (rc == 0) {
        std::cout << "test_supersite_cursor_parity: PASS\n";
        TEST_SUMMARY();
        return 0;
    }

    // Otherwise print the captured output for debugging and fail
    std::cerr << "=== Captured output from child (for debugging) ===\n";
    std::cerr << output << std::endl;
    std::cerr << "test_supersite_cursor_parity: FAIL - child test returned non-zero\n";
    return 1;
}
