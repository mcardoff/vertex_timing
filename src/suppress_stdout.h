#ifndef SUPPRESS_STDOUT_H
#define SUPPRESS_STDOUT_H

#include <cstdio>
#include <unistd.h>
#include <fcntl.h>

// RAII guard that redirects stdout to /dev/null for its lifetime.
// Used to suppress ROOT's "OBJ: TEfficiency..._clone" messages that are
// printed via printf (bypassing gErrorIgnoreLevel) when TEfficiency objects
// are constructed or drawn.
struct SuppressStdout {
  int saved_fd;
  SuppressStdout() {
    fflush(stdout);
    saved_fd = dup(STDOUT_FILENO);
    int devnull = open("/dev/null", O_WRONLY);
    dup2(devnull, STDOUT_FILENO);
    close(devnull);
  }
  ~SuppressStdout() {
    fflush(stdout);
    dup2(saved_fd, STDOUT_FILENO);
    close(saved_fd);
  }
};

#endif // SUPPRESS_STDOUT_H
