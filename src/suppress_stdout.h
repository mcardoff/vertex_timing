#ifndef SUPPRESS_STDOUT_H
#define SUPPRESS_STDOUT_H

// ---------------------------------------------------------------------------
// suppress_stdout.h
//   Provides SuppressStdout, an RAII guard that silences stdout for the
//   duration of its scope.  Needed because ROOT's TEfficiency printing
//   goes through printf (which bypasses gErrorIgnoreLevel) and cannot be
//   suppressed via any ROOT API.
// ---------------------------------------------------------------------------

#include <cstdio>
#include <unistd.h>
#include <fcntl.h>

// ---------------------------------------------------------------------------
// SuppressStdout
//   On construction: flushes stdout, duplicates the current fd, then
//   redirects STDOUT_FILENO to /dev/null.
//   On destruction: flushes again and restores the original fd.
//   Usage: instantiate as a block-scoped local variable around any code
//   that produces unwanted printf output.
// ---------------------------------------------------------------------------
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
