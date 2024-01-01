/// \file
/// \brief Memory allocation wrappers that exit on failure
///
/// Much Graphviz code is not in a position to gracefully handle failure of
/// dynamic memory allocation. The following wrappers provide a safe compromise
/// where allocation failure does not need to be handled, but simply causes
/// process exit. This is not ideal for external callers, but it is better than
/// memory corruption or confusing crashes.
///
/// Note that the wrappers also take a more comprehensive strategy of zeroing
/// newly allocated memory than `malloc`. This reduces the number of things
/// callers need to think about and has only a modest overhead.

#pragma once

#include <assert.h>
#include <cgraph/exit.h>
#include <cgraph/prisize_t.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline void *gv_calloc(size_t nmemb, size_t size) {

  if (nmemb > 0 && SIZE_MAX / nmemb < size) {
    fprintf(stderr,
            "integer overflow when trying to allocate "
            "%" PRISIZE_T " * %" PRISIZE_T " bytes\n",
            nmemb, size);
    graphviz_exit(EXIT_FAILURE);
  }

  void *p = calloc(nmemb, size);
  if (nmemb > 0 && size > 0 && p == NULL) {
    fprintf(stderr,
            "out of memory when trying to allocate %" PRISIZE_T " bytes\n",
            nmemb * size);
    graphviz_exit(EXIT_FAILURE);
  }

  return p;
}

static inline void *gv_alloc(size_t size) { return gv_calloc(1, size); }

static inline void *gv_realloc(void *ptr, size_t old_size, size_t new_size) {

  // make realloc with 0 size equivalent to free, even under C23 rules
  if (new_size == 0) {
    free(ptr);
    return NULL;
  }

  void *p = realloc(ptr, new_size);
  if (p == NULL) {
    fprintf(stderr,
            "out of memory when trying to allocate %" PRISIZE_T " bytes\n",
            new_size);
    graphviz_exit(EXIT_FAILURE);
  }

  // if this was an expansion, zero the new memory
  if (new_size > old_size) {
    memset((char *)p + old_size, 0, new_size - old_size);
  }

  return p;
}

static inline void *gv_recalloc(void *ptr, size_t old_nmemb, size_t new_nmemb,
                                size_t size) {

  assert(size > 0 && "attempt to allocate array of 0-sized elements");
  assert(old_nmemb < SIZE_MAX / size && "claimed previous extent is too large");

  // will multiplication overflow?
  if (new_nmemb > SIZE_MAX / size) {
    fprintf(stderr,
            "integer overflow when trying to allocate %" PRISIZE_T
            " * %" PRISIZE_T " bytes\n",
            new_nmemb, size);
    graphviz_exit(EXIT_FAILURE);
  }

  return gv_realloc(ptr, old_nmemb * size, new_nmemb * size);
}


static inline char *gv_strdup(const char *original) {

  char *copy = strdup(original);
  if (copy == NULL) {
    fprintf(stderr,
            "out of memory when trying to allocate %" PRISIZE_T " bytes\n",
            strlen(original) + 1);
    graphviz_exit(EXIT_FAILURE);
  }

  return copy;
}

static inline char *gv_strndup(const char *original, size_t length) {

  char *copy;

  copy = strndup(original, length);

  if (copy == NULL) {
    fprintf(stderr,
            "out of memory when trying to allocate %" PRISIZE_T " bytes\n",
            length + 1);
    graphviz_exit(EXIT_FAILURE);
  }

  return copy;
}
