#pragma once
#include <gperftools/tcmalloc.h>
#include <memory>
#include <cassert>

void * mem_alloc(size_t bytes) {
  void *p = tc_malloc(bytes);
  assert(p);
  return p;
}

void mem_free(void *p) {
  if (p) tc_free(p);
}

