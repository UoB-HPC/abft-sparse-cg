#include <cstdlib>
#include <cstring>
#include <iostream>

#include "CGContext.h"

std::list<CGContext::Entry> CGContext::context_list;

CGContext* CGContext::create(const char *target, const char *mode)
{
  // Find requested implementation and construct it
  for (auto entry : context_list)
  {
    if (!strcmp(entry.target, target) && !strcmp(entry.mode, mode))
    {
      return entry.constructor();
    }
  }

  std::cerr << std::endl
            << "No implementation found for " << target << "-" << mode
            << std::endl << std::endl;
  exit(1);
  return NULL;
}

void CGContext::list_contexts()
{
  std::cout << std::endl
            << "Registered contexts:"
            << std::endl;
  for (auto entry : context_list)
  {
    std::cout << "\t" << entry.target << "-" << entry.mode << std::endl;
  }
  std::cout << std::endl;
}
