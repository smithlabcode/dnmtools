/* Copyright (C) 2019-2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "dnmtools_utils.hpp"

#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

using std::string;
using std::copy;
using std::ostringstream;
using std::ostream_iterator;

auto
get_command_line(const int argc, const char **const argv) -> string {
  if (argc == 0) return string();
  std::ostringstream cmd;
  cmd << '"';
  copy(argv, argv + (argc - 1), ostream_iterator<const char *>(cmd, " "));
  cmd << argv[argc - 1] << '"';
  return cmd.str();
}
