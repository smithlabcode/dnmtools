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

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>

using std::copy;
using std::ostream_iterator;
using std::ostringstream;
using std::string;

auto
get_command_line(const int argc,
                 char *argv[]) -> std::string {  // NOLINT(*-c-arrays)
  if (argc == 0)
    return std::string{};
  std::ostringstream cmd;
  cmd << '"';
  // NOLINTBEGIN(*-pointer-arithmetic)
  copy(argv, argv + (argc - 1), ostream_iterator<const char *>(cmd, " "));
  cmd << argv[argc - 1] << '"';
  // NOLINTEND(*-pointer-arithmetic)
  return cmd.str();
}
