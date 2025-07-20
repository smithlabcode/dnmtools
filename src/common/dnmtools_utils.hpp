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

#ifndef DNMTOOLS_UTILS_HPP
#define DNMTOOLS_UTILS_HPP

#include <string>

auto
get_command_line(const int argc, char *argv[]) -> std::string;

#endif
