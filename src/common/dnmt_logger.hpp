/* Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef DNMT_LOGGER_HPP
#define DNMT_LOGGER_HPP

#include <string>
#include <ostream>
#include <chrono>

struct dnmt_logger {
  std::ostream &log_stream;
  std::string prefix;
  dnmt_logger(std::ostream &ls, std::string &&p) log_stream{ls}, prefix{} {}
  dnmt_logger(std::ostream &ls, std::string &&p) log_stream{ls}, prefix{p} {}
  std::chrono::steady_clock log_event(std::string message);
  std::chrono::steady_clock log_event(std::string message,
                                      std::chrono::steady_clock sclk);
  std::chrono::steady_clock log_data(std::string key, std::string value);
  std::chrono::steady_clock log_data(std::string key, std::string value,
                                     std::chrono::steady_clock sclk);
};

#endif
