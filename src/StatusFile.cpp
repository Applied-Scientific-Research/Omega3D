/*
 * StatusFile.cpp - Class to contain and compose a line-per-step status file
 *
 * (c)2019,21 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "StatusFile.h"

#ifdef _WIN32
  #include <ciso646>
#endif

#include <iostream>
#include <fstream>
#include <cassert>

// are we even using the status file?
bool
StatusFile::is_active() {
  return use_it;
}

// accept a filename, means that we should use it
void
StatusFile::set_filename(const std::string _fn) {
  assert(not _fn.empty() && "Filename is blank");
  use_it = true;
  fn = _fn;
}

std::string
StatusFile::get_filename() {
  return fn;
}

// begin writing a new data set to the file
void
StatusFile::reset_sim() {
  if (use_it) {
    num_lines = 0;
  }
}

// append a title and an int or float to the list to write
void
StatusFile::append_value(const std::string _name, const float _val) {
  names.push_back(_name);
  vals.push_back(_val);
}
void
StatusFile::append_value(const std::string _name, const int _val) {
  names.push_back(_name);
  vals.push_back(_val);
}

// append an int or float to the list to write
void
StatusFile::append_value(const float _val) {
  names.push_back("float");
  vals.push_back(_val);
}
void
StatusFile::append_value(const int _val) {
  names.push_back("int");
  vals.push_back(_val);
}

// write a line
void
StatusFile::write_line() {
  if (use_it) {
    std::ofstream outfile;
    assert(not fn.empty() && "Filename is blank");
    outfile.open(fn, std::ios::app);

    // is this a new data set?
    if (num_lines == 0) {
      // new data set needs white space
      if (num_sims > 0) outfile << std::endl;

      // write header text line
      if (format != csv) outfile << "# ";
      for (size_t i=0; i<names.size(); ++i) {
        outfile << names[i];
        if (i == names.size()-1) break;
        if (format == csv) outfile << ",";
        else outfile << " ";
      }
      outfile << std::endl;

      // keep count
      num_sims++;
    }

    // write data line
    for (size_t i=0; i<vals.size(); ++i) {
      auto& val = vals[i];
      std::visit([&outfile](const auto& v) { outfile << v; }, val);
      if (i == vals.size()-1) {
        outfile << std::endl;
      } else {
        if (format == csv) outfile << ",";
        else outfile << " ";
      }
    }
    num_lines++;
  }

  // empty the vector
  vals.clear();
}

// read from "runtime" json object
void
StatusFile::from_json(const nlohmann::json j) {

  if (j.find("statusFile") != j.end()) {
    std::string sfile = j["statusFile"];
    std::cout << "  status file name= " << sfile << std::endl;
    set_filename(sfile);
  }
  if (j.find("statusFormat") != j.end()) {
    std::string formatstr = j["statusFormat"];
    std::cout << "  status file format= " << formatstr << std::endl;
    if (formatstr == "csv") format = csv;
    else format = dat;
  }
}

// write to "runtime" json object
void
StatusFile::add_to_json(nlohmann::json& j) const {

  if (not fn.empty()) {
    j["statusFile"] = fn;

    // since we only have dat (default) and csv, this is easy
    if (format != dat) j["statusFormat"] = "csv";
  }
}

