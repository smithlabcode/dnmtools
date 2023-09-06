/* metagene (tsscpgplot): get data to plot methylation level around a TSS
 *
 * Copyright (C) 2010-2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <stdexcept>
#include <unordered_map>

#include "GenomicRegion.hpp"
#include "LevelsCounter.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::pair;
using std::runtime_error;
using std::sort;
using std::string;
using std::to_string;
using std::unordered_map;
using std::vector;

static MSite
tss_from_gene(const GenomicRegion &r) {
  MSite s;
  s.chrom = r.get_chrom();
  s.pos = (r.pos_strand() ? r.get_start() : r.get_end());
  s.strand = r.get_strand();
  return s;
}

static void
process_chrom(const uint32_t region_size, const vector<GenomicRegion> &genes,
              const pair<uint32_t, uint32_t> &bounds,
              const vector<MSite> &sites, vector<LevelsCounter> &levels) {
  const uint32_t twice_rs = 2 * region_size;

  auto comp = [](const MSite &a, const MSite &b) { return a.pos < b.pos; };

  for (auto i = bounds.first; i < bounds.second; ++i) {
    const MSite tss = tss_from_gene(genes[i]);
    MSite left = tss;
    left.pos = left.pos > region_size ? left.pos - region_size : 0;

    MSite right = tss;
    right.pos += region_size;

    const auto lim = upper_bound(cbegin(sites), cend(sites), right, comp);
    auto itr = lower_bound(cbegin(sites), cend(sites), left, comp);
    for (; itr != lim; ++itr) {
      const auto dist = itr->pos - left.pos;
      const auto base_idx = tss.strand == '+' ? dist : twice_rs - dist;
      levels[base_idx].update(*itr);
    }
  }
}

template<typename T> static void
collapse_bins(const uint32_t bin_size, vector<T> &v) {
  const uint32_t n_bins = std::ceil(static_cast<double>(v.size()) / bin_size);
  vector<T> vv(n_bins);
  for (auto i = 0u; i < v.size(); ++i) vv[i / bin_size] += v[i];
  v.swap(vv);
}

int
metagene(int argc, const char **argv) {
  constexpr auto description =
    "Compute the information needed for metagene plots of DNA methylation    \
     levels. The columns in the output correspond to the fields calculated   \
     globally by the `levels` and per-region by the `roi` command. Input     \
     for features is in BED format, and when present the 6th column is used  \
     to indicate strand. For features of non-zero width (where the 2nd and   \
     3rd columns are not identical) the negative strand will indicate that   \
     3rd column should be used. This means, for example, if the features are \
     genes, and the promoters are of interest, the strand will be used       \
     correctly.";

  try {
    string outfile;
    uint32_t region_size = 5000;
    bool verbose = false;
    bool show_progress = false;
    uint32_t bin_size = 50;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<features-bed> <counts>");
    opt_parse.add_opt("output", 'o', "output file (default: terminal)", false,
                      outfile);
    opt_parse.add_opt("size", 's', "analyze this size in both directions",
                      false, region_size);
    opt_parse.add_opt("bin", 'b', "bin size", false, bin_size);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string features_file_name = leftover_args.front();
    const string cpg_file_name = leftover_args.back();
    /**********************************************************************/

    if (verbose) cerr << "[loading feature annotations data]" << endl;
    vector<GenomicRegion> features;
    ReadBEDFile(features_file_name, features);
    sort(begin(features), end(features));
    if (verbose)
      cerr << "[number of features: " << features.size() << "]" << endl;

    // identify the start and end of ranges for each chromosome
    unordered_map<string, pair<uint32_t, uint32_t>> lookup;
    typedef decltype(lookup)::iterator lookup_itr;

    string chrom_name;
    auto prev_idx = 0u;
    for (auto i = 0u; i < features.size(); ++i)
      if (features[i].get_chrom() != chrom_name) {
        if (!chrom_name.empty()) lookup.insert({chrom_name, {prev_idx, i}});
        prev_idx = i;
        chrom_name = features[i].get_chrom();
      }
    lookup.insert({chrom_name, {prev_idx, features.size()}});
    if (verbose)
      cerr << "[number of chroms with features: " << lookup.size() << "]\n";

    auto pair_diff = [&lookup](const lookup_itr x) {
      return (x != end(lookup)) ? x->second.second - x->second.first : 0u;
    };

    vector<LevelsCounter> levels(2 * region_size);

    bamxx::bgzf_file cpgin(cpg_file_name, "r");
    if (!cpgin) throw runtime_error("failed to open file: " + cpg_file_name);

    uint32_t total_features = 0u;

    vector<MSite> sites;
    string line;
    chrom_name.clear();
    while (getline(cpgin, line)) {
      auto the_site = MSite(line);
      if (the_site.chrom != chrom_name) {
        if (!sites.empty()) {
          const auto bounds = lookup.find(chrom_name);
          if (bounds != end(lookup))
            process_chrom(region_size, features, bounds->second, sites, levels);
          const auto n_features = pair_diff(bounds);
          if (show_progress)
            cerr << "[sites=" << sites.size() << " features=" << n_features
                 << "]" << endl;
          total_features += n_features;
          sites.clear();
        }
        if (show_progress) cerr << "[processing: " << the_site.chrom << "]";
        chrom_name = the_site.chrom;
      }
      sites.push_back(std::move(the_site));
    }

    if (!sites.empty()) {
      const auto bounds = lookup.find(chrom_name);
      if (bounds != end(lookup))
        process_chrom(region_size, features, bounds->second, sites, levels);
      const auto n_features = pair_diff(bounds);
      if (show_progress)
        cerr << "[sites=" << sites.size() << " features=" << n_features << "]"
             << endl;
      total_features += n_features;
    }

    collapse_bins(bin_size, levels);

    if (verbose)
      cerr << "output columns:\n"
           << LevelsCounter::tostring_as_row_header() << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(of.is_open() ? of.rdbuf() : std::cout.rdbuf());

    for (auto i = 0u; i < levels.size(); ++i)
      out << i * bin_size << '\t' << levels[i].tostring_as_row() << endl;
  }
  catch (std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
