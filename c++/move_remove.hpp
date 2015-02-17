/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include <algorithm>
#include "./qmc_data.hpp"

namespace cthyb {

// Removal of C, C^dagger operator
class move_remove_c_cdag {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;
 int block_index, block_size;
 bool record_histograms;
 std::map<std::string, statistics::histogram_segment_bin> histos; // Analysis histograms
 double delta_tau;
 qmc_data::trace_t new_trace;
 time_pt tau1, tau2;

 public:
 
 //----------------------------------

 move_remove_c_cdag(int block_index, int block_size, qmc_data& data, mc_tools::random_generator& rng, bool record_histograms)
    : data(data),
      config(data.config),
      rng(rng),
      block_index(block_index),
      block_size(block_size),
      record_histograms(record_histograms) {
  if (record_histograms) {
   histos.insert({"remove_length_proposed", {0, config.beta(), 100, "histo_remove_length_proposed.dat"}});
   histos.insert({"remove_length_accepted", {0, config.beta(), 100, "histo_remove_length_accepted.dat"}});
  }
 }

 //----------------
 
 mc_weight_type attempt() {

 auto & proposed_data = data.proposed_data[block_index];
 proposed_data.active = false;

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
  std::cerr << "* Configuration before:" << std::endl << config;
  data.imp_trace.tree.graphviz(std::ofstream("tree_before"));
#endif

  auto& det = data.dets[block_index];

  // Pick up a couple of C, Cdagger to remove at random
  // Remove the operators from the traces
  int det_size = det.size();
  if (det_size == 0) return 0; // nothing to remove
  int num_c_dag = rng(det_size), num_c = rng(det_size);
#ifdef EXT_DEBUG
  std::cerr << "* Proposing to remove: ";
  std::cerr << num_c_dag << "-th Cdag(" << block_index << ",...), ";
  std::cerr << num_c << "-th C(" << block_index << ",...)" << std::endl;
#endif

  // now mark 2 nodes for deletion
  int rs1, rs2;
  std::tie(tau1,rs1) = data.imp_trace.try_delete(num_c, block_index, false);
  std::tie(tau2,rs2) = data.imp_trace.try_delete(num_c_dag, block_index, true);

  // record the length of the proposed removal
  delta_tau = double(tau2 - tau1);
  if (record_histograms) histos["remove_length_proposed"] << delta_tau;

  auto det_ratio = det.try_remove(num_c_dag, num_c);

  // proposition probability
  double t_ratio = std::pow(block_size * config.beta() / double(det_size), 2); // Size of the det before the try_delete!

  // For quick abandon 
  double random_number = rng.preview();
  if (random_number == 0.0) return 0;
  double p_yee = std::abs(det_ratio / t_ratio / data.trace);

  // recompute the trace
  if (data.use_proposed) new_trace = data.imp_trace.estimate(p_yee, 0.0);
  else new_trace = data.imp_trace.estimate(p_yee, random_number);

  if (new_trace == 0.0) {
#ifdef EXT_DEBUG
   std::cout << "trace == 0" << std::endl;
#endif
   return 0;
  }
  auto trace_ratio = new_trace / data.trace;
  if (!std::isfinite(trace_ratio)) TRIQS_RUNTIME_ERROR << "trace_ratio not finite" << new_trace << "  "<< data.trace<<"  "<< new_trace /data.trace ;

  mc_weight_type p = trace_ratio * det_ratio;
  mc_weight_type p_total = p / t_ratio;

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << trace_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p_total << std::endl;
#endif

  if (!std::isfinite(p)) TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p;
  if (!std::isfinite(p_total)) TRIQS_RUNTIME_ERROR << "p / t_ratio not finite p : " << p << " t_ratio :  "<< t_ratio;

  // Proposed move measurement
  if (data.use_proposed && det_size == 1 && std::abs(p_total) > 0.0) {
   proposed_data.active = true;
   proposed_data.tau1 = tau1;
   proposed_data.tau2 = tau2;
   proposed_data.inner1 = rs1;
   proposed_data.inner2 = rs2;
   proposed_data.acceptance = 1 - std::min(1.0,std::abs(p_total));
   proposed_data.delta = 1.0/det_ratio;
  }

  return p_total;
 }

 //----------------

 mc_weight_type accept() {

  // remove from the tree
  data.imp_trace.confirm_delete();

  // remove from the configuration
  config.erase(tau1);
  config.erase(tau2);

  // remove from the determinants
  data.dets[block_index].complete_operation();
  data.update_sign();
  data.trace = new_trace;
  if (record_histograms) histos["remove_length_accepted"] << delta_tau;

#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

  return data.current_sign / data.old_sign;
 }

 //----------------

 void reject() {
  // No Proposed move measurement
  data.proposed_data[block_index].active = false;

  data.imp_trace.cancel_delete();
#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

 }
};
}
