/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <phaser/phaser_header.h>

#include <io/haplotype_writer.h>
#include <io/graph_writer.h>

#include <fstream>

void phaser::write_files_and_finalise() {
    vrb.title("Finalization:");

    //step0: multi-threading
    if (options["thread"].as < int > () > 1) pthread_mutex_destroy(&mutex_workers);

    // Step 1: Final PBWT solve
    G.solve();

    // Step 2: Final multiallelic constraint enforcement (BEFORE H update)
    // This ensures G contains corrected genotypes before copying to H
    if (oneallele_enforcer.enabled() && multiallelic_map.size() > 0) {
        tac.clock();
        oneallele_enforcer.reset_epoch_stats();
        oneallele_enforcer.enforce(multiallelic_map, G, V, "final");
        
        const auto& epoch_stats = oneallele_enforcer.epoch_stats();
        if (epoch_stats.violations_found > 0 || epoch_stats.flips_applied > 0) {
            std::string mode_str;
            switch (oneallele_enforcer.mode()) {
                case shapeit5::modules::OneAlleleMode::TRANSITION: mode_str = "transition"; break;
                case shapeit5::modules::OneAlleleMode::MICRO: mode_str = "micro"; break;
                case shapeit5::modules::OneAlleleMode::MICRO_DONOR: mode_str = "micro-donor"; break;
            }
            vrb.bullet("Final multiallelic correction (" + mode_str + ") [violations=" + 
                      stb.str(epoch_stats.violations_found) + " / flipped=" + 
                      stb.str(epoch_stats.flips_applied) + "] (" + 
                      stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
        }
    }

    // Step 3: Now copy corrected G into H and transpose for output
    H.updateHaplotypes(G);
    H.transposeHaplotypes_H2V(false);

    //step1: writing best guess haplotypes in VCF/BCF file
    std::string oformat = options["output-format"].as < std::string > ();
    if (oformat == "graph")
        graph_writer(G, V).writeGraphs(options["bingraph"].as < std::string > ());
    else
        haplotype_writer(H, G, V, options["thread"].as < int > ()).writeHaplotypes(options["output"].as < std::string > (), oformat , options["input"].as < std::string > ());

    // Step 4: Report cumulative statistics
    if (enforce_oneallele) {
        const shapeit5::modules::OneAlleleStats& stats = oneallele_enforcer.stats();
        vrb.bullet("One-allele enforcement : positions=" + stb.str(stats.positions_checked) + " / violations=" + stb.str(stats.violations_found) + " / flips=" + stb.str(stats.flips_applied));
        if (!oneallele_stats_path.empty()) {
            std::ofstream ofs(oneallele_stats_path);
            if (!ofs.good()) {
                vrb.warning("Unable to write one-allele stats to [" + oneallele_stats_path + "]");
            } else {
                ofs << "positions_checked\t" << stats.positions_checked << "\n";
                ofs << "violations_found\t" << stats.violations_found << "\n";
                ofs << "flips_applied\t" << stats.flips_applied << "\n";
            }
        }
    }

    //step2: Measure overall running time
    vrb.bullet("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}