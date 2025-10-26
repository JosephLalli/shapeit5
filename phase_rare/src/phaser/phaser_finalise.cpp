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

void phaser::write_files_and_finalise() {
    vrb.title("Finalization:");

	//step0: multi-threading
	if (options["thread"].as < int > () > 1) pthread_mutex_destroy(&mutex_workers);

	//step1: writing best guess haplotypes in VCF/BCF file
	haplotype_writer writerH (H, G, V, options["thread"].as < int > ());
	writerH.setRegions(input_start, input_stop);
	writerH.writeHaplotypes(options["output"].as < std::string > (), options["input"].as < std::string > ());

    //step2: Measure overall running time
    vrb.bullet("Total running time = " + stb.str(tac.abs_time()) + " seconds");

    if (enforce_oneallele_rare) {
        vrb.bullet("One-allele (rare) : positions=" + stb.str(oneallele_rare_stats.positions_checked) +
                   " / sample-violations=" + stb.str(oneallele_rare_stats.sample_violations_found) +
                   " / extreme-violations=" + stb.str(oneallele_rare_stats.extreme_violations_found) +
                   " / flips=" + stb.str(oneallele_rare_stats.flips_applied));
        
        // Report enhanced metrics if available (non-basic modes)
        if (oneallele_rare_mode != OneAlleleRareMode::PP_BASIC && 
            (oneallele_rare_stats.sparse_donor_resolutions > 0 || oneallele_rare_stats.li_stephens_enhanced > 0)) {
            vrb.bullet("Enhanced resolution : donor-weighted=" + stb.str(oneallele_rare_stats.sparse_donor_resolutions) +
                       " / li-stephens=" + stb.str(oneallele_rare_stats.li_stephens_enhanced) +
                       " / genotype-changes=" + stb.str(oneallele_rare_stats.genotype_changes) +
                       " / phase-only=" + stb.str(oneallele_rare_stats.phase_only_changes));
        }
        if (!oneallele_rare_stats_path.empty()) {
            std::ofstream ofs(oneallele_rare_stats_path);
            if (!ofs.good()) {
                vrb.warning("Unable to write rare one-allele stats to [" + oneallele_rare_stats_path + "]");
            } else {
                // Base metrics (existing)
                ofs << "positions_checked\t" << oneallele_rare_stats.positions_checked << "\n";
                ofs << "sample_violations_found\t" << oneallele_rare_stats.sample_violations_found << "\n";
                ofs << "extreme_violations_found\t" << oneallele_rare_stats.extreme_violations_found << "\n";
                ofs << "flips_applied\t" << oneallele_rare_stats.flips_applied << "\n";
                
                // Enhanced metrics (new)
                ofs << "sparse_donor_resolutions\t" << oneallele_rare_stats.sparse_donor_resolutions << "\n";
                ofs << "pp_only_resolutions\t" << oneallele_rare_stats.pp_only_resolutions << "\n";
                ofs << "genotype_changes\t" << oneallele_rare_stats.genotype_changes << "\n";
                ofs << "phase_only_changes\t" << oneallele_rare_stats.phase_only_changes << "\n";
                ofs << "complex_enumeration_cases\t" << oneallele_rare_stats.complex_enumeration_cases << "\n";
                ofs << "li_stephens_enhanced\t" << oneallele_rare_stats.li_stephens_enhanced << "\n";
            }
        }
    }
}
