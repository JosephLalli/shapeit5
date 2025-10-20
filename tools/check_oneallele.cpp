#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include <htslib/vcf.h>
}

struct Counters {
    uint64_t positions = 0;
    uint64_t hap_violations = 0;
};

void evaluate_counts(const std::vector<int>& hap_counts, Counters& ctrs) {
    if (hap_counts.empty()) return;
    ctrs.positions++;
    const std::size_t samples = hap_counts.size() / 2;
    for (std::size_t i = 0; i < samples; ++i) {
        if (hap_counts[2 * i] > 1) ctrs.hap_violations++;
        if (hap_counts[2 * i + 1] > 1) ctrs.hap_violations++;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.bcf>" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string filename = argv[1];
    htsFile* fp = bcf_open(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Failed to open BCF: " << filename << std::endl;
        return EXIT_FAILURE;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Failed to read BCF header" << std::endl;
        bcf_close(fp);
        return EXIT_FAILURE;
    }

    const int nsamples = bcf_hdr_nsamples(hdr);
    if (nsamples <= 0) {
        std::cerr << "No samples found" << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return EXIT_FAILURE;
    }

    std::vector<int> hap_counts(static_cast<std::size_t>(nsamples) * 2, 0);
    Counters ctrs;

    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Failed to allocate BCF record" << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return EXIT_FAILURE;
    }

    int last_rid = -1;
    int32_t last_pos = -1;

    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    while (bcf_read(fp, hdr, rec) == 0) {
        if (rec->rid != last_rid || rec->pos != last_pos) {
            if (last_rid != -1) {
                evaluate_counts(hap_counts, ctrs);
            }
            std::fill(hap_counts.begin(), hap_counts.end(), 0);
            last_rid = rec->rid;
            last_pos = rec->pos;
        }

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt < 0) {
            std::cerr << "Warning: failed to get genotypes at pos " << rec->pos + 1 << std::endl;
            continue;
        }
        const int ploidy = ngt / nsamples;
        if (ploidy < 2) {
            std::cerr << "Warning: ploidy <2 at pos " << rec->pos + 1 << std::endl;
            continue;
        }

        for (int s = 0; s < nsamples; ++s) {
            int idx0 = s * ploidy;
            int idx1 = idx0 + 1;
            if (idx1 >= ngt) break;
            int allele0 = bcf_gt_allele(gt_arr[idx0]);
            int allele1 = bcf_gt_allele(gt_arr[idx1]);
            if (!bcf_gt_is_missing(gt_arr[idx0]) && allele0 > 0) hap_counts[2 * s] += 1;
            if (!bcf_gt_is_missing(gt_arr[idx1]) && allele1 > 0) hap_counts[2 * s + 1] += 1;
        }
    }

    if (last_rid != -1) {
        evaluate_counts(hap_counts, ctrs);
    }

    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    std::cout << "positions=" << ctrs.positions
              << "\thap_violations=" << ctrs.hap_violations << std::endl;

    return EXIT_SUCCESS;
}
