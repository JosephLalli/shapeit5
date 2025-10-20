#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include <htslib/vcf.h>
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.bcf> <output.bcf>" << std::endl;
        return EXIT_FAILURE;
    }

    const char* in_path = argv[1];
    const char* out_path = argv[2];

    htsFile* in = bcf_open(in_path, "rb");
    if (!in) {
        std::cerr << "Failed to open input" << std::endl;
        return EXIT_FAILURE;
    }
    bcf_hdr_t* hdr = bcf_hdr_read(in);
    if (!hdr) {
        std::cerr << "Failed to read header" << std::endl;
        bcf_close(in);
        return EXIT_FAILURE;
    }

    htsFile* out = bcf_open(out_path, "wb");
    if (!out) {
        std::cerr << "Failed to open output" << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(in);
        return EXIT_FAILURE;
    }
    if (bcf_hdr_write(out, hdr) != 0) {
        std::cerr << "Failed to write header" << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(in);
        bcf_close(out);
        return EXIT_FAILURE;
    }

    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    while (bcf_read(in, hdr, rec) == 0) {
        if (rec->d.id && std::strcmp(rec->d.id, "rs1_A_G") == 0) {
            int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
            if (ngt > 0) {
                const int nsamples = bcf_hdr_nsamples(hdr);
                const int ploidy = ngt / nsamples;
                if (ploidy >= 2) {
                    gt_arr[0] = bcf_gt_phased(0);
                    gt_arr[1] = bcf_gt_phased(bcf_gt_allele(gt_arr[1]));
                    bcf_update_genotypes(hdr, rec, gt_arr, ngt);
                }
            }
        }
        if (bcf_write(out, hdr, rec) != 0) {
            std::cerr << "Failed to write record" << std::endl;
            bcf_destroy(rec);
            free(gt_arr);
            bcf_hdr_destroy(hdr);
            bcf_close(in);
            bcf_close(out);
            return EXIT_FAILURE;
        }
    }

    free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    if (bcf_close(in) != 0 || bcf_close(out) != 0) {
        std::cerr << "Error closing files" << std::endl;
        return EXIT_FAILURE;
    }

    if (bcf_index_build3(out_path, NULL, 14, 1) != 0) {
        std::cerr << "Failed to build index" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
