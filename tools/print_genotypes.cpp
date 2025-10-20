#include <cstdlib>
#include <iostream>
#include <string>

extern "C" {
#include <htslib/vcf.h>
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <bcf>" << std::endl;
        return EXIT_FAILURE;
    }
    htsFile* fp = bcf_open(argv[1], "rb");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    bcf1_t* rec = bcf_init();
    int32_t* gt = nullptr; int ngt_arr = 0; int ns = bcf_hdr_nsamples(hdr);
    while (bcf_read(fp, hdr, rec) == 0) {
        std::cout << bcf_hdr_id2name(hdr, rec->rid) << ":" << rec->pos + 1 << "\t" << (rec->d.id ? rec->d.id : "") << "\t";
        int ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
        int ploidy = ngt / ns;
        std::cout << "S1=";
        for (int p = 0; p < ploidy; ++p) {
            int val = bcf_gt_allele(gt[p]);
            std::cout << val;
            if (p == 0) std::cout << "|";
        }
        std::cout << std::endl;
    }
    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return 0;
}
