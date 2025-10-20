#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

extern "C" {
#include <htslib/vcf.h>
}

void write_record(bcf_hdr_t* hdr, htsFile* fp, const char* chrom, int pos, const char* id,
                  const char* ref, const char* alt, int allele0, int allele1,
                  int nsamples, int ac_value) {
    bcf1_t* rec = bcf_init();
    if (!rec) {
        std::cerr << "Failed to allocate record" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    rec->rid = bcf_hdr_name2id(hdr, chrom);
    if (rec->rid < 0) {
        std::cerr << "Unknown contig: " << chrom << std::endl;
        std::exit(EXIT_FAILURE);
    }
    rec->pos = pos - 1;

    bcf_update_id(hdr, rec, id);
    std::string alleles = std::string(ref) + "," + alt;
    if (bcf_update_alleles_str(hdr, rec, alleles.c_str()) != 0) {
        std::cerr << "Failed to set alleles" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::vector<int32_t> g(static_cast<std::size_t>(nsamples) * 2);
    for (int s = 0; s < nsamples; ++s) {
        if (s == 0) {
            g[2 * s] = bcf_gt_phased(allele0);
            g[2 * s + 1] = bcf_gt_phased(allele1);
        } else {
            g[2 * s] = bcf_gt_phased(0);
            g[2 * s + 1] = bcf_gt_phased(0);
        }
    }
    if (bcf_update_genotypes(hdr, rec, g.data(), static_cast<int>(g.size())) != 0) {
        std::cerr << "Failed to set genotypes" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    int32_t ac = ac_value;
    int32_t an = nsamples * 2;
    if (bcf_update_info_int32(hdr, rec, "AC", &ac, 1) != 0) {
        std::cerr << "Failed to set AC" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (bcf_update_info_int32(hdr, rec, "AN", &an, 1) != 0) {
        std::cerr << "Failed to set AN" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (bcf_write(fp, hdr, rec) != 0) {
        std::cerr << "Failed to write record" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    bcf_destroy(rec);
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <output.bcf> [mode]" << std::endl;
        return EXIT_FAILURE;
    }

    const char* mode = (argc >= 3) ? argv[2] : "violation";
    bool fixed = std::string(mode) == "fixed";

    const int nsamples = 60;

    const char* out = argv[1];
    htsFile* fp = bcf_open(out, "wb");
    if (!fp) {
        std::cerr << "Failed to open output: " << out << std::endl;
        return EXIT_FAILURE;
    }

    bcf_hdr_t* hdr = bcf_hdr_init("w");
    bcf_hdr_append(hdr, "##fileformat=VCFv4.3");
    bcf_hdr_append(hdr, "##contig=<ID=1,length=1000000>");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ALT allele count\">");
    bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">");

    for (int i = 0; i < nsamples; ++i) {
        std::ostringstream name;
        name << "S" << (i + 1);
        if (bcf_hdr_add_sample(hdr, name.str().c_str()) != 0) {
            std::cerr << "Failed to add sample" << std::endl;
            return EXIT_FAILURE;
        }
    }
    bcf_hdr_add_sample(hdr, NULL);

    if (bcf_hdr_write(fp, hdr) != 0) {
        std::cerr << "Failed to write header" << std::endl;
        return EXIT_FAILURE;
    }

    write_record(hdr, fp, "1", 1000, "rs1_A_T", "A", "T", 1, 1, nsamples, 2);
    if (fixed) {
        write_record(hdr, fp, "1", 1000, "rs1_A_G", "A", "G", 0, 0, nsamples, 0);
    } else {
        write_record(hdr, fp, "1", 1000, "rs1_A_G", "A", "G", 1, 0, nsamples, 1);
    }
    write_record(hdr, fp, "1", 2000, "rs2_C_G", "C", "G", 0, 1, nsamples, 1);

    bcf_hdr_destroy(hdr);
    if (bcf_close(fp) != 0) {
        std::cerr << "Failed to close output" << std::endl;
        return EXIT_FAILURE;
    }

    if (bcf_index_build3(out, NULL, 14, 1) != 0) {
        std::cerr << "Failed to build index" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
