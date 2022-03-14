#include <models/genotype_checker.h>

genotype_checker::genotype_checker(haplotype_set & _H) : H(_H) {
	Errors = vector < vector < bool > > (H.IDXesti.size(), vector < bool > (H.n_variants, false));
}

genotype_checker::~genotype_checker() {
	Errors.clear();
}

void genotype_checker::check() {
	tac.clock();
	vrb.title("Check genotyping discordances");
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		for (int l = 0 ; l < H.n_variants ; l ++) {
			Errors[i][l] = (!H.Missing[H.IDXesti[i]][l]) && ((H.Htrue[2*H.IDXesti[i]+0][l] + H.Htrue[2*H.IDXesti[i]+1][l]) != (H.Hesti[2*H.IDXesti[i]+0][l] + H.Hesti[2*H.IDXesti[i]+1][l]));
		}
	}

	unsigned int n_genotyping_errors = 0;
	for (int i = 0 ; i < Errors.size() ; i++) for (int l = 0 ; l < Errors[i].size() ; l ++) n_genotyping_errors += Errors[i][l];
	vrb.bullet("#Genotyping errors = " + stb.str(n_genotyping_errors));
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}

void genotype_checker::writePerSample(string fout) {
	tac.clock();
	vrb.title("Writing genotyping discordances per sample in [" + fout + "]");
	output_file fdo (fout);
	for (int i = 0 ; i < H.IDXesti.size() ; i++) {
		int n_errors = 0, n_nmissing = 0;
		for (int l = 0 ; l < H.n_variants ; l ++) {
			n_errors += Errors[i][l];
			n_nmissing += !H.Missing[H.IDXesti[i]][l];
		}
		fdo << H.vecSamples[H.IDXesti[i]] << " " << n_errors << " " << n_nmissing << " " << stb.str(n_errors * 100.0f / n_nmissing, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}


void genotype_checker::writePerVariant(string fout) {
	tac.clock();
	vrb.title("Writing genotyping discordances per variant in [" + fout + "]");
	output_file fdo (fout);
	for (int l = 0 ; l < H.n_variants ; l ++) {
		int n_errors = 0, n_nmissing = 0;
		for (int i = 0 ; i < H.IDXesti.size() ; i++) {
			n_errors += Errors[i][l];
			n_nmissing += !H.Missing[H.IDXesti[i]][l];
		}
		fdo << H.RSIDs[l]  << " " << H.Positions[l] << " " << n_errors << " " << n_nmissing << " " << stb.str(n_errors * 100.0f / n_nmissing, 2) << endl;
	}
	fdo.close();
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
}
