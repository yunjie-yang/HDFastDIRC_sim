#include "plot_func.C"
void draw_DLLs()
{
	gStyle->SetOptStat(0);

	string infile_dir   = "/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/data_bggen";
	string outplot_dir  = "/media/sf_SharedFolderVM/FastDIRC_geometry/plots/recon";


	
	string infile_label;
	string plot_label;

	string pion_hist = "ll_diff_pion_tree";
	string kaon_hist = "ll_diff_kaon_tree";


	//string infile_name  = "hist_data_60785_000_015";
	//string plot_label_extra = "cut";

/*
	string infile_name  = "hist_data_beam_tree_23_all";
	string plot_label_extra = "Roman_tree";

	std::vector<std::vector<double>> PBins =  {
								{0.,1.5},
								{1.5, 2.5},
								{2.5, 3.5},
								{3.5, 4.5},
								{4.5, 12.}
						  };
	vector<int> rebins = {1000,400,1000,500,400};
	vector<int> ranges = {200,150,100,50,25};
*/
	string infile_name  = "hist_data_beam_tree_23_cut_bin";
	string plot_label_extra = "Roman_tree_bin";

	std::vector<std::vector<double>> PBins =  {
								{3.9, 4.15},
						  };
	vector<int> rebins = {500};
	vector<int> ranges = {60};


	TFile* infile = new TFile(Form("%s/%s.root",infile_dir.c_str(),infile_name.c_str()));
	string plot_title;
	for (size_t locPBin = 0; locPBin < PBins.size() ; locPBin++)
	{
		plot_title = Form("P: [%.2f, %.2f] GeV",PBins[locPBin][0],PBins[locPBin][1]);

		TH1F* hist_ll_diff_pion_tree = (TH1F*) infile -> Get(Form("ll_diff_pion_tree_PBin%d",int(locPBin)));
		TH1F* hist_ll_diff_kaon_tree = (TH1F*) infile -> Get(Form("ll_diff_kaon_tree_PBin%d",int(locPBin)));

		TH1F* hist_ll_diff_pion_sim_exact = (TH1F*) infile -> Get(Form("ll_diff_pion_sim_exact_PBin%d",int(locPBin)));
		TH1F* hist_ll_diff_kaon_sim_exact = (TH1F*) infile -> Get(Form("ll_diff_kaon_sim_exact_PBin%d",int(locPBin)));

		TCanvas* canv = new TCanvas("canv","canv",1000,400);
		canv->Divide(2,1);
		canv->cd(1);
		draw_compare_DLL(hist_ll_diff_pion_tree,hist_ll_diff_kaon_tree,rebins[locPBin],ranges[locPBin],plot_title);
		canv->cd(2);
		draw_compare_DLL(hist_ll_diff_pion_sim_exact,hist_ll_diff_kaon_sim_exact,rebins[locPBin],ranges[locPBin],plot_title + " (optimal)");

		canv->Print(Form("%s/dll_PBin%d_%s.pdf",outplot_dir.c_str(),int(locPBin),plot_label_extra.c_str()));
		delete canv;
	
	}




}
