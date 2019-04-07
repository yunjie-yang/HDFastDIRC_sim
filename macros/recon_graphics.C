void recon_graphics()
{
	gStyle->SetOptStat(0);

	string infile_dir   = "/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/data_bggen";
	string outplot_dir  = "/media/sf_SharedFolderVM/FastDIRC_geometry/plots/recon";

	int case_number = 7;

	int rebin = 400;
	double range = 25;
	
	string infile_label;
	string plot_label;
	string plot_title;

	string pion_hist = "ll_diff_pion_tree";
	string kaon_hist = "ll_diff_kaon_tree";

	switch (case_number)
	{
		case 1:
			infile_label = "data_g4_5000";
			plot_label = "data_g4_5000";
			plot_title = "data, p > 4 GeV (Nph/pdf = 20k)";
			break;
		case 2:
			infile_label = "data_g4_25000";
			plot_label = "data_g4_25000";
			plot_title = "data, p > 4 GeV (Nph/pdf = 100k)";
			break;
		case 3:
			infile_label = "data_l3_5000";
			plot_label = "data_l3_5000";
			plot_title = "data, p < 3 GeV (Nph/pdf = 20k)";
			rebin = 2000;
			range = 200;
			break;
		case 4:
			infile_label = "bggen_g4_5000";
			plot_label = "bggen_g4_5000";
			plot_title = "bggen, p > 4 GeV (Nph/pdf = 20k)";
			break;
		case 5:
			infile_label = "bggen_g4_25000";
			plot_label = "bggen_g4_25000";
			plot_title = "bggen, p > 4 GeV (Nph/pdf = 100k)";
			break;
		case 6:
			infile_label = "data_3_25000";
			plot_label = "data_3_25000_generate_hits";
			plot_title = "data, p:[2.8,3.2] GeV (Nph/pdf = 100k) with generated hits";
			range = 60;

			//infile_label = "data_4_25000";
			//plot_label = "data_4_25000_generate_hits";
			//plot_title = "data, p:[3.8,4.2] GeV (Nph/pdf = 100k) with generated hits";

			//infile_label = "data_g4_25000";
			//plot_label = "data_g4_25000_generate_hits";
			//plot_title = "data, p > 4 GeV (Nph/pdf = 100k) with generated hits";
			pion_hist = "ll_diff_pion";
			kaon_hist = "ll_diff_kaon";
			break;
		case 7:
			infile_label = "data_3_25000";
			plot_label = "data_3_25000";
			plot_title = "data, p: [2.8, 3.2] GeV (Nph/pdf = 100k)";
			range = 60;

			//infile_label = "data_4_25000";
			//plot_label = "data_4_25000";
			//plot_title = "data, p: [3.8, 4.2] GeV (Nph/pdf = 100k)";
			break;
	}


	string infile_name  = "hist_"+ infile_label;
	string outplot_name = "dll_" + plot_label;

	TFile* infile = new TFile(Form("%s/%s.root",infile_dir.c_str(),infile_name.c_str()));

	TH1F* dll_pion = (TH1F*) infile -> Get(pion_hist.c_str());
	TH1F* dll_kaon = (TH1F*) infile -> Get(kaon_hist.c_str());

	dll_pion -> SetLineColor(kRed);	
	dll_kaon -> SetLineColor(kBlue);

	dll_pion -> SetMarkerColor(kRed);	
	dll_kaon -> SetMarkerColor(kBlue);

	dll_pion -> SetMarkerColor(kRed);	
	dll_kaon -> SetMarkerColor(kBlue);

	dll_pion -> SetFillColorAlpha(kRed,0.5);	
	dll_kaon -> SetFillColorAlpha(kBlue,0.5);

	dll_pion -> Rebin(rebin);	
	dll_kaon -> Rebin(rebin);

	dll_pion -> Scale(1./dll_pion->GetMaximum());
	dll_kaon -> Scale(1./dll_kaon->GetMaximum());

	dll_pion -> SetTitle(plot_title.c_str());
	dll_pion -> GetXaxis() -> SetTitle("DLL");
	dll_pion -> GetYaxis() -> SetTitle("a.u.");
	dll_pion -> GetXaxis() -> SetRangeUser(-range,range);


	TLegend* leg = new TLegend(0.75,0.75,0.9,0.9);
	leg->AddEntry(dll_pion,"#pi");	
	leg->AddEntry(dll_kaon,"K");	

	TCanvas* canv = new TCanvas("canv","canv",600,400);
	dll_pion->Draw("HIST");
	dll_kaon->Draw("HIST SAME");
	leg->Draw();
	gPad->RedrawAxis();

	canv->Print(Form("%s/%s.pdf",outplot_dir.c_str(),outplot_name.c_str()));




}
