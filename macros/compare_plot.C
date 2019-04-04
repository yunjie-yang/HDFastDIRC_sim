void bwr_palette()
{
	const int NumPoints = 1000;
	static Int_t bwr_palette[NumPoints];
	static Bool_t initialized = kFALSE;

	Double_t Red[]    = {0., 1.0, 1.0};
	Double_t Green[]  = {0., 1.0, 0.0};
	Double_t Blue[]   = {1., 1.0, 0.0};
	Double_t Length[] = {0., .50, 1.0};
	if (!initialized)
	{
		Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue, NumPoints);
		for (int i=0;i<NumPoints;i++) bwr_palette[i] = FI+i;
		initialized = kTRUE;
		return;
	}
	gStyle->SetPalette(NumPoints,bwr_palette);

}
void draw_compare_hit_patterns(TFile* file_fastdirc, TFile* file_geant, string output_dir, string title_label)
{
	
	TH2F* hist_fastdirc = (TH2F*) file_fastdirc->Get("pion_dist_rowcol");
	TH2F* hist_geant    = (TH2F*) file_geant->Get("DIRC_truth/hTruthPixelHit_South");

	hist_fastdirc -> Scale(1./hist_fastdirc->GetEntries());
	hist_geant    -> Scale(1./hist_geant->GetEntries());

	TH2F* hist_ratio_fg    = (TH2F*) hist_geant->Clone();
	TH2F* hist_ratio_gf    = (TH2F*) hist_geant->Clone();
	for (int bin_x = 1; bin_x <= hist_ratio_fg->GetXaxis()->GetNbins(); bin_x++)
	for (int bin_y = 1; bin_y <= hist_ratio_fg->GetYaxis()->GetNbins(); bin_y++)
	{
		double content_fastdirc = hist_fastdirc->GetBinContent(bin_x,bin_y);
		double content_geant    = hist_geant->GetBinContent(bin_x,bin_y);
		double ratio_fg = 0.;
		double ratio_gf = 0.;
		if (content_fastdirc > 1e-5 && content_geant!=0.)
		{
			ratio_fg = content_fastdirc/content_geant;
		}
		if (content_geant > 1e-5 && content_fastdirc!=0.)
		{
			ratio_gf = content_geant/content_fastdirc;
		}
		hist_ratio_fg->SetBinContent(bin_x,bin_y,ratio_fg);
		hist_ratio_gf->SetBinContent(bin_x,bin_y,ratio_gf);
	}


	hist_ratio_fg->SetMaximum(2.);
	hist_ratio_fg->SetMinimum(0.);

	hist_ratio_gf->SetMaximum(2.);
	hist_ratio_gf->SetMinimum(0.);

	hist_ratio_fg->GetZaxis()->SetNdivisions(1010);
	hist_ratio_gf->GetZaxis()->SetNdivisions(1010);

	hist_fastdirc -> SetTitle(Form("FastDIRC %s",title_label.c_str()));
	hist_geant -> SetTitle(Form("HDGeant4 %s",title_label.c_str()));
	hist_ratio_fg -> SetTitle(Form("FastDIRC/HDGeant4 %s",title_label.c_str()));
	hist_ratio_gf -> SetTitle(Form("HDGeant4/FastDIRC %s",title_label.c_str()));

	TCanvas* canv = new TCanvas("canv","canv",2400,1600);
	canv->Divide(2,2);

	gStyle->SetPalette();

	canv->cd(1);
	gPad->SetRightMargin(0.15);
	hist_fastdirc -> Draw("colz");
	canv->cd(2);
	gPad->SetRightMargin(0.15);
	hist_geant -> Draw("colz");

	TExec *ex1 = new TExec("ex1","bwr_palette();");
	ex1->Draw();
	canv->cd(3);
	gPad->SetRightMargin(0.15);
	hist_ratio_fg -> Draw("colz");
	canv->cd(4);
	gPad->SetRightMargin(0.15);
	hist_ratio_gf -> Draw("colz");

	canv->Print(Form("%s/%s.pdf",output_dir.c_str(),title_label.c_str()));
	delete canv;

}

void draw_compare_timing(TFile* file_fastdirc, TFile* file_geant, string output_dir, string title_label)
{
	
	TH1F* hist_fastdirc      = (TH1F*) file_fastdirc->Get("pion_dist_t");
	TH1F* hist_geant         = (TH1F*) file_geant->Get("DIRC_truth/hTruthPixelHitTime_1D_South");
	TH1F* hist_geant_tfixed  = (TH1F*) file_geant->Get("DIRC_truth/hTruthPixelHitTime_1D_tfixed_South");


	hist_fastdirc        -> Scale(1./hist_fastdirc->GetMaximum());
	hist_geant           -> Scale(1./hist_geant->GetMaximum());
	hist_geant_tfixed    -> Scale(1./hist_geant_tfixed->GetMaximum());

	TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(hist_fastdirc,"FastDIRC");
	leg->AddEntry(hist_geant,"HDGeant4");

	TLegend* leg2 = new TLegend(0.5,0.7,0.9,0.9);
	leg2->AddEntry(hist_fastdirc,"FastDIRC");
	leg2->AddEntry(hist_geant_tfixed,"HDGeant4 (t_fixed)");

	hist_fastdirc      -> SetLineColor(kRed);
	hist_geant         -> SetLineColor(kBlue);
	hist_geant_tfixed  -> SetLineColor(kBlue+1);

	hist_fastdirc->SetTitle(Form("%s",title_label.c_str()));

	hist_fastdirc->GetXaxis()->SetRangeUser(0.,300.);
	hist_fastdirc->GetXaxis()->SetTitle("Hit Time (ns)");
	hist_fastdirc->GetYaxis()->SetTitle("a.u.");

	TCanvas* canv = new TCanvas("canv","canv",2000,1600);
	canv->Divide(2,2);
	
	canv->cd(1);
	hist_fastdirc->Draw("HIST");
	hist_geant->Draw("HIST SAME");
	leg->Draw("SAME");

	canv->cd(2);
	gPad->SetLogy();
	hist_fastdirc->Draw("HIST");
	hist_geant->Draw("HIST SAME");
	leg->Draw("SAME");

        canv->cd(3);
        hist_fastdirc->Draw("HIST");
        hist_geant_tfixed->Draw("HIST SAME");
        leg2->Draw("SAME");

        canv->cd(4);
        gPad->SetLogy();
        hist_fastdirc->Draw("HIST");
        hist_geant_tfixed->Draw("HIST SAME");
        leg2->Draw("SAME");

	canv->Print(Form("%s/timing_%s.pdf",output_dir.c_str(),title_label.c_str()));
	delete canv;

}

void draw_compare_hit_patterns_profile(TProfile2D* hist_fastdirc, TProfile2D* hist_geant, string output_dir, string plot_label, string plot_title)
{

	hist_fastdirc -> SetTitle(Form("FastDIRC: %s",plot_title.c_str()));
	hist_geant -> SetTitle(Form("HDGeant4: %s",plot_title.c_str()));

	TCanvas* canv = new TCanvas("canv","canv",2000,800);
	canv->Divide(2,1);

	gStyle->SetPalette();

	canv->cd(1);
	gPad->SetRightMargin(0.15);
	hist_fastdirc -> Draw("colz1");
	//hist_fastdirc -> Draw();
	canv->cd(2);
	gPad->SetRightMargin(0.15);
	hist_geant -> Draw("colz1");
	//hist_geant -> Draw();

	canv->Print(Form("%s/%s.pdf",output_dir.c_str(),plot_label.c_str()));
	delete canv;

}


void draw_compare_3D(TFile* file_fastdirc, TFile* file_geant, string output_dir, string title_label)
{
	
	TH3F* hist_fastdirc      = (TH3F*) file_fastdirc->Get("pion_dist_3D");
	TH3F* hist_geant         = (TH3F*) file_geant->Get("DIRC_truth/hTruthPixelHit3D_South");
	
	const int nsteps = 100;
	double step_size = 1;
	double start_point = 10.;
	for (int step_i = 0 ; step_i < nsteps ; step_i++ )
	{
		double step_str = start_point + step_i * step_size;
		double step_stp = step_str + step_size;
		string plot_dir = output_dir + "/" + title_label;

		string plot_title = Form("%.1f--%0.1f ns",step_str,step_stp);	
		string plot_label = Form("%d",step_i);	
	
		hist_fastdirc -> GetZaxis()->SetRangeUser(step_str,step_stp);
		hist_geant    -> GetZaxis()->SetRangeUser(step_str,step_stp);

		TProfile2D* hist2d_fastdirc = (TProfile2D*) hist_fastdirc -> Project3DProfile("yx");
		TProfile2D* hist2d_geant    = (TProfile2D*) hist_geant    -> Project3DProfile("yx");
		
		cout<<(plot_dir + "/" + plot_label).c_str()<<endl;
		draw_compare_hit_patterns_profile(hist2d_fastdirc,hist2d_geant,plot_dir,plot_label,plot_title);
	}

}

void draw_compare_hit_patterns_all(TFile* file_data, TFile* file_geant, string output_dir, string title_label)
{
	
	TH2F* hist_data     = (TH2F*) file_data->Get("hit_dist_rowcol_tree");
	TH2F* hist_geant    = (TH2F*) file_geant->Get("hit_dist_rowcol_tree");

	TH2F* hist_data_fastdirc   = (TH2F*) file_data->Get("pion_dist_rowcol");
	TH2F* hist_geant_fastdirc  = (TH2F*) file_geant->Get("pion_dist_rowcol");

/*
	hist_fastdirc -> Scale(1./hist_fastdirc->GetEntries());
	hist_geant    -> Scale(1./hist_geant->GetEntries());

	TH2F* hist_ratio_fg    = (TH2F*) hist_geant->Clone();
	TH2F* hist_ratio_gf    = (TH2F*) hist_geant->Clone();
	for (int bin_x = 1; bin_x <= hist_ratio_fg->GetXaxis()->GetNbins(); bin_x++)
	for (int bin_y = 1; bin_y <= hist_ratio_fg->GetYaxis()->GetNbins(); bin_y++)
	{
		double content_fastdirc = hist_fastdirc->GetBinContent(bin_x,bin_y);
		double content_geant    = hist_geant->GetBinContent(bin_x,bin_y);
		double ratio_fg = 0.;
		double ratio_gf = 0.;
		if (content_fastdirc > 1e-5 && content_geant!=0.)
		{
			ratio_fg = content_fastdirc/content_geant;
		}
		if (content_geant > 1e-5 && content_fastdirc!=0.)
		{
			ratio_gf = content_geant/content_fastdirc;
		}
		hist_ratio_fg->SetBinContent(bin_x,bin_y,ratio_fg);
		hist_ratio_gf->SetBinContent(bin_x,bin_y,ratio_gf);
	}


	hist_ratio_fg->SetMaximum(2.);
	hist_ratio_fg->SetMinimum(0.);

	hist_ratio_gf->SetMaximum(2.);
	hist_ratio_gf->SetMinimum(0.);

	hist_ratio_fg->GetZaxis()->SetNdivisions(1010);
	hist_ratio_gf->GetZaxis()->SetNdivisions(1010);
*/
	hist_data           -> SetTitle(Form("data %s",title_label.c_str()));
	hist_geant          -> SetTitle(Form("bggen %s",title_label.c_str()));
	hist_data_fastdirc  -> SetTitle(Form("data (FastDIRC) %s",title_label.c_str()));
	hist_geant_fastdirc -> SetTitle(Form("bggen (FastDIRC) %s",title_label.c_str()));
	//hist_ratio_fg -> SetTitle(Form("FastDIRC/HDGeant4 %s",title_label.c_str()));
	//hist_ratio_gf -> SetTitle(Form("HDGeant4/FastDIRC %s",title_label.c_str()));

	TCanvas* canv = new TCanvas("canv","canv",2400,1600);
	canv->Divide(2,2);

	gStyle->SetPalette();

	canv->cd(1);
	gPad->SetRightMargin(0.15);
	hist_data -> Draw("colz");
	canv->cd(2);
	gPad->SetRightMargin(0.15);
	hist_data_fastdirc -> Draw("colz");

	//TExec *ex1 = new TExec("ex1","bwr_palette();");
	//ex1->Draw();
	canv->cd(3);
	gPad->SetRightMargin(0.15);
	hist_geant -> Draw("colz");
	canv->cd(4);
	gPad->SetRightMargin(0.15);
	hist_geant_fastdirc -> Draw("colz");

	canv->Print(Form("%s/%s.pdf",output_dir.c_str(),title_label.c_str()));
	delete canv;

}

void compare_plot()
{
	gStyle->SetOptStat(0);
/*
	for (int i = 0 ; i < 12 ; i++)
	{
		TFile* file_fastdirc = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/OutOfBox/bar_%d/hist_bar_%d.root",i,i));
		//TFile* file_fastdirc = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/bar_%d/hist_bar_%d.root",i,i));
		TFile* file_geant    = new TFile(Form("/home/yunjiey/Documents/PerTrack/outputs/bar_%d/hd_root_bar_%d.root",i,i));

		draw_compare_hit_patterns(file_fastdirc,file_geant,Form("/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/OutOfBox"),Form("bar_%d",i));
	}
*/

/*
	vector<vector<int>> bar_theta_phi =    {
						{1,4,32},
						{3,2,45},
						{3,3,165},
						{5,1,-75},
						{6,7,-110},
						{13,5,141}
					       };

	for (size_t i = 0 ; i < bar_theta_phi.size(); i++)
	{
		char* label = Form("bar_%d_%d_%d",bar_theta_phi[i][0],bar_theta_phi[i][1],bar_theta_phi[i][2]);
		
		TFile* file_fastdirc = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/%s/hist_%s.root",label,label));
		TFile* file_geant    = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/hdgeant4_outputs/%s/hd_root.root",label));

		draw_compare_hit_patterns(file_fastdirc,file_geant,Form("/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/angles"),Form("%s",label));
		draw_compare_timing(file_fastdirc,file_geant,Form("/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/timing"),Form("%s",label));

	}
*/


/*
	vector<vector<double>> bar_theta_phi_x_y =    {
						//{1,4,32,12.47,-13.8358},
						//{14,4,40,-250,-68.8493},
						//{6,2.5,68.4,120.8,-32.9},
						{3,1.2,113.5,180.2,-25.1}
					       };

	for (size_t i = 0 ; i < bar_theta_phi_x_y.size(); i++)
	{
		string label = "bar_" + to_string(int(round(bar_theta_phi_x_y[i][0])))
				      + "_" + to_string(int(round(bar_theta_phi_x_y[i][1]))) 
				      + "_" + to_string(int(round(bar_theta_phi_x_y[i][2])))
				      + "_" + to_string(int(round(bar_theta_phi_x_y[i][3])))
				      + "_" + to_string(int(round(bar_theta_phi_x_y[i][4])));
		
		TFile* file_fastdirc = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/%s/hist.root",label.c_str()));
		TFile* file_geant    = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/hdgeant4_outputs/%s/hd_root.root",label.c_str()));

		string hit_pattern_output_dir = "/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/full_kinematics";
		draw_compare_hit_patterns(file_fastdirc,file_geant,hit_pattern_output_dir,label);

		string hit_timing_output_dir  = "/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/timing";
		draw_compare_timing(file_fastdirc,file_geant,hit_timing_output_dir,label);

	}
*/

	TFile* file_data  = new TFile("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/data_bggen/hist_data.root"); 
	TFile* file_bggen = new TFile("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/data_bggen/hist_bggen.root"); 
	string label = "x_0_5_bar3";
	string outputdir = "/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/data_bggen";
	draw_compare_hit_patterns_all(file_data,file_bggen,outputdir,label);

}
