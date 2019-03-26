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

void compare_plot()
{
	gStyle->SetOptStat(0);
	for (int i = 0 ; i < 12 ; i++)
	{
		TFile* file_fastdirc = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/bar_%d/hist_bar_%d.root",i,i));
		TFile* file_geant    = new TFile(Form("/home/yunjiey/Documents/PerTrack/outputs/bar_%d/hd_root_bar_%d.root",i,i));

		//TH2F* hist_fastdirc = (TH2F*) file_fastdirc->Get("pion_dist_rowcol_wide");
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

		hist_fastdirc -> SetTitle(Form("FastDIRC bar %d",i));
		hist_geant -> SetTitle(Form("HDGeant4 bar %d",i));
		hist_ratio_fg -> SetTitle(Form("FastDIRC/HDGeant4 bar %d",i));
		hist_ratio_gf -> SetTitle(Form("HDGeant4/FastDIRC bar %d",i));

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

		canv->Print(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/bar_%d.pdf",i));
		delete canv;
	}

}
