void compare_plot()
{
	gStyle->SetOptStat(0);
	for (int i = 0 ; i < 12 ; i++)
	{
		TCanvas* canv = new TCanvas("canv","canv",2000,800);
		TFile* file_fastdirc = new TFile(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/outputs/bar_%d/hist_bar_%d.root",i,i));
		TFile* file_geant    = new TFile(Form("/home/yunjiey/Documents/PerTrack/outputs/bar_%d/hd_root_bar_%d.root",i,i));

		//TH2F* hist_fastdirc = (TH2F*) file_fastdirc->Get("pion_dist_rowcol_wide");
		TH2F* hist_fastdirc = (TH2F*) file_fastdirc->Get("pion_dist_rowcol");
		TH2F* hist_geant    = (TH2F*) file_geant->Get("DIRC_truth/hTruthPixelHit_South");

		canv->Divide(2,1);
		canv->cd(1);
		hist_fastdirc -> SetTitle(Form("FastDIRC bar %d",i));
		hist_fastdirc -> Draw("colz");
		canv->cd(2);
		hist_geant -> SetTitle(Form("HDGeant4 bar %d",i));
		hist_geant -> Draw("colz");
		canv->Print(Form("/media/sf_SharedFolderVM/FastDIRC_geometry/plots/compare/bar_%d.pdf",i));
		delete canv;
	}

}
