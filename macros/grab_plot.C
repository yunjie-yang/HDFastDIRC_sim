void grab_plot()
{
	gStyle->SetOptStat(0);
	TCanvas* canv = new TCanvas("canv","canv",1200,800);
	canv->Divide(4,3);
	for (int i = 0 ; i < 12; i++)
	{
		TFile* inputfile = new TFile(Form("root_outputs/bar_%d.root",i));
		TH2F* hist = (TH2F*) inputfile->Get("pion_dist_rowcol");
		canv->cd(i+1);
		hist->SetTitle(Form("bar %d",i));
		hist->Draw("colz");
	}
	canv->Print("bars.pdf");





}
