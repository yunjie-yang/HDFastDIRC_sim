void draw_compare_DLL(TH1* hist_dll_pion, TH1* hist_dll_kaon, int rebin, double range, string plot_title) 
{

        hist_dll_pion -> SetLineColor(kBlue);
        hist_dll_kaon -> SetLineColor(kRed);

        hist_dll_pion -> SetMarkerColor(kBlue);
        hist_dll_kaon -> SetMarkerColor(kRed);

        hist_dll_pion -> SetMarkerColor(kBlue);
        hist_dll_kaon -> SetMarkerColor(kRed);

        hist_dll_pion -> SetFillColorAlpha(kBlue,0.5);
        hist_dll_kaon -> SetFillColorAlpha(kRed,0.5);

        hist_dll_pion -> Rebin(rebin);
        hist_dll_kaon -> Rebin(rebin);

        hist_dll_pion -> Scale(1./hist_dll_pion->GetMaximum());
        hist_dll_kaon -> Scale(1./hist_dll_kaon->GetMaximum());

        hist_dll_pion -> SetTitle(plot_title.c_str());
        hist_dll_pion -> GetXaxis() -> SetTitle("DLL");
        hist_dll_pion -> GetYaxis() -> SetTitle("a.u.");
        hist_dll_pion -> GetXaxis() -> SetRangeUser(-range,range);


	TF1 *ff_pion, *ff_kaon;
	double sep=0,m1=0,m2=0,s1=0,s2=0;

	//if(hist_dll_pion->GetEntries()>100){
	if(1)
	{
		hist_dll_pion->Fit("gaus","S");
    		ff_pion = hist_dll_pion->GetFunction("gaus");
    		ff_pion->SetLineColor(1);
		m1=ff_pion->GetParameter(1);
		s1=ff_pion->GetParameter(2);
	}
	//if(hist_dll_kaon->GetEntries()>100){
	if(1)
	{
		hist_dll_kaon->Fit("gaus","S");
		ff_kaon = hist_dll_kaon->GetFunction("gaus");
		ff_kaon->SetLineColor(1);
		m2=ff_kaon->GetParameter(1);
		s2=ff_kaon->GetParameter(2);
	}

	if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
	printf("m1  = %6.2f, s1 = %6.2f\n", m1, s1);
	printf("m2  = %6.2f, s2 = %6.2f\n", m2, s2);
	printf("sep = %6.2f\n", sep);


        TLegend* leg = new TLegend(0.725,0.675,0.9,0.9);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

        leg->AddEntry(hist_dll_pion,"#pi");
        leg->AddEntry(hist_dll_kaon,"K");
	leg->AddEntry((TObject*)0, Form("%2.2f #sigma",sep), "");

        hist_dll_pion->Draw("HIST");
        hist_dll_kaon->Draw("HIST SAME");
	ff_pion->Draw("SAME");
	ff_kaon->Draw("SAME");
        leg->Draw();
        gPad->RedrawAxis();
}
