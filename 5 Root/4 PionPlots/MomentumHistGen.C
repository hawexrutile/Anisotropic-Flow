//*Creates Histogram  for px py pz in the same pad
{
    TH1D *h1=new TH1D("h1","h1",40,-2,2);
    TH1D *h2=new TH1D("h2","h2", 40, -2, 2);
    TH1D *h3=new TH1D("h3","h3", 40, -2, 2);
    auto legend = new TLegend(0.1,0.7,0.48,0.9);

    int i;
    while ( i<listid.size()){
        h1->Fill(listpx[i]);
        h2->Fill(listpy[i]);
        h3->Fill(listpz[i]);
        i++;

    };
    h1->Draw();
    h1->SetLineColor(kBlue);
    h1->SetTitle("p Distribution; Px; Counts");
    
    h2->Draw("SAME");
    h2->SetLineColor(kRed);
    h2->SetTitle("p Distribution; Py; Counts");

    h3->Draw("SAME");
    h3->SetLineColor(kGreen);
    h3->SetTitle("p Distribution; Pz; Counts");

    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry("h1","px","l");
    legend->AddEntry("h2","py","l");
    legend->AddEntry("h3","pz","l");
    legend->Draw();
}