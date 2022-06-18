{

    c1= new TCanvas("c1","Triple Plot",1000,500);
    c1->Divide(3,1);

    TH1D *h1=new TH1D("h1","h1",40,-2,2);
    TH1D *h2=new TH1D("h2","h2", 40, -10, 6);
    TH1D *h3=new TH1D("h3","h3", 40, -2, 2);
    auto legend = new TLegend(0.1,0.7,0.48,0.9);

    ifstream pt;
    ifstream eta;
    ifstream phi;

    double valuept;
    double valueeta;
    double valuephi;
    pt.open("ofpt.txt");
    eta.open("ofeta.txt");
    phi.open("ofphi.txt");
    while(!pt.eof()){
        pt >> valuept;
        eta >> valueeta;
        phi >> valuephi;
        h1->Fill(valuept);
        h2->Fill(valueeta);
        h3->Fill(valuephi);

    };
    c1->cd(1); 
    h1->Draw();
    h1->SetLineColor(kBlue);
    h1->SetTitle("pt Distribution; pt; Counts");
    
    c1->cd(2); 
    h2->Draw();
    h2->SetLineColor(kRed);
    h2->SetTitle("eta Distribution; eta; Counts");

    c1->cd(3);
    h3->Draw();
    h3->SetLineColor(kGreen);
    h3->SetTitle("phi Distribution; phi; Counts");
}
