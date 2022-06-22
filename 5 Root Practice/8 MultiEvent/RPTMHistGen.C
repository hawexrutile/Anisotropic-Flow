{
    c1= new TCanvas("c1","Triple Plot",1000,500);
    c1->Divide(3,1);

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    
    int array[]={111,211};
    int m=0;
    for(int j : array){
        if (m==0){
            lister("ampt.dat",j);
            c1->cd(1);
            TH1D *h1=new TH1D("h1","h1",40,-2,2);
            for (double i:listpt){
                h1->Fill(i);
            };
            h1->Draw();
            h1->SetTitle("P_t Distribution;P_t;Counts");
            h1->SetLineColor(kRed);
            
            c1->cd(2);
            TH1D *h2=new TH1D("h2","h2", 40, -2, 2);
            for (double i:listeta){
                h2->Fill(i);
            };
            h2->Draw();
            h2->SetTitle("P_eta Distribution;P_eta;Counts");
            h2->SetLineColor(kRed);

            c1->cd(3);
            TH1D *h3=new TH1D("h3","h3", 40, -2, 2);
            for (double i:listphi){
                h3->Fill(i);
            };
            h3->Draw();
            h3->SetTitle("P_phi Distribution;P_phi;Counts");
            h3->SetLineColor(kRed);
            m++;
        }
        else{
            lister("ampt.dat",j);
            c1->cd(1);
            TH1D *h4=new TH1D("h4","h4",40,-2,2);
            for (double i:listpt){
                h4->Fill(i);
            };
            h4->Draw("SAME");
            h4->SetLineColor(kBlue);
            
            c1->cd(2);
            TH1D *h5=new TH1D("h5","h5", 40, -2, 2);
            for (double i:listeta){
                h5->Fill(i);
            };
            h5->Draw("SAME");
            h5->SetLineColor(kBlue);
            
            c1->cd(3);
            TH1D *h6=new TH1D("h6","h6", 40, -2, 2);
            for (double i:listphi){
                h6->Fill(i);
            };
            h6->Draw("SAME");
            h6->SetLineColor(kBlue);
            m++;

        };
    };
}
