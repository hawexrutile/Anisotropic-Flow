{
    c1= new TCanvas("c1","Triple Plot",1000,1000);
    c1->Divide(2,2);

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    
    int array[]={111,211,22};
    int m=0;
    for(int j : array){
        if (m==0||m==1){
            lister("ampt.dat",j);
            c1->cd(1);
            TH1D *h1=new TH1D("h1","h1",200,0,2);
            for (double i:listpt){
                h1->Fill(i);
            };
            h1->Draw();
            h1->SetTitle("P_t Distribution;P_t;Counts");
            h1->SetLineColor(kRed);
            
            c1->cd(2);
            TH1D *h2=new TH1D("h2","h2", 200, -10, 10);
            for (double i:listeta){
                h2->Fill(i);
            };
            h2->Draw();
            h2->SetTitle("Pseudo Rapidity Distribution;Pseudo Rapidity Eta;Counts");
            h2->SetLineColor(kRed);

            c1->cd(3);
            TH1D *h3=new TH1D("h3","h3", 200, -3.5, 3.5);
            for (double i:listphi){
                h3->Fill(i);
            };
             h3->Draw();
            h3->SetTitle("Azimuthal Distribution;Rapidity Phi;Counts");
            h3->SetLineColor(kRed);

            c1->cd(4);
            TH1D *h4=new TH1D("h4","h4", 200, -10, 10);
            for (double i:listyi){
                h4->Fill(i);
            };
            h4->Draw();
            h4->SetTitle("Rapidity Distribution;Rapidity y;Counts");
            h4->SetLineColor(kRed);
            m++;
        }
        else{
            lister("ampt.dat",j);
            c1->cd(1);
            TH1D *h5=new TH1D("h5","h5",200,0,2);
            for (double i:listpt){
                h5->Fill(i);
            };
            h5->Draw("SAME");
            h5->SetLineColor(kBlue);
            
            c1->cd(2);
            TH1D *h6=new TH1D("h6","h6", 200, -10, 10);
            for (double i:listeta){
                h6->Fill(i);
            };
            h6->Draw("SAME");
            h6->SetLineColor(kBlue);
            
            c1->cd(3);
            TH1D *h7=new TH1D("h7","h7", 200, -3.5, 3.5);
            for (double i:listphi){
                 h7->Fill(i);
            };
            h7->Draw("SAME");
            h7->SetLineColor(kBlue);
            m++;

            c1->cd(4);
            TH1D *h8=new TH1D("h8","h8", 200, -10, 10);
            for (double i:listyi){
                h8->Fill(i);
            };
            h8->Draw("SAME");
            h8->SetLineColor(kBlue);
            m++;

        };
    };
}
