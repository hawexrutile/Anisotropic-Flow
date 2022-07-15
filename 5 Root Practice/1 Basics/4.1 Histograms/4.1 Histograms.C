{
    h = new TH1D("h","First Histogram",40,0,10);  //TH(histogram)1(Dimension)D(Doubles) Name Title NoofBins Begin and end
                                                  //No of bins dhould be square root of the number of entries  
    for (int i=0; i<100000;i++){
    h->Fill(gRandom->Gaus(5,1));
    }
    h->SetMinimum(0);
    h->Draw();
/*  h->Fill(0);
    h->Fill(2);
    h->Fill(2.5);
    h->Fill(3.5);
    h->Fill(2.5);
    h->Fill(0);
    h->Fill(4.506);
    h->Fill(5);
    h->Fill(1.11);*/
    std::cout<< "Counts in each bin\n";
    for (int i=1; i<=h->GetNbinsX(); i++){
        cout<<h->GetBinCenter(i) <<"\t"<< h->GetBinContent(i)<<"\n";  //tab is compulsory
    }
}