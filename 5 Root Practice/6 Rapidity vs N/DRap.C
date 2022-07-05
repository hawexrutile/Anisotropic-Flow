{  //Code for plotting d(N)/d(N\y) vs y;  rapidity is denoted by y
    lister("ampt.dat",111);
    int BinCount=50;
    int GraphLowerLim=-10;
    int GraphUpperLim=10;
    
    TH1D *h1= new TH1D("h1","rapidity distribution",BinCount,GraphLowerLim,GraphUpperLim);  //!be carefull here; we need binwidth here such that its monotonously increasing or decreasing(not at gaussian midpoint). so not sqrt(total no of trails)
    double listdndy[50];              //aray for containing dy/dn--------y axis
    double listycent[50];              //array for contayning centres of each bin (rapidity)-----x axis
    for (double i :listyi){               //*Filling the Histogram
        h1->Fill(i);                    //list phi is vector containing rapidity entries
    };
    for (int j=0 ;j<BinCount;j++){            //*Taking counts from each bin, subtracting adjacent bin heights and dividing by bin width
        
        double l= h1->GetBinContent(j);
        double u= h1->GetBinContent(j+1);
//        listdndy[j]=(((5-(-5)/40))/(u-l));
        // if (u!=l){
         listdndy[j]=(u-l)/((((GraphUpperLim*1.0)-(GraphLowerLim*1.0))/(BinCount*1.0)));
         listycent[j]=(h1->GetBinCenter(j));
        //  cout<<listdndy[j]<<"      "<<listycent[j]<<"\t"<< j<<"\n";
        // };
        
    };
    TGraph *g2= new TGraph(BinCount,listycent,listdndy);
    g2-> SetMarkerStyle(22);
    g2->Draw("AP");
    g2->SetTitle("dN/dy vs y;Rapidty(y);dN/dy")
};