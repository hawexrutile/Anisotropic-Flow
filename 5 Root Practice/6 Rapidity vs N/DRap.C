{  //Code for plotting d(rapidity)/d(N) vs rapidity;  rapidity is denoted by y
    
    TH1D *h1= new TH1D("h1","rapidity distribution",40,-5,5);
    double listdydn[40];              //aray for containing dy/dn--------y axis
    double listycent[40];              //array for contayning centres of each bin (rapidity)-----x axis
    for (double i :listphi){               //*Filling the Histogram
        h1->Fill(i);                    //list phi is vector containing rapidity entries
    }
    for (int j=0 ;j<40;j++){            //*Taking counts from each bin, subtracting adjacent bin heights and dividing by bin width
        
        double l= h1->GetBinContent(j);
        double u= h1->GetBinContent(j+1);
//        listdydn[j]=(((5-(-5)/40))/(u-l));
        if (u!=l){
         listdydn[j]=(((5.-(-5.)/40.))/(u-l));
         listycent[j]=(h1->GetBinCenter(j));
        cout<<listdydn[j]<<"      "<<listycent[j]<<"\t"<< j<<"\n";
        }
        
    }
    TGraph *g2= new TGraph(40,listdydn,listycent);
    g2-> SetMarkerStyle(21);
    g2->Draw("AP");
};