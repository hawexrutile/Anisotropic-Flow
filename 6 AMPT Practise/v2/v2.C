{
    double v2=0;
    double IntegNum=0;
    double IntegDen=0;
    int BinCount=45;                //!Always has to be even, as  we are gonna divide by 2
    int HalfBinCount=BinCount/2;
    double GraphLowerLim=-1.6;
    double GraphUpperLim=1.6;
    vector<double> CosArray;
    int j=0;
    int lpos1=0;             //array for contayning centres of each bin (Phi)-----x axis
    int lpos2=0;
    TH1D *PhiHist=new TH1D("PhiHist","N vs Phi for pion",BinCount,GraphLowerLim,GraphUpperLim);
    TH1D *CosHist=new TH1D("CosHist","Cos(Phi-Psi) vs Phi for pion",BinCount,GraphLowerLim,GraphUpperLim);
    while(j<(listtrail.size()-1)){                     
        double SinSum=0; 
        double CosSum=0;



        while ((lpos1>=listtrail[j]) && (lpos1<listtrail[j+1])){   //This while loop is to calculate Psi2
            SinSum=SinSum +sin(2*listphi[lpos1]);        //summation of cos and sin phi
            CosSum=CosSum +cos(2*listphi[lpos1]);      
            lpos1++;
        };
        double psi2=(1/2)*atan(SinSum/CosSum);         //for each event we calcaulate Psi2


        while ((lpos2>=listtrail[j]) && (lpos2<listtrail[j+1])){  //This while loop is for filling the phi and cos(Phi-Psi)histogram for all the events
            CosArray.push_back(cos(listphi[lpos2]-psi2));               //The calculated Psi2 is subtracted from the phis of each particle for the given event and filled into a histogrma
            PhiHist->Fill(listphi[lpos2]);                        //
            lpos2++;
        };
        j++;
    };
    double ListdNdPhi[45];              //aray for containing dN/dPhi--------y axis
    double ListPhiCentre[45];          //!Cant put BinCount Variable here cuz array size has to be pre defined
    for (int j=0 ;j<HalfBinCount;j++){            //*Taking counts from each bin, subtracting adjacent bin heights and dividing by bin width
        double l= PhiHist->GetBinContent(HalfBinCount+j);//l=lower;u=upper;w=width
        double u= PhiHist->GetBinContent(HalfBinCount+j+1);//We are using half bin count cuz the range of integration is from 0-2Pi and not -Pi to Pi.ie fro the middle of the plot.We do this as the function in the integrand is odd.
        double w= PhiHist->GetBinWidth(j);

        
        double ContentCos= CosHist->GetBinContent(HalfBinCount+j);//l=lower;u=upper;w=width
        IntegNum=IntegNum+(ContentCos*(((u-l)*1.0)/(w*1.0))*w);
        IntegDen=IntegDen+(((u-l)*1.0)/(w*1.0))*w;


        ListdNdPhi[j]=((u-l)*1.0)/(w*1.0);
        ListPhiCentre[j]=(PhiHist->GetBinCenter(j));

    };
    v2=IntegNum/IntegDen;
    TGraph *PhiDistHist= new TGraph(BinCount,ListPhiCentre,ListdNdPhi);
    PhiDistHist-> SetMarkerStyle(22);
    PhiDistHist->Draw("AP");
    PhiDistHist->SetTitle("dN/dPhi vs Phi for pion");
    // PhiHist->Draw();
    cout<<v2;

};