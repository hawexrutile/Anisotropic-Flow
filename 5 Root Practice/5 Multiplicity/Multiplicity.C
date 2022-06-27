{   c1=new TCanvas("c1","Multiplicity",1000,500);
    c1->Divide(2,1);
    int j=0;                        //j is used to add  to pinter to read 2 consecutive points from the vector
    int l=0;
    lister("ampt.dat",99);
    TH1D *pion= new TH1D("pion","Pion Multiplicity",10,3200,4000);
    TH1D *photon= new TH1D("photon","Photon Multiplicity",8,300,500);

    while(int (j<listtrail.size()-1)){
        vector<int> temp={};
        while ((l>=listtrail[j]) && (l<listtrail[j+1])){
            temp.push_back(listid[l]);
            l++;
        };
        map<int,int> lol=pcounter(temp);    
        //  cout<<lol.at(111)+lol.at(211)<<"\n";
        //  cout<<lol.at(22)<<"\n";
        pion->Fill(lol.at(111)+lol.at(211));
        photon->Fill(lol.at(22));
        j++;
    };
    c1->cd(1);
    pion->Draw();
    pion->SetLineColor(kBlue);
    pion->SetTitle("Pion Multiplicity;Counts(N);Multiplicity");
    c1->cd(2);
    photon->Draw();
    photon->SetLineColor(kRed);
    pion->SetTitle("Photon Multiplicity;Counts(N);Multiplicity");
}
