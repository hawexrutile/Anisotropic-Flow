{
    int i;
    TH1D *proton= new TH1D("proton","protonhist",10,100,120);
    while (i<100){
        for (auto x : lol){
            if (x.first==2212){
                proton->Fill(x.second);
            }
        }
        i++;
    }
    proton->Draw();
}
