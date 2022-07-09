
void ErrorBar(){
    TCanvas *c1=new TCanvas();
    TGraphErrors *gr= new TGraphErrors();
    ifstream file;
    file.open("TrialData.txt",ios::in);

    double x, y, z ,ex ,ey;
    int n=0;
    while(!file.eof()){
        file >> x>>y>>ex>>ey;

        n=gr->GetN();
        gr->SetPoint(n,x,y);
        gr->SetPointError(n,ex,ey);
    };
    

    gr->Draw("A*");


};