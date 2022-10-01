 {
        g1 = new TGraph();
    double x,y;
    for(int i=0;i<20; i++){       //int is imp
        x=0.5*i;
        y=2*x+5+gRandom->Gaus(0,2);//Not Gauss
        g1-> SetPoint(i,x,y);   //Hack for adding points
    };
    g1->SetMarkerStyle(22);
    g1-> Draw("AP");

    fit= new TF1("fit","[0]*x+[1]",0,5);  //[0] and [1] ard the fir params respectively; Chi(greak) means square of distance
//  fit= new TF1("fit","pol4",0,5);      // pol4 means  polinomila of degree 4
//  fit= new TF1("fit","gaus",0,5);    

    fit->SetParaNames("Slope","y-intercept");
    g1->Fit(fit);
 }