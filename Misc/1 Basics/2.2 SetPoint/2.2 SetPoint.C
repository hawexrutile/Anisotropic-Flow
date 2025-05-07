//Add points on the go
{
    g1 = new TGraph();
    for(int i=0;i<15; i++){       //int is imp
        g1-> SetPoint(i,i,i*i);   //Hack for adding points
    };
    g1-> Draw("AP");

    g2 = new TF1("g2","x^2",0,15);//for comparison
    g2-> Draw("SAME");

}