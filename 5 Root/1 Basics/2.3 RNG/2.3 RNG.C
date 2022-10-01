//Random No generator
{
    std::cout << gRandom->Uniform(10)<<"\n";   //Random to the given nO
    std::cout << gRandom->Uniform(1,2)<<"\n";  //Random between 1,2
    std::cout << gRandom->Gaus()<<"\n";        //Default is "mean 0" and "width 1"

    g1 = new TGraph();
    for(int i=0;i<100000; i++){       //int is imp
        g1-> SetPoint(i,gRandom->Gaus(5,1),gRandom->Gaus(0,2));   //Hack for adding points
    };
    g1-> Draw("AP");

}