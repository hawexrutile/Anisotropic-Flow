//Double Plots

{
    c1= new TCanvas("c1","Double Plot",1000,50);        //
    c1->Divide(2,1);

    c1->cd(1); 
    f1= new TF1("f1","x^3-sin(x)",-2,2); 
    f1->Draw();
    f1->SetLineColor(kBlue);
    f1->SetTitle("MyFunction;A;B");

    c1->cd(2);
    f2= new TF1("f2","sin(x)",-2,2); 
    f2->Draw();                         //No "SAME" arg 
    f2->SetLineColor(kRed);
    f2->SetTitle("MyFunction2;A2;B2");
}
