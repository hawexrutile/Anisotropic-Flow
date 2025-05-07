//Practise 1: ploting simple equations
{
    f1= new TF1("f1","x^3-sin(x)",-2,2); 
    f1->Draw();
    f1->SetLineColor(kBlue);
    f1->SetTitle("MyFunction;A;B");

    f2= new TF1("f2","sin(x)",-2,2); 
    f2->Draw("SAME");
    f2->SetLineColor(kRed);
    f2->SetTitle("MyFunction2;A2;B2");
}
