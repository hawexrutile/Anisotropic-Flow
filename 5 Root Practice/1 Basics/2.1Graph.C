//Calling Data and plotting it

{
    g1 = new TGraph("TrialData.txt");
    g1-> SetMarkerStyle(21);
    g1->Draw("APL");                      //A calculates the Axis; P for only daa points; L for straight line; C for curved etc
}