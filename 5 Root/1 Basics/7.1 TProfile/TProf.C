//TProfile;s Fill  takes in two values for a particular x the avg of all the coresponding y values is ploted
//TProfile is a 2D histogram unlike TH2D takes binning in only the x direstion

#include <TProfile.h>
#include <TRandom2.h>
void TProf(){
    TProfile *avgplot =new TProfile("avgplot", "Plots averge values", 100,0,10,"S"); //!"S" is used inside the constructor to make the errorbars the size of the standard-deviation. and 100 is the number of bins(Its kinda like a histogram+TGraph+Error)
    TRandom2 *rand = new TRandom2();
    int i=0;

    while (i<1000){
        avgplot->Fill(2,rand->Rndm());
        i++;

    }
    avgplot->Draw();


};