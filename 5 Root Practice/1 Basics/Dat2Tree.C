void Dat2Tree(){
    double x, y;               //TO STORE EACH OF THE VARIABLE IN EACH LINE
    ifstream dat;              
    dat.open("TrialData.txt",ios::in);
    TFile *amptroot=new TFile("amptroot.root","create");       //CREATES A ROOT FILE
    TTree *tree=new TTree("tree","tree");             //Creates a Tree
    tree->Branch("x",&x,"x/D");           //Creates branches for the tree
    tree->Branch("y",&y,"y/D");           //"name",&pointer,"name/DataType"

    while(!dat.eof()){                //write data in each branch of the tree
        dat>>x>>y;
        tree->Fill();            //Fills the current x , y values in the respective branch of the tree
    };
    amptroot->Write();               //Stores the tree in the amptroot.root file
    amptroot->Close();                //close the root and dat file
    dat.close();




};