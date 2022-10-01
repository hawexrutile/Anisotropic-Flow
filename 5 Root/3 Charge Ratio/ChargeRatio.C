//*Particle charge identifier
{
    pdgDb = TDatabasePDG::Instance();
    Double_t q;
    Double_t c=0;
    Double_t n=0;
    for (int x : listid){
        TParticlePDG* pdgP  = pdgDb->GetParticle(x);
        Double_t q = pdgP->Charge() / 3;
        if (q==0){
            n++;
        }
        else {
            c++;
        }
    }

    cout << c/(c+n)<<" \t charged particle ratio\n";
    cout << n/(c+n)<<" \t neutral particle ratio\n";

    
}
//fix both 