{
    double v2=0;
    double sum=0;
    int j=0;
    int lpos1=0;             //array for contayning centres of each bin (Phi)-----x axis
    int lpos2=0;
    int Tot=0;
    while(j<(listtrail.size()-1)){                     
        double SinSum=0; 
        double CosSum=0;
        while ((lpos1>=listtrail[j]) && (lpos1<listtrail[j+1])){   //This while loop is to calculate Psi2
            SinSum=SinSum +sin(2*listphi[lpos1]);        //summation of cos and sin phi
            CosSum=CosSum +cos(2*listphi[lpos1]);      
            lpos1++;
        };
        double psi2=(1/2)*atan(SinSum/CosSum);         //for each event we calcaulate Psi2
        while ((lpos2>=listtrail[j]) && (lpos2<listtrail[j+1])){  //This while loop is for filling the phi and cos(Phi-Psi)histogram for all the events
            if (listphi[lpos2]>=0){
            sum=sum+cos(listphi[lpos2]-psi2);               //The calculated Psi2 is subtracted from the phis of each particle for the given event and filled into a histogrma  
            Tot++;
            };
            lpos2++;
        };
        j++;
    };
    v2=sum/Tot;
    cout<<v2;

};