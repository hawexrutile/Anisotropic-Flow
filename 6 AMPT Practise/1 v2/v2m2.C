//This is a method to plot v2 vs pt

{
    double Pmin=*min_element(listpt.begin(),listpt.end());     //To create bins where v2 values of respective bins are stored
    double Pmax=*max_element(listpt.begin(),listpt.end());
    int bin=100;                                                //set bin value here; the number of points you want to plot
    double bw=(Pmax-Pmin)/(1.0*bin);
    int i=0;
    double Bmin=0;
    double Bmax=0;
    double binarray[100];                    //! total no of bins has to be set here(dont use the bin variable cuz its an array)
    double v2array[100];                     //! same here
    while (i<bin){               //for each pt-bin we run the loop on all the trails
        double v2=0;
        double sum=0;
        int j=0;
        int lpos1=0;             //looping int fot calculationg Psi
        int lpos2=0;             //looping int for calculating tot
        int Tot=0;               //Variable to calculate v2
        Bmin=Pmin+i*(bw);
        Bmax=Pmin+(i+1)*(bw);
        double bc=Pmin+((2*i+1)/2)*bw;    //bin centre
        while(j<(listtrail.size()-1)){   //for each event
            double SinSum=0; 
            double CosSum=0;
            while ((lpos1>=listtrail[j]) && (lpos1<listtrail[j+1])){   //This while loop is to calculate Psi2. Here listtrail is a list of no events in each trail
                    SinSum=SinSum +sin(2*listphi[lpos1]);        //summation of cos and sin phi
                    CosSum=CosSum +cos(2*listphi[lpos1]);        //?should Psi2 be for all the trails in the event or for the particular particle or for that particular ptbin
                    lpos1++;

            };
            double psi2=(1/2)*atan(SinSum/CosSum);         //for each event we calcaulate Psi2
            while ((lpos2>=listtrail[j]) && (lpos2<listtrail[j+1])){  //for that particular event
                if ((listpt[lpos2]>=Bmin)&&(listpt[lpos2]<Bmax)){     //if the ith particle's pt belongs to the bin
                    if (listphi[lpos2]>=0){                           //if phi is greater than 0 ; cuz we are averaging from 0-2pi
                        sum=sum+cos(listphi[lpos2]-psi2);               
                        Tot++;
                    };
                };
                lpos2++;
            };
            j++;
        };

        v2=sum/Tot;
        binarray[i]=bc;           //These arrays are used for plotting
        v2array[i]=v2;
        i++;
    };
    // for (double i:v2array ){
    //     cout<<i;
    // };
    TGraph *g2= new TGraph(bin,binarray,v2array);
    g2-> SetMarkerStyle(22);
    g2->Draw("AP");
    g2->SetTitle("v2 vs pt ;Pt;v2")

};