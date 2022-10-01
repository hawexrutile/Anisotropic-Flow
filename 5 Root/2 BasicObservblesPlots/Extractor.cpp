//*Extractor can be used to extract individual colums from a ampt dat file to respective txt files.

#include <iostream>         //General Stream
#include <fstream>          //Opening and Writing on files
#include <sstream>          //for selecting 
#include <string>
#include <cmath>
#include <limits>

using namespace std;

int main() {
    double sum = 0;
    string x;
    ifstream inFile;        //initialises write objects
    ofstream ofpx;
    ofstream ofpy;
    ofstream ofpz;
    ofstream ofpt;
    ofstream ofeta;
    ofstream ofphi;
    
    
    inFile.open("ampt.dat");   //opens a file for writing dat
    ofpx.open("ofpx.txt");   
    ofpy.open("ofpy.txt");   
    ofpz.open("ofpz.txt");   
    ofpt.open("ofpt.txt");   
    ofeta.open("ofeta.txt");   
    ofphi.open("ofphi.txt");   
                             //fstream is used here to save the dat file into its object called InFile
    if (!inFile) {           
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
    while (getline(inFile, x)) {    
        istringstream ss(x);               //storing string in a "istringstream" object called ss
        string id;string px;string py;string pz;string m;string x;string y;string z;string pt;string eta;string phi;string t;
        int iid;double dpx;double dpy;double dpz;double dm;double dx;double dy;double dz;double dpt;double deta;double dphi;double dt;
        ss >> id;                         //String till first white space is witten in the first row of the respective txt file
        ss >> px;                         //String from first whire space till second white space
        ss >> py;                         //likewise
        ss >> pz;                         //Kochirakoso
        ss >> m;
        ss >> x;
        ss >> y;
        ss >> z;
        ss >> t;

        dpx= stod(px);           //Converting string to double
        dpy= stod(py);           //phi here is rapidity
        dpz=  stod(pz);   
        dpt=  sqrt(stod(px) * stod(px) + stod(py) * stod(py));   
        dphi=  0.5 * log(abs((stod(m)+stod(pz))/(stod(m)-stod(pz))));   

        if (abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))!=0 &&
            to_string(log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))))!="inf"){
            deta= log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz)))));
            ofeta << deta << "\n";   
        }
/*        if (abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))==0){
            eta=  -1 * numeric_limits<double>::infinity();
        }
        else {
            eta=log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))); 
        }
*/ //?i TRIED TO TAKE CARE OF INF" BUT IT DINT WORK PLAS LOOK INTO IT
        ofpx << dpx << "\n";           //Converting string to double
        ofpy << dpy << "\n";           //phi here is rapidity
        ofpz << dpz << "\n";   
        ofpt << dpt << "\n";   
        ofphi << dphi << "\n";   
        
    }
    
    inFile.close();
    ofpx.close();
    ofpy.close();
    ofpz.close();
    ofpt.close();
    ofeta.close();
    ofphi.close();

}



