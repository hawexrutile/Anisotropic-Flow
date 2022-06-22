#include <iostream>         //General Stream
#include <fstream>          //Opening and Writing on files
#include <sstream>          //for selecting 
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>         //for sorting
#include <process.h>
#include <limits>            //for using inf;didt workout though

using namespace std;
    vector<int> listid;                                     //*outside cuz need to call em later for ither functions
    vector<double> listpx;
    vector<double> listpy;
    vector<double> listpz;
    vector<double> listpt;
    vector<double> listeta;
    vector<double> listphi;

void extractor(string datfilename){                         //*ANCHOR Extractor can be used to extract individual colums from a ampt dat file to respective txt files.
    double sum = 0;
    string x;
    ifstream inFile;        //initialises write objects
    ofstream ofpx;
    ofstream ofpy;
    ofstream ofpz;
    ofstream ofpt;
    ofstream ofeta;
    ofstream ofphi;
    
    
    inFile.open(datfilename);   //opens a file for writing dat
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

void lister(string datfilename,int pid){                    //*ANCHOR Lister is used to get vector output of all the observables and ID 
    double sum = 0;                                         //TODO give pid default of 99 if no input
    string x;                                               //! pid is 99 for all particles and respective pid for others
    ifstream inFile;        //initialises write objects

    
    inFile.open(datfilename);
    if (!inFile) {           
        cout << "Unable to open file";
        exit(1);                           //! terminate with error
    }
    int n=0;
    int i=0;
    while (getline(inFile, x)) { 
        string eventno; string testno; string trailno; string b; string Npart1; string Npart2; string Npart1_el; string Npart1_inel; string Npart2_el; string Npart2_inel;

        istringstream ss(x); //storing string in a "istringstream" object called ss              
        string id; string m; string px;string py;string pz; string x;string y;string z;string pt; string eta;string phi;string t;
        
        
        if (n==i){
            ss>> eventno;
            ss>> testno;
            ss>>trailno;
            i=i+stoi(trailno);

        }
        else {
            ss >> id;                         //String till first white space is witten in the first row of the respective txt file
            ss >> px;                         //String from first white space till second white space
            ss >> py;                         //likewise
            ss >> pz;                         //Kochirakoso
            ss >> m;
            ss >> x;
            ss >> y;
            ss >> z;
            ss >> t;
            if (pid==99){
                if (abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))!=0 && 
                to_string(log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))))!="inf"){
                listeta.push_back(log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))) );   
                }   
                listid.push_back(stoi(id));
                listpx.push_back(stod(px));           //Converting string to double
                listpy.push_back(stod(py) );           //phi here is rapidity
                listpz.push_back(stod(pz) );   
                listpt.push_back(sqrt(stod(px) * stod(px) + stod(py) * stod(py)) );   
                listphi.push_back(0.5 * log(abs((stod(m)+stod(pz))/(stod(m)-stod(pz)))) );    
            n++;
            }
            else if  (stoi(id)==pid){
                listid.push_back(stoi(id));
                listpx.push_back(stod(px));           //Converting string to double
                listpy.push_back(stod(py) );           //phi here is rapidity
                listpz.push_back(stod(pz) );   
                listpt.push_back(sqrt(stod(px) * stod(px) + stod(py) * stod(py)) );   
                listphi.push_back(0.5 * log(abs((stod(m)+stod(pz))/(stod(m)-stod(pz)))) );   
                if (abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))!=0 && 
                    to_string(log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))))!="inf"){
                    listeta.push_back(log(abs((stod(px)+stod(py))/(stod(px)+stod(py)+(2*stod(pz))))) );   
                };
            n++;   
            };    
        };
        
    };
}

/*
template <typename S>                  //!without this cout wont be able to list vectors
ostream& operator<<(ostream& os,       //TODO (I deleteed cout for reading vectors) by mistake, remake it
                    const vector<S>& vector)
{
    // Printing all the elements
    // using <<
    for (auto element : vector) {
        os << element << " ";
    }
    return os;
}
*/

map<int,int> pcounter(vector<int> idlist){                  //*ANCHOR for counting particle and maps them to the respective id's
    map<int, int> idmap;
    vector<int> sortedid=idlist;                    //initialized here, but sorted later
    int first;
    int n;
    sort(sortedid.begin(), sortedid.end());         // sort function only intakes pointer
    while (sortedid.size()){
        first=sortedid[0];
        n=0;
        int size=sortedid.size();                   // ! fixing size here to be used in while, keeping it directly in while loop causes recursive errors
        while(n<size && sortedid[0]==first){
                n++;
                sortedid.erase(sortedid.begin());
        }
        idmap.insert({first, n});
    }
    return idmap;        
};

//void amptrun(){                   //!Complete it
//    _execl("C:\\Windows\\System32\\wsl.exe")
//}




int main(){
    int sum;
//    extractor("ampt.dat");
    lister("ampt.dat",99);
    // map<int,int> lol=pcounter(listid);
    // cout << "\nKEY   ELEMENT\n";
    // for(auto x: lol){
	// cout << x.first << "-->" <<x.second <<endl;
    // sum=sum+x.second;
    // };
    // cout<<"no of elems:\t"<< sum<<"\n";
    // return 0;
    int g=0;    
     for (double l : listphi){
        g++;
     }
     cout<< g <<"\t";
};



