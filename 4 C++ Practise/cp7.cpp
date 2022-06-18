#include <iostream>         //General Stream
#include <fstream>          //Opening and Writing on files
#include <sstream>          //for selecting 
#include <string>

using namespace std;

int main() {
    double sum = 0;
    string x;
    ifstream inFile;
    
    inFile.open("ampt.dat");                              //fstream is used here to save the dat file into its object called InFile
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    
    while (getline(inFile, x)) {    
        istringstream ss(x);                              //storing string in a "istringstream" object called ss
        string id;
        string px;
        string py;
        string pz;
        string m;
        string x;
        string y;
        string z;
        string t;
        
        ss >> id;                         //String till first white space
        ss >> px;                         //String from first whire space till second white space
        ss >> py;                         //likewise
        ss >> pz;                         //Kochirakoso
        ss >> m;
        ss >> x;
        ss >> y;
        ss >> z;
        ss >> t;

        cout << stod(m) << "\n";           //Converting string to double
    }
    
    inFile.close();

}



