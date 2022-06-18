#include <iostream>
int main() {


  double tempf;
  double tempc;
  std::cout<<"enter the temperature in Fahrenheit \n";
  std::cin>>tempf;
  
  
  
  // Ask the user
  
  
  tempc = (tempf - 32) / 1.8;
  
  std::cout << "The temp is " << tempc << " degrees Celsius.\n";
  
}