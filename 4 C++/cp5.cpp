#include <iostream>

void name_x_times(std::string name, int x) {
    while (x > 0) {
        std::cout << name;
        x--;

    }
}

int main() {

    std::string my_name = "Alan";
    int some_number = 5; 
    name_x_times(my_name, some_number);

}