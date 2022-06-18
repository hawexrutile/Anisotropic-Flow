#include <iostream>

// Define needs_water() here:

std::string needs_water(int days, bool is_succulent) {
    if (is_succulent && days < 13) {
        return "Don't water the plant!";
    }
    else if (is_succulent && days > 12) {
        return"Go ahead and give the plant a little water.";
    }
    else {
        return"Time to water the plant.";
    }
}

int main() {

    std::cout << needs_water(10, false) << "\n";

}