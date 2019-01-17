#include "psass.h"


int main(int argc, char *argv[]) {

    std::cout << "PSASS started." << std::endl;

    Psass psass(argc, argv);
    psass.run();

    std::cout << "\nPSASS ended successfully." << std::endl;
    return 0;
}
