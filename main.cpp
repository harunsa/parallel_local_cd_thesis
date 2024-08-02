#include "../../dyncomm/DynCommDriver.cpp"
using namespace Spectra;
using namespace std;
#include <bulk/bulk.hpp>
#include <bulk/backends/thread/thread.hpp>


int main() {
    DynCommDriver::multiscale_driver("synthV1000T100.txt");
    return 0;
}

