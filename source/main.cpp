#include "forman/formangradient.h"

using namespace std;

int main(int argc, char* argv[])
{

    Timer time;

    //reading the input
    FormanGradient grad = FormanGradient(argc,argv);

    time.start();
    grad.computeFormanGradient(true);
    time.stop();
    cout << "Forman gradient computed in " << time.getElapsedTime() << " seconds" << endl;


    time.start();
    grad.computePersistentHomology();
    time.stop();
    cout << "Persistent homology computed in " << time.getElapsedTime() << " seconds" << endl;

    return 0;
}

