// Minimalistic GRBMaker for semiparametric regression
// Based on work by Jean-Baptiste Sauvan and Josh Bendavid

#include <iostream>
#include <string>
#include <time.h>

#include "SemiparametricGBRMaker.h"

int main(int argc, char** argv)
{
    if(argc!=2)
    {
        std::cout << "Usage: regression.exe configurationFile\n";
        return 1;
    }

    std::cout << "INFO: Creating new regression\n";

    // measure start time
    time_t tStart, tEnd;
    time(&tStart);

    SemiparametricGBRMaker* regMaker = new SemiparametricGBRMaker();
    std::string parameterFile(argv[1]);
    regMaker->readParameterFile( parameterFile );
    regMaker->run();
    regMaker->close();

    // measure end time
    time(&tEnd);
    double t = difftime(tEnd, tStart);
    int hours = int(t)/3600;
    int min = (int(t)/60)%60;
    std::cout << "INFO: Elapsed time = " << t << " s\n";
    std::cout << "                   = " << hours << " h " << min << " min\n"; 

}