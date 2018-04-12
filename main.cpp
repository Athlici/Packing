#include <iostream>
//#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <tuple>

using std::vector,std::tuple,std::get,std::string;

#include "IpIpoptApplication.hpp"
#include "Struct.cpp"
#include "Transform.cpp"
#include "PhiObj.cpp"
#include "PhiFunc.cpp"
#include "Objective.cpp"
#include "dNLP.cpp"
#include "Helpers.cpp"

int main(int argc, char** argv) {

    PhiPolygon P = regularPolygon(3);

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetNumericValue("tol", 1e-9);
//    app->Options()->SetStringValue("mu_strategy", "adaptive");
//    app->Options()->SetIntegerValue("max_iter", 100);

    ApplicationReturnStatus status = app->Initialize();

    vector<double> res = boundCircMod(app,P);

    for(int i=0;i<3;i++)
        std::cout<< res[i] << "\n";
}
