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
#include "Helpers.cpp"
#include "Objective.cpp"
#include "dNLP.cpp"

int main(int argc, char** argv) {

//Set up geometry of objects
    PhiCircCompl C = PhiCircCompl(circle(point(0,0),1));
    PhiPolygon P = regularPolygon(3);

    Scale f = Scale(0);
    Translate g = Translate(1);

    PhiFunc* phi = phiFunc(C,f,P,g);

    double x[] = {0.75,0,0};
    while(phi->eval(x)<0){
        x[0]*=2;
    }
    Objective* obj = new FirstVar();
    vector<var> vars = {var(0,2e19),var(-2e19,2e19),var(-2e19,2e19)};
    SmartPtr<TNLP> boundcirc = new dNLP(obj,vars,phi->getIneqs(x),x);


    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetNumericValue("tol", 1e-9);
//    app->Options()->SetStringValue("mu_strategy", "adaptive");
//    app->Options()->SetStringValue("derivative_test", "second-order");
//    app->Options()->SetIntegerValue("max_iter", 2);

    ApplicationReturnStatus status = app->Initialize();

    status = app->OptimizeTNLP(boundcirc);

//Calculate PhiFunctions
//    PhiFunc(A,B);

//Find initial configuration

//Iterate until stationary
//  Set up IPOPT-Modell
//  Run IPOPT-Modell to new configuration

}
