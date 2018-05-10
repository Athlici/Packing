#include <iostream>
//#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <tuple>
#include <eigen3/Eigen/Dense>
#include <coin/IpIpoptApplication.hpp>

using std::vector,std::tuple,std::get,std::string;
using Eigen::Matrix2d,Eigen::Vector2d,Eigen::RowVector2d;

#include "Struct.cpp"
#include "Transform.cpp"
#include "PhiFunc.cpp"
#include "PhiObj.cpp"
#include "Objective.cpp"
#include "dNLP.cpp"
#include "Helpers.cpp"

int main(int argc, char** argv) {

//    PhiCompObj* P = regularPolygon(3);
//    double s2 = 1/sqrt(2);
//    PhiCompObj* P = new PhiCircSeg(point(-s2,s2),point(s2,s2),circle(point(0,0),1));

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetNumericValue("tol", 1e-9);
    app->Options()->SetStringValue("linear_solver", "ma57");
//    app->Options()->SetIntegerValue("max_iter", 100);

    ApplicationReturnStatus status = app->Initialize();

//    circle res = boundCircMod(app,P);
//    std::cout<< res.r << "\n";
//    std::cout<< res.p.x << "\n";
//    std::cout<< res.p.y << "\n";

    int n = 5;
    vector<circle> tmp(n);
    for(int i=0;i<n;i++) tmp[i] = circle(point(0,0),1);
    vector<point> pts = circlePack(tmp);
//    for(int i=0;i<n;i++) std::cout << "(" << pts[i].x << "," << pts[i].y << ")\n";

    Objective* obj = new FirstVar();
    vector<var> vars(3*n+1);
    vars[0] = var(0,2e19);
    for(int i=0;i<n;i++){
        vars[3*i+1] = var(-2e19,2e19);
        vars[3*i+2] = var(-2e19,2e19);
        vars[3*i+3] = var(-2*M_PI,2*M_PI);
    }

    Scale f = Scale(0);
    PhiInfObj* C = new PhiCircCompl(circle(point(0,0),1));

    vector<RotTrans> rt(n);
    for(int i=0;i<n;i++)
        rt[i] = RotTrans(3*i+1);
    PhiCompObj* P = regularPolygon(5);

    vector<PhiFunc*> foo(n*(n+1)/2);
    for(int i=0;i<n;i++)
        foo[i] = phiFunc(C,f,P,rt[i]);
    int ind = n;
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            foo[ind++] = phiFunc(P,rt[i],P,rt[j]);

    PhiFunc* phi = new PhiFuncNode(true,foo);

    double x[3*n+1];
    x[0]=1;
    for(int i=0;i<n;i++){
        x[3*i+1] = pts[i].x;
        x[3*i+2] = pts[i].y;
        x[3*i+3] = 0;
    }

    while(phi->eval(x)<0)
        x[0]*=2;

    bool newseg;
    do{
        SmartPtr<dNLP> nlp = new dNLP(obj,vars,phi->getIneqs(x),x);
        app->OptimizeTNLP(nlp);
        newseg = phi->getIneqs(x) != phi->getIneqs(nlp->res);
        for(int i=0;i<3*n+1;i++) x[i] = nlp->res[i];
    }while(newseg);

    std::cout << "{";
    for(int i=0;i<3*n+1;i++)
        std::cout << x[i] << ",";
    std::cout << "}\n";

}
