//#include <stdint.h>
#include <tuple>
#include <math.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <coin/IpIpoptApplication.hpp>
#include "tinyxml2.h"

using std::vector,std::tuple,std::get,std::string;
using Eigen::Matrix2d,Eigen::Vector2d,Eigen::RowVector2d;
using Eigen::Matrix3d,Eigen::Vector3d,Eigen::RowVector3d;
using namespace tinyxml2;

#include "Struct.cpp"
#include "Transform.cpp"
#include "PhiFunc.cpp"
#include "PhiObj.cpp"
#include "Objective.cpp"
#include "dNLP.cpp"
#include "Helpers.cpp"
#include "Containers.cpp"

double randperm=20,randorient=5;

int main(int argc, char** argv) {

    //Ipopt initilization
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", 2);
    app->Options()->SetIntegerValue("max_iter", 1000);
//    app->Options()->SetNumericValue("tol", 1e-6);
    app->Options()->SetStringValue("linear_solver", "ma57");
//    app->Options()->SetStringValue("linear_solver", "mumps");
    app->Options()->SetStringValue("accept_every_trial_step", "yes"); //if it's stupid and it works, is it really stupid?
//    app->Options()->SetStringValue("derivative_test", "second-order");
    ApplicationReturnStatus status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << "Error: Couldn't initialize IPOPT!\n";
        return (int) status;
    }

    std::default_random_engine randgen(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> pidist(-M_PI,M_PI);

    //read objects from xml input file
    XMLInterface xml;
    if(xml.load(argv[1])){
        std::cout << "Error: Couldn't load input file!\n";
        return 1;
    }
    Model* model = xml.parse();

    int n = model->objs.size();

    //calculate bounding circles for all objects
    vector<circle> bc(n);
    for(int i=0;i<n;i++){
        bc[i] = boundCircMod(app,model->objs[i]);
        model->objs[i]->move(bc[i].p); //move to bounding circle to origin
    }

    //create distance functions for all objects
    PhiFunc* phi = model->createPhiFunc();

    //circle packing from random object permutation
    vector<vector<double>> initials(randperm*randorient);
    for(int k=0;k<randperm;k++){
        vector<int> perm(n);
        vector<double> bcp(n);
        for(int i=0;i<n;i++) perm[i] = i;
        shuffle(perm.begin(),perm.end(),randgen);
        for(int i=0;i<n;i++) bcp[perm[i]] = bc[i].r;
        vector<point> pts(n),ptsp = circlePack(bcp);
        for(int i=0;i<n;i++) pts[i] = ptsp[perm[i]];

        for(int l=0;l<randorient;l++){
            //create feasible initial parameters
            vector<double> x(3*n+1);
            x[0]=1;
            for(int i=0;i<n;i++){
                double t=pidist(randgen),s=sin(t),c=cos(t);
                x[3*i+1] = pts[i].x;
                x[3*i+2] = pts[i].y;
                x[3*i+3] = t;
            }
            for(int i=0;i<10&&phi->eval(x.data())<-0.001;i++)
                x[0]*=2;
            initials[randorient*k+l] = x;
        }
    }

    //iterative solver
    Objective* f = model->f;
    vector<double> bestsol = initials[0];
    for(int i=0;i<randperm*randorient;i++){
        double* x = initials[i].data();
        bool notstalled;
        do{
            double prior = f->eval(n,x);
            if(phi->eval(x) < 0.01){
                x[0]*=1.03;
                for(int j=0;j<n;j++){
                    x[3*j+1]*=1.03;
                    x[3*j+2]*=1.03;
                }
            }
            SmartPtr<dNLP> nlp = new dNLP(f,model->vars,phi->getIneqs(x),x);
            status = app->OptimizeTNLP(nlp);
            double fv=f->eval(n,nlp->res);
            if(status == Solve_Succeeded && fv < prior){
                if(fv < f->eval(n,bestsol.data()) && phi->eval(nlp->res) > -0.001){
                    for(int j=0;j<3*n+1;j++)
                        bestsol[j] = nlp->res[j];
                    std::cout << fv << "\n";
                }
                if(phi->getIneqs(x) != phi->getIneqs(nlp->res)){
                    notstalled = true;
                    for(int j=0;j<3*n+1;j++)
                        x[j] = nlp->res[j];
                } else
                    notstalled = false;
            } else 
                notstalled = false;
        }while(notstalled);
    }

//    vector<string> res = phi->print(x);
//    for(int i=0;i<res.size();i++)
//      std::cout<<res[i]<<"\n";

    //write solution into XML file
    xml.write(bestsol,bc);
}
