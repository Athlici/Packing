#include <iostream>
//#include <stdint.h>
#include <math.h>
#include <string>
#include <vector>
#include <tuple>
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

int main(int argc, char** argv) {

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetNumericValue("tol", 1e-9);
    app->Options()->SetStringValue("linear_solver", "ma57");
//    app->Options()->SetStringValue("derivative_test", "second-order");

    ApplicationReturnStatus status = app->Initialize();

    vector<var> vars = {var(0,2e19)};
    Objective* f = new FirstVar();
    PhiInfObj* C;
    Scale sc = Scale(0);
    XMLNode* box;
    vector<PhiCompObj*> objs;
    vector<RotTrans> rt;
    vector<XMLElement*> objsXML;

    XMLDocument doc;
    doc.LoadFile(argv[1]);
    XMLNode* root = doc.FirstChild();
    for(XMLNode* node=root->FirstChild();node!=nullptr;node=node->NextSibling()){
        if(!strcmp(node->Value(),"Container")){
            box = node;
            XMLElement* cont = box->FirstChildElement();
            if(!strcmp(cont->Value(),"CircCompl"))
                C = new PhiCircCompl(toCircle(cont));
        }if(!strcmp(node->Value(),"Objects")){
            for(XMLElement* obj=node->FirstChildElement();obj!=nullptr;obj=obj->NextSiblingElement()){
                rt.push_back(RotTrans(vars.size()));
                vars.push_back(var(-2e19,2e19));
                vars.push_back(var(-2e19,2e19));
                vars.push_back(var(-2*M_PI,2*M_PI));
                objsXML.push_back(obj);
                if(!strcmp(obj->Name(),"Polygon")){
                    vector<point> points;
                    for(XMLElement* p=obj->FirstChildElement("Point");p!=nullptr;p=p->NextSiblingElement("Point"))
                        points.push_back(toPoint(p));
                    objs.push_back(new PhiPolygon(points));
                }if(!strcmp(obj->Name(),"CircSeg")){ //TODO: Check orientation, calculate cx,cy
                    XMLElement* p0=obj->FirstChildElement("Point");
                    XMLElement* p1=p0->NextSiblingElement("Point");
                    XMLElement* c=obj->FirstChildElement("Circle");
                    objs.push_back(new PhiCircSeg(toPoint(p0),toPoint(p1),toCircle(c)));
                }if(!strcmp(obj->Name(),"Hat")){
                    XMLElement* p0=obj->FirstChildElement("Point");
                    XMLElement* p1=p0->NextSiblingElement("Point");
                    XMLElement* p2=p1->NextSiblingElement("Point");
                    XMLElement* c=obj->FirstChildElement("Circle");
                    objs.push_back(new PhiHat(toPoint(p0),toPoint(p1),toPoint(p2),toCircle(c,-1)));
                }
            }
        }
    }

    int n = objs.size();
    vector<circle> bc(n);
    for(int i=0;i<n;i++) bc[i] = boundCircMod(app,objs[i]);

    vector<point> pts = circlePack(bc);
    vector<PhiFunc*> comp(n*(n+1)/2);
    for(int i=0;i<n;i++)
        comp[i] = phiFunc(C,sc,objs[i],rt[i]);
    int ind = n;
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            comp[ind++] = phiFunc(objs[i],rt[i],objs[j],rt[j]);
    PhiFunc* phi = new PhiFuncNode(true,comp);

    double x[3*n+1];
    x[0]=1;
    for(int i=0;i<n;i++){
        x[3*i+1] = pts[i].x+bc[i].p.x;
        x[3*i+2] = pts[i].y+bc[i].p.y;
        x[3*i+3] = 0;
    }
    for(int i=0;i<10&&phi->eval(x)<-0.01;i++){
        x[0]*=2;
    }   //TODO: Contingency

    bool newseg;
    do{
        SmartPtr<dNLP> nlp = new dNLP(f,vars,phi->getIneqs(x),x);
        app->OptimizeTNLP(nlp);
//        newseg = phi->getIneqs(x) != phi->getIneqs(nlp->res);
        newseg = nlp->res[0] < x[0];
        for(int i=0;i<3*n+1;i++) x[i] = nlp->res[i];
    }while(newseg);

    vector<string> res = phi->print(x);
    for(int i=0;i<res.size();i++)
      std::cout<<res[i]<<"\n";

    XMLElement* tmp = doc.NewElement("Solution");
    tmp->SetAttribute("r",x[0]);
    box->InsertEndChild(tmp);
    for(int i=0;i<n;i++){
        tmp = doc.NewElement("Solution");
        tmp->SetAttribute("x",x[3*i+1]);
        tmp->SetAttribute("y",x[3*i+2]);
        tmp->SetAttribute("phi",x[3*i+3]);
        objsXML[i]->InsertEndChild(tmp);
    }
    doc.SaveFile("out.xml");
}
