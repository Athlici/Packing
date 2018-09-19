#include <iostream>
#include <math.h>

#include "tinyxml2.h"

using namespace tinyxml2;
using std::to_string;

int main(int argc, char** argv) {

    XMLDocument xmlDoc;
    XMLNode* pBox = xmlDoc.NewElement("Container");
    xmlDoc.InsertFirstChild(pBox);
    XMLNode* pObjs = xmlDoc.NewElement("Objects");
    xmlDoc.InsertEndChild(pObjs);

    XMLElement* cc = xmlDoc.NewElement("CircCompl");
    cc->SetAttribute("x",0.0);
    cc->SetAttribute("y",0.0);
    cc->SetAttribute("r",1.0);
    pBox->InsertFirstChild(cc);

    int n=5;

    for(int i=0;i<n;i++){
        double phi=2*M_PI/n;
        XMLElement* pol = xmlDoc.NewElement("Polygon");
        for(int j=0;j<5;j++){
            XMLElement* point = xmlDoc.NewElement("Point");
            point->SetAttribute("x",to_string(cos(phi*j)).c_str());
            point->SetAttribute("y",to_string(sin(phi*j)).c_str());
            pol->InsertEndChild(point);
        }
        pObjs->InsertEndChild(pol);
    }

//    double s2 = 1/sqrt(2);
//    PhiCompObj* P = new PhiCircSeg(point(-s2,s2),point(s2,s2),circle(point(0,0),1));

    xmlDoc.SaveFile("test.xml");
}
