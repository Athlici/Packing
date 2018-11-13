class Model{
  public:
    vector<var> vars = {var(0,2e19)};
    Objective* f = new FirstVar();
    PhiInfObj* C;
    Scale sc = Scale(0);
    vector<PhiCompObj*> objs;
    vector<RotTrans> rt;

    Model() {};

    PhiFunc* createPhiFunc(){
        int n = objs.size();
        vector<PhiFunc*> comp(n*(n+1)/2);
        for(int i=0;i<n;i++)
            comp[i] = phiFunc(C,sc,objs[i],rt[i]);
        int ind = n;
        for(int i=0;i<n;i++)
            for(int j=i+1;j<n;j++)
                comp[ind++] = phiFunc(objs[i],rt[i],objs[j],rt[j]);
        return new PhiFuncNode(true,comp);
    }
};

class XMLInterface{
    XMLDocument doc;
    XMLNode* box;
    vector<XMLElement*> objsXML;

    PhiCompObj* parseObj(XMLElement* obj){
        if(!strcmp(obj->Name(),"Union")){
            vector<PhiCompObj*> nodes;
            for(XMLElement* n=obj->FirstChildElement();n!=nullptr;n=n->NextSiblingElement())
                nodes.push_back(parseObj(n));
            return new PhiCompNode(nodes);
        }
        if(!strcmp(obj->Name(),"Polygon")){
            vector<point> points;
            for(XMLElement* p=obj->FirstChildElement("Point");p!=nullptr;p=p->NextSiblingElement("Point"))
                points.push_back(toPoint(p));
            return new PhiPolygon(points);
        }if(!strcmp(obj->Name(),"CircSeg")){ //TODO: Check orientation, calculate cx,cy
            XMLElement* p0=obj->FirstChildElement("Point");
            XMLElement* p1=p0->NextSiblingElement("Point");
            XMLElement* c=obj->FirstChildElement("Circle");
            return new PhiCircSeg(toPoint(p0),toPoint(p1),toCircle(c));
        }if(!strcmp(obj->Name(),"Hat")){
            XMLElement* p0=obj->FirstChildElement("Point");
            XMLElement* p1=p0->NextSiblingElement("Point");
            XMLElement* p2=p1->NextSiblingElement("Point");
            XMLElement* c=obj->FirstChildElement("Circle");
            return new PhiHat(toPoint(p0),toPoint(p1),toPoint(p2),toCircle(c,-1));
        }
        return 0;
    }

  public:
    XMLInterface() {};

    bool load(char* path){
        return doc.LoadFile(path);
    }

    Model* parse(){
        Model* model = new Model();
        for(XMLNode* node=doc.FirstChild()->FirstChild();node!=nullptr;node=node->NextSibling()){
            if(!strcmp(node->Value(),"Container")){
                box = node;
                XMLElement* cont = box->FirstChildElement();
                if(!strcmp(cont->Value(),"CircCompl"))
                    model->C = new PhiCircCompl(toCircle(cont));
            }if(!strcmp(node->Value(),"Objects")){
                for(XMLElement* obj=node->FirstChildElement();obj!=nullptr;obj=obj->NextSiblingElement()){
                    objsXML.push_back(obj);
                    model->rt.push_back(RotTrans(model->vars.size()));
                    model->vars.push_back(var(-2e19,2e19));
                    model->vars.push_back(var(-2e19,2e19));
                    model->vars.push_back(var(-2*M_PI,2*M_PI));
                    model->objs.push_back(parseObj(obj));
                }
            }
        }
        return model;
    }

    bool write(vector<double> sol,vector<circle> bc){
        XMLElement* tmp = doc.NewElement("Solution");
        tmp->SetAttribute("r",sol[0]);
        box->InsertEndChild(tmp);
        for(int i=0;i<bc.size();i++){
            double s=sin(sol[3*i+3]),c=cos(sol[3*i+3]),x=bc[i].p.x,y=bc[i].p.y;
            tmp = doc.NewElement("Solution");
            tmp->SetAttribute("x",sol[3*i+1]+c*x-s*y);
            tmp->SetAttribute("y",sol[3*i+2]+s*x+c*y);
            tmp->SetAttribute("phi",sol[3*i+3]);
            objsXML[i]->InsertEndChild(tmp);
        }
        return doc.SaveFile("out.xml");
    }
};
