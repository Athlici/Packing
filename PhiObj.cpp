class PhiCircCompl;

class PhiCompObj{
    public:
        virtual PhiFunc* phiFunc(RotTrans f, PhiCompObj* B, RotTrans g) = 0;
        virtual PhiFunc* phiFuncCC(Translate g, PhiCircCompl* C, Scale f) = 0;
};

class PhiInfObj{
    public:
        virtual PhiFunc* phiFunc(Scale f, PhiCompObj* B, Translate g) = 0;
};

class PhiCircCompl: public PhiInfObj{
    public:
        circle c;
        PhiCircCompl(circle ci) : c(ci) {};

        PhiFunc* phiFunc(Scale f, PhiCompObj* B, Translate g){        
            return B->phiFuncCC(g, this, f);
        }
};

//class PhiLine: public PhiCompObj{};

class PhiCompNode: PhiCompObj{
        vector<PhiCompObj*> nodes;
    public:
        PhiCompNode(vector<PhiCompObj*> ni) : nodes(ni) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* B, RotTrans g){
            vector<PhiFunc*> res(nodes.size());
            for(int i=0;i<nodes.size();i++)
                res[i] = nodes[i]->phiFunc(f,B,g);
            return new PhiFuncNode(true,res);
        }

        PhiFunc* phiFuncCC(Translate g,PhiCircCompl* C, Scale f){
            vector<PhiFunc*> res(nodes.size());
            for(int i=0;i<nodes.size();i++)
                res[i] = nodes[i]->phiFuncCC(g,C,f);
            return new PhiFuncNode(true,res);
        }
};

class PhiPolygon: public PhiCompObj{
    public:
        vector<point> p;
        PhiPolygon(vector<point> pi) : p(pi) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* B, RotTrans g){}

        PhiFunc* phiFuncCC(Translate g, PhiCircCompl* C, Scale f){
            vector<PhiFunc*> comp(p.size());
            for(int i=0;i<p.size();i++)
                comp[i] = new PhiFuncScCcTrCl(C->c,f,p[i],g);
            return new PhiFuncNode(true,comp);
        }
};

class PhiCircSeg: public PhiCompObj{
    public:
        point p0,p1;
        circle pc;
        PhiCircSeg(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* B, RotTrans g){}

        PhiFunc* phiFuncCC(Translate g, PhiCircCompl* C, Scale f){
            return new PhiFuncNode(true,{
                new PhiFuncScCcTrCl(C->c,f,p0,g),
                new PhiFuncScCcTrCl(C->c,f,p1,g),
                new PhiFuncNode(false,{
                    new PhiFuncScCcTrCl(C->c,f,pc,g),
                    new PhiFuncHScCcTrCs(g,pc,p0, 1),
                    new PhiFuncHScCcTrCs(g,pc,p1,-1)})});
        }
};

class PhiHat: public PhiCompObj{
    public:
        point p0,p1,p2;
        circle pc;
        PhiHat(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* B, RotTrans g){}

        PhiFunc* phiFuncCC(Translate g, PhiCircCompl* C, Scale f){
            return new PhiFuncNode(true,{
                new PhiFuncScCcTrCl(C->c,f,p0,g),
                new PhiFuncScCcTrCl(C->c,f,p1,g),
                new PhiFuncScCcTrCl(C->c,f,p2,g)});
        }
};

PhiFunc* phiFunc(PhiInfObj* A, Scale f, PhiCompObj* B, Translate g){
    return A->phiFunc(f,B,g);
}

PhiFunc* phiFunc(PhiCompObj* A, RotTrans f, PhiCompObj* B, RotTrans g){
    return A->phiFunc(f,B,g);
}
