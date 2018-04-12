class PhiCompObj{
    public:
        
};

class PhiCircCompl: public PhiCompObj{
    public:
        circle c;
        PhiCircCompl(circle ci) : c(ci) {};
};

class PhiLine: public PhiCompObj{
};

class PhiPolygon: public PhiCompObj{
    public:
        vector<point> p;
        PhiPolygon(vector<point> pi) : p(pi) {};
};

class PhiCircSeg: public PhiCompObj{
    public:
        point p0,p1;
        circle pc;
        PhiCircSeg(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci) {};
};

class PhiHat: public PhiCompObj{
    public:
        point p0,p1;
        circle pc;
        PhiHat(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci) {};
};

//PhiFunc phiFunc(PhiCompObj A, Transform f, PhiCompObj B, Transform g){
//}
