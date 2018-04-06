//class PhiCompObjPrim{};

class PhiCircCompl{
    public:
        circle c;
        PhiCircCompl(circle ci) : c(ci) {};
};

class PhiPolygon{
    public:
        vector<point> p;
        PhiPolygon(vector<point> pi) : p(pi) {};
};

class PhiHat{
    public:
        point p0,p1,pc;
        double r;
};
