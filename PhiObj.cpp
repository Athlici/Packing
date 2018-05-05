class PhiCircCompl;
class PhiPolygon;
class PhiCircSeg;
class PhiHat;

//TODO: Can the symmetric definitions be handled with templates?

class PhiCompObj{
    public:
        virtual PhiFunc* phiFunc(RotTrans f, PhiCompObj* O, RotTrans g) = 0;
        virtual PhiFunc* phiFunc(RotTrans f, PhiPolygon* P, RotTrans g) = 0;
        virtual PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g) = 0;
        virtual PhiFunc* phiFunc(RotTrans f, PhiHat*     H, RotTrans g) = 0;
        virtual PhiFunc* phiFunc(RotTrans g, PhiCircCompl* Cc, Scale f) = 0;
};

class PhiInfObj{
    public:
        virtual PhiFunc* phiFunc(Scale f, PhiCompObj* B, RotTrans g) = 0;
};

class PhiCircCompl: public PhiInfObj{
    public:
        circle c;
        PhiCircCompl(circle ci) : c(ci) {};

        PhiFunc* phiFunc(Scale f, PhiCompObj* B, RotTrans g){
            return B->phiFunc(g, this, f);
        }
};

//class PhiLine: public PhiCompObj{};

class PhiCompNode: PhiCompObj{
        vector<PhiCompObj*> nodes;
    public:
        PhiCompNode(vector<PhiCompObj*> ni) : nodes(ni) {};

        template <class T> PhiFunc* phiFunc(RotTrans f,T* O, RotTrans g){
            vector<PhiFunc*> res(nodes.size());
            for(int i=0;i<nodes.size();i++)
                res[i] = nodes[i]->phiFunc(f,O,g);
            return new PhiFuncNode(true,res);
        }
};

class PhiPolygon: public PhiCompObj{
    public:
        vector<point> p;
        PhiPolygon(vector<point> pi) : p(pi) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiPolygon* Q, RotTrans g){
            int n = p.size(),m = Q->p.size();
            vector<PhiFunc*> comp(n+m);
            for(int i=0;i<n;i++){
                vector<PhiFunc*> tmp(m);
                for(int j=0;j<m;j++)
                    tmp[j] = new PhiFuncLnRtClRt(p[(i+1)%n],p[i],f,Q->p[j],g);
                    comp[i] = new PhiFuncNode(true,tmp);
            }
            for(int i=0;i<m;i++){
                vector<PhiFunc*> tmp(n);
                for(int j=0;j<n;j++)
                    tmp[j] = new PhiFuncLnRtClRt(Q->p[(i+1)%m],Q->p[i],g,p[j],f);
                comp[n+i] = new PhiFuncNode(true,tmp);
            }
            return new PhiFuncNode(false,comp);
        }

        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){}

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g){}

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            vector<PhiFunc*> comp(p.size());
            for(int i=0;i<p.size();i++)
                comp[i] = new PhiFuncScCcRtCl(C->c,f,p[i],g);
            return new PhiFuncNode(true,comp);
        }
};

class PhiCircSeg: public PhiCompObj{
    public:
        point p0,p1;
        circle pc;
        PhiCircSeg(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiPolygon* P, RotTrans g){
            return P->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){}

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g){}

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            return new PhiFuncNode(true,{
                new PhiFuncScCcRtCl(C->c,f,p0,g),
                new PhiFuncScCcRtCl(C->c,f,p1,g),
                new PhiFuncNode(false,{
                    new PhiFuncScCcRtCl(C->c,f,pc,g),
                    new PhiFuncHScCcRtCs(g,pc,p0, 1),
                    new PhiFuncHScCcRtCs(g,pc,p1,-1)})});
        }
};

class PhiHat: public PhiCompObj{
    public:
        point p0,p1,p2;
        circle pc;
        PhiHat(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiPolygon* P, RotTrans g){
            return P->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){
            return C->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g){}

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            return new PhiFuncNode(true,{
                new PhiFuncScCcRtCl(C->c,f,p0,g),
                new PhiFuncScCcRtCl(C->c,f,p1,g),
                new PhiFuncScCcRtCl(C->c,f,p2,g)});
        }
};

PhiFunc* phiFunc(PhiInfObj* A, Scale f, PhiCompObj* B, RotTrans g){
    return A->phiFunc(f,B,g);
}

PhiFunc* phiFunc(PhiCompObj* A, RotTrans f, PhiCompObj* B, RotTrans g){
    return A->phiFunc(f,B,g);
}
