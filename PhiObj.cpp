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
        tuple<point,point> moveLine(point p0,point p1,double r){
            double dx = p0.x-p1.x, dy = p0.y-p1.y;
            double c = r/sqrt(dx*dx+dy*dy);
            return tuple<point,point>(point(p0.x-dy*c,p0.y+dx*c),
                                      point(p1.x-dy*c,p1.y+dy*c));
        }

    public:
        vector<point> p;
        int n;

        PhiPolygon(vector<point> pi) : p(pi),n(pi.size()) {};

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f,point q,RotTrans g){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = new PhiFuncLnRtClRt(p[i],p[(i+1)%n],f,q,g);
            return new PhiFuncNode(false,comp);
        }

        PhiFunc* phiFunc(RotTrans f,circle c,RotTrans g){
            vector<PhiFunc*> comp(2*n);
            vector<point>    exts(2*n);
            for(int i=0;i<n;i++){
                tuple<point,point> tmp = moveLine(p[i],p[(i+1)%n],c.r);
                exts[2*i  ] = get<0>(tmp);
                exts[2*i+1] = get<1>(tmp);
            }
            for(int i=0;i<n;i++){
                comp[2*i  ] = new PhiFuncLnRtClRt(exts[2*i],exts[2*i+1],f,c.p,g);
                comp[2*i+1] = new PhiFuncNode(true,{
                    new PhiFuncLnRtClRt(exts[(n+2*i-1)%n],exts[2*i],f,c.p,g),
                    new PhiFuncClRtClRt(c,g,p[i],f)});
            }
            return new PhiFuncNode(false,comp);
        }

        PhiFunc* phiFunc(RotTrans f,point p1,point p2,RotTrans g){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = new PhiFuncLnRtClRt(p1,p2,g,p[i],f);
            return new PhiFuncNode(true,comp);
        }

        PhiFunc* phiFunc(RotTrans f, PhiPolygon* Q, RotTrans g){
            int m = Q->n;
            vector<PhiFunc*> comp(n+m);
            for(int i=0;i<n;i++)
                comp[i] = Q->phiFunc(g,p[i],p[(i+1)%n],f);
            for(int i=0;i<m;i++)
                comp[n+i] = phiFunc(f,Q->p[i],Q->p[(i+1)%m],g);
            return new PhiFuncNode(false,comp);
        }

        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g);
    
        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g);

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = new PhiFuncScCcRtCl(C->c,f,p[i],g);
            return new PhiFuncNode(true,comp);
        }
};

class PhiCircSeg: public PhiCompObj{
        static point circCompletion(point p0,point p1,point pc){
            double x0=p0.x, y0=p0.y, x1=p1.x, y1=p1.y, xc=pc.x, yc=pc.y;
            double c = ((x0-xc)*(x1-xc)+(y0-yc)*(y1-yc))/
                       (xc*(y0-y1)+yc*(x1-x0)+x0*y1-x1*y0);
            return point(xc-x0-x1+c*(y0-y1),yc-y0-y1+c*(x0-x1));
        }
    public:
        point p0,p1;
        circle pc;
        PhiPolygon P;

        PhiCircSeg(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci),
            P(PhiPolygon({p0,p1,this->circCompletion(p0i,p1i,pci.p)})) {};

        template <class T> PhiFunc* phiFunc(RotTrans f,T* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){
            return new PhiFuncNode(false,{
                P.phiFunc(f,C->pc,g),C->P.phiFunc(g,pc,f),
                P.phiFunc(f,&C->P,g),new PhiFuncClRtClRt(pc,f,C->pc,g)});
        }

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g);

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

PhiFunc* PhiPolygon::phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){
    return new PhiFuncNode(false,{C->P.phiFunc(g,this,f),this->phiFunc(f,C->pc,g)});
}

class PhiHat: public PhiCompObj{
    public:
        point p0,p1,p2;
        circle pc;
        PhiPolygon P;

        PhiHat(point p0i,point p1i,point p2i,circle pci) :
            p0(p0i),p1(p1i),p2(p2i),pc(pci),P(PhiPolygon({p0i,p1i,p2i})) {};

        template <class T> PhiFunc* phiFunc(RotTrans f,T* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFuncG(RotTrans f, PhiPolygon* Q, RotTrans g){
            int m = Q->n;
            vector<PhiFunc*> comp(2*m+2);    //P1,P2 not in Q, if q over l then q in C
            for(int i=0;i<m;i++)
                comp[i] = new PhiFuncNode(false,{
                    new PhiFuncLnRtClRt(p0,p1,f,Q->p[i],g),
                    new PhiFuncClRtClRt(pc,f,Q->p[i],g,-1)});
            comp[m  ] = Q->phiFunc(g,p0,f);
            comp[m+1] = Q->phiFunc(g,p1,f);
            return new PhiFuncNode(true,comp);
        }

        PhiFunc* phiFuncG(RotTrans f, PhiCircSeg* D, RotTrans g){
            vector<PhiFunc*> comp(5);
//            return new PhiFuncNode(false,{})
        }

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g){}

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            return new PhiFuncNode(true,{
                new PhiFuncScCcRtCl(C->c,f,p0,g),
                new PhiFuncScCcRtCl(C->c,f,p1,g),
                new PhiFuncScCcRtCl(C->c,f,p2,g)});
        }
};

PhiFunc* PhiPolygon::phiFunc(RotTrans f, PhiHat* H, RotTrans g){
    return new PhiFuncNode(false,{H->P.phiFunc(g,this,f),H->phiFuncG(g,this,f)});
}

PhiFunc* PhiCircSeg::phiFunc(RotTrans f, PhiHat* H, RotTrans g){
    return new PhiFuncNode(false,{H->P.phiFunc(g,this,f),H->phiFuncG(g,this,f)});
}

PhiFunc* phiFunc(PhiInfObj* A, Scale f, PhiCompObj* B, RotTrans g){
    return A->phiFunc(f,B,g);
}

PhiFunc* phiFunc(PhiCompObj* A, RotTrans f, PhiCompObj* B, RotTrans g){
    return A->phiFunc(f,B,g);
}
