class PhiCircCompl;     //Circle Complement
class PhiLineCompl;     //Line Complement
class PhiPolygon;       //Polygon
class PhiCircSeg;       //Circle Segment
class PhiHat;           //Hat (Circle Complement Segment)

class PhiCompObj{       //Composite Objects (includes primitives)
    public:
        //every object implements a distance function to all other objects
        virtual PhiFunc* phiFunc(RotTrans, PhiCompObj*, RotTrans) = 0;
        virtual PhiFunc* phiFunc(RotTrans, PhiPolygon*, RotTrans) = 0;
        virtual PhiFunc* phiFunc(RotTrans, PhiCircSeg*, RotTrans) = 0;
        virtual PhiFunc* phiFunc(RotTrans, PhiHat*    , RotTrans) = 0;
        virtual PhiFunc* phiFunc(RotTrans, PhiCircCompl*,  Scale) = 0;
        virtual PhiFunc* phiFunc(RotTrans, PhiLineCompl*,  Scale) = 0;

        //helper function to virtually move composite objects
        virtual void move(point) = 0;
};
 
class PhiInfObj{        //Infinite Objects
    public:
        virtual PhiFunc* phiFunc(Scale f, PhiCompObj* B, RotTrans g) = 0;
};

class PhiCircCompl: public PhiInfObj{
    public:
        circle c;
        PhiCircCompl(circle ci) : c(ci) {}

        PhiFunc* phiFunc(Scale f, PhiCompObj* B, RotTrans g){
            return B->phiFunc(g, this, f);
        }
};

class PhiLineCompl: public PhiInfObj{
    public:
        point p0,p1;
        PhiLineCompl(point p0i, point p1i) : p0(p0i),p1(p1i) {}

        PhiFunc* phiFunc(Scale f, PhiCompObj* B, RotTrans g){
            return B->phiFunc(g, this, f);
        }
};

//class PhiLine: public PhiCompObj{};

//Union of Composite Objects
class PhiCompNode: public PhiCompObj{
        vector<PhiCompObj*> nodes;
    public:
        PhiCompNode(vector<PhiCompObj*> ni) : nodes(ni) {}

        template <class S,class T> PhiFunc* phiFuncM(RotTrans f,S* O, T g){
            vector<PhiFunc*> res(nodes.size());
            for(int i=0;i<nodes.size();i++)
                res[i] = nodes[i]->phiFunc(f,O,g);
            return new PhiFuncNode(true,res);
        }
        
        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O,RotTrans g){
            return phiFuncM(f,O,g);
        }
        
        PhiFunc* phiFunc(RotTrans f,PhiPolygon* P,RotTrans g){
            return phiFuncM(f,P,g);
        }
        
        PhiFunc* phiFunc(RotTrans f,PhiCircSeg* D,RotTrans g){
            return phiFuncM(f,D,g);
        }
        
        PhiFunc* phiFunc(RotTrans f,PhiHat* H,RotTrans g){
            return phiFuncM(f,H,g);
        }
        
        PhiFunc* phiFunc(RotTrans f,PhiCircCompl* C,Scale g){
            return phiFuncM(f,C,g);
        }

        PhiFunc* phiFunc(RotTrans f, PhiLineCompl* L, Scale g){
            return phiFuncM(f,L,g);
        }

        void move(point p){
            for(PhiCompObj* obj : nodes)
                obj->move(p);
        }
};

//Polygon (assumed to be convex)
class PhiPolygon: public PhiCompObj{
    public:
        vector<point> p;    //Counter clockwise enumeration of points
        int n;

        PhiPolygon(vector<point> pi) : p(pi),n(pi.size()) {}

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        //a point is outside iff it is to the right of one side
        PhiFunc* phiFunc(RotTrans f,point q,RotTrans g){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = PhiFunc::phiFunc(p[i],p[(i+1)%n],f,q,g);
            return new PhiFuncNode(false,comp);
        }

        //a circle is outside iff its center is outside the minkowsky sum 
        PhiFunc* phiFunc(RotTrans f,circle c,RotTrans g){
            vector<PhiFunc*> comp(2*n);
            vector<point>    exts(2*n);
            for(int i=0;i<n;i++){
                tuple<point,point> tmp = moveLine(p[i],p[(i+1)%n],c.r);
                exts[2*i  ] = get<0>(tmp);
                exts[2*i+1] = get<1>(tmp);
            }
            for(int i=0;i<n;i++){
                comp[2*i  ] = PhiFunc::phiFunc(exts[2*i],exts[2*i+1],f,c.p,g);
                comp[2*i+1] = new PhiFuncNode(true,{
                    PhiFunc::phiFunc(exts[(2*n+2*i-1)%(2*n)],exts[2*i],f,c.p,g),
//                    PhiFunc::phiFunc(exts[2*i],exts[(n+2*i-1)%n],f,c.p,g),
                    PhiFunc::phiFunc(c,g,p[i],f)});
            }
            return new PhiFuncNode(false,comp);
        }

        PhiFunc* phiFunc(RotTrans f,point p1,point p2,RotTrans g){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = PhiFunc::phiFunc(p1,p2,g,p[i],f);
            return new PhiFuncNode(true,comp);
        }

        //a convex and concave polygon don't intersect iff no edges are contained
        PhiFunc* phiFunc(RotTrans f, PhiPolygon* Q, RotTrans g){
            int m = Q->n;
            vector<PhiFunc*> comp(n+m);
            for(int i=0;i<n;i++)
                comp[i] = Q->phiFunc(g,p[i],p[(i+1)%n],f);
            for(int i=0;i<m;i++)
                comp[n+i] = phiFunc(f,Q->p[i],Q->p[(i+1)%m],g);
            return new PhiFuncNode(false,comp);
        }

        //implemented in the complementary classes
        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g);
        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g);

        //no intersection iff no point is contained
        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = PhiFunc::phiFunc(C->c,f,p[i],g);
            return new PhiFuncNode(true,comp);
        }

        //no intersection iff no point is contained
        PhiFunc* phiFunc(RotTrans g, PhiLineCompl* L, Scale f){
            vector<PhiFunc*> comp(n);
            for(int i=0;i<n;i++)
                comp[i] = PhiFunc::phiFunc(L->p0,L->p1,f,p[i],g);
            return new PhiFuncNode(true,comp);
        }

        void move(point pd){
            for(int i=0;i<n;i++)
                p[i].move(pd);
        }
};

class PhiCircSeg: public PhiCompObj{
        //find the point completing the bounding triangle of the circle segment
        static point circCompletion(point p0,point p1,point pc){
            double x0=p0.x, y0=p0.y, x1=p1.x, y1=p1.y, xc=pc.x, yc=pc.y;
            double c = ((x0-xc)*(x1-xc)+(y0-yc)*(y1-yc))/
                       (xc*(y0-y1)+x0*(y1-yc)+x1*(yc-y0));
            return point(x0+x1-xc+c*(y0-y1),y0+y1-yc+c*(x1-x0));
        }
    public:
        point p0,p1;    //Eliminate as well?
        circle pc;
        PhiPolygon P;

        PhiCircSeg(point p0i,point p1i,circle pci) : p0(p0i),p1(p1i),pc(pci),
            P(PhiPolygon({p0,p1,this->circCompletion(p0i,p1i,pci.p)})) {
        }

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f,PhiPolygon* P, RotTrans g){
            return P->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){
            return new PhiFuncNode(false,{
                P.phiFunc(f,C->pc,g),C->P.phiFunc(g,pc,f),
                P.phiFunc(f,&C->P,g),PhiFunc::phiFunc(pc,f,C->pc,g)});
        }

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g);

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, RotTrans f){
            vector<PhiFunc*> comp = {PhiFunc::phiFunc(C->c,f,p0,g),
                                     PhiFunc::phiFunc(C->c,f,p1,g)};
            if(C->c.r+pc.r>0)
                comp.push_back(new PhiFuncNode(false,{
                      PhiFunc::phiFunc(C->c,f,pc,g),
                      PhiFunc::phiFunc(C->c,f,pc,p0,g, 1),
                      PhiFunc::phiFunc(C->c,f,pc,p1,g,-1)}));
            return new PhiFuncNode(true,comp);
        }

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            return new PhiFuncNode(true,{
                PhiFunc::phiFunc(C->c,f,p0,g),
                PhiFunc::phiFunc(C->c,f,p1,g),
                new PhiFuncNode(false,{ //TODO add C->pc.r+pc.r>=0
                    PhiFunc::phiFunc(C->c,f,pc,g),
                    PhiFunc::phiFunc(C->c,f,pc,p0,g,-1),
                    PhiFunc::phiFunc(C->c,f,pc,p1,g, 1)})});
        }

        PhiFunc* phiFunc(RotTrans g, PhiLineCompl* L, Scale f){
            return new PhiFuncNode(true,{
                PhiFunc::phiFunc(L->p0,L->p1,f,p0,g),
                PhiFunc::phiFunc(L->p0,L->p1,f,p1,g),
                new PhiFuncNode(false,{
                    PhiFunc::phiFunc(L->p0,L->p1,f,pc,g),
                    PhiFunc::phiFunc(L->p0,L->p1,f,P.p[2],g)})});
        }

        void move(point pd){
            p0.move(pd);p1.move(pd);pc.move(pd);
            P.move(pd);
        }
};

PhiFunc* PhiPolygon::phiFunc(RotTrans f, PhiCircSeg* C, RotTrans g){
    return new PhiFuncNode(false,{C->P.phiFunc(g,this,f),this->phiFunc(f,C->pc,g)});
}

class PhiHat: public PhiCompObj{
    public:
        point p0,p1,p2;     //Get rid of this?
        circle pc;          //Radius needs to be negative
        PhiPolygon P;

        PhiHat(point p0i,point p1i,point p2i,circle ci) :
            p0(p0i),p1(p1i),p2(p2i),pc(ci),P(PhiPolygon({p0i,p1i,p2i})) {}

        PhiHat(point p0i,point p1i,point p2i,double r) :
            p0(p0i),p1(p1i),p2(p2i),P(PhiPolygon({p0i,p1i,p2i})) {
            double dx = p0.x-p1.x,dy = p0.y-p1.y;
            double s = sqrt(4*r*r/(dx*dx+dy*dy)-1);
            pc = circle(point((p0.x+p1.x+dy*s)/2,(p0.y+p1.y-dx*s)/2),r);
        }

        PhiFunc* phiFunc(RotTrans f,PhiCompObj* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f,PhiPolygon* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFunc(RotTrans f,PhiCircSeg* O, RotTrans g){
            return O->phiFunc(g, this, f);
        }

        PhiFunc* phiFuncG(RotTrans f, PhiPolygon* Q, RotTrans g){
            int m = Q->n;
            vector<PhiFunc*> comp(m+2);    //P1,P2 not in Q, if q over l then q in C
            for(int i=0;i<m;i++)
                comp[i] = new PhiFuncNode(false,{
                    PhiFunc::phiFunc(p0,p1,f,Q->p[i],g),
                    PhiFunc::phiFunc(pc,f,Q->p[i],g)});
            comp[m  ] = Q->phiFunc(g,p0,f);
            comp[m+1] = Q->phiFunc(g,p1,f);
            return new PhiFuncNode(true,comp);
        }

        PhiFunc* phiFuncG(RotTrans f, circle c, RotTrans g){
            return new PhiFuncNode(false,{
                PhiFunc::phiFunc(pc,f,c,g),
                PhiFunc::phiFunc(p0,p1,f,c,g),
                new PhiFuncNode(true,{
            PhiFunc::phiFunc(get<0>(moveLine(p0,p1,c.r)),get<0>(moveLine(p0,p2,c.r)),f,c.p,g),
            PhiFunc::phiFunc(get<1>(moveLine(p2,p1,c.r)),get<1>(moveLine(p0,p1,c.r)),f,c.p,g),
            PhiFunc::phiFunc(c,g,p0,f),PhiFunc::phiFunc(c,g,p1,f)})});
        }

        PhiFunc* phiFuncG(RotTrans f, PhiCircSeg* D, RotTrans g){
            return new PhiFuncNode(false,{
                D->P.phiFunc(g,p0,p1,f),
                D->phiFunc(g,new PhiCircCompl(pc),f),
                phiFuncG(f,D->pc,g),
                new PhiFuncNode(true,{
                    D->phiFunc(g,new PhiPolygon({p0,p2,point(2*p0.x-p1.x,2*p0.y-p1.y)}),f),
                    PhiFunc::phiFunc(pc,f,D->p1,g),
                    PhiFunc::phiFunc(D->p0,D->p1,g,p1,f),
                    PhiFunc::phiFunc(D->p1,D->p0,g,p0,f)}),
                new PhiFuncNode(true,{
                    D->phiFunc(g,new PhiPolygon({p1,point(2*p1.x-p0.x,2*p1.y-p0.y),p2}),f),
                    PhiFunc::phiFunc(pc,f,D->p0,g),
                    PhiFunc::phiFunc(D->p0,D->p1,g,p0,f),
                    PhiFunc::phiFunc(D->p1,D->p0,g,p1,f)})});
        }

        PhiFunc* phiFunc(RotTrans f, PhiHat* H, RotTrans g){
            return new PhiFuncNode(false,{
                P.phiFunc(f,H,g),phiFuncG(f,&H->P,g),
                new PhiFuncNode(true,{
                        PhiFunc::phiFunc(p2,p0,f,H->p1,g),
                        PhiFunc::phiFunc(H->p2,H->p0,g,p1,f),
                        PhiFunc::phiFunc(pc,f,H->p0,g),
                        PhiFunc::phiFunc(H->pc,g,p0,f)}),
                new PhiFuncNode(true,{
                        PhiFunc::phiFunc(p1,p2,f,H->p0,g),
                        PhiFunc::phiFunc(H->p1,H->p2,g,p0,f),
                        PhiFunc::phiFunc(pc,f,H->p1,g),
                        PhiFunc::phiFunc(H->pc,g,p1,f)})});
        }

        PhiFunc* phiFunc(RotTrans g, PhiCircCompl* C, Scale f){
            return P.phiFunc(g,C,f);
        }

        PhiFunc* phiFunc(RotTrans g, PhiLineCompl* L, Scale f){
            return P.phiFunc(g,L,f);
        }

        void move(point pd){
            p0.move(pd);p1.move(pd);p2.move(pd);pc.move(pd);
            P.move(pd);
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
