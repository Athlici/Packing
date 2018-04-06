class PhiFuncPrim;

class PhiFunc{
    public:
        virtual double eval(const double* x) = 0;
        virtual vector<PhiFuncPrim*> getIneqs(const double* x) = 0;
        virtual vector<int> getIndices(const double* x) = 0;
};

class PhiFuncPrim: public PhiFunc{
    public:
        virtual void getD1(const double* x,double* res) = 0;
        virtual vector<int> getD1ind() = 0;
        virtual void getD2(const double* x,double* res) = 0;
        virtual vector<tuple<int,int>> getD2ind() = 0;
        vector<PhiFuncPrim*> getIneqs (const double* x){return {this};}
        vector<int> getIndices (const double* x){return {};}
};

class PhiFuncNode: public PhiFunc{
    private:
        bool sense; //0 -> max,1 -> min
        vector<PhiFunc*> nodes;

        int argmax(const double* x){
            int ind = 0;
            double val = nodes[0]->eval(x);
            for(int i=1;i<nodes.size();i++){
                double tmp = nodes[i]->eval(x);
                if(tmp>val){
                    val = tmp;
                    ind = i;
                }
            }
            return ind;
        }

    public:
        PhiFuncNode(bool s,vector<PhiFunc*> n) : sense(s),nodes(n) {};

        double eval(const double* x){
            double res = nodes[0]->eval(x);
            for(int i=1;i<nodes.size();i++){
                double tmp = nodes[i]->eval(x);
                if(sense ^ (tmp>res))
                    res=tmp;
            }
            return res;
        }

        vector<PhiFuncPrim*> getIneqs(const double* x){
            if(sense){
                vector<PhiFuncPrim*> res = {};
                for(int i=0;i<nodes.size();i++){
                    vector<PhiFuncPrim*> tmp = nodes[i]->getIneqs(x);
                    res.insert(res.end(),tmp.begin(),tmp.end());
                }
                return res;
            } else
                return nodes[argmax(x)]->getIneqs(x);
        }

        vector<int> getIndices(const double* x){
            if(sense){
                vector<int> res = {};
                for(int i=0;i<nodes.size();i++){
                    vector<int> tmp = nodes[i]->getIndices(x);
                    res.insert(res.end(),tmp.begin(),tmp.end());
                }
                return res;
            } else {
                int ind = argmax(x);
                vector<int> res = nodes[ind]->getIndices(x);
                res.push_back(ind);
                return res;
            }
        }
};

class PhiFuncScCcTrPt : public PhiFuncPrim{ //TODO: scale cc.p
    double dx,dy,c,c2;
    int i,j;
  public:
    PhiFuncScCcTrPt(circle cc,Scale f,point p,Translate g){
        dx = cc.p.x-p.x;
        dy = cc.p.y-p.y;
        c = cc.r;
        c2= 2*c*c;
        i = f.ind;
        j = g.ind;
    }

    double eval(const double* x){
        double r=c*x[i],a=x[j]-dx,b=x[j+1]-dy;
        return r*r-a*a-b*b;
    }

    vector<int> getD1ind(){
        return {i,j,j+1};
    }

    void getD1(const double* x,double* res){
        res[0] = c2*x[i];
        res[1] = 2*(dx-x[j]);
        res[2] = 2*(dy-x[j+1]);
    }

    vector<tuple<int,int>> getD2ind(){
        return {{i,i},{j,j},{j+1,j+1}};
    }

    void getD2(const double* x,double* res){
        res[0] = c2; res[1] = -2; res[2] = -2;
    }
};

class PhiFuncLnRTPtRT : public PhiFuncPrim{
  public:
    PhiFuncLnRTPtRT(point l1, point l2, RotTrans f, point p, RotTrans g){
    };
    
    double eval(const double* x){
        return 0;
    }

    vector<int> getD1ind(){
        return {};
    }

    void getD1(const double* x,double* res){}

    vector<tuple<int,int>> getD2ind(){
        return {};
    }

    void getD2(const double* x,double* res){}
};

PhiFunc* phiFunc(PhiCircCompl C,Scale f,PhiPolygon P,Translate g){
    vector<PhiFunc*> comp(P.p.size());
    for(int i=0;i<P.p.size();i++)
        comp[i] = new PhiFuncScCcTrPt(C.c,f,P.p[i],g);
    return new PhiFuncNode(true,comp);
}

PhiFunc* phiFunc(PhiPolygon P,RotTrans f,PhiPolygon Q,RotTrans g){
    int n = P.p.size(),m = Q.p.size();
    vector<PhiFunc*> comp(n+m);
    for(int i=0;i<n;i++){
        vector<PhiFunc*> tmp(m);
        for(int j=0;j<m;j++)
            tmp[j] = new PhiFuncLnRTPtRT(P.p[i],P.p[(i+1)%n],f,Q.p[j],g);
        comp[i] = new PhiFuncNode(true,tmp);
    }
    for(int i=0;i<m;i++){
        vector<PhiFunc*> tmp(n);
        for(int j=0;j<n;j++)
            tmp[j] = new PhiFuncLnRTPtRT(Q.p[i],Q.p[(i+1)%m],g,P.p[j],f);
        comp[n+i] = new PhiFuncNode(true,tmp);
    }
    return new PhiFuncNode(false,comp);
}
