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

class PhiFuncScCcRtCl : public PhiFuncPrim{ //Only unit-circle around (0,0)
    double px,py,rc;
    double d1,d2;
    int i,j;
  public:
    PhiFuncScCcRtCl(circle cc,Scale f,circle c,RotTrans g){
        px = c.p.x; py = c.p.y; rc = c.r;
        i = f.ind; j = g.ind;
    }

    PhiFuncScCcRtCl(circle cc,Scale f,point p,RotTrans g){
        //PhiFuncScCcRtCl(cc,f,circle(p,0),g);
        px = p.x; py = p.y; rc = 0;
        i = f.ind; j = g.ind;
    }

    double eval(const double* x){
        double s = sin(x[j+2]), c = cos(x[j+2]);
        double dr = x[i]  -rc;
        double dx = x[j]  +px*c+py*s;
        double dy = x[j+1]-px*s+py*c;
        return dr*dr-dx*dx-dy*dy;
    }

    vector<int> getD1ind(){
        return {i,j,j+1,j+2};
    }

    void getD1(const double* x,double* res){
        double s = sin(x[j+2]), c = cos(x[j+2]);
        res[0] = 2*(x[i]-rc);
        res[1] =-2*(x[j]  +px*c+py*s);
        res[2] =-2*(x[j+1]-px*s+py*c);
        res[3] = 2*(x[j]*(px*s-py*c)+x[j+1]*(px*c+py*s));
    }

    vector<tuple<int,int>> getD2ind(){
        return {{i,i},{j,j},{j+1,j+1},{j,j+2},{j+1,j+2},{j+2,j+2}};
    }

    void getD2(const double* x,double* res){
        double s = sin(x[j+2]), c = cos(x[j+2]);
        res[0] = 2; res[1] = -2; res[2] = -2;
        res[3] = 2*(px*s-py*c);
        res[4] = 2*(px*c+py*s);
        res[5] = 2*(x[j]*(px*c+py*s)+x[j+1]*(py*c-px*s));
    }
};

class PhiFuncHScCcRtCs : public PhiFuncPrim{
    double dx,dy,det;
    int i;
  public:
    PhiFuncHScCcRtCs(RotTrans f, circle c, point p, double s){
      dx = s*(c.p.y-p.y);
      dy = s*(p.x-c.p.x);
      det= s*(p.x*c.p.y-p.y*c.p.x);
      i = f.ind;
    };

    double eval(const double* x){
        double s = sin(x[i+2]), c = cos(x[i+2]);
        return x[i]*(dx*c+dy*s)+x[i+1]*(-dx*s+dy*c)+det;
    }

    vector<int> getD1ind(){
        return {i,i+1,i+2};
    }

    void getD1(const double* x,double* res){
        double s = sin(x[i+2]), c = cos(x[i+2]);
        res[0] =  dx*c+dy*s;
        res[1] = -dx*s+dy*c;
        res[2] = x[i]*(dy*c-dx*s)-x[i+1]*(dx*c+dy*s);
    }

    vector<tuple<int,int>> getD2ind(){
        return {{i,i+2},{i+1,i+2},{i+2,i+2}};
    }

    void getD2(const double* x,double* res){
        double s = sin(x[i+2]), c = cos(x[i+2]);
        res[0] = -dx*s+dy*c;
        res[1] = -dx*c-dy*s;
        res[2] = -x[i]*(dy*s+dx*c)+x[i+1]*(dx*s-dy*c);
    }
};

class PhiFuncLnRtClRt : public PhiFuncPrim{
    Matrix<double,2,2> A,B,T;
    double c;
    int i,j;
  public:
    PhiFuncLnRtClRt(point l1, point l2, RotTrans f, point p, RotTrans g){
        double dx=l1.x-l2.x,dy=l1.y-l2.y,n=1/sqrt(dx*dx+dy*dy);
        dx *= n; dy *= n;
        A << dx*p.y-dy*p.x, -dy*p.y-dx*p.x, dx*p.x+dy*p.y, dx*p.y-dy*p.x;
        B << dx, -dy, dy, dx;
        c = n*(l2.x*l1.y-l1.x*l2.y);
        i = f.ind; j = g.ind;
    };
    
    double eval(const double* x){
        Vector2d    f1(sin(x[i+2]),cos(x[i+2]));
        RowVector2d f2(sin(x[j+2]),cos(x[j+2])),y(x[j]-x[i],x[j+1]-x[i+1]);
        return (f2*A+y*B)*f1+c;
    }

    vector<int> getD1ind(){
        return {i,i+1,i+2,j,j+1,j+2};
    }

    void getD1(const double* x,double* res){
        Vector2d    f1(sin(x[i+2]),cos(x[i+2])),f1d(f1[1],-f1[0]);
        RowVector2d f2(sin(x[j+2]),cos(x[j+2])),f2d(f2[1],-f2[0]);
        RowVector2d y(x[j]-x[i],x[j+1]-x[i+1]);
        Vector2d    Bf1 = B*f1;
        double tmp[] = {-Bf1[0],-Bf1[1],(f2*A+y*B)*f1d,Bf1[0],Bf1[1],f2d*A*f1};
        for(int i=0;i<6;i++) res[i] = tmp[i];
    }

    vector<tuple<int,int>> getD2ind(){
        if(i<j)
            return {{i,i+2},{i+1,i+2},{i+2,i+2},{j+2,j+2},{i+2,j},{i+2,j+1},{i+2,j+2}};
        else
            return {{i,i+2},{i+1,i+2},{i+2,i+2},{j+2,j+2},{j,i+2},{j+1,i+2},{j+2,i+2}};
    }

    void getD2(const double* x,double* res){
        Vector2d    f1(sin(x[i+2]),cos(x[i+2])),f1d(f1[1],-f1[0]);
        RowVector2d f2(sin(x[j+2]),cos(x[j+2])),f2d(f2[1],-f2[0]);
        RowVector2d y(x[j]-x[i],x[j+1]-x[i+1]);
        Vector2d    Bf1d = B*f1d;
        RowVector2d f2A=f2*A;
        double tmp[] {-Bf1d[0],-Bf1d[1],-f2A*f1,-(f2A+y*B)*f1,Bf1d[0],Bf1d[1],f2d*A*f1d};
        for(int i=0;i<7;i++) res[i] = tmp[i];
    }
};

class PhiFuncClRtClRt : public PhiFuncPrim{
  public:
    PhiFuncClRtClRt(point l1, point l2, RotTrans f, point p, RotTrans g){
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
