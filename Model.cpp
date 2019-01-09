class Model{
  public:
    vector<var> vars;
    Objective* f;
    vector<PhiInfObj*> C;
    Scale sc = 0;
    vector<PhiCompObj*> objs;
    vector<RotTrans> rt;

    Model() {};

    PhiFunc* createPhiFunc(){
        int n = objs.size(),m = C.size();
        vector<PhiFunc*> comp(n*(n-1)/2+n*m);
        for(int i=0;i<n;i++)
            for(int j=0;j<m;j++)
                comp[i*m+j] = phiFunc(C[j],sc,objs[i],rt[i]);
        int ind = n*m;
        for(int i=0;i<n;i++)
            for(int j=i+1;j<n;j++)
                comp[ind++] = phiFunc(objs[i],rt[i],objs[j],rt[j]);
        return new PhiFuncNode(true,comp);
    }
};
