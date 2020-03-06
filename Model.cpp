//class to store the complete problem
class Model{
  public:
    vector<var> vars;           //variables
    Objective* f;               //objective funtion
    vector<PhiInfObj*> C;       //container object (union of primitive objects)
    Scale sc = 0;               //container modification function TODO: make more general
    vector<PhiCompObj*> objs;   //list of objects that shall be fitted
    vector<RotTrans> rt;        //list of their transformations

    Model() {};

    //generate the no-fit constraints
    PhiFunc* createPhiFunc(){
        int n = objs.size(),m = C.size();
        vector<PhiFunc*> comp(n*(n-1)/2+n*m);
        //must lie inside the container
        for(int i=0;i<n;i++)
            for(int j=0;j<m;j++)
                comp[i*m+j] = phiFunc(C[j],sc,objs[i],rt[i]);
        int ind = n*m;
        //must not lie inside each other
        for(int i=0;i<n;i++)
            for(int j=i+1;j<n;j++)
                comp[ind++] = phiFunc(objs[i],rt[i],objs[j],rt[j]);
        return new PhiFuncNode(true,comp);
    }
};
