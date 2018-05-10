class Transform{
//    public:
//        virtual PhiFunc* phiFunc(PhiCompObj* P, PhiCompObj* Q, Transform* f) = 0;
};

class Scale: public Transform{
  public:
    int ind;
    Scale(int i) : ind(i) {}
};

//class Translate: public Transform{
//  public:
//    int ind;
//    Translate(int i) : ind(i) {}
//};

class RotTrans: public Transform{
  public:
    int ind;
    RotTrans(){};
    RotTrans(int i) : ind(i) {}
};

//PhiFunc* phiFunc()
