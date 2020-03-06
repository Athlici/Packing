class Transform{
//    public:
//        virtual PhiFunc* phiFunc(PhiCompObj* P, PhiCompObj* Q, Transform* f) = 0;
};

//scaling function, only for the container
class Scale: public Transform{
  public:
    int ind;
    Scale(){};
    Scale(int i) : ind(i) {}
};

//class Translate: public Transform{    <- imitated by RotTrans with third variable fixed 0
//  public:
//    int ind;
//    Translate(int i) : ind(i) {}
//};

//rotate and then translate
class RotTrans: public Transform{
  public:
    int ind;
    RotTrans(){};
    RotTrans(int i) : ind(i) {}
};
