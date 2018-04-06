class Transform{};

class Scale: public Transform{
  public:
    int ind;
    Scale(int i) : ind(i) {}
};

class Translate: public Transform{
  public:
    int ind;
    Translate(int i) : ind(i) {}
};

class RotTrans: public Transform{
  public:
    int ind;
    RotTrans(int i) : ind(i) {}
};
