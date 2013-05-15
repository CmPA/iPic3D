#ifndef injIFields
#define injIFields

struct injInfoFields{
    double*** ExITemp;
    double*** EyITemp;
    double*** EzITemp;
    double*** BxITemp;
    double*** ByITemp;
    double*** BzITemp;
    int Nxsize_store;
    int Nysize_store;
    
    injInfoFields(int Nxsize,int Nysize, int Nzsize);
    ~injInfoFields();
};


struct injInfoParticles
{
    double*** VxITemp;
    double*** VyITemp;
    double*** VzITemp;
    
    double*** VthxITemp;
    double*** VthyITemp;
    double*** VthzITemp;
    
    int*** Npcelx_array;
    int*** Npcely_array;
    int*** Npcelz_array;
    
    double*** QITemp; // weight of injected particles: make it variable, keeping the amount of injected particles.
    
    double*** RoITemp;
    int Nxsize_store;
    int Nysize_store;
    
    injInfoParticles(int Nxsize,int Nysize, int Nzsize);
    ~injInfoParticles();
};
/*
inline injInfoFields::injInfoFields(int Nxsize,int Nysize, int Nzsize)
{        ExITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    EyITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    EzITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    BxITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    ByITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    BzITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    Nxsize_store=Nxsize;
    Nysize_store=Nysize;
};
 inline injInfoFields::~injInfoFields()
{
        delArr3(ExITemp,Nxsize_store,Nysize_store);
        delArr3(EyITemp,Nxsize_store,Nysize_store);
        delArr3(EzITemp,Nxsize_store,Nysize_store);
        delArr3(BxITemp,Nxsize_store,Nysize_store);
        delArr3(ByITemp,Nxsize_store,Nysize_store);
        delArr3(BzITemp,Nxsize_store,Nysize_store);
};

inline injInfoParticles::injInfoParticles(int Nxsize,int Nysize, int Nzsize)
{
    VxITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    VyITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    VzITemp=newArr3(double,Nxsize,Nysize,Nzsize);

    VthxITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    VthyITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    VthzITemp=newArr3(double,Nxsize,Nysize,Nzsize);

    Npcelx_array=newArr3(int,Nxsize,Nysize,Nzsize);
    Npcely_array=newArr3(int,Nxsize,Nysize,Nzsize);
    Npcelz_array=newArr3(int,Nxsize,Nysize,Nzsize);

    QITemp=newArr3(double,Nxsize,Nysize,Nzsize);
    RoITemp=newArr3(double,Nxsize,Nysize,Nzsize);

    Nxsize_store=Nxsize;
    Nysize_store=Nysize;
};

 inline injInfoParticles::~injInfoParticles()
{
        delArr3(VxITemp,Nxsize_store,Nysize_store);
        delArr3(VyITemp,Nxsize_store,Nysize_store);

        delArr3(VzITemp,Nxsize_store,Nysize_store);

        delArr3(VthxITemp,Nxsize_store,Nysize_store);
        delArr3(VthyITemp,Nxsize_store,Nysize_store);
        delArr3(VthzITemp,Nxsize_store,Nysize_store);

        delArr3(Npcelx_array,Nxsize_store,Nysize_store);
        delArr3(Npcely_array,Nxsize_store,Nysize_store);
        delArr3(Npcelz_array,Nxsize_store,Nysize_store);

        delArr3(QITemp,Nxsize_store,Nysize_store);
        delArr3(RoITemp,Nxsize_store,Nysize_store);
};
*/
#endif
