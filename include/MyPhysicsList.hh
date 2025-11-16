
#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class MyPhysicsList : public G4VModularPhysicsList
{
public:
    MyPhysicsList();        // 构造函数
    ~MyPhysicsList() override = default; // 析构函数

protected:
    void SetCuts() override;  // 粒子 cut
};

#endif
