//====================== include/ActionInitialization.hh ======================
#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;
class RunAction;
class SteppingAction;

class ActionInitialization : public G4VUserActionInitialization
{
public:
  explicit ActionInitialization(DetectorConstruction* det);
  ~ActionInitialization() override;

  void BuildForMaster() const override;
  void Build() const override;

private:
  DetectorConstruction* fDetector; // geometriye erişmek için
};

#endif
