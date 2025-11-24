#ifndef STEPPING_ACTION_HH
#define STEPPING_ACTION_HH

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;
class RunAction;

// Sadece RunAction pointer'ı tutuyoruz
class SteppingAction : public G4UserSteppingAction
{
public:
  explicit SteppingAction(RunAction* runAction);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step* step) override;

private:
  RunAction* fRunAction;  // hiçbir şekilde DetectorConstruction tutmuyoruz
};

#endif
