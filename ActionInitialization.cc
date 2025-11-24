//====================== src/ActionInitialization.cc ======================
#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* det)
  : G4VUserActionInitialization(),
    fDetector(det)
{
}

ActionInitialization::~ActionInitialization() = default;

void ActionInitialization::BuildForMaster() const
{
  // Sadece run-level istatistikler vs.
  auto runAction = new RunAction(fDetector);
  SetUserAction(runAction);
}

void ActionInitialization::Build() const
{
  // 1) Primary kaynak
  auto primGen = new PrimaryGeneratorAction();
  SetUserAction(primGen);

  // 2) RunAction
  auto runAction = new RunAction(fDetector);
  SetUserAction(runAction);

  // 3) SteppingAction – SADECE RunAction pointer’ı alıyor
  auto stepping = new SteppingAction(runAction);
  SetUserAction(stepping);
}
