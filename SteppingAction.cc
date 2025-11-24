#include "SteppingAction.hh"
#include "RunAction.hh"

#include "G4Step.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(RunAction* runAction)
  : G4UserSteppingAction(),
    fRunAction(runAction)
{
}

SteppingAction::~SteppingAction() = default;

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fRunAction) {
    return;
  }

  // Bu adımda depolanan enerji
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.) {
    return;
  }

  // Derinlik hesabı için pre-step noktasının pozisyonu
  const auto* prePoint = step->GetPreStepPoint();
  const G4ThreeVector pos = prePoint->GetPosition();

  // Z ekseni boyunca ilerleyen bir demet varsayıyoruz.
  // Derinlik: z=0 giriş düzlemi, pozitif z derinlik.
  G4double depth_mm = pos.z() / mm;
  if (depth_mm < 0.) {
    // Absorber öncesi adımları sayma
    return;
  }

  // MeV cinsine çevir
  G4double edep_MeV = edep / MeV;

  // RunAction içindeki derinlik doz histogramına ekle
  fRunAction->AddDepthDose(depth_mm, edep_MeV);
}
