//====================== include/RunAction.hh ======================
#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class G4Run;
class DetectorConstruction;

class RunAction : public G4UserRunAction
{
public:
  explicit RunAction(DetectorConstruction* det);
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  // SteppingAction'dan çağrılacak
  void AddDepthDose(G4double depth_mm, G4double edep_MeV);

  // İstersen makrodan derinlik gridini değiştirebilirsin
  void SetDepthRange(G4double zmin_mm, G4double zmax_mm, G4int nbins);

  // Veri var mı?
  bool HasDepthDoseData() const;

  // Skalar parametreler
  G4double GetMu_cm()       const { return fMu_cm; }
  G4double GetHvl_mm()      const { return fHvl_mm; }
  G4double GetTransFrac()   const { return fTransFrac; }
  G4double GetTotalDoseGy() const { return fTotalDose_Gy; }
  G4double GetDoseRMSGy()   const { return fDoseRMS_Gy; }

  // Derinlik–doz dağılımı
  const std::vector<G4double>& GetDepthPositions()  const { return fDepthPos_mm;    }
  const std::vector<G4double>& GetDepthDoseBins()   const { return fDepthDose_Gy;   }
  const std::vector<G4double>& GetDepthDoseErrors() const { return fDepthDoseErr_Gy;}
  G4double GetDepthBinWidth() const { return fDepthBinWidth_mm; }

private:
  void ComputeAttenuationParameters();

  DetectorConstruction*   fDetector;

  // Grid
  G4double                 fDepthMin_mm;
  G4double                 fDepthMax_mm;
  G4int                    fDepthNBins;
  G4double                 fDepthBinWidth_mm;

  std::vector<G4double>    fDepthPos_mm;
  std::vector<G4double>    fDepthDose_Gy;      // bin başına Gy
  std::vector<G4double>    fDepthDoseErr_Gy;   // basit RMS tahmini

  // Attenuation parametreleri
  G4double                 fMu_cm;
  G4double                 fHvl_mm;
  G4double                 fTransFrac;
  G4double                 fTotalDose_Gy;
  G4double                 fDoseRMS_Gy;

  G4double                 fMassDensity_gcm3;
};

#endif
