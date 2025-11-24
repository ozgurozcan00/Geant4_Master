//====================== src/RunAction.cc ======================
#include "RunAction.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

RunAction::RunAction(DetectorConstruction* det)
  : G4UserRunAction(),
    fDetector(det),
    fDepthMin_mm(0.0),
    fDepthMax_mm(100.0),
    fDepthNBins(100),
    fDepthBinWidth_mm(0.0),
    fMu_cm(0.0),
    fHvl_mm(0.0),
    fTransFrac(0.0),
    fTotalDose_Gy(0.0),
    fDoseRMS_Gy(0.0),
    fMassDensity_gcm3(1.0)   // basit varsayılan
{
  SetDepthRange(fDepthMin_mm, fDepthMax_mm, fDepthNBins);
}

void RunAction::SetDepthRange(G4double zmin_mm, G4double zmax_mm, G4int nbins)
{
  if (zmax_mm <= zmin_mm) {
    zmin_mm = 0.0;
    zmax_mm = 100.0;
  }
  if (nbins <= 0) nbins = 100;

  fDepthMin_mm      = zmin_mm;
  fDepthMax_mm      = zmax_mm;
  fDepthNBins       = nbins;
  fDepthBinWidth_mm = (fDepthMax_mm - fDepthMin_mm) / (G4double)fDepthNBins;

  fDepthPos_mm.resize(fDepthNBins);
  fDepthDose_Gy.assign(fDepthNBins, 0.0);
  fDepthDoseErr_Gy.assign(fDepthNBins, 0.0);

  for (G4int i = 0; i < fDepthNBins; ++i) {
    fDepthPos_mm[i] = fDepthMin_mm + (i + 0.5) * fDepthBinWidth_mm;
  }
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  // Her run başında histogramı sıfırla
  fDepthDose_Gy.assign(fDepthNBins, 0.0);
  fDepthDoseErr_Gy.assign(fDepthNBins, 0.0);

  fMu_cm        = 0.0;
  fHvl_mm       = 0.0;
  fTransFrac    = 0.0;
  fTotalDose_Gy = 0.0;
  fDoseRMS_Gy   = 0.0;

  // Yoğunluk bilgisini detector'dan almayı deneyebilirsin,
  // şimdilik yoksa 1 g/cm3 kalsın.
  // Ör: fMassDensity_gcm3 = det->GetMassDensity_gcm3();
}

void RunAction::AddDepthDose(G4double depth_mm, G4double edep_MeV)
{
  if (fDepthNBins <= 0) return;
  if (depth_mm < fDepthMin_mm || depth_mm >= fDepthMax_mm) return;

  G4int bin = (G4int)((depth_mm - fDepthMin_mm) / fDepthBinWidth_mm);
  if (bin < 0) bin = 0;
  if (bin >= fDepthNBins) bin = fDepthNBins - 1;

  // Basit doz hesabı: 1 cm^2 alan varsay, her bin için
  // hacim = alan * dz
  G4double dz_cm   = (fDepthBinWidth_mm / 10.0);     // mm -> cm
  G4double volume  = 1.0 * dz_cm;                    // cm^3
  G4double mass_kg = (fMassDensity_gcm3 * volume) / 1000.0; // g -> kg

  if (mass_kg <= 0.) mass_kg = 1.e-6;

  const G4double MeV_to_J = 1.602176634e-13;         // J/MeV
  G4double edep_J  = edep_MeV * MeV_to_J;
  G4double dose_Gy = edep_J / mass_kg;               // Gy = J/kg

  fDepthDose_Gy[bin] += dose_Gy;
}

bool RunAction::HasDepthDoseData() const
{
  for (auto v : fDepthDose_Gy) {
    if (v > 0.) return true;
  }
  return false;
}

void RunAction::ComputeAttenuationParameters()
{
  if (!HasDepthDoseData()) {
    fMu_cm        = 0.0;
    fHvl_mm       = 0.0;
    fTransFrac    = 0.0;
    fTotalDose_Gy = 0.0;
    fDoseRMS_Gy   = 0.0;
    return;
  }

  // Toplam doz ve RMS
  G4double sum = 0.0;
  G4double sum2 = 0.0;
  for (auto d : fDepthDose_Gy) {
    sum  += d;
    sum2 += d * d;
  }
  fTotalDose_Gy = sum;
  G4double mean = sum / (G4double)fDepthNBins;
  G4double var  = sum2 / (G4double)fDepthNBins - mean * mean;
  if (var < 0.) var = 0.;
  fDoseRMS_Gy = std::sqrt(var);

  // Basit attenuation fit'i: I(z) ~ I0 * exp(-mu z)
  // I0 olarak ilk bin'i al.
  G4double I0 = fDepthDose_Gy.front();
  if (I0 <= 0.) {
    fMu_cm     = 0.0;
    fHvl_mm    = 0.0;
    fTransFrac = 0.0;
    return;
  }

  // HVL: I(z) = I0/2 olduğu nokta
  G4double target = 0.5 * I0;
  G4double hvl_mm = 0.0;
  bool found = false;
  for (G4int i = 1; i < fDepthNBins; ++i) {
    G4double z0 = fDepthPos_mm[i-1];
    G4double z1 = fDepthPos_mm[i];
    G4double y0 = fDepthDose_Gy[i-1];
    G4double y1 = fDepthDose_Gy[i];

    if ((y0 >= target && y1 <= target) ||
        (y0 <= target && y1 >= target)) {
      // lineer interpolasyon
      G4double t = (target - y0) / (y1 - y0);
      if (t < 0.) t = 0.;
      if (t > 1.) t = 1.;
      hvl_mm = z0 + t * (z1 - z0);
      found = true;
      break;
    }
  }

  if (found && hvl_mm > 0.) {
    G4double hvl_cm = hvl_mm / 10.0;
    fMu_cm  = std::log(2.0) / hvl_cm;
    fHvl_mm = hvl_mm;
  } else {
    fMu_cm  = 0.0;
    fHvl_mm = 0.0;
  }

  // Top kalınlık için transmisyon
  G4double total_thickness_cm =
      (fDepthMax_mm - fDepthMin_mm) / 10.0;
  if (fMu_cm > 0.0 && total_thickness_cm > 0.0)
    fTransFrac = std::exp(-fMu_cm * total_thickness_cm);
  else
    fTransFrac = 0.0;

  // Çok kaba bir hata tahmini: sqrt(N) tipi
  for (G4int i = 0; i < fDepthNBins; ++i) {
    G4double d = fDepthDose_Gy[i];
    if (d > 0.)
      fDepthDoseErr_Gy[i] = std::sqrt(d); // boyutsal olarak tutarlı değil ama grafikte iş görür
    else
      fDepthDoseErr_Gy[i] = 0.0;
  }
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  ComputeAttenuationParameters();

  G4cout << "==== Run summary (RunAction) ====" << G4endl;
  G4cout << "  Depth bins : " << fDepthNBins
         << "  [" << fDepthMin_mm << " , " << fDepthMax_mm << "] mm"
         << G4endl;
  G4cout << "  Total dose ≈ " << fTotalDose_Gy / gray << " Gy" << G4endl;
  if (fMu_cm > 0.) {
    G4cout << "  mu        : " << fMu_cm   << " 1/cm" << G4endl;
    G4cout << "  HVL       : " << fHvl_mm  << " mm"   << G4endl;
    G4cout << "  TransFrac : " << fTransFrac << G4endl;
  }
}
