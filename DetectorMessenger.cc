#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include <G4UnitsTable.hh>
#include <G4Exception.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>

#include <sstream>

DetectorMessenger::DetectorMessenger(DetectorConstruction* det)
  : G4UImessenger(),
    fDetector(det),
    fShieldDir(nullptr),
    fGeomDir(nullptr),
    fGDMLDir(nullptr),
    fSampleCmd(nullptr),
    fThicknessCmd(nullptr),
    fSizeXYCmd(nullptr),
    fShapeCmd(nullptr),
    fGDMLFileCmd(nullptr)
{
  // Ana /shield/ directory
  fShieldDir = new G4UIdirectory("/shield/");
  fShieldDir->SetGuidance("Shielding geometry control.");

  // Geometri alt dizini
  fGeomDir = new G4UIdirectory("/shield/geom/");
  fGeomDir->SetGuidance("Geometry parameters for shielding.");

  // GDML alt dizini
  fGDMLDir = new G4UIdirectory("/shield/gdml/");
  fGDMLDir->SetGuidance("GDML import control.");

  // /shield/geom/setSample
  fSampleCmd = new G4UIcmdWithAnInteger("/shield/geom/setSample", this);
  fSampleCmd->SetGuidance("Set sample index (0..9).");
  fSampleCmd->SetParameterName("sampleIndex", false);
  fSampleCmd->SetRange("sampleIndex>=0 && sampleIndex<=9");

  // /shield/geom/setThickness
  fThicknessCmd =
      new G4UIcmdWithADoubleAndUnit("/shield/geom/setThickness", this);
  fThicknessCmd->SetGuidance("Set shield thickness (full thickness).");
  fThicknessCmd->SetParameterName("thickness", false);
  fThicknessCmd->SetUnitCategory("Length");
  fThicknessCmd->SetDefaultUnit("mm");

  // /shield/geom/setSizeXY  halfX unitX halfY unitY
  fSizeXYCmd = new G4UIcommand("/shield/geom/setSizeXY", this);
  fSizeXYCmd->SetGuidance("Set absorber half sizes in X and Y.");
  {
    auto* pX = new G4UIparameter("halfX", 'd', false);
    auto* uX = new G4UIparameter("unitX", 's', false);
    auto* pY = new G4UIparameter("halfY", 'd', false);
    auto* uY = new G4UIparameter("unitY", 's', false);
    fSizeXYCmd->SetParameter(pX);
    fSizeXYCmd->SetParameter(uX);
    fSizeXYCmd->SetParameter(pY);
    fSizeXYCmd->SetParameter(uY);
  }

  // /shield/geom/setShape
  fShapeCmd = new G4UIcmdWithAString("/shield/geom/setShape", this);
  fShapeCmd->SetGuidance("Set absorber shape: box, cylinder, sphere.");
  fShapeCmd->SetParameterName("shape", false);

  // /shield/gdml/setFile
  fGDMLFileCmd = new G4UIcmdWithAString("/shield/gdml/setFile", this);
  fGDMLFileCmd->SetGuidance("Set GDML geometry file path.");
  fGDMLFileCmd->SetParameterName("filename", false);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fSampleCmd;
  delete fThicknessCmd;
  delete fSizeXYCmd;
  delete fShapeCmd;
  delete fGDMLFileCmd;

  delete fGeomDir;
  delete fGDMLDir;
  delete fShieldDir;
}

//---------------------------------------------------------
//  Komutları DetectorConstruction'a yönlendir
//---------------------------------------------------------
void DetectorMessenger::SetNewValue(G4UIcommand* command,
                                    G4String newValue)
{
  if (command == fSampleCmd) {
    fDetector->SetSampleIndex(fSampleCmd->GetNewIntValue(newValue));
  }
  else if (command == fThicknessCmd) {
    G4double val = fThicknessCmd->GetNewDoubleValue(newValue);
    fDetector->SetShieldThickness(val);
  }
  else if (command == fSizeXYCmd) {
    // Beklenen format:
    //   /shield/geom/setSizeXY  25 mm  25 mm
    std::istringstream iss(newValue);
    G4double vX = 0., vY = 0.;
    G4String uX, uY;
    iss >> vX >> uX >> vY >> uY;
    if (!iss) {
      G4ExceptionDescription ed;
      ed << "Bad syntax for /shield/geom/setSizeXY.\n"
         << "Expected: halfX unitX halfY unitY (e.g. \"25 mm 25 mm\")";
      G4Exception("DetectorMessenger::SetNewValue",
                  "ShieldSizeXY", JustWarning, ed);
      return;
    }

    G4double valX = vX * G4UnitDefinition::GetValueOf(uX);
    G4double valY = vY * G4UnitDefinition::GetValueOf(uY);
    fDetector->SetSizeXY(valX, valY);
  }
  else if (command == fShapeCmd) {
    fDetector->SetShape(newValue);
  }
  else if (command == fGDMLFileCmd) {
    fDetector->SetGDMLFile(newValue);
    G4cout << "[DetectorMessenger] GDML file set to: "
           << newValue << G4endl;
  }
}
