#ifndef DetectorMessenger_hh
#define DetectorMessenger_hh 1

#include <G4UImessenger.hh>
#include <G4UIcommand.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithAString.hh>
#include <G4String.hh>

class DetectorConstruction;

/**
 *  /shield/ ve altındaki tüm geometri/GDML komutlarını yöneten messenger
 *
 *  - /shield/geom/setSample    (int 0..9)
 *  - /shield/geom/setThickness (double length)
 *  - /shield/geom/setSizeXY    (halfX unitX halfY unitY)
 *  - /shield/geom/setShape     (box|cylinder|sphere)
 *  - /shield/gdml/setFile      (string path)
 */
class DetectorMessenger : public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction* det);
  virtual ~DetectorMessenger() override;

  virtual void SetNewValue(G4UIcommand* command,
                           G4String newValue) override;

private:
  DetectorConstruction*   fDetector;

  G4UIdirectory*          fShieldDir;
  G4UIdirectory*          fGeomDir;
  G4UIdirectory*          fGDMLDir;

  G4UIcmdWithAnInteger*       fSampleCmd;
  G4UIcmdWithADoubleAndUnit*  fThicknessCmd;
  G4UIcommand*                fSizeXYCmd;
  G4UIcmdWithAString*         fShapeCmd;
  G4UIcmdWithAString*         fGDMLFileCmd;
};

#endif
