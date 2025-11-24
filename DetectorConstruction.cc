#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4GeometryManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4GDMLParser.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4RunManager.hh>
#include <G4Exception.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>

#include <algorithm>  // std::max

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fShieldThickness(10.0 * mm),
    fHalfX(50.0 * mm),
    fHalfY(50.0 * mm),
    fSampleIndex(0),
    fShape("box"),
    fGDMLFile(""),
    fWorldPhys(nullptr),
    fMessenger(nullptr)
{
  // UI komutlarını sağlayan messenger
  fMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger;
}

//---------------------------------------------------------
//  Şekil set etme (box / cylinder / sphere)
//---------------------------------------------------------
void DetectorConstruction::SetShape(const G4String& name)
{
  G4String lower = name;
  lower.toLower();

  if (lower == "box" || lower == "slab" || lower == "rectangle") {
    fShape = "box";
  }
  else if (lower == "cylinder" || lower == "tube" || lower == "tubs") {
    fShape = "cylinder";
  }
  else if (lower == "sphere" || lower == "ball") {
    fShape = "sphere";
  }
  else {
    G4cout << "[DetectorConstruction] Unknown shape '" << name
           << "'. Using 'box' as default." << G4endl;
    fShape = "box";
  }
}

//---------------------------------------------------------
//  İç slab geometrisi (fallback + default)
//---------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::ConstructInternalGeometry()
{
  // Eski geometriyi temizle (ReinitializeGeometry için güvenli)
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  auto nist = G4NistManager::Instance();

  // --- WORLD ---
  G4double worldSize = 2.0 * m;
  auto worldSolid = new G4Box("WorldSolid",
                              0.5 * worldSize,
                              0.5 * worldSize,
                              0.5 * worldSize);

  auto worldMat = nist->FindOrBuildMaterial("G4_AIR");
  auto worldLV  = new G4LogicalVolume(worldSolid, worldMat, "WorldLV");

  fWorldPhys = new G4PVPlacement(
      nullptr,
      G4ThreeVector(),
      worldLV,
      "WorldPV",
      nullptr,
      false,
      0,
      true);

  auto worldVis = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9, 0.02));
  worldVis->SetForceWireframe(true);
  worldLV->SetVisAttributes(worldVis);

  // --- ABSORBER / SHIELD ---

  // Parametreleri güvence altına al
  G4double hx = (fHalfX > 0.) ? fHalfX : 50.0 * mm;
  G4double hy = (fHalfY > 0.) ? fHalfY : 50.0 * mm;
  G4double hz = (fShieldThickness > 0.) ? 0.5 * fShieldThickness : 5.0 * mm;

  // Malzeme: şimdilik örnek olarak Al; sample index'e göre
  // ileride Al5083, Al+WC, Al+B4C karışımlarını burada tanımlarsın.
  G4Material* shieldMat = nist->FindOrBuildMaterial("G4_Al");

  // Şekle göre solid seç
  G4VSolid* shieldSolid = nullptr;

  if (fShape == "sphere") {
    G4double r = std::max(std::max(hx, hy), hz);
    shieldSolid = new G4Sphere("ShieldSolid",
                               0.0 * mm, r,
                               0.0 * deg, 360.0 * deg,
                               0.0 * deg, 180.0 * deg);
  }
  else if (fShape == "cylinder") {
    G4double radius = std::max(hx, hy);
    shieldSolid = new G4Tubs("ShieldSolid",
                             0.0 * mm, radius,
                             hz,
                             0.0 * deg, 360.0 * deg);
  }
  else {
    // box / slab
    shieldSolid = new G4Box("ShieldSolid", hx, hy, hz);
  }

  auto shieldLV = new G4LogicalVolume(shieldSolid, shieldMat, "ShieldLV");

  new G4PVPlacement(
      nullptr,
      G4ThreeVector(0, 0, 0),
      shieldLV,
      "ShieldPV",
      worldLV,
      false,
      0,
      true);

  auto visShield = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.4));
  visShield->SetForceSolid(true);
  shieldLV->SetVisAttributes(visShield);

  G4cout
      << "[DetectorConstruction] Internal geometry built:\n"
      << "  shape     = " << fShape << "\n"
      << "  thickness = " << fShieldThickness / mm << " mm\n"
      << "  halfX     = " << hx / mm << " mm\n"
      << "  halfY     = " << hy / mm << " mm\n"
      << "  sampleIdx = " << fSampleIndex << G4endl;

  return fWorldPhys;
}

//---------------------------------------------------------
//  Construct: GDML varsa önce onu, yoksa slab
//---------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // GDML yolu boş değilse önce GDML world'ünü dene
  if (!fGDMLFile.empty()) {
    G4GDMLParser parser;

    G4cout << "[DetectorConstruction] Reading GDML: "
           << fGDMLFile << G4endl;

    parser.Read(fGDMLFile, false);

    // <setup> ismine dokunmadan, default setup'tan world al
    G4VPhysicalVolume* world = parser.GetWorldVolume();

    // Geriye dönük uyum için "Default" ismini de dene
    if (!world) {
      G4cout << "[DetectorConstruction] parser.GetWorldVolume() returned null, "
             << "trying with name \"Default\" ..." << G4endl;
      world = parser.GetWorldVolume("Default");
    }

    if (!world) {
      G4ExceptionDescription ed;
      ed << "GDML file " << fGDMLFile
         << " does not define a valid world (no <setup>/<world>).\n"
         << "Falling back to internal slab geometry.";
      G4Exception("DetectorConstruction::Construct",
                  "GDML_NO_WORLD", JustWarning, ed);

      // fallback: iç slab
      return ConstructInternalGeometry();
    }

    fWorldPhys = world;

    G4cout << "[DetectorConstruction] Using world from GDML: "
           << fGDMLFile << G4endl;

    return fWorldPhys;
  }

  // GDML yoksa direkt slab
  return ConstructInternalGeometry();
}

//---------------------------------------------------------
//  Örnek sample ismi (UI veya overlay label için)
//---------------------------------------------------------
std::string DetectorConstruction::GetSampleName() const
{
  switch (fSampleIndex) {
    case 0: return "Al5083";
    case 1: return "Al5083 + 5% WC";
    case 2: return "Al5083 + 10% WC";
    case 3: return "Al5083 + 15% WC";
    case 4: return "Al5083 + 5% B4C";
    case 5: return "Al5083 + 10% B4C";
    case 6: return "Al5083 + 15% B4C";
    case 7: return "Al5083 + 5% Hybrid";
    case 8: return "Al5083 + 10% Hybrid";
    case 9: return "Al5083 + 15% Hybrid";
    default: return "UnknownSample";
  }
}
