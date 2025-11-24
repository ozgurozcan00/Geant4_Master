#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh 1

#include <G4VUserDetectorConstruction.hh>
#include <G4SystemOfUnits.hh>
#include <G4String.hh>

#include <string>

class G4VPhysicalVolume;
class DetectorMessenger;

/**
 *  Slab tabanlı zırh geometrisi + opsiyonel GDML import
 *
 *  - Eğer SetGDMLFile() ile bir GDML yolu verilmişse:
 *      -> Construct() içinde önce GDML world okunur.
 *      -> GDML içinde geçerli bir <setup>/<world> yoksa
 *         iç slab geometrisine geri düşer.
 *  - GDML boş ise doğrudan iç slab geometrisi kurulur.
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction() override;

  // Ana Construct
  virtual G4VPhysicalVolume* Construct() override;

  // SD/field kullanmıyorsak boş bırakıyoruz
  virtual void ConstructSDandField() override {}

  // --- Parametre set/get (UI komutlarıyla uyumlu) ---

  // Kalınlık (tam kalınlık, mm biriminde saklanıyor)
  void     SetShieldThickness(G4double val) { fShieldThickness = val; }
  G4double GetShieldThickness() const       { return fShieldThickness; }

  // Numune ID (0..9)
  void  SetSampleIndex(G4int i) { fSampleIndex = i; }
  G4int GetSampleIndex() const  { return fSampleIndex; }

  // Slab yarım boyutları (X,Y)
  void SetSizeXY(G4double halfX, G4double halfY)
  {
    fHalfX = halfX;
    fHalfY = halfY;
  }

  G4double GetHalfX() const { return fHalfX; }
  G4double GetHalfY() const { return fHalfY; }

  // 2D scoring map için alias fonksiyonlar
  G4double GetScoringHalfX() const { return fHalfX; }
  G4double GetScoringHalfY() const { return fHalfY; }

  // Şekil seçimi (box / cylinder / sphere)
  void      SetShape(const G4String& name);
  G4String  GetShape() const { return fShape; }

  // GDML dosya yolu
  void                SetGDMLFile(const std::string& path) { fGDMLFile = path; }
  const std::string&  GetGDMLFile() const                  { return fGDMLFile; }

  // UI / overlay label için
  std::string GetSampleName() const;

private:
  // İç slab geometrisini kuran fonksiyon
  G4VPhysicalVolume* ConstructInternalGeometry();

private:
  // Parametreler (mm cinsinden)
  G4double fShieldThickness;   // Tam kalınlık (mm)
  G4double fHalfX;             // Absorber yarım boy X (mm)
  G4double fHalfY;             // Absorber yarım boy Y (mm)
  G4int    fSampleIndex;       // 0..9

  // Geometri şekli: "box", "cylinder", "sphere"
  G4String fShape;

  // GDML yolu (boş ise kullanılmaz)
  std::string fGDMLFile;

  // Dünya hacmi
  G4VPhysicalVolume* fWorldPhys;

  // UI messenger (/shield/geom/..., /shield/gdml/...)
  DetectorMessenger* fMessenger;
};

#endif
