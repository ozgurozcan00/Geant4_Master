#include <QApplication>
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QFileDialog>
#include <QMessageBox>
#include <QGroupBox>
#include <QFormLayout>
#include <QTabWidget>
#include <QPainter>
#include <QVector>
#include <QCheckBox>
#include <QProgressBar>
#include <QTextStream>
#include <QFile>
#include <QFileInfo>
#include <QPoint>

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

// Geant4 Başlıkları
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4GDMLParser.hh"
#include "G4Material.hh"
#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4MaterialTable.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

// ====================================================================
// GLOBAL VERİ VE SİNYAL YÖNETİMİ
// ====================================================================
struct HitPoint {
    double x, y; // mm
};

struct SimResult {
    std::vector<double> energies; // MeV
    std::vector<double> angles;   // derece
    std::vector<HitPoint> hits;   // mm
    int    totalParticlesFired = 0;
    int    particlesDetected    = 0;
    double thicknessUsed        = 0.0; // cm
};
SimResult gSimResult;

// Yalnız primer parçacıkları say? (secondary hariç)
bool gCountOnlyPrimary = false;

// İlerleme çubuğunu güncellemek için callback arayüzü
class IProgressUpdate {
public:
    virtual ~IProgressUpdate() = default;
    virtual void UpdateProgress(int current, int total) = 0;
};
IProgressUpdate* gProgressHandler = nullptr;

class CustomDetector;
class CustomGenerator;
CustomDetector* gDetector   = nullptr;
CustomGenerator* gGenerator = nullptr;

// ====================================================================
// GRAFİK WIDGET'I (SPEKTRUM, ISI HARİTASI, AÇISAL DAĞILIM)
// ====================================================================
class ResultWidget : public QWidget {
    Q_OBJECT
public:
    enum Mode { Spectrum, Heatmap, Angular };

    ResultWidget(QWidget* parent = nullptr)
        : QWidget(parent), mMode(Spectrum), mLogScale(false), mMaxEnergy(1.0)
    {
        setBackgroundRole(QPalette::Base);
        setAutoFillBackground(true);
        setMinimumHeight(450);
    }

    void setData(const SimResult& res, double maxE) {
        mRes       = res;
        mMaxEnergy = (maxE > 0.0 ? maxE : 1.0);
        update();
    }

    void setMode(Mode m)        { mMode = m; update(); }
    void setLogScale(bool log)  { mLogScale = log; update(); }

    void saveImage() {
        QString fileName = QFileDialog::getSaveFileName(
            this, "Grafiği Kaydet", "analiz_grafigi.png", "PNG Files (*.png)");
        if (!fileName.isEmpty()) {
            QPixmap pixmap(this->size());
            this->render(&pixmap);
            pixmap.save(fileName);
        }
    }

protected:
    void paintEvent(QPaintEvent*) override {
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing);

        int w = width();
        int h = height();
        int margin = 60;
        QRect graphRect(margin, margin, w - 2 * margin, h - 2 * margin);

        painter.fillRect(rect(), Qt::white);
        painter.setPen(Qt::black);
        painter.drawRect(graphRect);

        if (mRes.particlesDetected == 0 && mRes.totalParticlesFired > 0) {
            painter.drawText(rect(), Qt::AlignCenter,
                             "Analiz Tamamlandı: Zırh %100 Başarılı\n(Hiçbir Radyasyon Geçmedi)");
            return;
        }
        if (mRes.totalParticlesFired == 0) {
            painter.drawText(rect(), Qt::AlignCenter,
                             "Simülasyon Bekleniyor...");
            return;
        }

        if (mMode == Spectrum) {
            drawHistogram(painter, graphRect, mRes.energies, mMaxEnergy,
                          "Enerji Spektrumu", "Enerji (MeV)");
        } else if (mMode == Angular) {
            drawHistogram(painter, graphRect, mRes.angles, 90.0,
                          "Saçılma Açısı Dağılımı", "Açı (Derece)");
        } else {
            drawHeatmap(painter, graphRect);
        }
    }

private:
    SimResult mRes;
    double   mMaxEnergy;
    Mode     mMode;
    bool     mLogScale;

    // Generic histogram çizici
    void drawHistogram(QPainter& p, QRect r,
                       const std::vector<double>& data,
                       double maxVal,
                       const QString& title,
                       const QString& xLabel)
    {
        p.setPen(Qt::black);
        p.drawText(QPoint(width()/2 - 120, 30), title);

        if (data.empty()) {
            p.drawText(r, Qt::AlignCenter, "Veri yok");
            return;
        }

        int numBins = 100;
        QVector<int> bins(numBins, 0);
        double binWidth = maxVal / numBins;
        if (binWidth <= 0) binWidth = 1.0;

        int maxCount = 0;
        for (double val : data) {
            if (val < 0) continue;
            int idx = static_cast<int>(val / binWidth);
            if (idx >= 0 && idx < numBins) {
                bins[idx]++;
                if (bins[idx] > maxCount) maxCount = bins[idx];
            }
        }
        if (maxCount == 0) maxCount = 1;

        double barWidth = static_cast<double>(r.width()) / numBins;
        p.setBrush(QColor(50, 120, 220, 150));
        p.setPen(Qt::NoPen);

        for (int i = 0; i < numBins; ++i) {
            double val = bins[i];
            double hRatio;
            if (mLogScale) {
                if (val < 1) val = 0.1;
                hRatio = std::log10(val) /
                         std::log10(maxCount > 1 ? maxCount : 10.0);
                if (hRatio < 0) hRatio = 0;
            } else {
                hRatio = val / maxCount;
            }
            double barH = hRatio * r.height();
            p.drawRect(QRectF(r.left() + i * barWidth,
                              r.bottom() - barH,
                              barWidth,
                              barH));
        }

        p.setPen(Qt::black);
        p.drawText(r.bottomLeft() + QPoint(0, 20), "0");
        p.drawText(r.bottomRight() + QPoint(-60, 20),
                   QString::number(maxVal, 'f', 2));
        p.drawText(QPoint(r.center().x() - 40, r.bottom() + 35), xLabel);

        QString yLabel = mLogScale ? "Log(N)" : QString("N = %1").arg(maxCount);
        p.drawText(r.topLeft() + QPoint(-55, 10), yLabel);
    }

    void drawHeatmap(QPainter& p, QRect r) {
        p.setPen(Qt::black);
        p.drawText(QPoint(width()/2 - 120, 30), "Radyasyon Çıkış Haritası");

        if (mRes.hits.empty()) {
            p.drawText(r, Qt::AlignCenter, "Veri yok");
            return;
        }

        int gridSz = 50;
        QVector<QVector<int>> grid(gridSz, QVector<int>(gridSz, 0));

        // Veri tabanlı bounding box
        double minX = 1e30, maxX = -1e30;
        double minY = 1e30, maxY = -1e30;
        for (const auto& hit : mRes.hits) {
            if (hit.x < minX) minX = hit.x;
            if (hit.x > maxX) maxX = hit.x;
            if (hit.y < minY) minY = hit.y;
            if (hit.y > maxY) maxY = hit.y;
        }
        if (minX >= maxX) { minX -= 10.0; maxX += 10.0; }
        if (minY >= maxY) { minY -= 10.0; maxY += 10.0; }

        int maxHits = 0;
        for (const auto& hit : mRes.hits) {
            int x = static_cast<int>(((hit.x - minX) / (maxX - minX)) * gridSz);
            int y = static_cast<int>(((hit.y - minY) / (maxY - minY)) * gridSz);
            if (x < 0) x = 0;
            if (x >= gridSz) x = gridSz - 1;
            if (y < 0) y = 0;
            if (y >= gridSz) y = gridSz - 1;
            grid[x][y]++;
            if (grid[x][y] > maxHits) maxHits = grid[x][y];
        }
        if (maxHits == 0) maxHits = 1;

        double cellW = static_cast<double>(r.width()) / gridSz;
        double cellH = static_cast<double>(r.height()) / gridSz;

        for (int i = 0; i < gridSz; ++i) {
            for (int j = 0; j < gridSz; ++j) {
                if (grid[i][j] > 0) {
                    double intensity = static_cast<double>(grid[i][j]) / maxHits;
                    int red  = static_cast<int>(255 * intensity);
                    int blue = static_cast<int>(255 * (1.0 - intensity));
                    p.fillRect(QRectF(r.left() + i * cellW,
                                      r.bottom() - (j + 1) * cellH,
                                      cellW, cellH),
                               QColor(red, 0, blue));
                }
            }
        }

        p.setPen(Qt::black);
        p.setBrush(Qt::NoBrush);
        p.drawRect(r);
        p.drawText(QPoint(r.left(), r.bottom() + 25),
                   "XY Düzlemi (Dedektör Düzlemi)");
    }
};

// ====================================================================
// GEANT4 MALZEME OLUŞTURMA
// ====================================================================

G4Material* CreateCustomMaterial(const G4String& name,
                                 const G4String& formula,
                                 double density)
{
    auto* nist = G4NistManager::Instance();

    // Eğer G4_ ile başlayan hazır bir malzeme kodu ise direkt NIST'ten al
    if (formula.find("G4_") != std::string::npos) {
        G4Material* m = nist->FindOrBuildMaterial(formula);
        if (m) return m;
    }

    if (auto* existing = G4Material::GetMaterial(name, false)) {
        return existing;
    }

    // Basit stokiyometri: "H 2 C 1 O 1" gibi bekliyoruz
    std::vector<G4String> tokens;
    std::istringstream iss(formula);
    G4String token;
    while (iss >> token) {
        tokens.push_back(token);
    }

    if (tokens.size() < 2 || tokens.size() % 2 != 0) {
        G4cerr << "Uyarı: Malzeme formülü hatalı: " << formula << G4endl;
        return nist->FindOrBuildMaterial("G4_AIR");
    }

    G4int nComponents = static_cast<G4int>(tokens.size() / 2);
    auto* mat = new G4Material(name, density, nComponents);

    for (size_t i = 0; i < tokens.size(); i += 2) {
        auto* element = nist->FindOrBuildElement(tokens[i]);
        if (!element) {
            G4cerr << "Uyarı: Element bulunamadı: " << tokens[i] << G4endl;
            continue;
        }
        G4int natoms = std::stoi(tokens[i+1]);
        mat->AddElement(element, natoms);
    }
    return mat;
}

// ====================================================================
// GEANT4 SINIFLARI
// ====================================================================

class CustomDetector : public G4VUserDetectorConstruction {
public:
    CustomDetector()
        : fScoringVolume(nullptr),
          fCadFileName(""),
          fUseCAD(false),
          fMatFormula("G4_Pb"),
          fDensity(11.34 * g/cm3),
          fThickness(2.5 * cm),
          fDetDistance(20.0 * cm),
          fDetSize(20.0 * cm)
    {}

    void SetCADFile(const G4String& path) {
        fCadFileName = path;
        fUseCAD      = !path.empty();
    }

    void ClearCAD() {
        fCadFileName.clear();
        fUseCAD = false;
    }

    void SetMaterial(const G4String& formula, double dens) {
        fMatFormula = formula;
        fDensity    = dens;
    }

    void SetThickness(double t) {
        fThickness = t;
    }

    void SetDetectorParams(double distCm, double sizeCm) {
        if (distCm <= 0.0) distCm = 20.0;
        if (sizeCm <= 0.0) sizeCm = 20.0;
        fDetDistance = distCm * cm;
        fDetSize     = sizeCm * cm;
    }

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    G4VPhysicalVolume* Construct() override {
        auto* nist = G4NistManager::Instance();
        auto* worldMat  = nist->FindOrBuildMaterial("G4_AIR");
        auto* shieldMat = CreateCustomMaterial("CustomShield", fMatFormula, fDensity);

        auto* solidWorld = new G4Box("World", 0.5*m, 0.5*m, 0.5*m);
        auto* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
        auto* physWorld  = new G4PVPlacement(
            nullptr, G4ThreeVector(), logicWorld,
            "World", nullptr, false, 0, true);

        // Zırh
        if (fUseCAD && !fCadFileName.empty()) {
            G4GDMLParser parser;
            parser.Read(fCadFileName, false);

            G4LogicalVolume* cadLogic = parser.GetVolume("Shield");
            if (!cadLogic) {
                G4cout << "Uyarı: GDML içinde 'Shield' bulunamadı, world volume kullanılıyor." << G4endl;
                cadLogic = parser.GetWorldVolume()->GetLogicalVolume();
            }

            if (cadLogic) {
                cadLogic->SetMaterial(shieldMat);
                new G4PVPlacement(nullptr, G4ThreeVector(0,0,0),
                                  cadLogic, "CAD_Shield",
                                  logicWorld, false, 0, true);
            }
        } else {
            auto* solidShield =
                new G4Box("Shield", 10*cm, 10*cm, fThickness/2.0);
            auto* logicShield =
                new G4LogicalVolume(solidShield, shieldMat, "Shield");
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,0),
                              logicShield, "Shield",
                              logicWorld, false, 0, true);
        }

        // Dedektör düzlemi (20x20 cm, ayarlanabilir, çok ince)
        G4double detHalfSize = fDetSize / 2.0;
        G4double detHalfTh   = 0.5 * mm;
        G4double detZ        = 0.0;

        if (fUseCAD) {
            // CAD'de kalınlık bilinmediği için doğrudan dünya merkezinden mesafe
            detZ = fDetDistance;
        } else {
            detZ = fThickness / 2.0 + fDetDistance;
        }

        auto* solidDet =
            new G4Box("DetectorPlane", detHalfSize, detHalfSize, detHalfTh);
        auto* logicDet =
            new G4LogicalVolume(solidDet, worldMat, "DetectorPlaneLV");
        new G4PVPlacement(nullptr, G4ThreeVector(0,0,detZ),
                          logicDet, "DetectorPlane",
                          logicWorld, false, 0, true);

        fScoringVolume = logicDet;

        return physWorld;
    }

private:
    G4LogicalVolume* fScoringVolume;
    G4String         fCadFileName;
    G4bool           fUseCAD;
    G4String         fMatFormula;
    double           fDensity;
    double           fThickness;

    double           fDetDistance; // world zırh arkasındaki mesafe
    double           fDetSize;     // dedektör boyutu (kenar)
};

class CustomGenerator : public G4VUserPrimaryGeneratorAction {
public:
    enum SourceMode {
        kPencil  = 0,
        kFan     = 1,
        kIsoHalf = 2,
        kIsoFull = 3
    };

    CustomGenerator()
        : mMode(kPencil),
          mSourceDistance(20.0 * cm),
          mSourceRadius(0.0 * cm),
          mFanAngle(10.0 * deg)
    {
        fParticleGun = new G4ParticleGun(1);
        SetParticle("gamma");
        SetEnergy(0.662); // MeV
    }

    ~CustomGenerator() override {
        delete fParticleGun;
    }

    void SetParticle(const G4String& name) {
        auto* p = G4ParticleTable::GetParticleTable()->FindParticle(name);
        if (p) {
            fParticleGun->SetParticleDefinition(p);
        } else {
            G4cerr << "Uyarı: Parçacık bulunamadı: " << name << G4endl;
        }
    }

    void SetEnergy(double energyMeV) {
        if (energyMeV <= 0.0) energyMeV = 0.1;
        fParticleGun->SetParticleEnergy(energyMeV * MeV);
    }

    void SetSourceMode(int idx) {
        if (idx < 0 || idx > 3) idx = 0;
        mMode = static_cast<SourceMode>(idx);
    }

    void SetSourceDistance(double distCm) {
        if (distCm <= 0.0) distCm = 20.0;
        mSourceDistance = distCm * cm;
    }

    void SetSourceRadius(double rCm) {
        if (rCm < 0.0) rCm = 0.0;
        mSourceRadius = rCm * cm;
    }

    void SetFanAngle(double angDeg) {
        if (angDeg <= 0.0) angDeg = 10.0;
        mFanAngle = angDeg * deg;
    }

    void GeneratePrimaries(G4Event* anEvent) override {
        // Temel kaynak noktası: z = -mSourceDistance
        G4ThreeVector pos(0., 0., -mSourceDistance);
        G4ThreeVector dir(0., 0., 1.);

        // Kaynak spot yarıçapı (disk)
        if (mSourceRadius > 0.0) {
            double u   = G4UniformRand();
            double r   = mSourceRadius * std::sqrt(u);
            double phi = 2.0 * CLHEP::pi * G4UniformRand();
            pos.setX(r * std::cos(phi));
            pos.setY(r * std::sin(phi));
        }

        switch (mMode) {
        case kPencil:
            dir = G4ThreeVector(0., 0., 1.);
            break;

        case kFan:
        {
            double theta = mFanAngle * G4UniformRand(); // [0, fan]
            double phi   = 2.0 * CLHEP::pi * G4UniformRand();
            dir = G4ThreeVector(std::sin(theta) * std::cos(phi),
                                std::sin(theta) * std::sin(phi),
                                std::cos(theta));
            if (dir.z() <= 0) dir.setZ(std::fabs(dir.z()));
        }
            break;

        case kIsoHalf: // +Z yarı uzay
        {
            double u     = G4UniformRand();   // cosθ in [0,1]
            double cosTh = u;
            double sinTh = std::sqrt(1.0 - u*u);
            double phi   = 2.0 * CLHEP::pi * G4UniformRand();
            dir = G4ThreeVector(sinTh * std::cos(phi),
                                sinTh * std::sin(phi),
                                cosTh);
        }
            break;

        case kIsoFull: // 4π
        default:
        {
            double u     = 2.0 * G4UniformRand() - 1.0; // cosθ in [-1,1]
            double cosTh = u;
            double sinTh = std::sqrt(1.0 - u*u);
            double phi   = 2.0 * CLHEP::pi * G4UniformRand();
            dir = G4ThreeVector(sinTh * std::cos(phi),
                                sinTh * std::sin(phi),
                                cosTh);
        }
            break;
        }

        fParticleGun->SetParticlePosition(pos);
        fParticleGun->SetParticleMomentumDirection(dir);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

private:
    G4ParticleGun* fParticleGun;
    SourceMode     mMode;
    G4double       mSourceDistance;
    G4double       mSourceRadius;
    G4double       mFanAngle;
};

class CustomEventAction : public G4UserEventAction {
public:
    void BeginOfEventAction(const G4Event* evt) override {
        G4int eventID = evt->GetEventID();
        if (gProgressHandler && eventID % 100 == 0) {
            auto* rm = G4RunManager::GetRunManager();
            G4int nEvents = rm->GetNumberOfEventsToBeProcessed();
            gProgressHandler->UpdateProgress(eventID, nEvents);
        }
    }
};

class CustomSteppingAction : public G4UserSteppingAction {
public:
    explicit CustomSteppingAction(const CustomDetector* det)
        : fDet(det) {}

    void UserSteppingAction(const G4Step* step) override {
        auto* scoring = fDet->GetScoringVolume();
        if (!scoring) return;

        auto preVol  = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
        auto postVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        if (!preVol || !postVol) return;

        if (gCountOnlyPrimary && step->GetTrack()->GetParentID() != 0) {
            return;
        }

        // Dedektör düzlemine giriş anını yakala
        if (preVol->GetLogicalVolume() != scoring &&
            postVol->GetLogicalVolume() == scoring)
        {
            G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();

            double energyMeV =
                step->GetPostStepPoint()->GetKineticEnergy() / MeV;

            G4ThreeVector direction =
                step->GetTrack()->GetMomentumDirection();
            double angleRad =
                direction.angle(G4ThreeVector(0,0,1));
            double angleDeg = angleRad * 180.0 / CLHEP::pi;

            gSimResult.energies.push_back(energyMeV);
            gSimResult.angles.push_back(angleDeg);
            gSimResult.hits.push_back({pos.x(), pos.y()});
            gSimResult.particlesDetected++;
        }
    }

private:
    const CustomDetector* fDet;
};

class CustomRunAction : public G4UserRunAction {
public:
    void BeginOfRunAction(const G4Run*) override {
        gSimResult.energies.clear();
        gSimResult.angles.clear();
        gSimResult.hits.clear();
        gSimResult.particlesDetected    = 0;
        gSimResult.totalParticlesFired  = 0;
        gSimResult.thicknessUsed        = 0.0;
    }

    void EndOfRunAction(const G4Run*) override {}
};

class CustomActionInit : public G4VUserActionInitialization {
public:
    explicit CustomActionInit(CustomDetector* det)
        : fDet(det) {}

    void BuildForMaster() const override {
        SetUserAction(new CustomRunAction);
    }

    void Build() const override {
        auto* gen = new CustomGenerator;
        SetUserAction(gen);
        gGenerator = gen;

        SetUserAction(new CustomRunAction);
        SetUserAction(new CustomEventAction);
        SetUserAction(new CustomSteppingAction(fDet));
    }

private:
    CustomDetector* fDet;
};

// ====================================================================
// PROFESYONEL QT5 ARAYÜZÜ
// ====================================================================
class MainWindow : public QWidget, public IProgressUpdate {
    Q_OBJECT

public:
    explicit MainWindow(G4RunManager* runMgr)
        : QWidget(nullptr),
          fRunManager(runMgr),
          mUsingCAD(false)
    {
        setWindowTitle("Profesyonel Zırhlama Analiz Laboratuvarı");
        resize(1200, 800);
        gProgressHandler = this;

        setStyleSheet(
            "QGroupBox { font-weight: bold; border: 1px solid #ccc; "
            "border-radius: 6px; margin-top: 10px; background-color: #f9f9f9; } "
            "QGroupBox::title { subcontrol-origin: margin; "
            "subcontrol-position: top left; padding: 0 5px; "
            "background-color: #f9f9f9; } "
            "QPushButton { border-radius: 4px; padding: 6px; }"
        );

        auto* mainLayout = new QHBoxLayout(this);

        // --- SOL PANEL ---
        QWidget* leftPanel = new QWidget();
        leftPanel->setFixedWidth(420);
        auto* leftLayout = new QVBoxLayout(leftPanel);

        // 1. Zırh Özellikleri
        auto* grpGeo = new QGroupBox("1. Zırh Özellikleri", this);
        auto* layGeo = new QFormLayout(grpGeo);

        cmbPresets = new QComboBox(this);
        cmbPresets->addItems({
            "-- Özel Malzeme --",
            "Kurşun (Pb)",
            "Beton (Standart)",
            "Su",
            "Demir (Fe)",
            "Alüminyum",
            "Borlu Polietilen"
        });
        connect(cmbPresets,
                QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &MainWindow::onPresetChanged);

        txtFormula   = new QLineEdit("G4_Pb", this);
        txtDensity   = new QLineEdit("11.34", this);
        txtThickness = new QLineEdit("2.5", this);

        btnLoadCad = new QPushButton("CAD Geometrisi Yükle (.gdml)", this);
        connect(btnLoadCad, &QPushButton::clicked,
                this, &MainWindow::loadCadFile);

        layGeo->addRow("Hazır Malzeme:",   cmbPresets);
        layGeo->addRow("Kimyasal Formül:", txtFormula);
        layGeo->addRow("Yoğunluk (g/cm³):", txtDensity);
        layGeo->addRow("Kalınlık (cm):",    txtThickness);
        layGeo->addRow(btnLoadCad);

        // 2. Radyasyon Kaynağı
        auto* grpRad = new QGroupBox("2. Radyasyon Kaynağı", this);
        auto* layRad = new QFormLayout(grpRad);

        cmbParticle = new QComboBox(this);
        cmbParticle->addItems({"gamma", "e-", "alpha", "neutron", "proton"});

        txtEnergy = new QLineEdit("0.662", this); // MeV
        txtEvents = new QLineEdit("10000", this);

        cmbSourceGeom = new QComboBox(this);
        cmbSourceGeom->addItems({
            "Pencil Beam",
            "Fan Beam",
            "İzotropik (yarım uzay)",
            "İzotropik (4π)"
        });

        txtSourceDist   = new QLineEdit("20.0", this); // cm
        txtSourceRadius = new QLineEdit("0.0", this);  // cm
        txtFanAngle     = new QLineEdit("10.0", this); // derece

        chkPrimaryOnly = new QCheckBox("Sadece primer parçacıkları say", this);
        connect(chkPrimaryOnly, &QCheckBox::toggled,
                this, &MainWindow::updateCountingMode);

        layRad->addRow("Parçacık Tipi:",       cmbParticle);
        layRad->addRow("Enerji (MeV):",        txtEnergy);
        layRad->addRow("Olay Sayısı:",         txtEvents);
        layRad->addRow("Kaynak Geometrisi:",   cmbSourceGeom);
        layRad->addRow("Kaynak Z Mesafesi (cm):", txtSourceDist);
        layRad->addRow("Kaynak Spot Yarıçapı (cm):", txtSourceRadius);
        layRad->addRow("Fan Açısı (derece):", txtFanAngle);
        layRad->addRow(chkPrimaryOnly);

        // 3. Kalınlık Taraması
        auto* grpScan = new QGroupBox("3. Kalınlık Taraması (Opsiyonel)", this);
        auto* layScan = new QFormLayout(grpScan);

        txtScanMin    = new QLineEdit("0.5", this);
        txtScanMax    = new QLineEdit("10.0", this);
        txtScanStep   = new QLineEdit("0.5", this);
        txtScanEvents = new QLineEdit("5000", this);

        btnRunScan = new QPushButton("Kalınlık Taramasını Başlat", this);
        btnRunScan->setStyleSheet(
            "background-color: #4CAF50; color: white; font-weight: bold;");
        connect(btnRunScan, &QPushButton::clicked,
                this, &MainWindow::runThicknessScan);

        layScan->addRow("Min Kalınlık (cm):", txtScanMin);
        layScan->addRow("Max Kalınlık (cm):", txtScanMax);
        layScan->addRow("Adım (cm):",         txtScanStep);
        layScan->addRow("Olay/Sweep:",        txtScanEvents);
        layScan->addWidget(btnRunScan);

        // 4. Dedektör Düzlemi
        auto* grpDet = new QGroupBox("4. Dedektör Düzlemi", this);
        auto* layDet = new QFormLayout(grpDet);

        txtDetDistance = new QLineEdit("20.0", this); // cm
        txtDetSize     = new QLineEdit("20.0", this); // cm

        layDet->addRow("Zırh arkasındaki mesafe (cm):", txtDetDistance);
        layDet->addRow("Dedektör boyutu (kenar, cm):",  txtDetSize);

        // Ana Run butonu
        btnRun = new QPushButton("SİMÜLASYONU BAŞLAT", this);
        btnRun->setStyleSheet(
            "background-color: #2196F3; color: white; font-weight: bold; "
            "padding: 15px; font-size: 14px;");
        connect(btnRun, &QPushButton::clicked,
                this, &MainWindow::runSimulation);

        progressBar = new QProgressBar(this);
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignCenter);

        leftLayout->addWidget(grpGeo);
        leftLayout->addWidget(grpRad);
        leftLayout->addWidget(grpScan);
        leftLayout->addWidget(grpDet);
        leftLayout->addStretch();
        leftLayout->addWidget(progressBar);
        leftLayout->addWidget(btnRun);

        // --- SAĞ PANEL ---
        auto* tabs = new QTabWidget(this);

        // Sekme 1: Grafikler
        QWidget* tabGraph = new QWidget();
        auto* layGraph = new QVBoxLayout(tabGraph);

        auto* toolBar = new QHBoxLayout();
        chkLogScale = new QCheckBox("Logaritmik Ölçek", this);
        connect(chkLogScale, &QCheckBox::toggled,
                this, &MainWindow::updateGraphSettings);

        cmbGraphMode = new QComboBox(this);
        cmbGraphMode->addItem("Enerji Spektrumu (Energy Spectrum)");
        cmbGraphMode->addItem("Isı Haritası (Heatmap)");
        cmbGraphMode->addItem("Açısal Dağılım (Angular Distribution)");
        connect(cmbGraphMode,
                QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &MainWindow::updateGraphSettings);

        QPushButton* btnSaveImg = new QPushButton("Resim Kaydet", this);
        connect(btnSaveImg, &QPushButton::clicked,
                this, &MainWindow::savePlot);

        QPushButton* btnSaveCsv = new QPushButton("CSV Olarak Kaydet", this);
        connect(btnSaveCsv, &QPushButton::clicked,
                this, &MainWindow::saveCsv);

        toolBar->addWidget(cmbGraphMode);
        toolBar->addWidget(chkLogScale);
        toolBar->addStretch();
        toolBar->addWidget(btnSaveImg);
        toolBar->addWidget(btnSaveCsv);

        resultWidget = new ResultWidget(this);
        layGraph->addLayout(toolBar);
        layGraph->addWidget(resultWidget);

        // Sekme 2: Detaylı Rapor
        QWidget* tabCalc = new QWidget();
        auto* layCalc = new QVBoxLayout(tabCalc);
        txtReport = new QLabel("Analiz sonuçları burada görünecek...", this);
        txtReport->setAlignment(Qt::AlignTop | Qt::AlignLeft);
        txtReport->setStyleSheet(
            "font-family: Monospace; font-size: 13px; "
            "background-color: white; padding: 15px; border: 1px solid #ddd;");
        txtReport->setWordWrap(true);
        txtReport->setTextFormat(Qt::PlainText);
        layCalc->addWidget(txtReport);

        tabs->addTab(tabGraph, "Görsel Analiz");
        tabs->addTab(tabCalc, "Fiziksel Rapor & HVL");

        mainLayout->addWidget(leftPanel);
        mainLayout->addWidget(tabs);
    }

    // IProgressUpdate
    void UpdateProgress(int current, int total) override {
        if (total <= 0) total = 1;
        progressBar->setMaximum(total);
        progressBar->setValue(current);
        QApplication::processEvents();
    }

public slots:
    void onPresetChanged(int index) {
        if (index == 1) { setMat("G4_Pb",        "11.34"); }
        else if (index == 2) { setMat("G4_CONCRETE", "2.3");  }
        else if (index == 3) { setMat("G4_WATER",    "1.0");  }
        else if (index == 4) { setMat("G4_Fe",       "7.87"); }
        else if (index == 5) { setMat("G4_Al",       "2.70"); }
        else if (index == 6) { setMat("H 2 C 1",     "0.93"); }
    }

    void setMat(const QString& f, const QString& d) {
        txtFormula->setText(f);
        txtDensity->setText(d);
    }

    void updateCountingMode(bool checked) {
        gCountOnlyPrimary = checked;
    }

    void loadCadFile() {
        QString fileName = QFileDialog::getOpenFileName(
            this, "CAD Yükle", "", "GDML (*.gdml)");
        if (!fileName.isEmpty()) {
            QFileInfo fi(fileName);
            btnLoadCad->setText("Yüklendi: " + fi.fileName());
            if (gDetector) {
                gDetector->SetCADFile(fileName.toStdString());
            }
            txtThickness->setEnabled(false);
            mUsingCAD = true;
        }
    }

    void applySourceSettings() {
        if (!gGenerator) return;

        gGenerator->SetSourceMode(cmbSourceGeom->currentIndex());
        double dist  = txtSourceDist->text().toDouble();
        double rad   = txtSourceRadius->text().toDouble();
        double fAng  = txtFanAngle->text().toDouble();

        gGenerator->SetSourceDistance(dist);
        gGenerator->SetSourceRadius(rad);
        gGenerator->SetFanAngle(fAng);

        gGenerator->SetParticle(cmbParticle->currentText().toStdString());
        double E = txtEnergy->text().toDouble();
        if (E <= 0.0) E = 0.1;
        gGenerator->SetEnergy(E);
    }

    void applyDetectorSettings() {
        if (!gDetector) return;

        double dDist = txtDetDistance->text().toDouble();
        double dSize = txtDetSize->text().toDouble();
        gDetector->SetDetectorParams(dDist, dSize);
    }

    void applyMaterialSettings(bool allowThickness) {
        if (!gDetector) return;

        QString formulaStr = txtFormula->text();
        double densVal = txtDensity->text().toDouble();
        if (densVal <= 0.0) densVal = 1.0;

        gDetector->SetMaterial(formulaStr.toStdString(),
                               densVal * g/cm3);

        if (!mUsingCAD && allowThickness) {
            double thick = txtThickness->text().toDouble();
            if (thick <= 0.0) thick = 0.1;
            gDetector->SetThickness(thick * cm);
        }
    }

    void runSimulation() {
        if (!fRunManager || !gDetector) {
            QMessageBox::critical(this, "Hata",
                                  "RunManager veya Dedektör tanımlı değil!");
            return;
        }
        if (!gGenerator) {
            QMessageBox::critical(this, "Hata",
                                  "Primary Generator henüz oluşturulmadı.\n"
                                  "Programı kapatıp yeniden deneyin.");
            return;
        }

        btnRun->setEnabled(false);
        btnRunScan->setEnabled(false);
        btnRun->setText("Hesaplanıyor...");
        progressBar->setValue(0);

        applyMaterialSettings(true);
        applyDetectorSettings();
        applySourceSettings();

        fRunManager->ReinitializeGeometry();

        int nEvents = txtEvents->text().toInt();
        if (nEvents <= 0) nEvents = 1000;
        gSimResult.totalParticlesFired = nEvents;
        gSimResult.thicknessUsed       =
            (mUsingCAD ? 0.0 : txtThickness->text().toDouble());

        fRunManager->Initialize();
        fRunManager->BeamOn(nEvents);

        UpdateProgress(nEvents, nEvents);
        updateGraphSettings();
        generateReport();

        btnRun->setText("SİMÜLASYONU BAŞLAT");
        btnRun->setEnabled(true);
        btnRunScan->setEnabled(true);
        QMessageBox::information(this, "Başarılı", "Simülasyon tamamlandı.");
    }

    void runThicknessScan() {
        if (!fRunManager || !gDetector || !gGenerator) {
            QMessageBox::critical(this, "Hata",
                                  "RunManager / Dedektör / Generator hazır değil.");
            return;
        }
        if (mUsingCAD) {
            QMessageBox::warning(this, "Uyarı",
                                 "CAD geometrisi ile kalınlık taraması (analitik HVL) "
                                 "anlamlı değil.\n"
                                 "Önce CAD kullanımını devre dışı bırakın.");
            return;
        }

        double tMin  = txtScanMin->text().toDouble();
        double tMax  = txtScanMax->text().toDouble();
        double tStep = txtScanStep->text().toDouble();
        int    nEv   = txtScanEvents->text().toInt();

        if (tMin <= 0.0 || tMax <= tMin || tStep <= 0.0) {
            QMessageBox::warning(this, "Uyarı",
                                 "Kalınlık taraması için geçersiz aralık/adım.");
            return;
        }
        if (nEv <= 0) nEv = 1000;

        applyMaterialSettings(false);   // kalınlığı döngü içinde ayarlayacağız
        applyDetectorSettings();
        applySourceSettings();

        struct ScanPoint {
            double t;
            double transmission;
            double mu;
            double hvl;
            double tvl;
        };
        std::vector<ScanPoint> scanData;

        btnRun->setEnabled(false);
        btnRunScan->setEnabled(false);
        progressBar->setValue(0);
        progressBar->setMaximum(100);

        int nSteps = static_cast<int>(std::floor((tMax - tMin)/tStep)) + 1;
        if (nSteps <= 0) nSteps = 1;

        for (int i = 0; i < nSteps; ++i) {
            double t = tMin + i * tStep;
            if (t > tMax + 1e-9) break;

            gDetector->SetThickness(t * cm);
            fRunManager->ReinitializeGeometry();

            gSimResult.totalParticlesFired = nEv;
            gSimResult.thicknessUsed       = t;

            fRunManager->Initialize();
            fRunManager->BeamOn(nEv);

            double I0 = nEv;
            double I  = gSimResult.particlesDetected;
            double tr = (I0 > 0.0 ? I / I0 : 0.0);
            double mu = 0.0, hvl = 0.0, tvl = 0.0;

            if (tr > 0.0 && tr < 1.0) {
                mu  = -std::log(tr) / t;
                hvl = std::log(2.0) / mu;
                tvl = std::log(10.0) / mu;
            }

            scanData.push_back({t, tr, mu, hvl, tvl});

            int prog = static_cast<int>( (100.0 * (i+1)) / nSteps );
            progressBar->setValue(prog);
            QApplication::processEvents();
        }

        // CSV kaydet
        QString fileName = QFileDialog::getSaveFileName(
            this, "Kalınlık Taramasını CSV Olarak Kaydet",
            "thickness_scan.csv", "CSV Files (*.csv)");
        if (!fileName.isEmpty()) {
            QFile file(fileName);
            if (file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
                QTextStream s(&file);
                s << "Thickness(cm),Transmission,mu(1/cm),HVL(cm),TVL(cm)\n";
                for (const auto& sp : scanData) {
                    s << sp.t << "," << sp.transmission << ","
                      << sp.mu << "," << sp.hvl << "," << sp.tvl << "\n";
                }
                file.close();
            }
        }

        // Rapor içine özet ekle
        std::stringstream ss;
        ss << "\n\n=== KALINLIK TARAMASI ÖZETİ ===\n";
        ss << "Taramada kullanılan aralık: " << tMin << " - " << tMax
           << " cm, adım: " << tStep << " cm, olay/sweep: " << nEv << "\n";
        if (!scanData.empty()) {
            ss << "Örnek noktalar:\n";
            for (size_t i = 0; i < scanData.size() && i < 5; ++i) {
                ss << "  t = " << scanData[i].t << " cm, "
                   << "T = " << scanData[i].transmission
                   << ", mu = " << scanData[i].mu << " 1/cm\n";
            }
            if (scanData.size() > 5) ss << "  ... (" << scanData.size()-5
                                        << " nokta daha)\n";
        }
        QString old = txtReport->text();
        txtReport->setText(old + QString::fromStdString(ss.str()));

        btnRun->setEnabled(true);
        btnRunScan->setEnabled(true);
        QMessageBox::information(this, "Tamam",
                                 "Kalınlık taraması tamamlandı. CSV dosyasını inceleyebilirsiniz.");
    }

    void updateGraphSettings() {
        resultWidget->setLogScale(chkLogScale->isChecked());

        if (cmbGraphMode->currentIndex() == 0)
            resultWidget->setMode(ResultWidget::Spectrum);
        else if (cmbGraphMode->currentIndex() == 1)
            resultWidget->setMode(ResultWidget::Heatmap);
        else
            resultWidget->setMode(ResultWidget::Angular);

        double maxE = txtEnergy->text().toDouble();
        resultWidget->setData(gSimResult, maxE);
    }

    void savePlot() {
        resultWidget->saveImage();
    }

    void saveCsv() {
        if (gSimResult.energies.empty()) {
            QMessageBox::warning(this, "Uyarı",
                                 "Kaydedilecek veri yok.");
            return;
        }

        QString fileName = QFileDialog::getSaveFileName(
            this, "CSV Kaydet", "analiz_sonuclari.csv",
            "CSV Files (*.csv)");
        if (fileName.isEmpty()) return;

        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
            QTextStream stream(&file);
            stream << "Enerji(MeV),Aci(Deg),PozisyonX(mm),PozisyonY(mm)\n";
            for (size_t i = 0; i < gSimResult.energies.size(); ++i) {
                stream << gSimResult.energies[i] << ","
                       << gSimResult.angles[i] << ","
                       << gSimResult.hits[i].x << ","
                       << gSimResult.hits[i].y << "\n";
            }
            file.close();
        }
    }

    void generateReport() {
        double I0 = gSimResult.totalParticlesFired;
        double I  = gSimResult.particlesDetected;
        double x  = gSimResult.thicknessUsed; // cm

        double transmission = (I0 > 0.0 ? I / I0 : 0.0);

        // Ortalama enerji ve açı
        double meanE = 0.0, rmsE = 0.0;
        if (!gSimResult.energies.empty()) {
            for (double e : gSimResult.energies) meanE += e;
            meanE /= gSimResult.energies.size();
            for (double e : gSimResult.energies) rmsE += (e-meanE)*(e-meanE);
            rmsE = std::sqrt(rmsE / gSimResult.energies.size());
        }

        double meanAng = 0.0, rmsAng = 0.0;
        if (!gSimResult.angles.empty()) {
            for (double a : gSimResult.angles) meanAng += a;
            meanAng /= gSimResult.angles.size();
            for (double a : gSimResult.angles) rmsAng += (a-meanAng)*(a-meanAng);
            rmsAng = std::sqrt(rmsAng / gSimResult.angles.size());
        }

        std::stringstream ss;
        ss << "=== DETAYLI ZIRHLAMA ANALİZ RAPORU ===\n\n";
        ss << "PARAMETRELER:\n";
        ss << "  - Malzeme: " << txtFormula->text().toStdString()
           << " (" << txtDensity->text().toStdString() << " g/cm^3)\n";
        if (!mUsingCAD) {
            ss << "  - Kalınlık: " << x << " cm\n";
        } else {
            ss << "  - Geometri: CAD (GDML) tabanlı, etkin kalınlık analitik olarak hesaplanmadı.\n";
        }
        ss << "  - Parçacık: " << cmbParticle->currentText().toStdString()
           << " (E = " << txtEnergy->text().toStdString() << " MeV)\n";
        ss << "  - Kaynak geometrisi: "
           << cmbSourceGeom->currentText().toStdString() << "\n";
        ss << "  - Sayım modu: "
           << (gCountOnlyPrimary ? "Yalnız primer parçacıklar" : "Tüm parçacıklar")
           << "\n\n";

        ss << "İSTATİSTİKLER:\n";
        ss << "  - Toplam Parçacık (I0): " << static_cast<int>(I0) << "\n";
        ss << "  - Geçen Parçacık (I)  : " << static_cast<int>(I) << "\n";
        ss << "  - Sönümleme Oranı     : "
           << std::fixed << std::setprecision(4)
           << (1.0 - transmission) * 100.0 << " %\n";
        ss << "  - İletim Oranı        : "
           << (transmission * 100.0) << " %\n\n";

        ss << "ENERJİ VE AÇI İSTATİSTİKLERİ (GEÇEN PARÇACIKLAR):\n";
        if (!gSimResult.energies.empty()) {
            ss << "  - Ortalama Enerji     : " << meanE << " MeV\n";
            ss << "  - RMS Enerji Genişliği: " << rmsE << " MeV\n";
        } else {
            ss << "  - Enerji istatistiği: Veri yok.\n";
        }
        if (!gSimResult.angles.empty()) {
            ss << "  - Ortalama Açı        : " << meanAng << " derece\n";
            ss << "  - RMS Açı Genişliği   : " << rmsAng << " derece\n";
        } else {
            ss << "  - Açı istatistiği: Veri yok.\n";
        }
        ss << "\n";

        if (!mUsingCAD) {
            if (transmission > 0.0 && transmission < 1.0 && x > 0.0) {
                double mu  = -std::log(transmission) / x;   // 1/cm
                double hvl = std::log(2.0) / mu;            // cm
                double tvl = std::log(10.0) / mu;           // cm

                ss << "HESAPLANAN FİZİKSEL DEĞERLER:\n";
                ss << "  - Lineer Zayıflatma Katsayısı (μ): "
                   << mu << " 1/cm\n";
                ss << "  - Yarı Değer Kalınlığı (HVL)     : "
                   << hvl << " cm\n";
                ss << "  - 1/10 Değer Kalınlığı (TVL)     : "
                   << tvl << " cm\n";
                ss << "    (Not: Bu değerler simülasyon istatistiğine dayalıdır)\n";
            } else if (transmission == 0.0 && I0 > 0.0) {
                ss << "SONUÇ: TAM SÖNÜMLEME.\n";
                ss << "Radyasyon zırhı geçemedi. "
                   << "Daha hassas bir HVL hesabı için olay sayısını artırmak gerekir.\n";
            } else {
                ss << "SONUÇ: İstatistik yetersiz veya iletim oranı 1.0'a çok yakın.\n";
            }
        } else {
            ss << "HVL / TVL: CAD geometrisi için doğrudan kalınlık parametresine dayalı "
                  "analitik HVL hesaplanmıyor.\n"
                  "İstenirse belirli yönlerde efektif kalınlık hesabı için ek geometri analizi\n"
                  "modülü eklenebilir.\n";
        }

        txtReport->setText(QString::fromStdString(ss.str()));
    }

private:
    G4RunManager*  fRunManager;
    QLineEdit*     txtFormula;
    QLineEdit*     txtDensity;
    QLineEdit*     txtThickness;
    QLineEdit*     txtEnergy;
    QLineEdit*     txtEvents;

    QLineEdit*     txtScanMin;
    QLineEdit*     txtScanMax;
    QLineEdit*     txtScanStep;
    QLineEdit*     txtScanEvents;

    QLineEdit*     txtDetDistance;
    QLineEdit*     txtDetSize;

    QLineEdit*     txtSourceDist;
    QLineEdit*     txtSourceRadius;
    QLineEdit*     txtFanAngle;

    QPushButton*   btnRun;
    QPushButton*   btnRunScan;
    QPushButton*   btnLoadCad;
    QComboBox*     cmbPresets;
    QComboBox*     cmbParticle;
    QComboBox*     cmbGraphMode;
    QComboBox*     cmbSourceGeom;
    QCheckBox*     chkLogScale;
    QCheckBox*     chkPrimaryOnly;
    QProgressBar*  progressBar;
    ResultWidget*  resultWidget;
    QLabel*        txtReport;

    bool           mUsingCAD;
};

// ====================================================================
// main
// ====================================================================
int main(int argc, char** argv) {
    QApplication app(argc, argv);

    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

    // Dedektör
    gDetector = new CustomDetector();
    runManager->SetUserInitialization(gDetector);

    // Fizik Listesi: Shielding (gamma + nötron zırhlama için uygun referans list)
    G4PhysListFactory physFactory;
    G4VUserPhysicsList* physicsList =
        physFactory.GetReferencePhysList("Shielding");
    runManager->SetUserInitialization(physicsList);

    // Aksiyonlar
    runManager->SetUserInitialization(new CustomActionInit(gDetector));

    MainWindow window(runManager);
    window.show();

    int ret = app.exec();

    delete runManager;
    return ret;
}

#include "shielding_qt_sim.moc"
