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

// ====================================================================
// GLOBAL VERİ VE SİNYAL YÖNETİMİ
// ====================================================================
struct HitPoint {
    double x, y; // mm cinsinden
};

struct SimResult {
    std::vector<double> energies; // MeV cinsinden
    std::vector<double> angles;   // derece
    std::vector<HitPoint> hits;   // mm
    int    totalParticlesFired = 0;
    int    particlesDetected    = 0;
    double thicknessUsed        = 0.0; // cm cinsinden (rapor için)
};
SimResult gSimResult;

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

    // Generic Histogram Çizici (Enerji veya Açı için)
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

        // Eksen etiketleri
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

        // 20 cm x 20 cm'lik alana projeksiyon (mm cinsinden 200 x 200)
        double sizeXY = 20.0 * cm; // burada cm -> 10*mm, pos.x pos.y mm cinsinden
        int maxHits = 0;

        for (const auto& hit : mRes.hits) {
            int x = static_cast<int>((hit.x + sizeXY/2.0) / sizeXY * gridSz);
            int y = static_cast<int>((hit.y + sizeXY/2.0) / sizeXY * gridSz);
            if (x >= 0 && x < gridSz && y >= 0 && y < gridSz) {
                grid[x][y]++;
                if (grid[x][y] > maxHits) maxHits = grid[x][y];
            }
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
                   "XY Düzlemi (Zırhın Arkası)");
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

    // Aynı isimli malzeme daha önce oluşturulmuşsa tekrar yaratma
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
          fThickness(2.5 * cm)
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
        fThickness = t; // fUseCAD bayrağını değiştirmiyoruz
    }

    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    G4VPhysicalVolume* Construct() override {
        auto* nist = G4NistManager::Instance();
        auto* worldMat  = nist->FindOrBuildMaterial("G4_AIR");
        auto* shieldMat = CreateCustomMaterial("CustomShield", fMatFormula, fDensity);

        // Dünya
        auto* solidWorld = new G4Box("World", 0.5*m, 0.5*m, 0.5*m);
        auto* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
        auto* physWorld  = new G4PVPlacement(
            nullptr, G4ThreeVector(), logicWorld,
            "World", nullptr, false, 0, true);

        // CAD mi, basit kutu mu?
        fScoringVolume = nullptr;

        if (fUseCAD && !fCadFileName.empty()) {
            G4GDMLParser parser;
            parser.Read(fCadFileName, false);

            // GDML içinde "Shield" isminde bir volume varsa onu al, yoksa world volume'yi kullan
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
                fScoringVolume = cadLogic;
            }
        } else {
            // Basit dikdörtgen zırh
            auto* solidShield =
                new G4Box("Shield", 10*cm, 10*cm, fThickness/2.0);
            auto* logicShield =
                new G4LogicalVolume(solidShield, shieldMat, "Shield");
            new G4PVPlacement(nullptr, G4ThreeVector(0,0,0),
                              logicShield, "Shield",
                              logicWorld, false, 0, true);
            fScoringVolume = logicShield;
        }

        return physWorld;
    }

private:
    G4LogicalVolume* fScoringVolume;
    G4String         fCadFileName;
    G4bool           fUseCAD;
    G4String         fMatFormula;
    double           fDensity;
    double           fThickness;
};

class CustomGenerator : public G4VUserPrimaryGeneratorAction {
public:
    CustomGenerator()
    {
        fParticleGun = new G4ParticleGun(1);
        SetParticle("gamma");
        SetEnergy(0.662); // MeV
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -20.0*cm));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
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

    // energyMeV: doğrudan MeV cinsinden giriliyor
    void SetEnergy(double energyMeV) {
        if (energyMeV <= 0.0) energyMeV = 0.1;
        fParticleGun->SetParticleEnergy(energyMeV * MeV);
    }

    void GeneratePrimaries(G4Event* anEvent) override {
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

private:
    G4ParticleGun* fParticleGun;
};

// --- Event Action (İlerleme Çubuğu İçin) ---
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

        // Zırh hacminden Dünya hacmine geçiş anını yakala
        if (preVol->GetLogicalVolume() == scoring &&
            postVol->GetLogicalVolume()->GetName() == "World")
        {
            G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
            if (pos.z() > 0) {
                // Enerji: MeV cinsine çevir
                double energyMeV =
                    step->GetTrack()->GetKineticEnergy() / MeV;

                // Açı Hesabı (Z ekseni ile)
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

    void EndOfRunAction(const G4Run*) override {
        // İstersen burada global istatistikleri yazdırabilirsin
    }
};

class CustomActionInit : public G4VUserActionInitialization {
public:
    explicit CustomActionInit(CustomDetector* det)
        : fDet(det) {}

    void BuildForMaster() const override {
        SetUserAction(new CustomRunAction);
    }

    void Build() const override {
        // Generator oluşturup global pointer'a atıyoruz
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
        resize(1100, 750);
        gProgressHandler = this; // callback bağla

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
        leftPanel->setFixedWidth(350);
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

        layGeo->addRow("Hazır Malzeme:",  cmbPresets);
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

        layRad->addRow("Parçacık Tipi:", cmbParticle);
        layRad->addRow("Enerji (MeV):",  txtEnergy);
        layRad->addRow("Olay Sayısı:",   txtEvents);

        // Kontrol Butonları
        btnRun = new QPushButton("SİMÜLASYONU BAŞLAT", this);
        btnRun->setStyleSheet(
            "background-color: #2196F3; color: white; font-weight: bold; "
            "padding: 15px; font-size: 14px;"
        );
        connect(btnRun, &QPushButton::clicked,
                this, &MainWindow::runSimulation);

        progressBar = new QProgressBar(this);
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignCenter);

        leftLayout->addWidget(grpGeo);
        leftLayout->addWidget(grpRad);
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

    // IProgressUpdate Arayüzü Implementasyonu
    void UpdateProgress(int current, int total) override {
        if (total <= 0) total = 1;
        progressBar->setMaximum(total);
        progressBar->setValue(current);
        QApplication::processEvents(); // Arayüz donmasın
    }

public slots:
    void onPresetChanged(int index) {
        // Hazır malzemeler
        if (index == 1) { setMat("G4_Pb",        "11.34"); } // Kurşun
        else if (index == 2) { setMat("G4_CONCRETE", "2.3"); }  // Beton
        else if (index == 3) { setMat("G4_WATER",    "1.0"); }  // Su
        else if (index == 4) { setMat("G4_Fe",       "7.87"); } // Demir
        else if (index == 5) { setMat("G4_Al",       "2.70"); } // Alüminyum
        else if (index == 6) { setMat("H 2 C 1",     "0.93"); } // Basit PE
    }

    void setMat(const QString& f, const QString& d) {
        txtFormula->setText(f);
        txtDensity->setText(d);
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

    void runSimulation() {
        if (!fRunManager || !gDetector) {
            QMessageBox::critical(this, "Hata",
                                  "RunManager veya Dedektör tanımlı değil!");
            return;
        }
        if (!gGenerator) {
            QMessageBox::critical(this, "Hata",
                                  "Primary Generator henüz oluşturulmadı."
                                  "\nProgramı kapatıp yeniden deneyin.");
            return;
        }

        btnRun->setEnabled(false);
        btnRun->setText("Hesaplanıyor...");
        progressBar->setValue(0);

        // 1. Geometri ayarları
        QString formulaStr = txtFormula->text();
        double densVal = txtDensity->text().toDouble();
        if (densVal <= 0.0) densVal = 1.0;

        gDetector->SetMaterial(formulaStr.toStdString(),
                               densVal * g/cm3);

        double thick = txtThickness->text().toDouble();
        if (!mUsingCAD) {
            if (thick <= 0.0) thick = 0.1;
            gDetector->SetThickness(thick * cm);
        }

        fRunManager->ReinitializeGeometry();

        // 2. Parçacık ayarları
        if (gGenerator) {
            gGenerator->SetParticle(cmbParticle->currentText().toStdString());
            double E = txtEnergy->text().toDouble();
            if (E <= 0.0) E = 0.1;
            gGenerator->SetEnergy(E); // MeV
        }

        // 3. Olay sayısı
        int nEvents = txtEvents->text().toInt();
        if (nEvents <= 0) nEvents = 1000;
        gSimResult.totalParticlesFired = nEvents;
        gSimResult.thicknessUsed       = mUsingCAD ? 0.0 : thick;

        fRunManager->Initialize();
        fRunManager->BeamOn(nEvents);

        // 4. Sonuçlar ve UI
        UpdateProgress(nEvents, nEvents);
        updateGraphSettings();
        generateReport();

        btnRun->setText("SİMÜLASYONU BAŞLAT");
        btnRun->setEnabled(true);
        QMessageBox::information(this, "Başarılı", "Simülasyon tamamlandı.");
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
           << " (E = " << txtEnergy->text().toStdString() << " MeV)\n\n";

        ss << "İSTATİSTİKLER:\n";
        ss << "  - Toplam Parçacık (I0): " << static_cast<int>(I0) << "\n";
        ss << "  - Geçen Parçacık (I)  : " << static_cast<int>(I) << "\n";
        ss << "  - Sönümleme Oranı     : "
           << std::fixed << std::setprecision(4)
           << (1.0 - transmission) * 100.0 << " %\n";
        ss << "  - İletim Oranı        : "
           << (transmission * 100.0) << " %\n\n";

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
    QPushButton*   btnRun;
    QPushButton*   btnLoadCad;
    QComboBox*     cmbPresets;
    QComboBox*     cmbParticle;
    QComboBox*     cmbGraphMode;
    QCheckBox*     chkLogScale;
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

    // Qt ile aynı thread'de güvenli olması için SerialOnly kullanıyoruz
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

    // 1. Dedektör
    gDetector = new CustomDetector();
    runManager->SetUserInitialization(gDetector);

    // 2. Fizik Listesi (Shielding)
    G4PhysListFactory physFactory;
    G4VUserPhysicsList* physicsList =
        physFactory.GetReferencePhysList("Shielding");
    runManager->SetUserInitialization(physicsList);

    // 3. Aksiyonlar
    runManager->SetUserInitialization(new CustomActionInit(gDetector));

    MainWindow window(runManager);
    window.show();

    int ret = app.exec();

    delete runManager;
    return ret;
}

#include "shielding_qt_sim.moc"
