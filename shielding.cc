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
#include <QProgressBar> // YENİ
#include <QTextStream>  // YENİ (CSV Kaydı için)
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

// Geant4 Başlıkları
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh" // YENİ
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4AnalysisManager.hh"
#include "G4GDMLParser.hh"
#include "G4PhysListFactory.hh" // YENİ: Hazır Fizik Listeleri için

// ====================================================================
// GLOBAL VERİ VE SİNYAL YÖNETİMİ
// ====================================================================
struct HitPoint {
    double x, y;
};

struct SimResult {
    std::vector<double> energies;
    std::vector<double> angles; // YENİ: Saçılma açıları (derece)
    std::vector<HitPoint> hits;
    int totalParticlesFired = 0;
    int particlesDetected = 0;
    double thicknessUsed = 0;
};
SimResult gSimResult;

// İlerleme çubuğunu güncellemek için callback arayüzü
class IProgressUpdate {
public:
    virtual void UpdateProgress(int current, int total) = 0;
};
IProgressUpdate* gProgressHandler = nullptr;

class CustomDetector;
class CustomGenerator;
CustomDetector* gDetector = nullptr;
CustomGenerator* gGenerator = nullptr;

// ====================================================================
// GRAFİK WIDGET'I (SPEKTRUM, ISI HARİTASI, AÇISAL DAĞILIM)
// ====================================================================
class ResultWidget : public QWidget {
public:
    enum Mode { Spectrum, Heatmap, Angular };

    ResultWidget(QWidget* parent = nullptr) : QWidget(parent), mMode(Spectrum), mLogScale(false) {
        setBackgroundRole(QPalette::Base);
        setAutoFillBackground(true);
        setMinimumHeight(450);
    }

    void setData(const SimResult& res, double maxE) {
        mRes = res;
        mMaxEnergy = maxE;
        update();
    }

    void setMode(Mode m) { mMode = m; update(); }
    void setLogScale(bool log) { mLogScale = log; update(); }
    
    void saveImage() {
        QString fileName = QFileDialog::getSaveFileName(this, "Grafiği Kaydet", "analiz_grafigi.png", "PNG Files (*.png)");
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
            painter.drawText(rect().center(), "Analiz Tamamlandı: Zırh %100 Başarılı (Hiçbir Radyasyon Geçmedi)");
            return;
        }
        if (mRes.totalParticlesFired == 0) {
            painter.drawText(rect().center(), "Simülasyon Bekleniyor...");
            return;
        }

        if (mMode == Spectrum) drawHistogram(painter, graphRect, mRes.energies, mMaxEnergy, "Enerji Spektrumu", "Enerji (MeV)");
        else if (mMode == Angular) drawHistogram(painter, graphRect, mRes.angles, 90.0, "Saçılma Açısı Dağılımı", "Açı (Derece)");
        else drawHeatmap(painter, graphRect);
    }

private:
    SimResult mRes;
    double mMaxEnergy = 1.0;
    Mode mMode;
    bool mLogScale;

    // Generic Histogram Çizici (Enerji veya Açı için)
    void drawHistogram(QPainter& p, QRect r, const std::vector<double>& data, double maxVal, QString title, QString xLabel) {
        p.drawText(width()/2 - 100, 30, title);
        
        if (data.empty()) return;

        int numBins = 100;
        QVector<int> bins(numBins, 0);
        double binWidth = maxVal / numBins;
        int maxCount = 0;

        for (double val : data) {
            int idx = static_cast<int>(val / binWidth);
            if (idx >= 0 && idx < numBins) bins[idx]++;
        }
        for (int c : bins) if (c > maxCount) maxCount = c;
        if (maxCount == 0) maxCount = 1;

        double barWidth = (double)r.width() / numBins;
        p.setBrush(QColor(50, 120, 220, 150));
        p.setPen(Qt::NoPen);

        for (int i = 0; i < numBins; ++i) {
            double val = bins[i];
            double hRatio;
            if (mLogScale) {
                if (val < 1) val = 0.1; 
                hRatio = std::log10(val) / std::log10(maxCount > 1 ? maxCount : 10);
                if (hRatio < 0) hRatio = 0;
            } else {
                hRatio = val / maxCount;
            }
            double barH = hRatio * r.height();
            p.drawRect(QRectF(r.left() + i * barWidth, r.bottom() - barH, barWidth, barH));
        }

        // Eksenler
        p.setPen(Qt::black);
        p.drawText(r.bottomLeft() + QPoint(0, 20), "0");
        p.drawText(r.bottomRight() + QPoint(-40, 20), QString::number(maxVal, 'f', 1));
        p.drawText(r.center().x(), r.bottom() + 35, xLabel);
        p.drawText(r.topLeft() + QPoint(-45, 10), mLogScale ? "Log" : QString::number(maxCount));
    }

    void drawHeatmap(QPainter& p, QRect r) {
        p.drawText(width()/2 - 100, 30, "Radyasyon Çıkış Haritası");
        int gridSz = 50;
        QVector<QVector<int>> grid(gridSz, QVector<int>(gridSz, 0));
        double sizeXY = 20.0 * cm; 
        int maxHits = 0;

        for (const auto& hit : mRes.hits) {
            int x = static_cast<int>((hit.x + sizeXY/2) / sizeXY * gridSz);
            int y = static_cast<int>((hit.y + sizeXY/2) / sizeXY * gridSz);
            if (x >= 0 && x < gridSz && y >= 0 && y < gridSz) {
                grid[x][y]++;
                if (grid[x][y] > maxHits) maxHits = grid[x][y];
            }
        }

        double cellW = (double)r.width() / gridSz;
        double cellH = (double)r.height() / gridSz;

        for (int i = 0; i < gridSz; ++i) {
            for (int j = 0; j < gridSz; ++j) {
                if (grid[i][j] > 0) {
                    double intensity = (double)grid[i][j] / maxHits;
                    int red = static_cast<int>(255 * intensity);
                    int blue = static_cast<int>(255 * (1.0 - intensity));
                    p.fillRect(QRectF(r.left() + i * cellW, r.bottom() - (j + 1) * cellH, cellW, cellH), QColor(red, 0, blue));
                }
            }
        }
        p.setPen(Qt::black);
        p.setBrush(Qt::NoBrush);
        p.drawRect(r);
        p.drawText(r.bottom() + 20, "XY Düzlemi (Zırhın Arkası)");
    }
};

// ====================================================================
// GEANT4 SINIFLARI
// ====================================================================

G4Material* CreateCustomMaterial(const G4String& name, const G4String& formula, double density) {
    G4NistManager* nist = G4NistManager::Instance();
    if (formula.contains("G4_")) return nist->FindOrBuildMaterial(formula);

    std::vector<G4String> tokens;
    std::istringstream iss(formula);
    G4String token;
    while (iss >> token) tokens.push_back(token);

    if (tokens.size() % 2 != 0) return nist->FindOrBuildMaterial("G4_AIR");

    G4Material* mat = new G4Material(name, density, tokens.size() / 2);
    for (size_t i = 0; i < tokens.size(); i += 2) {
        G4Element* element = nist->FindOrBuildElement(tokens[i]);
        if (element) mat->AddElement(element, std::stoi(tokens[i+1]));
    }
    return mat;
}

class CustomDetector : public G4VUserDetectorConstruction {
public:
    CustomDetector() 
        : fScoringVolume(nullptr), fCadFileName(""), fUseCAD(false), 
          fMatFormula("G4_Pb"), fDensity(11.34*g/cm3), fThickness(2.5*cm) {}
    
    void SetCADFile(G4String path) { fCadFileName = path; fUseCAD = true; }
    void SetMaterial(G4String formula, double dens) { fMatFormula = formula; fDensity = dens; }
    void SetThickness(double t) { fThickness = t; fUseCAD = false; }
    
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume* Construct() {
        G4NistManager* nist = G4NistManager::Instance();
        G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
        G4Material* shieldMat = CreateCustomMaterial("CustomShield", fMatFormula, fDensity);

        G4Box* solidWorld = new G4Box("World", 0.5*m, 0.5*m, 0.5*m);
        G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, worldMat, "World");
        G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);

        if (fUseCAD && fCadFileName != "") {
            G4GDMLParser parser;
            parser.Read(fCadFileName, false);
            G4LogicalVolume* cadLogic = parser.GetVolume("Shield");
            if (!cadLogic) cadLogic = parser.GetWorldVolume()->GetLogicalVolume();
            if (cadLogic) {
                cadLogic->SetMaterial(shieldMat);
                new G4PVPlacement(0, G4ThreeVector(0,0,0), cadLogic, "CAD_Shield", logicWorld, false, 0, true);
                fScoringVolume = cadLogic;
            }
        } else {
            G4Box* solidShield = new G4Box("Shield", 10*cm, 10*cm, fThickness/2.0);
            G4LogicalVolume* logicShield = new G4LogicalVolume(solidShield, shieldMat, "Shield");
            new G4PVPlacement(0, G4ThreeVector(0,0,0), logicShield, "Shield", logicWorld, false, 0, true);
            fScoringVolume = logicShield;
        }
        return physWorld;
    }

private:
    G4LogicalVolume* fScoringVolume;
    G4String fCadFileName;
    G4bool fUseCAD;
    G4String fMatFormula;
    double fDensity;
    double fThickness;
};

class CustomGenerator : public G4VUserPrimaryGeneratorAction {
public:
    CustomGenerator() {
        fParticleGun = new G4ParticleGun(1);
        SetParticle("gamma");
        fParticleGun->SetParticleEnergy(0.662 * MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -20.0*cm));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }
    ~CustomGenerator() { delete fParticleGun; }
    
    void SetParticle(G4String name) {
        G4ParticleDefinition* p = G4ParticleTable::GetParticleTable()->FindParticle(name);
        if (p) fParticleGun->SetParticleDefinition(p);
    }
    void SetEnergy(double energyMeV) { fParticleGun->SetParticleEnergy(energyMeV * MeV); }
    void GeneratePrimaries(G4Event* anEvent) { fParticleGun->GeneratePrimaryVertex(anEvent); }

private:
    G4ParticleGun* fParticleGun;
};

// --- Event Action (İlerleme Çubuğu İçin) ---
class CustomEventAction : public G4UserEventAction {
public:
    void BeginOfEventAction(const G4Event* evt) override {
        int eventID = evt->GetEventID();
        // Her 100 olayda bir UI güncellemesi yap (Sık yaparsak yavaşlar)
        if (gProgressHandler && eventID % 100 == 0) {
            int nEvents = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
            gProgressHandler->UpdateProgress(eventID, nEvents);
        }
    }
};

class CustomSteppingAction : public G4UserSteppingAction {
public:
    CustomSteppingAction(const CustomDetector* det) : fDet(det) {}
    void UserSteppingAction(const G4Step* step) override {
        if (!fDet->GetScoringVolume()) return;
        auto vol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
        auto nextVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        
        if (vol && nextVol && vol->GetLogicalVolume() == fDet->GetScoringVolume() && 
            nextVol->GetLogicalVolume()->GetName() == "World") {
             G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
             if (pos.z() > 0) {
                double energy = step->GetTrack()->GetKineticEnergy();
                
                // Açı Hesaplama: Z ekseni ile yapılan açı
                G4ThreeVector direction = step->GetTrack()->GetMomentumDirection();
                double angle = direction.angle(G4ThreeVector(0,0,1)) * 180.0 / 3.14159;

                gSimResult.energies.push_back(energy);
                gSimResult.angles.push_back(angle);
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
        gSimResult.particlesDetected = 0;
        gSimResult.totalParticlesFired = 0;
    }
    void EndOfRunAction(const G4Run*) override { }
};

class CustomActionInit : public G4VUserActionInitialization {
public:
    CustomActionInit(CustomDetector* det) : fDet(det) {}
    void BuildForMaster() const override { SetUserAction(new CustomRunAction); }
    void Build() const override {
        SetUserAction(new CustomGenerator);
        SetUserAction(new CustomRunAction);
        SetUserAction(new CustomEventAction); // Progress Bar için eklendi
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
    MainWindow(G4RunManager* runMgr) : fRunManager(runMgr) {
        setWindowTitle("Profesyonel Zırhlama Analiz Laboratuvarı");
        resize(1100, 750);
        gProgressHandler = this; // Callback'i bağla

        setStyleSheet(
            "QGroupBox { font-weight: bold; border: 1px solid #ccc; border-radius: 6px; margin-top: 10px; background-color: #f9f9f9; } "
            "QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; background-color: #f9f9f9; } "
            "QPushButton { border-radius: 4px; padding: 6px; }"
        );

        auto mainLayout = new QHBoxLayout(this);

        // --- SOL PANEL ---
        QWidget* leftPanel = new QWidget();
        leftPanel->setFixedWidth(350);
        auto leftLayout = new QVBoxLayout(leftPanel);

        // 1. Zırh Özellikleri
        auto grpGeo = new QGroupBox("1. Zırh Özellikleri", this);
        auto layGeo = new QFormLayout(grpGeo);
        
        cmbPresets = new QComboBox(this);
        cmbPresets->addItems({"-- Özel Malzeme --", "Kurşun (Pb)", "Beton (Standart)", "Su", "Demir (Fe)", "Alüminyum", "Borlu Polietilen"});
        connect(cmbPresets, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onPresetChanged);
        
        txtFormula = new QLineEdit("G4_Pb", this);
        txtDensity = new QLineEdit("11.34", this);
        txtThickness = new QLineEdit("2.5", this);
        btnLoadCad = new QPushButton("CAD Geometrisi Yükle (.gdml)", this);
        connect(btnLoadCad, &QPushButton::clicked, this, &MainWindow::loadCadFile);

        layGeo->addRow("Hazır Malzeme:", cmbPresets);
        layGeo->addRow("Kimyasal Formül:", txtFormula);
        layGeo->addRow("Yoğunluk (g/cm3):", txtDensity);
        layGeo->addRow("Kalınlık (cm):", txtThickness);
        layGeo->addRow(btnLoadCad);

        // 2. Radyasyon Kaynağı
        auto grpRad = new QGroupBox("2. Radyasyon Kaynağı", this);
        auto layRad = new QFormLayout(grpRad);
        cmbParticle = new QComboBox(this);
        cmbParticle->addItems({"gamma", "e-", "alpha", "neutron", "proton"});
        txtEnergy = new QLineEdit("0.662", this); 
        txtEvents = new QLineEdit("10000", this);
        layRad->addRow("Parçacık Tipi:", cmbParticle);
        layRad->addRow("Enerji (MeV):", txtEnergy);
        layRad->addRow("Olay Sayısı:", txtEvents);

        // Kontrol Butonları
        btnRun = new QPushButton("SİMÜLASYONU BAŞLAT", this);
        btnRun->setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; padding: 15px; font-size: 14px;");
        connect(btnRun, &QPushButton::clicked, this, &MainWindow::runSimulation);
        
        progressBar = new QProgressBar(this);
        progressBar->setValue(0);
        progressBar->setAlignment(Qt::AlignCenter);

        leftLayout->addWidget(grpGeo);
        leftLayout->addWidget(grpRad);
        leftLayout->addStretch();
        leftLayout->addWidget(progressBar);
        leftLayout->addWidget(btnRun);

        // --- SAĞ PANEL ---
        auto tabs = new QTabWidget(this);
        
        // Sekme 1: Grafikler
        QWidget* tabGraph = new QWidget();
        auto layGraph = new QVBoxLayout(tabGraph);
        
        auto toolBar = new QHBoxLayout();
        chkLogScale = new QCheckBox("Logaritmik Ölçek", this);
        connect(chkLogScale, &QCheckBox::toggled, this, &MainWindow::updateGraphSettings);
        
        cmbGraphMode = new QComboBox(this);
        cmbGraphMode->addItem("Enerji Spektrumu (Energy Spectrum)");
        cmbGraphMode->addItem("Isı Haritası (Heatmap)");
        cmbGraphMode->addItem("Açısal Dağılım (Angular Distribution)");
        connect(cmbGraphMode, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::updateGraphSettings);
        
        QPushButton* btnSaveImg = new QPushButton("Resim Kaydet", this);
        connect(btnSaveImg, &QPushButton::clicked, this, &MainWindow::savePlot);

        QPushButton* btnSaveCsv = new QPushButton("CSV Olarak Kaydet", this);
        connect(btnSaveCsv, &QPushButton::clicked, this, &MainWindow::saveCsv);

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
        auto layCalc = new QVBoxLayout(tabCalc);
        txtReport = new QLabel("Analiz sonuçları burada görünecek...", this);
        txtReport->setAlignment(Qt::AlignTop);
        txtReport->setStyleSheet("font-family: Monospace; font-size: 13px; background-color: white; padding: 15px; border: 1px solid #ddd;");
        txtReport->setWordWrap(true);
        layCalc->addWidget(txtReport);

        tabs->addTab(tabGraph, "Görsel Analiz");
        tabs->addTab(tabCalc, "Fiziksel Rapor & HVL");

        mainLayout->addWidget(leftPanel);
        mainLayout->addWidget(tabs);
    }

    // IProgressUpdate Arayüzü Implementasyonu
    void UpdateProgress(int current, int total) override {
        progressBar->setMaximum(total);
        progressBar->setValue(current);
        QApplication::processEvents(); // Arayüzün donmasını engeller
    }

public slots:
    void onPresetChanged(int index) {
        // Hazır malzemeler
        if (index == 1) { setMat("G4_Pb", "11.34"); } // Kurşun
        else if (index == 2) { setMat("G4_CONCRETE", "2.3"); } // Beton
        else if (index == 3) { setMat("G4_WATER", "1.0"); } // Su
        else if (index == 4) { setMat("G4_Fe", "7.87"); } // Demir
        else if (index == 5) { setMat("G4_Al", "2.70"); } // Alüminyum
        else if (index == 6) { setMat("H 2 C 1", "0.93"); } // Basitleştirilmiş Polietilen (Kabaca)
    }

    void setMat(QString f, QString d) {
        txtFormula->setText(f);
        txtDensity->setText(d);
    }

    void loadCadFile() {
        QString fileName = QFileDialog::getOpenFileName(this, "CAD Yükle", "", "GDML (*.gdml)");
        if (!fileName.isEmpty()) {
            btnLoadCad->setText("Yüklendi: " + QFileInfo(fileName).fileName());
            gDetector->SetCADFile(fileName.toStdString());
            txtThickness->setEnabled(false); 
        }
    }

    void runSimulation() {
        btnRun->setEnabled(false);
        btnRun->setText("Hesaplanıyor...");
        progressBar->setValue(0);
        
        // 1. Geometri
        double thick = txtThickness->text().toDouble();
        gDetector->SetThickness(thick * cm);
        gDetector->SetMaterial(txtFormula->text().toStdString(), txtDensity->text().toDouble() * g/cm3);
        fRunManager->ReinitializeGeometry(); 

        // 2. Parçacık
        if (gGenerator) {
            gGenerator->SetParticle(cmbParticle->currentText().toStdString());
            gGenerator->SetEnergy(txtEnergy->text().toDouble());
        }

        // 3. Çalıştır
        int nEvents = txtEvents->text().toInt();
        gSimResult.totalParticlesFired = nEvents;
        gSimResult.thicknessUsed = thick;

        fRunManager->Initialize();
        fRunManager->BeamOn(nEvents);

        // 4. Sonuçlar
        UpdateProgress(nEvents, nEvents);
        updateGraphSettings();
        generateReport();

        btnRun->setText("SİMÜLASYONU BAŞLAT");
        btnRun->setEnabled(true);
        QMessageBox::information(this, "Başarılı", "Simülasyon tamamlandı.");
    }

    void updateGraphSettings() {
        resultWidget->setLogScale(chkLogScale->isChecked());
        if (cmbGraphMode->currentIndex() == 0) resultWidget->setMode(ResultWidget::Spectrum);
        else if (cmbGraphMode->currentIndex() == 1) resultWidget->setMode(ResultWidget::Heatmap);
        else resultWidget->setMode(ResultWidget::Angular);
        
        double maxE = txtEnergy->text().toDouble();
        resultWidget->setData(gSimResult, maxE);
    }

    void savePlot() { resultWidget->saveImage(); }

    void saveCsv() {
        QString fileName = QFileDialog::getSaveFileName(this, "CSV Kaydet", "analiz_sonuclari.csv", "CSV Files (*.csv)");
        if (fileName.isEmpty()) return;
        
        QFile file(fileName);
        if (file.open(QIODevice::WriteOnly)) {
            QTextStream stream(&file);
            stream << "Enerji(MeV),Aci(Derece),PozisyonX(mm),PozisyonY(mm)\n";
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
        double I = gSimResult.particlesDetected;
        double x = gSimResult.thicknessUsed; 
        double transmission = (I0 > 0) ? (I / I0) : 0;
        
        std::stringstream ss;
        ss << "=== DETAYLI ZIRHLAMA ANALİZ RAPORU ===\n\n";
        ss << "PARAMETRELER:\n";
        ss << "  - Malzeme: " << txtFormula->text().toStdString() << " (" << txtDensity->text().toStdString() << " g/cm3)\n";
        ss << "  - Kalınlık: " << x << " cm\n";
        ss << "  - Parçacık: " << cmbParticle->currentText().toStdString() << " (" << txtEnergy->text().toStdString() << " MeV)\n\n";
        
        ss << "İSTATİSTİKLER:\n";
        ss << "  - Toplam Parçacık (I0): " << (int)I0 << "\n";
        ss << "  - Geçen Parçacık (I)  : " << (int)I << "\n";
        ss << "  - Sönümleme Oranı     : %" << std::fixed << std::setprecision(4) << (1.0 - transmission) * 100.0 << "\n";
        ss << "  - İletim Oranı        : %" << (transmission * 100.0) << "\n\n";

        if (transmission > 0 && transmission < 1) {
            double mu = -std::log(transmission) / x;
            double hvl = std::log(2) / mu;
            double tvl = std::log(10) / mu;
            
            ss << "HESAPLANAN FİZİKSEL DEĞERLER:\n";
            ss << "  - Lineer Zayıflatma Katsayısı (μ): " << mu << " /cm\n";
            ss << "  - Yarı Değer Kalınlığı (HVL)     : " << hvl << " cm\n";
            ss << "  - 1/10 Değer Kalınlığı (TVL)     : " << tvl << " cm\n";
            ss << "    (Not: Bu değerler simülasyon istatistiğine dayalıdır)\n";
        } else if (transmission == 0) {
            ss << "SONUÇ: TAM SÖNÜMLEME.\n";
            ss << "Radyasyon zırhı geçemedi. Daha hassas ölçüm için olay sayısını artırın.\n";
        }

        txtReport->setText(QString::fromStdString(ss.str()));
    }

private:
    G4RunManager* fRunManager;
    QLineEdit *txtFormula, *txtDensity, *txtThickness, *txtEnergy, *txtEvents;
    QPushButton *btnRun, *btnLoadCad;
    QComboBox *cmbPresets, *cmbParticle, *cmbGraphMode;
    QCheckBox *chkLogScale;
    QProgressBar *progressBar;
    ResultWidget* resultWidget;
    QLabel *txtReport;
};

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    
    // Geant4 Run Manager
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    
    // 1. Dedektör
    gDetector = new CustomDetector();
    runManager->SetUserInitialization(gDetector);
    
    // 2. Fizik Listesi (PRO DEĞİŞİKLİK)
    // "Shielding" fizik listesi Geant4 fabrikasından yüklenir.
    // Bu liste, nükleer koruma hesapları için optimize edilmiştir.
    G4PhysListFactory physFactory;
    G4VUserPhysicsList* physicsList = physFactory.GetReferencePhysList("Shielding");
    runManager->SetUserInitialization(physicsList);

    // 3. Eylemler
    runManager->SetUserInitialization(new CustomActionInit(gDetector));

    MainWindow window(runManager);
    window.show();

    int ret = app.exec();
    delete runManager;
    return ret;
}

#include "shielding_qt_sim.moc"
