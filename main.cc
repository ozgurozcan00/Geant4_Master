//====================== main.cc ======================
#define QT_NO_KEYWORDS

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include <G4RunManagerFactory.hh>
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#include <G4VisExecutive.hh>
#include <G4PhysListFactory.hh>
#include <G4SystemOfUnits.hh>

#include <sstream>
#include <string>
#include <vector>

//---------------- Yardımcı: UI komutu çalıştır ----------------
static G4int ApplyCommandChecked(const G4String& cmd)
{
  auto uiMgr = G4UImanager::GetUIpointer();
  G4int code = uiMgr->ApplyCommand(cmd);
  if (code != 0) {
    G4cout << "[UI CMD] \"" << cmd << "\" -> return code = "
           << code << G4endl;
  }
  return code;
}

#ifdef G4UI_USE_QT

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QLineEdit>
#include <QPushButton>
#include <QMessageBox>
#include <QDialog>
#include <QPainter>
#include <QPaintEvent>
#include <QPen>
#include <QPolygon>
#include <QPoint>
#include <QString>
#include <QColor>
#include <QFileDialog>
#include <QTextStream>
#include <QCheckBox>
#include <QTabWidget>

//---------------- Renk skalası ----------------
static QColor HeatColor(double t)
{
  if (t < 0.) t = 0.;
  if (t > 1.) t = 1.;
  int r = (int)(255 * t);
  int g = 0;
  int b = (int)(255 * (1.0 - t));
  return QColor(r, g, b);
}

//======================================================================
// Depth–dose grafiği
//======================================================================
class DepthDoseDialog : public QDialog
{
public:
  DepthDoseDialog(const std::vector<G4double>& bins,
                  const std::vector<G4double>& errors,
                  G4double binWidthMm,
                  G4double thicknessMm,
                  G4double doseGy,
                  G4double doseRMSGy,
                  G4double mu_cm,
                  G4double hvl_mm,
                  G4double transFrac,
                  QWidget* parent = nullptr)
    : QDialog(parent),
      m_bins(bins),
      m_errors(errors),
      m_binWidth(binWidthMm),
      m_thicknessMm(thicknessMm),
      m_doseGy(doseGy),
      m_doseRMSGy(doseRMSGy),
      m_mu_cm(mu_cm),
      m_hvl_mm(hvl_mm),
      m_transFrac(transFrac)
  {
    setWindowTitle("Depth–dose profile (with statistical errors)");
    resize(750, 480);
  }

protected:
  void paintEvent(QPaintEvent* /*event*/) override
  {
    QPainter p(this);
    p.fillRect(rect(), Qt::white);

    if (m_bins.empty() || m_binWidth <= 0.) {
      p.setPen(Qt::black);
      p.drawText(rect().center(), "No depth–dose data");
      return;
    }

    int n = static_cast<int>(m_bins.size());
    if ((int)m_errors.size() != n) {
      m_errors.clear();
      m_errors.resize(n, 0.0);
    }

    std::vector<double> cum(n, 0.0);
    double sum = 0.;
    double maxY = 0.;
    for (int i = 0; i < n; ++i) {
      sum += m_bins[i];
      cum[i] = sum;
      double yTop = m_bins[i] + (i < (int)m_errors.size() ? m_errors[i] : 0.0);
      if (yTop > maxY) maxY = yTop;
    }
    double maxCum = (n > 0) ? cum.back() : 0.0;
    if (maxY   <= 0.) maxY   = 1.0;
    if (maxCum <= 0.) maxCum = 1.0;

    const int marginLeft   = 70;
    const int marginRight  = 220;
    const int marginTop    = 20;
    const int marginBottom = 50;

    QRect plotRect(marginLeft,
                   marginTop,
                   width() - marginLeft - marginRight,
                   height() - marginTop - marginBottom);

    p.setPen(Qt::black);
    p.drawRect(plotRect);

    // eksen yazıları
    p.drawText(plotRect.center().x(), height() - 10, "depth (mm)");

    p.save();
    p.translate(20, plotRect.center().y());
    p.rotate(-90);
    p.drawText(0, 0, "Dose (Gy)");
    p.restore();

    double maxDepth = m_binWidth * n; // mm

    // X ekseni ticks
    p.setPen(Qt::black);
    int xTicks = 5;
    for (int i = 0; i <= xTicks; ++i) {
      double frac = (double)i / xTicks;
      double xVal = frac * maxDepth;
      int x = plotRect.left() + (int)(frac * plotRect.width());
      int y1 = plotRect.bottom();
      int y2 = y1 + 5;
      p.drawLine(x, y1, x, y2);

      QString label = QString::number(xVal, 'f', 1);
      int tx = x - 15;
      int ty = y2 + 15;
      p.drawText(tx, ty, label);
    }

    // Y ekseni ticks
    int yTicks = 4;
    for (int i = 0; i <= yTicks; ++i) {
      double frac = (double)i / yTicks;
      double yVal = frac * maxY;
      int y = plotRect.bottom() - (int)(frac * plotRect.height());
      int x1 = plotRect.left() - 5;
      int x2 = plotRect.left();
      p.drawLine(x1, y, x2, y);

      QString label = QString::number(yVal, 'g', 3);
      int tx = x1 - 50;
      int ty = y + 5;
      p.drawText(tx, ty, label);
    }

    p.setRenderHint(QPainter::Antialiasing, true);

    // Depth–dose eğrisi
    if (m_bins.size() >= 2) {
      QPolygon poly;
      poly.reserve(n);

      p.setPen(QPen(Qt::blue, 2));

      for (int i = 0; i < n; ++i) {
        double x_data = (i + 0.5) * m_binWidth; // mm
        double y_data = m_bins[i];              // Gy

        double x_norm = x_data / maxDepth;
        double y_norm = y_data / maxY;

        int x = plotRect.left() +
                static_cast<int>(x_norm * plotRect.width());
        int y = plotRect.bottom() -
                static_cast<int>(y_norm * plotRect.height());

        poly << QPoint(x, y);
      }
      p.drawPolyline(poly);
    }

    // Hata çubukları
    p.setPen(QPen(Qt::darkGray, 1));
    for (int i = 0; i < n; ++i) {
      double err = m_errors[i];
      if (err <= 0.) continue;

      double x_data = (i + 0.5) * m_binWidth;
      double y_data = m_bins[i];

      double y_low  = y_data - err;
      double y_high = y_data + err;
      if (y_low < 0.) y_low = 0.;

      double x_norm = x_data / maxDepth;
      double yLow_n = y_low  / maxY;
      double yHigh_n= y_high / maxY;

      int x = plotRect.left() +
              static_cast<int>(x_norm * plotRect.width());
      int y1 = plotRect.bottom() -
               static_cast<int>(yHigh_n * plotRect.height());
      int y2 = plotRect.bottom() -
               static_cast<int>(yLow_n  * plotRect.height());

      p.drawLine(x, y1, x, y2);
      p.drawLine(x - 3, y1, x + 3, y1);
      p.drawLine(x - 3, y2, x + 3, y2);
    }

    // Kümülatif eğri (normalize)
    if (m_bins.size() >= 2) {
      QPolygon polyCum;
      polyCum.reserve(n);

      p.setPen(QPen(Qt::red, 2, Qt::DashLine));

      std::vector<double> cum2(n, 0.0);
      double sum2 = 0.;
      for (int i = 0; i < n; ++i) {
        sum2 += m_bins[i];
        cum2[i] = sum2;
      }
      double maxCum2 = (n > 0) ? cum2.back() : 1.0;
      if (maxCum2 <= 0.) maxCum2 = 1.0;

      for (int i = 0; i < n; ++i) {
        double x_data = (i + 0.5) * m_binWidth;
        double y_data = cum2[i];

        double x_norm = x_data / maxDepth;
        double y_norm = y_data / maxCum2;

        int x = plotRect.left() +
                static_cast<int>(x_norm * plotRect.width());
        int y = plotRect.bottom() -
                static_cast<int>(y_norm * plotRect.height());

        polyCum << QPoint(x, y);
      }
      p.drawPolyline(polyCum);
    }

    // Sağ panel – özet bilgiler
    QRect infoRect(plotRect.right() + 10,
                   plotRect.top(),
                   marginRight - 20,
                   plotRect.height());

    p.setPen(Qt::black);
    p.drawRect(infoRect.adjusted(0, 0, -1, -1));

    int yText = infoRect.top() + 20;
    int xText = infoRect.left() + 10;

    auto drawLine = [&](const QString& text) {
      p.drawText(xText, yText, text);
      yText += 20;
    };

    drawLine(QString("Thickness : %1 mm").arg(m_thicknessMm, 0, 'f', 2));
    drawLine(QString("Dose      : %1 Gy").arg(m_doseGy, 0, 'g', 4));
    drawLine(QString("Dose RMS  : %1 Gy").arg(m_doseRMSGy, 0, 'g', 4));
    drawLine(QString("µ         : %1 1/cm").arg(m_mu_cm, 0, 'g', 4));
    drawLine(QString("HVL       : %1 mm").arg(m_hvl_mm, 0, 'g', 4));
    drawLine(QString("Trans.    : %1 %").arg(100.0 * m_transFrac, 0, 'f', 2));
  }

private:
  std::vector<G4double> m_bins;
  mutable std::vector<G4double> m_errors;
  G4double              m_binWidth;
  G4double              m_thicknessMm;
  G4double              m_doseGy;
  G4double              m_doseRMSGy;
  G4double              m_mu_cm;
  G4double              m_hvl_mm;
  G4double              m_transFrac;
};

#endif // G4UI_USE_QT

//---------------- main ----------------
int main(int argc, char** argv)
{
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);

  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  G4PhysListFactory physFactory;
  auto physicsList = physFactory.GetReferencePhysList("FTFP_BERT");
  runManager->SetUserInitialization(physicsList);

  runManager->SetUserInitialization(new ActionInitialization(detConstruction));

  // geometriyi kur
  runManager->Initialize();

  // vis
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  auto UImanager = G4UImanager::GetUIpointer();

  if (argc == 1) {
    auto ui = new G4UIExecutive(argc, argv, "Qt");

#ifdef G4UI_USE_QT
    //-------------------- Shielding Control paneli -------------------
    QWidget* panel = new QWidget();
    panel->setWindowTitle("Shielding Control");
    panel->resize(420, 720);

    auto mainLayout = new QVBoxLayout(panel);
    QTabWidget* tabs = new QTabWidget(panel);

    // Tab 1: Geometry + GDML
    QWidget* tabGeom = new QWidget(tabs);
    auto tabGeomLayout = new QVBoxLayout(tabGeom);

    auto sampleGroup  = new QGroupBox("Sample / Material", tabGeom);
    auto sampleLayout = new QVBoxLayout(sampleGroup);
    auto sampleLabel  = new QLabel("Sample ID (0..9):", sampleGroup);
    auto sampleCombo  = new QComboBox(sampleGroup);

    sampleCombo->addItem("0 - Al5083", 0);
    sampleCombo->addItem("1 - Al5083 + 5% WC", 1);
    sampleCombo->addItem("2 - Al5083 + 10% WC", 2);
    sampleCombo->addItem("3 - Al5083 + 15% WC", 3);
    sampleCombo->addItem("4 - Al5083 + 5% B4C", 4);
    sampleCombo->addItem("5 - Al5083 + 10% B4C", 5);
    sampleCombo->addItem("6 - Al5083 + 15% B4C", 6);
    sampleCombo->addItem("7 - Al5083 + 5% Hybrid", 7);
    sampleCombo->addItem("8 - Al5083 + 10% Hybrid", 8);
    sampleCombo->addItem("9 - Al5083 + 15% Hybrid", 9);

    sampleLayout->addWidget(sampleLabel);
    sampleLayout->addWidget(sampleCombo);
    sampleGroup->setLayout(sampleLayout);

    auto geomGroup  = new QGroupBox("Geometry", tabGeom);
    auto geomLayout = new QVBoxLayout(geomGroup);
    auto thickLabel = new QLabel("Thickness (mm):", geomGroup);
    auto thickSpin  = new QDoubleSpinBox(geomGroup);
    thickSpin->setRange(0.1, 500.0);
    thickSpin->setDecimals(2);
    thickSpin->setSingleStep(0.5);
    thickSpin->setValue(detConstruction->GetShieldThickness() / mm);

    auto halfXLabel = new QLabel("Absorber half-size X (mm):", geomGroup);
    auto halfXSpin  = new QDoubleSpinBox(geomGroup);
    halfXSpin->setRange(1.0, 1000.0);
    halfXSpin->setDecimals(1);
    halfXSpin->setSingleStep(1.0);
    halfXSpin->setValue(50.0);

    auto halfYLabel = new QLabel("Absorber half-size Y (mm):", geomGroup);
    auto halfYSpin  = new QDoubleSpinBox(geomGroup);
    halfYSpin->setRange(1.0, 1000.0);
    halfYSpin->setDecimals(1);
    halfYSpin->setSingleStep(1.0);
    halfYSpin->setValue(50.0);

    geomLayout->addWidget(thickLabel);
    geomLayout->addWidget(thickSpin);
    geomLayout->addWidget(halfXLabel);
    geomLayout->addWidget(halfXSpin);
    geomLayout->addWidget(halfYLabel);
    geomLayout->addWidget(halfYSpin);
    geomGroup->setLayout(geomLayout);

    // GDML import
    auto gdmlGroup  = new QGroupBox("GDML Geometry import", tabGeom);
    auto gdmlLayout = new QVBoxLayout(gdmlGroup);
    auto gdmlPathEdit   = new QLineEdit(gdmlGroup);
    gdmlPathEdit->setPlaceholderText("No GDML file selected");
    auto gdmlBrowseButton = new QPushButton("Browse GDML...", gdmlGroup);
    auto gdmlLoadButton   = new QPushButton("Load GDML geometry", gdmlGroup);

    gdmlLayout->addWidget(gdmlPathEdit);
    gdmlLayout->addWidget(gdmlBrowseButton);
    gdmlLayout->addWidget(gdmlLoadButton);
    gdmlGroup->setLayout(gdmlLayout);

    tabGeomLayout->addWidget(sampleGroup);
    tabGeomLayout->addWidget(geomGroup);
    tabGeomLayout->addWidget(gdmlGroup);
    tabGeomLayout->addStretch(1);
    tabGeom->setLayout(tabGeomLayout);

    // Tab 2: Source & Visualization
    QWidget* tabSourceVis = new QWidget(tabs);
    auto tabSourceVisLayout = new QVBoxLayout(tabSourceVis);

    auto sourceGroup  = new QGroupBox("Source (Primary)", tabSourceVis);
    auto sourceLayout = new QVBoxLayout(sourceGroup);
    auto particleLabel = new QLabel("Particle:", sourceGroup);
    auto particleCombo = new QComboBox(sourceGroup);
    particleCombo->addItem("gamma");
    particleCombo->addItem("e-");
    auto energyLabel = new QLabel("Energy (MeV):", sourceGroup);
    auto energySpin  = new QDoubleSpinBox(sourceGroup);
    energySpin->setRange(0.01, 100.0);
    energySpin->setDecimals(3);
    energySpin->setSingleStep(0.1);
    energySpin->setValue(1.0);
    sourceLayout->addWidget(particleLabel);
    sourceLayout->addWidget(particleCombo);
    sourceLayout->addWidget(energyLabel);
    sourceLayout->addWidget(energySpin);
    sourceGroup->setLayout(sourceLayout);

    auto visGroup  = new QGroupBox("Visualization", tabSourceVis);
    auto visLayout = new QVBoxLayout(visGroup);
    auto showTracksCheck = new QCheckBox("Show trajectories", visGroup);
    showTracksCheck->setChecked(true);
    auto accumulateCheck = new QCheckBox("Accumulate events in viewer", visGroup);
    accumulateCheck->setChecked(true);
    auto trackVerbLabel  = new QLabel("Tracking verbose (0–3):", visGroup);
    auto trackVerbSpin   = new QSpinBox(visGroup);
    trackVerbSpin->setRange(0, 3);
    trackVerbSpin->setValue(0);
    visLayout->addWidget(showTracksCheck);
    visLayout->addWidget(accumulateCheck);
    visLayout->addWidget(trackVerbLabel);
    visLayout->addWidget(trackVerbSpin);
    visGroup->setLayout(visLayout);

    tabSourceVisLayout->addWidget(sourceGroup);
    tabSourceVisLayout->addWidget(visGroup);
    tabSourceVisLayout->addStretch(1);
    tabSourceVis->setLayout(tabSourceVisLayout);

    // Tab 3: Run & Analysis
    QWidget* tabRun = new QWidget(tabs);
    auto tabRunLayout = new QVBoxLayout(tabRun);

    auto runGroup  = new QGroupBox("Run & Analysis", tabRun);
    auto runLayout = new QVBoxLayout(runGroup);

    auto nEventsLabel = new QLabel("Events (beamOn):", runGroup);
    auto nEventsEdit  = new QLineEdit(runGroup);
    nEventsEdit->setText("200000");

    auto applyGeomButton   = new QPushButton("Apply geometry", runGroup);
    auto applySrcButton    = new QPushButton("Apply source", runGroup);
    auto openViewerButton  = new QPushButton("Open viewer / draw geometry", runGroup);
    auto runButton         = new QPushButton("Run (beamOn)", runGroup);
    auto plotButton        = new QPushButton("Show depth–dose plot", runGroup);

    runLayout->addWidget(nEventsLabel);
    runLayout->addWidget(nEventsEdit);
    runLayout->addWidget(applyGeomButton);
    runLayout->addWidget(applySrcButton);
    runLayout->addWidget(openViewerButton);
    runLayout->addWidget(runButton);
    runLayout->addWidget(plotButton);
    runGroup->setLayout(runLayout);

    tabRunLayout->addWidget(runGroup);
    tabRunLayout->addStretch(1);
    tabRun->setLayout(tabRunLayout);

    tabs->addTab(tabGeom,      "Geometry");
    tabs->addTab(tabSourceVis, "Source/Vis");
    tabs->addTab(tabRun,       "Run/Analysis");

    mainLayout->addWidget(tabs);
    panel->setLayout(mainLayout);
    panel->show();

    //---------------- GDML browse / load ----------------
    QObject::connect(gdmlBrowseButton, &QPushButton::clicked,
                     [gdmlPathEdit]() {
      QString fn = QFileDialog::getOpenFileName(nullptr,
                     "Select GDML file",
                     QString(),
                     "GDML Files (*.gdml);;All Files (*)");
      if (!fn.isEmpty())
        gdmlPathEdit->setText(fn);
    });

    QObject::connect(gdmlLoadButton, &QPushButton::clicked,
                     [gdmlPathEdit, detConstruction]() {
      QString path = gdmlPathEdit->text();
      if (path.isEmpty()) {
        QMessageBox::warning(nullptr, "GDML",
                             "Önce bir GDML dosyası seç.");
        return;
      }

      detConstruction->SetGDMLFile(path.toStdString());
      auto* rm = G4RunManager::GetRunManager();
      if (rm) {
        rm->GeometryHasBeenModified();
        rm->ReinitializeGeometry();
      }
      ApplyCommandChecked("/run/initialize");

      ApplyCommandChecked("/vis/drawVolume");
      ApplyCommandChecked("/vis/viewer/refresh");

      QMessageBox::information(nullptr, "GDML",
                               "GDML geometry yüklendi.");
    });

    //---------------- Geometry apply ----------------
    QObject::connect(applyGeomButton, &QPushButton::clicked,
                     [sampleCombo, thickSpin, halfXSpin, halfYSpin]() {
      int sampleId = sampleCombo->currentData().toInt();
      double thick = thickSpin->value();
      double halfX = halfXSpin->value();
      double halfY = halfYSpin->value();

      std::ostringstream cmdSample;
      cmdSample << "/shield/geom/setSample " << sampleId;
      std::ostringstream cmdThick;
      cmdThick << "/shield/geom/setThickness " << thick << " mm";
      std::ostringstream cmdSizeXY;
      cmdSizeXY << "/shield/geom/setSizeXY "
                << halfX << " mm "
                << halfY << " mm";

      ApplyCommandChecked(cmdSample.str().c_str());
      ApplyCommandChecked(cmdThick.str().c_str());
      ApplyCommandChecked(cmdSizeXY.str().c_str());

      auto* rm = G4RunManager::GetRunManager();
      if (rm) {
        rm->GeometryHasBeenModified();
        rm->ReinitializeGeometry();
      }
      ApplyCommandChecked("/run/initialize");

      ApplyCommandChecked("/vis/drawVolume");
      ApplyCommandChecked("/vis/viewer/refresh");
    });

    //---------------- Source apply ----------------
    QObject::connect(applySrcButton, &QPushButton::clicked,
                     [particleCombo, energySpin]() {
      std::string particle = particleCombo->currentText().toStdString();
      double      energy   = energySpin->value();

      std::ostringstream cmdPart;
      cmdPart << "/gun/particle " << particle;
      std::ostringstream cmdE;
      cmdE << "/gun/energy " << energy << " MeV";

      ApplyCommandChecked(cmdPart.str().c_str());
      ApplyCommandChecked(cmdE.str().c_str());
    });

    //---------------- Viewer open ----------------
    QObject::connect(openViewerButton, &QPushButton::clicked,
                     []() {
      ApplyCommandChecked("/run/initialize");
      ApplyCommandChecked("/vis/open OGL 600x600-0+0");
      ApplyCommandChecked("/vis/drawVolume");
      ApplyCommandChecked("/vis/scene/add/axes 0 0 0 10 cm");
      ApplyCommandChecked("/vis/viewer/set/autoRefresh true");
      ApplyCommandChecked("/vis/viewer/refresh");
    });

    //---------------- Run (beamOn) ----------------
    QObject::connect(runButton, &QPushButton::clicked,
                     [nEventsEdit, showTracksCheck, accumulateCheck, trackVerbSpin]() {
      bool ok = false;
      int nEv = nEventsEdit->text().toInt(&ok);
      if (!ok || nEv <= 0) nEv = 10000;

      int verbose = trackVerbSpin->value();
      {
        std::ostringstream cmd;
        cmd << "/tracking/verbose " << verbose;
        ApplyCommandChecked(cmd.str().c_str());
      }

      if (showTracksCheck->isChecked()) {
        ApplyCommandChecked("/vis/scene/add/trajectories smooth");
        ApplyCommandChecked("/vis/modeling/trajectories/create/drawByCharge");
        ApplyCommandChecked("/vis/modeling/trajectories/selectByCharge all");
      }

      if (accumulateCheck->isChecked()) {
        ApplyCommandChecked("/vis/scene/endOfEventAction accumulate");
      } else {
        ApplyCommandChecked("/vis/scene/endOfEventAction refresh");
      }

      std::ostringstream cmdRun;
      cmdRun << "/run/beamOn " << nEv;
      ApplyCommandChecked(cmdRun.str().c_str());

      ApplyCommandChecked("/vis/viewer/refresh");
    });

    //---------------- Depth–dose grafiği ----------------
    QObject::connect(plotButton, &QPushButton::clicked,
                     [detConstruction]() {
      auto* rm = G4RunManager::GetRunManager();
      const auto* runAction =
          dynamic_cast<const RunAction*>(rm->GetUserRunAction());
      if (!runAction || !runAction->HasDepthDoseData()) {
        QMessageBox::warning(nullptr, "Depth–dose plot",
                             "Henüz depth–dose verisi yok.\n"
                             "Önce bir run (beamOn) çalıştır.");
        return;
      }

      const auto& bins   = runAction->GetDepthDoseBins();
      const auto& errors = runAction->GetDepthDoseErrors();
      double      binWidth  = runAction->GetDepthBinWidth();

      double thicknessMm = detConstruction->GetShieldThickness() / mm;
      double doseGy      = runAction->GetTotalDoseGy();
      double doseRMSGy   = runAction->GetDoseRMSGy();
      double mu_cm       = runAction->GetMu_cm();
      double hvl_mm      = runAction->GetHvl_mm();
      double transFrac   = runAction->GetTransFrac();

      DepthDoseDialog dlg(bins, errors,
                          binWidth,
                          thicknessMm,
                          doseGy,
                          doseRMSGy,
                          mu_cm,
                          hvl_mm,
                          transFrac);
      dlg.exec();
    });

#endif // G4UI_USE_QT

    ui->SessionStart();
    delete ui;
  }
  else {
    G4String macroFile = argv[1];
    G4String command   = "/control/execute " + macroFile;
    ApplyCommandChecked(command);
  }

  delete visManager;
  delete runManager;

  return 0;
}
