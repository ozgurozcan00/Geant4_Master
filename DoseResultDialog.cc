#include "DoseResultsDialog.hh"

#include "RunAction.hh"  // RunAction pointer'ı için

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QTableWidget>
#include <QHeaderView>
#include <QPushButton>
#include <QFileDialog>
#include <QTextStream>

#include <cmath>

//----------------------------------------------------------
// Kurucular
//----------------------------------------------------------

DoseResultsDialog::DoseResultsDialog(QWidget* parent)
  : QDialog(parent)
{
  buildUi();
  fillSummaryLabels();
  fillTable();
}

DoseResultsDialog::DoseResultsDialog(double mu_cm,
                                     double hvl_mm,
                                     double tvl_mm,
                                     double totalDoseGy,
                                     double transFrac,
                                     QWidget* parent)
  : QDialog(parent),
    m_mu_cm(mu_cm),
    m_hvl_mm(hvl_mm),
    m_tvl_mm(tvl_mm),
    m_totalDoseGy(totalDoseGy),
    m_transFrac(transFrac)
{
  buildUi();
  fillSummaryLabels();
  fillTable();
}

DoseResultsDialog::DoseResultsDialog(const std::vector<G4double>& depthMm,
                                     const std::vector<G4double>& doseGy,
                                     double mu_cm,
                                     double hvl_mm,
                                     double tvl_mm,
                                     double totalDoseGy,
                                     double transFrac,
                                     QWidget* parent)
  : QDialog(parent),
    m_depthMm(depthMm),
    m_doseGy(doseGy),
    m_mu_cm(mu_cm),
    m_hvl_mm(hvl_mm),
    m_tvl_mm(tvl_mm),
    m_totalDoseGy(totalDoseGy),
    m_transFrac(transFrac)
{
  buildUi();
  fillSummaryLabels();
  fillTable();
}

/// RunAction’dan direkt alan kurucu
DoseResultsDialog::DoseResultsDialog(const RunAction* runAction,
                                     QWidget* parent)
  : QDialog(parent)
{
  if (runAction) {
    // Derinlik–doz verilerini çek
    m_depthMm = runAction->GetDepthPositions();
    m_doseGy  = runAction->GetDepthDoseGy();

    m_mu_cm       = runAction->GetMu_cm();
    m_hvl_mm      = runAction->GetHvl_mm();
    m_totalDoseGy = runAction->GetTotalDoseGy();
    m_transFrac   = runAction->GetTransFrac();

    // TVL'yi µ'dan hesapla (1/cm -> mm)
    if (m_mu_cm > 0.0) {
      // HVL = ln(2)/µ, TVL = ln(10)/µ ; µ [1/cm] -> mm için ×10
      m_tvl_mm = std::log(10.0) / m_mu_cm * 10.0;
      if (m_hvl_mm <= 0.0) {
        m_hvl_mm = std::log(2.0) / m_mu_cm * 10.0;
      }
    }
  }

  buildUi();
  fillSummaryLabels();
  fillTable();
}

DoseResultsDialog::~DoseResultsDialog() = default;

//----------------------------------------------------------
// Dışarıdan set
//----------------------------------------------------------

void DoseResultsDialog::SetResults(const std::vector<G4double>& depthMm,
                                   const std::vector<G4double>& doseGy,
                                   double mu_cm,
                                   double hvl_mm,
                                   double tvl_mm,
                                   double totalDoseGy,
                                   double transFrac)
{
  m_depthMm     = depthMm;
  m_doseGy      = doseGy;
  m_mu_cm       = mu_cm;
  m_hvl_mm      = hvl_mm;
  m_tvl_mm      = tvl_mm;
  m_totalDoseGy = totalDoseGy;
  m_transFrac   = transFrac;

  fillSummaryLabels();
  fillTable();
}

//----------------------------------------------------------
// UI kurulum
//----------------------------------------------------------

void DoseResultsDialog::buildUi()
{
  setWindowTitle("Dose / HVL / TVL Results");
  resize(700, 500);

  auto* mainLayout = new QVBoxLayout(this);

  // Üstte özet alanı
  auto* summaryLayout = new QVBoxLayout();
  m_labelMu        = new QLabel(this);
  m_labelHvl       = new QLabel(this);
  m_labelTvl       = new QLabel(this);
  m_labelTotalDose = new QLabel(this);
  m_labelTransFrac = new QLabel(this);

  summaryLayout->addWidget(m_labelMu);
  summaryLayout->addWidget(m_labelHvl);
  summaryLayout->addWidget(m_labelTvl);
  summaryLayout->addWidget(m_labelTotalDose);
  summaryLayout->addWidget(m_labelTransFrac);

  auto* summaryBox = new QWidget(this);
  summaryBox->setLayout(summaryLayout);

  // Ortada derinlik–doz tablosu
  m_table = new QTableWidget(this);
  m_table->setColumnCount(2);
  QStringList headers;
  headers << "Depth (mm)" << "Dose (Gy)";
  m_table->setHorizontalHeaderLabels(headers);
  m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  m_table->verticalHeader()->setVisible(false);

  // Altta butonlar
  auto* btnLayout = new QHBoxLayout();
  btnLayout->addStretch(1);
  m_btnSaveCsv = new QPushButton("Export to CSV", this);
  auto* btnClose = new QPushButton("Close", this);

  btnLayout->addWidget(m_btnSaveCsv);
  btnLayout->addWidget(btnClose);

  mainLayout->addWidget(summaryBox);
  mainLayout->addWidget(m_table);
  mainLayout->addLayout(btnLayout);

  setLayout(mainLayout);

  // Sinyal/slot
  connect(m_btnSaveCsv, &QPushButton::clicked,
          this, &DoseResultsDialog::onSaveToCsv);
  connect(btnClose, &QPushButton::clicked,
          this, &DoseResultsDialog::accept);
}

//----------------------------------------------------------
// Özet label doldurma
//----------------------------------------------------------

void DoseResultsDialog::fillSummaryLabels()
{
  if (!m_labelMu) return;

  m_labelMu->setText(
      QString("µ (linear attenuation) : %1 1/cm")
      .arg(m_mu_cm, 0, 'g', 5));

  m_labelHvl->setText(
      QString("HVL (1/2 value layer)  : %1 mm")
      .arg(m_hvl_mm, 0, 'g', 5));

  m_labelTvl->setText(
      QString("TVL (1/10 value layer) : %1 mm")
      .arg(m_tvl_mm, 0, 'g', 5));

  m_labelTotalDose->setText(
      QString("Total dose in scoring  : %1 Gy")
      .arg(m_totalDoseGy, 0, 'g', 5));

  m_labelTransFrac->setText(
      QString("Transmission fraction  : %1 %%")
      .arg(100.0 * m_transFrac, 0, 'f', 2));
}

//----------------------------------------------------------
// Tablo doldurma
//----------------------------------------------------------

void DoseResultsDialog::fillTable()
{
  if (!m_table) return;

  int n = std::min(m_depthMm.size(), m_doseGy.size());
  m_table->setRowCount(n);

  for (int i = 0; i < n; ++i) {
    auto* itemDepth = new QTableWidgetItem(
        QString::number(m_depthMm[i], 'f', 3));
    auto* itemDose  = new QTableWidgetItem(
        QString::number(m_doseGy[i], 'g', 5));

    m_table->setItem(i, 0, itemDepth);
    m_table->setItem(i, 1, itemDose);
  }
}

//----------------------------------------------------------
// CSV Export
//----------------------------------------------------------

void DoseResultsDialog::onSaveToCsv()
{
  QString fileName =
      QFileDialog::getSaveFileName(this,
        "Save dose results to CSV",
        QString(),
        "CSV Files (*.csv);;All Files (*)");

  if (fileName.isEmpty())
    return;

  QFile file(fileName);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    return;
  }

  QTextStream out(&file);
  out << "Depth_mm,Dose_Gy\n";

  int n = std::min(m_depthMm.size(), m_doseGy.size());
  for (int i = 0; i < n; ++i) {
    out << m_depthMm[i] << "," << m_doseGy[i] << "\n";
  }
}
