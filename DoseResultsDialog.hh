#ifndef DoseResultsDialog_hh
#define DoseResultsDialog_hh 1

#include <QDialog>
#include <vector>

#include "globals.hh"  // G4double için

class QLabel;
class QTableWidget;
class QPushButton;
class RunAction;  // forward declaration

/// HVL–TVL–dose özetini ve (varsa) depth–dose verisini gösteren Qt diyaloğu
class DoseResultsDialog : public QDialog
{
  Q_OBJECT

public:
  /// Sadece sayısal sonuçlarla kullanılan basit kurucu
  explicit DoseResultsDialog(double mu_cm,
                             double hvl_mm,
                             double tvl_mm,
                             double totalDoseGy,
                             double transFrac,
                             QWidget* parent = nullptr);

  /// Derinlik–doz profilini de göstermek istersen, bu kurucuyu kullan
  DoseResultsDialog(const std::vector<G4double>& depthMm,
                    const std::vector<G4double>& doseGy,
                    double mu_cm,
                    double hvl_mm,
                    double tvl_mm,
                    double totalDoseGy,
                    double transFrac,
                    QWidget* parent = nullptr);

  /// En genel, boş kurucu
  explicit DoseResultsDialog(QWidget* parent = nullptr);

  /// RunAction’dan direkt okuyan kurucu
  explicit DoseResultsDialog(const RunAction* runAction,
                             QWidget* parent = nullptr);

  ~DoseResultsDialog() override;

  /// Sonradan set etmek istersen
  void SetResults(const std::vector<G4double>& depthMm,
                  const std::vector<G4double>& doseGy,
                  double mu_cm,
                  double hvl_mm,
                  double tvl_mm,
                  double totalDoseGy,
                  double transFrac);

private slots:
  void onSaveToCsv();

private:
  void buildUi();
  void fillSummaryLabels();
  void fillTable();

private:
  // Veriler
  std::vector<G4double> m_depthMm;
  std::vector<G4double> m_doseGy;

  double m_mu_cm       = 0.0;
  double m_hvl_mm      = 0.0;
  double m_tvl_mm      = 0.0;
  double m_totalDoseGy = 0.0;
  double m_transFrac   = 0.0;

  // Qt bileşenleri
  QLabel*        m_labelMu        = nullptr;
  QLabel*        m_labelHvl       = nullptr;
  QLabel*        m_labelTvl       = nullptr;
  QLabel*        m_labelTotalDose = nullptr;
  QLabel*        m_labelTransFrac = nullptr;
  QTableWidget*  m_table          = nullptr;
  QPushButton*   m_btnSaveCsv     = nullptr;
};

#endif // DoseResultsDialog_hh
