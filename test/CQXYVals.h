#include <QWidget>
#include <CXYVals.h>

class CQXYValsTest;

class CQXYValsCanvas : public QWidget {
  Q_OBJECT

 public:
  CQXYValsCanvas(CQXYValsTest *text);

 private:
  void paintEvent(QPaintEvent *) override;

  void mousePressEvent  (QMouseEvent *me) override;
  void mouseMoveEvent   (QMouseEvent *me) override;
  void mouseReleaseEvent(QMouseEvent *me) override;

  void keyPressEvent(QKeyEvent *) override;

  void updatePolygons();
  void getPolygons();

  QPolygon toQPolygon(const CXYVals::Polygon &polygon) const;

 private:
  typedef std::vector<QPolygon> Polygons;

  struct ValueData {
    CXYValsInside           xyvals;
    CXYValsInside::Polygons ipolygons;
    CXYValsInside::Polygons opolygons;
    Polygons                iqpolygons;
    Polygons                oqpolygons;
  };

  CQXYValsTest* test_;
  ValueData     valueData_;
  ValueData     valueData1_;
  bool          pressed_ { false };
  QPoint        pressPos_;
  QPoint        releasePos_;
};

class CQXYValsTest : public QWidget {
  Q_OBJECT

 public:
  CQXYValsTest();

 private:
  CQXYValsCanvas *canvas_;
};
