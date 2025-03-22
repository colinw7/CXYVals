#include <CQXYVals.h>
#include <CXYVals.h>

#include <QApplication>
#include <QPainter>
#include <QLabel>
#include <QHBoxLayout>
#include <QMouseEvent>
#include <QKeyEvent>

int
main(int argc, char **argv)
{
  QApplication app(argc, argv);

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);

    if (arg == "-test") {
      bool rc = CXYValsInside::unitTest();

      exit(rc);
    }
  }

  //---

  CQXYValsTest test;

  test.resize(800, 800);

  test.show();

  return app.exec();
}

CQXYValsTest::
CQXYValsTest()
{
  setObjectName("test");

  QHBoxLayout *layout = new QHBoxLayout(this);
  layout->setMargin(0); layout->setSpacing(0);

  canvas_ = new CQXYValsCanvas(this);

  layout->addWidget(canvas_);
}

//-----------

CQXYValsCanvas::
CQXYValsCanvas(CQXYValsTest *test) :
 QWidget(nullptr), test_(test)
{
  setObjectName("canvas");

  setFocusPolicy(Qt::StrongFocus);

  //setMouseTracking(true);
}

void
CQXYValsCanvas::
paintEvent(QPaintEvent *)
{
  static QColor colors[] = {
    // blue
    QColor(0x31,0x82,0xBD),
    QColor(0x6B,0xAE,0xD6),
    QColor(0x9E,0xCA,0xE1),
    QColor(0xC6,0xDB,0xEF),

    // orange
    QColor(0xE6,0x55,0x0D),
    QColor(0xFD,0x8D,0x3C),
    QColor(0xFD,0xAE,0x6B),
    QColor(0xFD,0xD0,0xA2),

    // green
    QColor(0x31,0xA3,0x54),
    QColor(0x74,0xC4,0x76),
    QColor(0xA1,0xD9,0x9B),
    QColor(0xC7,0xE9,0xC0),

    // purple
    QColor(0x75,0x6B,0xB1),
    QColor(0x9E,0x9A,0xC8),
    QColor(0xBC,0xBD,0xDC),
    QColor(0xDA,0xDA,0xEB),

    // gray
    QColor(0x63,0x63,0x63),
    QColor(0x96,0x96,0x96),
    QColor(0xBD,0xBD,0xBD),
    QColor(0xD9,0xD9,0xD9),
  };

  static int numColors = 20;

  //---

  QPainter painter(this);

  painter.fillRect(rect(), Qt::white);

  if (pressed_) {
    int x1 = std::min(pressPos_.x(), releasePos_.x());
    int y1 = std::min(pressPos_.y(), releasePos_.y());
    int x2 = std::max(pressPos_.x(), releasePos_.x());
    int y2 = std::max(pressPos_.y(), releasePos_.y());

    int w = x2 - x1;
    int h = y2 - y1;

    if (w > 0 && h > 0) {
      QRect r(x1, y1, w, h);

      painter.setPen(Qt::red);

      painter.drawRect(r);
    }
  }

  //---

  const std::vector<double> &xv = valueData_.xyvals.xvals();
  const std::vector<double> &yv = valueData_.xyvals.yvals();

  painter.setPen(Qt::black);

  for (uint iy = 1; iy < yv.size(); ++iy) {
    int y1 = yv[iy - 1];
    int y2 = yv[iy    ];

    int h = y2 - y1;

    for (uint ix = 1; ix < xv.size(); ++ix) {
      int x1 = xv[ix - 1];
      int x2 = xv[ix    ];

      int w = x2 - x1;

      QRect r(x1, y1, w, h);

      CXYValsInside::InsideValue val = valueData_.xyvals.insideVal(ix - 1, iy - 1);

      if (val > 0) {
        QColor fc = colors[(val - 1) % numColors];

        painter.fillRect(r, fc);
      }

      painter.setPen(Qt::black);

      painter.drawRect(r);
    }
  }

  //---

  for (uint i = 0; i < valueData_.oqpolygons.size(); ++i) {
    QPen pen(Qt::green);
    pen.setWidth(3);

    painter.setPen(pen);

    painter.drawPolygon(valueData_.oqpolygons[i]);
  }

  QFontMetrics fm(font());

  for (uint i = 0; i < valueData_.ipolygons.size(); ++i) {
    const CXYVals::Polygon &polygon = valueData_.ipolygons[i];

    double                     area = valueData_.xyvals.polygonArea (polygon);
    CXYValsInside::InsideValue val  = valueData_.xyvals.polygonValue(polygon);

    double xc, yc;

    polygon.centroid(xc, yc);

    QString text = QString("%1 (%2)").arg(val).arg(area);

    int tw = fm.horizontalAdvance(text);

    painter.setPen(Qt::black);

    painter.drawText(xc - tw/2, yc, text);
  }
}

void
CQXYValsCanvas::
mousePressEvent(QMouseEvent *me)
{
  pressed_    = true;
  pressPos_   = me->pos();
  releasePos_ = pressPos_;

  update();
}

void
CQXYValsCanvas::
mouseMoveEvent(QMouseEvent *me)
{
  releasePos_ = me->pos();

  update();
}

void
CQXYValsCanvas::
mouseReleaseEvent(QMouseEvent *me)
{
  pressed_    = false;
  releasePos_ = me->pos();

  double x1 = pressPos_  .x();
  double y1 = pressPos_  .y();
  double x2 = releasePos_.x();
  double y2 = releasePos_.y();

  //---

  std::vector<double> x { x1, x2, x2, x1};
  std::vector<double> y { y1, y1, y2, y2};

  CXYVals::Polygon poly(x, y);

  valueData_.ipolygons .push_back(poly);
  valueData_.iqpolygons.push_back(toQPolygon(poly));

  //---

  updatePolygons();

  update();
}

void
CQXYValsCanvas::
updatePolygons()
{
  valueData_.xyvals.initValues(valueData_.ipolygons);

  valueData_.opolygons .clear();
  valueData_.oqpolygons.clear();

  if (valueData_.xyvals.getPolygons(valueData_.opolygons)) {
    for (uint j = 0; j < valueData_.opolygons.size(); ++j) {
      const CXYValsInside::Polygon &polygon = valueData_.opolygons[j];

      valueData_.oqpolygons.push_back(toQPolygon(polygon));
    }
  }

  valueData1_ = valueData_;
}

void
CQXYValsCanvas::
getPolygons()
{
  valueData_.ipolygons .clear();
  valueData_.iqpolygons.clear();

  if (valueData_.xyvals.getPolygons(valueData_.ipolygons)) {
    for (uint j = 0; j < valueData_.ipolygons.size(); ++j) {
      const CXYValsInside::Polygon &polygon = valueData_.ipolygons[j];

      valueData_.iqpolygons.push_back(toQPolygon(polygon));
    }
  }

  valueData_.opolygons  = valueData_.ipolygons;
  valueData_.oqpolygons = valueData_.iqpolygons;
}

QPolygon
CQXYValsCanvas::
toQPolygon(const CXYVals::Polygon &polygon) const
{
  QPolygon poly;

  for (int i = 0; i < polygon.size(); ++i)
    poly.push_back(QPoint(polygon.x[i], polygon.y[i]));

  poly.push_back(QPoint(polygon.x[0], polygon.y[0]));

  return poly;
}

void
CQXYValsCanvas::
keyPressEvent(QKeyEvent *ke)
{
  if      (ke->key() == Qt::Key_Delete) {
    // clear grid
    valueData_.ipolygons .clear();
    valueData_.iqpolygons.clear();

    updatePolygons();

    update();
  }
  else if (ke->key() == Qt::Key_C) {
    valueData1_ = valueData_;

    // combine connected cells
    valueData_.xyvals.combineInside();

    valueData_.xyvals.setPolygonInsideValues();

    getPolygons();

    update();
  }
  else if (ke->key() == Qt::Key_F) {
    valueData1_ = valueData_;

    // fill
    valueData_.xyvals.fill();

    getPolygons();

    update();
  }
  else if (ke->key() == Qt::Key_H) {
    valueData1_ = valueData_;

    // fill disconnected rows
    valueData_.xyvals.fillDisconnectedRows();

    valueData_.xyvals.setPolygonInsideValues();

    getPolygons();

    update();
  }
  else if (ke->key() == Qt::Key_V) {
    valueData1_ = valueData_;

    // fill disconnected columns
    valueData_.xyvals.fillDisconnectedColumns();

    valueData_.xyvals.setPolygonInsideValues();

    getPolygons();

    update();
  }
  else if (ke->key() == Qt::Key_S) {
    // display smallest
    CXYValsInside::InsideValue i = valueData_.xyvals.findSmallest();

    std::cerr << i << "(" << valueData_.xyvals.valueArea(i) << ")" << std::endl;
  }
  else if (ke->key() == Qt::Key_P) {
    // print
    valueData_.xyvals.print(std::cerr); std::cerr << std::endl;
  }
  else if (ke->key() == Qt::Key_R) {
    // replace input polygons with extracted polygons
    valueData_.ipolygons  = valueData_.opolygons;
    valueData_.iqpolygons = valueData_.oqpolygons;

    updatePolygons();

    update();
  }
  else if (ke->key() == Qt::Key_U) {
    valueData_ = valueData1_;

    update();
  }
}
