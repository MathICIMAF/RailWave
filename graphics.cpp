#include "graphics.h"

graphics::graphics(QObject *parent) :
    QObject(parent)
{
}

graphics::graphics(QList<QString> materialInfo, FMatrix pF, FMatrix pVf, bool modesShape, QList<int> pIllRows,QCustomPlot *cPlot):QObject()
{
    X = pF;
    Y = pVf;
    illRows = pIllRows;

    measureX = "KHz";
    measureY = "mm/us";
    xAxisLabel = "Frequency";
    yAxisLabel = "Phase Velocity";
    X.multi_escalar_this(1e-3);
    Y.multi_escalar_this(1e-3);
    GetMaxMinCoords(X, Y);
    this->customPlot = cPlot;
    InitCustomPlot("", autoScaleXMax, autoScaleYMax);
    InitCurvesData(pF, pVf, NULL, -1);

}

void graphics::GetMaxMinCoords(FMatrix X, FMatrix Y)
{
    QList<double> sortedX, sortedY;
    for(int i = 0;i < X.rows; i++)
    {
        for(int j = 0; j < X.cols; j++)
        {
            if(!sortedX.contains(X.GetElem(i,j)))
                InsertInOrder(X.GetElem(i,j), sortedX);
            if(!sortedY.contains(Y.GetElem(i,j)))
                InsertInOrder(Y.GetElem(i,j), sortedY);
        }
    }
    autoScaleXMin = sortedX.at(0);
    autoScaleYMin = sortedY.at(0);
    autoScaleXMax = sortedX.at(sortedX.size()-1);
    autoScaleYMax = sortedY.at(floor((sortedY.size()-1)*0.9));
}

void graphics::InsertInOrder(double elem, QList<double> &list)
{
    for(int i = 0; i < list.size(); i++)
    {
        if(list[i] > elem)
        {
            list.insert(i,elem);
            return;
        }
    }
    list.push_back(elem);
}

void graphics::InitCustomPlot(QString title, int maxCoordX, int maxCoordY)
{
    //InitAxis
    customPlot->setInteractions(QCustomPlot::iRangeDrag | QCustomPlot::iRangeZoom | QCustomPlot::iSelectAxes | QCustomPlot::iSelectPlottables);
    customPlot->setRangeDrag(Qt::Horizontal|Qt::Vertical);
    customPlot->setRangeZoom(Qt::Horizontal|Qt::Vertical);
    customPlot->setupFullAxesBox();
    customPlot->setTitle(title);
    customPlot->xAxis->setLabel(xAxisLabel + " " + measureX);
    customPlot->yAxis->setLabel(yAxisLabel + " " + measureY);
    customPlot->xAxis->setRange(0, maxCoordX);
    customPlot->yAxis->setRange(0, maxCoordY);

    //Init Title
    QFont font;
    font.setPixelSize(12);
    customPlot->setTitleFont(font);

    //Init legend
    customPlot->legend->setVisible(false);

    //Signals+Slots
    // connect slot that ties some axis selections together (especially opposite axes):
    connect(customPlot, SIGNAL(selectionChangedByUser()), this, SLOT(selectionChanged()));
    connect(customPlot,SIGNAL(plottableClick(QCPAbstractPlottable*,QMouseEvent*)), this, SLOT(NewSelected(QCPAbstractPlottable*))); //?
    // connect slots that takes care that when an axis is selected, only that direction can be dragged and zoomed:
    connect(customPlot, SIGNAL(mousePress(QMouseEvent*)), this, SLOT(mousePress(QMouseEvent*)));
    connect(customPlot, SIGNAL(mouseWheel(QWheelEvent*)), this, SLOT(mouseWheel(QWheelEvent*)));
    // make bottom and left axes transfer their ranges to top and right axes:
    connect(customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->xAxis2, SLOT(setRange(QCPRange)));
    connect(customPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->yAxis2, SLOT(setRange(QCPRange)));
    // connect slot that shows a message in the status bar when a graph is clicked:
    connect(customPlot, SIGNAL(plottableClick(QCPAbstractPlottable*,QMouseEvent*)), this, SLOT(graphClicked(QCPAbstractPlottable*)));
    connect(customPlot, SIGNAL(mouseMove(QMouseEvent*)), this, SLOT(graphMouseMove(QMouseEvent*)));
}

void graphics::selectionChanged()
{
  /*
   normally, axis base line, axis tick labels and axis labels are selectable separately, but we want
   the user only to be able to select the axis as a whole, so we tie the selected states of the tick labels
   and the axis base line together. However, the axis label shall be selectable individually.

   The selection state of the left and right axes shall be synchronized as well as the state of the
   bottom and top axes.

   Further, we want to synchronize the selection of the graphs with the selection state of the respective
   legend item belonging to that graph. So the user can select a graph by either clicking on the graph itself
   or on its legend item.
  */

  // make top and bottom axes be selected synchronously, and handle axis and tick labels as one selectable object:
  if (customPlot->xAxis->selected().testFlag(QCPAxis::spAxis) || customPlot->xAxis->selected().testFlag(QCPAxis::spTickLabels) ||
      customPlot->xAxis2->selected().testFlag(QCPAxis::spAxis) || customPlot->xAxis2->selected().testFlag(QCPAxis::spTickLabels))
  {
    customPlot->xAxis2->setSelected(QCPAxis::spAxis|QCPAxis::spTickLabels);
    customPlot->xAxis->setSelected(QCPAxis::spAxis|QCPAxis::spTickLabels);
  }
  // make left and right axes be selected synchronously, and handle axis and tick labels as one selectable object:
  if (customPlot->yAxis->selected().testFlag(QCPAxis::spAxis) || customPlot->yAxis->selected().testFlag(QCPAxis::spTickLabels) ||
      customPlot->yAxis2->selected().testFlag(QCPAxis::spAxis) || customPlot->yAxis2->selected().testFlag(QCPAxis::spTickLabels))
  {
    customPlot->yAxis2->setSelected(QCPAxis::spAxis|QCPAxis::spTickLabels);
    customPlot->yAxis->setSelected(QCPAxis::spAxis|QCPAxis::spTickLabels);
  }
}

void graphics::NewSelected(QCPAbstractPlottable *plottable)
{
    QCPGraph *selected = (QCPGraph*)plottable;
    if(plottable->name().size() >0)
    {
        fil = selected->row;
        col = selected->col;
        int type = selected->type;
        selected->addToLegend();
        QString title = xAxisLabel + " = " + plottable->name().split(',')[0] +
                        " " + measureX + ", " + yAxisLabel + " = " +
                        plottable->name().split(',')[1] + " " + measureY;
        customPlot->setTitle(title);
        emit(SelectionChanged(fil,col,type));
    }
}

void graphics::mousePress(QMouseEvent* event)
{
  // if an axis is selected, only allow the direction of that axis to be dragged
  // if no axis is selected, both directions may be dragged

  if (customPlot->xAxis->selected().testFlag(QCPAxis::spAxis))
    customPlot->setRangeDrag(customPlot->xAxis->orientation());
  else if (customPlot->yAxis->selected().testFlag(QCPAxis::spAxis))
    customPlot->setRangeDrag(customPlot->yAxis->orientation());
  else
    customPlot->setRangeDrag(Qt::Horizontal|Qt::Vertical);

  if(customPlot->plottableAt(event->posF(),true) != NULL)
      NewSelected(customPlot->plottableAt(event->posF(),true));//,event->posF()
  else
  {
      customPlot->setTitle("");
      emit(SelectionChanged(-1,-1,-1));
  }
}

void graphics::mouseWheel(QWheelEvent* event)
{
  // if an axis is selected, only allow the direction of that axis to be zoomed
  // if no axis is selected, both directions may be zoomed

  if (customPlot->xAxis->selected().testFlag(QCPAxis::spAxis))
    customPlot->setRangeZoom(customPlot->xAxis->orientation());
  else if (customPlot->yAxis->selected().testFlag(QCPAxis::spAxis))
    customPlot->setRangeZoom(customPlot->yAxis->orientation());
  else
    customPlot->setRangeZoom(Qt::Horizontal|Qt::Vertical);
}

void graphics::graphClicked(QCPAbstractPlottable *item)
{
    NewSelected(item);
}

void graphics::graphMouseMove(QMouseEvent *event)
{
    if(customPlot->plottableAt(event->posF(),true) != NULL /*&& ui->actionFollow_motion->isChecked()*/ )
        NewSelected(customPlot->plottableAt(event->posF(),true));
    else
    {
        customPlot->setTitle("");
    }
}

void graphics::InitCurvesData(FMatrix X, FMatrix Y, QColor *color, int curveClasification)
{
    //if(ui->actionShow_Curves->isChecked())
        ShowCurvesData(X,Y,color,curveClasification);
    //else ShowGraphs(X,Y,color,curveClasification);
}

void graphics::ShowCurvesData(FMatrix X, FMatrix Y, QColor* color, int curveClasification)
{
    QPen graphPen;
    graphPen.setWidthF(2);
    QList<QColor> colors;
    if(color == NULL)
        color = new QColor(25,25,112);
    // colors = GetColors(Y.rows);
    for (int i=0; i<Y.cols; i++)
    {
        QCPGraphCurve *newCurve = new QCPGraphCurve(customPlot->xAxis,customPlot->yAxis);
        customPlot->addPlottable(newCurve);

        for(int j = 0 ; j < Y.rows;j++)
        {
            double x = X.GetElem(j,i);
            double y = Y.GetElem(j,i);
            EigenData eigen(x, y, j, i, curveClasification,
                            QString::number(x)+", " + QString::number(y));
            newCurve->addData(eigen);
            if(colors.size() != 0)
                graphPen.setColor(colors[i]);
            else graphPen.setColor(*color);
            newCurve->setPen(graphPen);
            newCurve->setSelectable(false);
        }
    }
    customPlot->replot();
}

void graphics::ShowGraphs(FMatrix X, FMatrix Y, QColor *color, int curveClasification)
{
    double x, y;
    QPen graphPen;
    graphPen.setWidthF(2);
    QList<QColor> colors;
    if(color == NULL)
        color = new QColor(25,25,112);
     colors = GetColors(Y.cols);
    for (int i=0; i<Y.cols; i++)
    {
        int pepe = 0;
        for(int j = 0 ; j < Y.rows;j++)
        {
            x = X.GetElem(j,i);
            y = Y.GetElem(j,i);
            QString name = QString::number(x)+", " + QString::number(y);
            customPlot->addGraph();
            customPlot->graph()->setName(name);
            customPlot->graph()->setScatterStyle(QCP::ssDot);
            customPlot->graph()->setScatterSize(4);
            customPlot->graph()->addData(x, y);
            if(colors.size() != 0)
                graphPen.setColor(colors[i]);
            else graphPen.setColor(*color);
//            if(illRows.contains(j))
//                graphPen.setColor(*red);
            customPlot->graph()->setPen(graphPen);
            customPlot->graph()->removeFromLegend();
            customPlot->graph()->row = j;
            customPlot->graph()->col = i;
            customPlot->graph()->type = curveClasification;
        }
    }
    customPlot->replot();
}

void graphics::ShowCurves(bool showC){
    customPlot->clearPlottables();
    customPlot->replot();
    if(showC && X.cols > 0)
        ShowCurvesData(X,Y,NULL,-1);
    else if(!showC && X.cols > 0)
        ShowGraphs(X,Y,NULL,-1);
    customPlot->replot();
}

void graphics::ChangeMeasures(QString newMeasure, double inc, QString axis)
{
    customPlot->clearPlottables();
    if(axis == "X")
    {
        measureX = newMeasure;
        autoScaleIncX *= inc;
    }
    else
    {
        measureY = newMeasure;
        autoScaleIncY *= inc;
    }

    if(X.cols > 0 && Y.cols > 0)
    {
        if(axis == "X")
        {
            X = X.multi_escalar(inc);
            InitCurvesData(X, Y, NULL, -1);
        }
        else
        {
            Y = Y.multi_escalar(inc);
            InitCurvesData(X, Y, NULL, -1);
        }
    }
    if(axis == "X")
          customPlot->xAxis->setRange(0, customPlot->xAxis->range().upper*inc);
    else
        customPlot->yAxis->setRange(0, customPlot->yAxis->range().upper*inc);

    customPlot->xAxis->setLabel(xAxisLabel + " " + measureX);
    customPlot->yAxis->setLabel(yAxisLabel + " " + measureY);
    customPlot->replot();
}

QList<QColor> graphics::GetColors(int elemsCount)
{
    QList<QColor> result;
    for(int i = 0 ; i < elemsCount; i ++)
    {
        int r = rand()%245+10;
        int g = rand()%245+10;
        int b = rand()%245+10;
        result.push_back(QColor(r,g,b));
    }
    return result;
}
