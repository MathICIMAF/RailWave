#include "datawindow.h"
#include "ui_datawindow.h"

datawindow::datawindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::datawindow)
{
    ui->setupUi(this);
}
datawindow::datawindow(int prefinemetSteps, double prefinementThresh, double ri, double pthickness,
                       QList<double> px, QList<double> py, QList<QList<int> > ptriangles,
                       QList<QList<int> > pBoundary, double pshearVel, double plongVel, double pdensity,
                           double pmaxWaveNumber, double pstep, int pcurvesNumber, bool psymmetry,
                       ShapeType pstructure, ProblemType problemType)
{
    ui = new Ui::datawindow();
    ui->setupUi(this);

    structure = pstructure;
    symmetries = psymmetry;
    x = px;
    y = py;
    triangles = ptriangles;
    boundaryIndex = pBoundary;
    refinementSteps = prefinemetSteps;
    refinementThresh = prefinementThresh;

    connect(&worker,SIGNAL(posChanged(int)),this,SLOT(setProgressBarPos(int)));
    connect(&worker,SIGNAL(finished()),this,SLOT(finish()));
    connect(&worker,SIGNAL(maxIteration(int)),ui->progressBar,SLOT(setMaximum(int)));
    worker.problemType = problemType;
    worker.maxWavesNumber = pmaxWaveNumber;
    worker.step = pstep;
    worker.autovalNumber = pcurvesNumber;
    worker.longVel = plongVel;//young
    worker.shearVel = pshearVel;//poisson
    worker.density = pdensity;

    ir = ri;
    thickness = pthickness;

    /*if(symmetries)
    {
        if(structure == PipeEnum || structure == BarEnum)
        {
            QList<Layer*> layrs;
            Material mtr("",pdensity,pshearVel,plongVel);
            if(structure == PipeEnum)
            {
                PipeLayer* pipe = new PipeLayer(ir,ir+thickness,mtr,true);
                layrs.push_back(pipe);
            }
            else
            {
                BarLayer* bar = new BarLayer(thickness,mtr,true);
                layrs.push_back(bar);
            }
            Geometry *window = new Geometry(layrs);
            window->GetXYTrianglesBoundaryRelation(symmetricalX, symmetricalY, symmetricalTriangles,
                                                   xZeroBoundaryIndices, yZeroBoundaryIndices);
            ir = window->InteriorRadius();
            thickness = window->Thickness();

            delete(window);
        }
        else
        {
            //else es rail. Cargar la malla de 1/2 del rail
        }
    }
    */
    Compute();
}

datawindow::~datawindow()
{
    if(worker.isRunning())
        worker.terminate();
    delete ui;
}

void  datawindow::Compute()
{
    emit(CurvesToShow(false, this));
    if(symmetries)
    {
        worker.symmetry = true;
        worker.symmetricalX = symmetricalX;
        worker.symmetricalY = symmetricalY;
        worker.symmetricalTriangles = symmetricalTriangles;
        worker.xZeroBoundaryIndices = xZeroBoundaryIndices;
        worker.yZeroBoundaryIndices = yZeroBoundaryIndices;
    }
    else
    {
        worker.symmetry = false;
        worker.x = x;
        worker.y = y;
        worker.triangles = triangles;
    }

    worker.start();
}
QList<QList<GVector2> > datawindow::GetBoundary()
{
    QList<QList<GVector2> > result;
    for(int i = 0; i<boundaryIndex.count(); i++)
    {
        result.append(QList<GVector2>());
        for(int j = 0; j < boundaryIndex[i].count(); j++)
            result[i].push_back(GVector2(x.at(boundaryIndex[i][j]), y.at(boundaryIndex[i][j])));
    }
    return result;
}
QList<QList<QList<GVector3> > >datawindow::GetDisplacements(int boundaries, int count,int sections)
{
    QList<QList<QList<GVector3> > >result;
    srand ( time(NULL) );

    for(int k = 0; k < boundaries; k++)
    {
        result.append(QList<QList<GVector3> >());
        for(int i = 0; i < sections; i++)
        {
            QList<GVector3> iesimo;
            for(int j=0;j<count;j++)
                iesimo.append(GVector3(0,0,0));
            result[k].append(iesimo);
        }
    }
    return result;
}

void datawindow::finish()
{
  NXIT = worker.NXIh;
  illRows = worker.illConditionatedRowsChis;
  boundary.clear();
  for(int i = 0; i < boundaryIndex.count(); i++)
  {
      boundary.append(QList<GVector2>());
      for(int j = 0; j < boundaryIndex[i].count(); j++)
          boundary[i].append(GVector2(x.at(boundaryIndex[i][j]),y.at(boundaryIndex[i][j])));
  }
  if(symmetries)
  {
      lF = worker.lF;
      lVf = worker.lVf;
      lVg = worker.lVg;
      tF = worker.tF;
      tVf = worker.tVf;
      tVg = worker.tVg;
      fF = worker.fF;
      fVf = worker.fVf;
      fVg = worker.fVg;
      longitudinalEigens = worker.longitudinalEigens;
      torsionalEigens = worker.torsionalEigens;
      flexuralEigens = worker.flexionalEigens;}
  else
  {
      f = worker.fM;
      vf = worker.vf;
      vg = worker.vg;
      eigenValVects = worker.eigenValsVects;
  }
  secTime = worker.secTime;
  emit(CurvesToShow(true, this));
}
void datawindow::on_cancelBtn_clicked()
{
    this->close();
}
void datawindow::closeEvent(QCloseEvent *event)
{
    if(worker.isRunning())
    {
        if(QMessageBox::question(this, "ATENTION!!!!",
                              "If you quit now, you will cancel the computation proccess. Are you sure you want to quit?",
                              QMessageBox::Yes, QMessageBox::No)
           != QMessageBox::No)
            emit(Cancelled(this));
        else event->ignore();
    }
}
void datawindow::setProgressBarPos(int pos)
{
    ui->progressBar->setValue(pos);
    if(this->isHidden())
        emit(SetProgressBar(pos));
}
