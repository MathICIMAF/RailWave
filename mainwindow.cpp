#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include "Surfaces.h"
#include <QLabel>
#include <QMap>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    materialInfo.begin();
    railLoaded = false;
    graphloaded = false;
    animated = false;
    //ui->secInitial->pixmap()->scaledToWidth(ui->secInitial->width());

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionOpen_Rail_triggered()
{    
   if(railLoaded){
       QMessageBox msgBox;
       msgBox.setText("You must to close current rail!");
       msgBox.exec();
       return;
   }
   QString fileName = QFileDialog::getOpenFileName(this,tr("Open Rail"), ".",tr("Rail files (*.rl)"));
   QFile file(fileName);
   if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
       return;
   QMap<QString,double> dic;
   QTextStream text(&file);
   QList<double> values;   
   while(!text.atEnd()){
       QString line = text.readLine();

       QStringList list = line.split(':',QString::SkipEmptyParts);
       QString left = list[0];       
       left.replace(" ","");
       QString right = list[1];
       right.replace(" ","");
       double temp;
       if(left[0].isLower()){
           temp = toDoubleAngle(right);
           dic[left] = temp;
       }
       else
           temp = toDoubleLength(right);
           dic[left] = temp;
       //values.append(((QString)list[1]).toDouble());
   }
   pol = new RailPolygon(dic["HD"],dic["FD"],dic["BD"],dic["HW"],dic["W"],dic["BW"],dic["R1H"],dic["R2H"],dic["RHB"],dic["alpha"],dic["beta"],dic["theta"]);
   //pol = new RailPolygon(values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9],values[10],values[11]);
   pol->InitialPolygon();
   pol->SubdividePolygon(5);
   raild = new RailDesigned(pol,28,28);
   mesh = new Mesh(raild);
   raild->target_level = 2;
   mesh->adapt_refine(0);
   front = pol->Polygon();
   radius = pol->GetRadius();
   ShowRail(pol->Polygon(),mesh,pol->GetRadius());
   WriteOFFMesh(mesh,"Rail");
   //ui->rail3D->ShowShape(mesh,ui->secNumbSp->value(),ui->spacingSp->value());
   /*Dictionary d;
   d.Fill(text.readAll());
   QList<PointS> initialData = CreatePoints(d);
   if(initialData.count() == 0)
   {
       QMessageBox msgBox;
       msgBox.setText("Selected file has not the correct format.");
       msgBox.exec();
   }
   else{
       QList<PointS> curva1 = DelimitCurves(1, initialData);
       QList<PointS> curva2 = DelimitCurves(2, initialData);
       QList<PointS> curva3 = DelimitCurves(3, initialData);
       QList<PointS> curva4 = DelimitCurves(4, initialData);
       QList<PointS> curva5 = DelimitCurves(5, initialData);
       QList<float> ws = WReader();


       result = AproximateSpline(curva1, curva2, curva3, curva4, curva5, ws);
       ui->railSec->drawControl(initialData);
       this->rail = new Rail(result,28,28);
       this->mesh = new Mesh(rail);
       //rail->target_level = 1;
       //mesh->adapt_refine(0);
       DrawRail(result);// ui->railSec->drawBorder(result);
       ui->railSec->drawMesh(mesh);
       ui->rail3D->ShowShape(mesh,ui->secNumbSp->value(),ui->spacingSp->value());      
       ui->animatePB->setEnabled(false);      
       ui->railInitial->hide();
       ui->curvesInit->hide();
       ui->secInitial->hide();
   }*/
}


void MainWindow::DrawRail(QList<PointS> res){
    for(int i = 0; i < res.size(); i++){
        QPointF p(res[i]._x,res[i]._y);
        front.append(p);
    }
    ui->railSec->drawBorder(front);
}

void MainWindow::on_actionClose_triggered()
{
    this->close();
}

void MainWindow::on_zoomIn_clicked()
{    
    if(railLoaded)
        ui->rail3D->ZoomIn();
}

void MainWindow::on_zoomOut_clicked()
{
    if(railLoaded)
        ui->rail3D->ZoomOut();
}

void MainWindow::on_actionClose_Rail_triggered()
{
    if(railLoaded && !animated){
        ui->railSec->CleanControl();
        ui->rail3D->CleanControl();
        railLoaded = false;
        graphloaded = false;
        ui->animatePB->setEnabled(false);
        ui->actionShow_Curves->setEnabled(false);        
        delete rail;
        ui->dispCurvs->clearPlottables();
        ui->dispCurvs->replot();        
    }
    else if(animated){
        QMessageBox msgBox;
        msgBox.setText("You must to stop animation!");
        msgBox.exec();
    }
}

void MainWindow::on_computeButt_clicked()
{
    if(railLoaded){
        if(animated){
            QMessageBox msgBox;
            msgBox.setText("You must to stop animation!");
            msgBox.exec();
            return;
        }
        ui->graphPropGB->setEnabled(false);
        ui->curvesGB->setEnabled(false);
        ui->matPropGB->setEnabled(false);
        ui->rail3dGB->setEnabled(false);
        ui->sectionGB->setEnabled(false);
        ui->meshRefGB->setEnabled(false);
        ui->matPropGB->setEnabled(false);
        ui->solTypeGB->setEnabled(false);
        ui->menuDisplay->setEnabled(false);
        ui->actionShow_Curves->setEnabled(false);
        ui->menuFile->setEnabled(false);
        ui->animatePB->setEnabled(false);
        maxChi = ui->maxWaveSp->value();

        //El paso se calcula para obtener 100 iteraciones
        //del calculo del problema generalizado de valores propios
        int maxNumeroIter = 100;
        step = maxChi/maxNumeroIter;//ui->stepSp->value();

        shearVel = ui->shearVelText->text().toInt();
        longVel = ui->longVelocText->text().toInt();
        density = ui->densityText->text().toInt();
        materialName = ui->matNameCBox->itemText(0);
        curvesNo = ui->curvNumbSp->value();
        refinementSteps = 0;
        refinementThresh = 0;
        structureType = PipeEnum;
        materialInfo.clear();
        materialInfo.append("Rail");
        materialInfo.append("Free Boundary");
        materialInfo.append(materialName);
        materialInfo.append("");
        materialInfo.append("Density = " + QString::number(density) + " kg/m3");
        materialInfo.append("Long.Vel = " + QString::number(longVel) + " m/s");
        materialInfo.append("Shear.Vel = " + QString::number(shearVel) + " m/s");
        if(ui->rbLinear->isChecked())
            probType = Linear;
        else
            probType = Quadratic;
        datawindow *dataWin = new datawindow(refinementSteps,refinementThresh,-1,-1,ui->rail3D->x,ui->rail3D->y,
                                             ui->rail3D->triangles,ui->rail3D->boundaryIndices,shearVel,longVel,
                                             density,maxChi,step,curvesNo,false,structureType,probType);
        connect(dataWin, SIGNAL(Cancelled(datawindow*)),this, SLOT(Cancelled(datawindow*)));
        connect(dataWin,SIGNAL(CurvesToShow(bool,datawindow*)), this, SLOT(ShowDisplayMenu(bool, datawindow*)));
        dataWin->show();        
    }
}

void MainWindow::on_secNumbSp_valueChanged(int value)
{    
    if(railLoaded)
        ui->rail3D->SetSections(value);
}

void MainWindow::on_spacingSp_valueChanged(double value)
{
    if(railLoaded)
        ui->rail3D->SetSpacing(value);
}

void MainWindow::Cancelled(datawindow *sender)
{
    if(railLoaded)
        delete sender;
    if(!graphloaded){
        ui->menuDisplay->setEnabled(false);
        ui->actionShow_Curves->setEnabled(false);
    }
    ui->graphPropGB->setEnabled(true);
    ui->curvesGB->setEnabled(true);
    ui->matPropGB->setEnabled(true);
    ui->rail3dGB->setEnabled(true);
    ui->sectionGB->setEnabled(true);
    ui->meshRefGB->setEnabled(true);
    ui->matPropGB->setEnabled(true);
    ui->solTypeGB->setEnabled(true);
    ui->menuFile->setEnabled(true);
}

void MainWindow::ShowDisplayMenu(bool val, datawindow *sender){
    ui->graphPropGB->setEnabled(true);
    ui->curvesGB->setEnabled(true);
    ui->matPropGB->setEnabled(true);
    ui->rail3dGB->setEnabled(true);
    ui->sectionGB->setEnabled(true);
    ui->meshRefGB->setEnabled(true);
    ui->matPropGB->setEnabled(true);
    ui->solTypeGB->setEnabled(true);
    ui->menuDisplay->setEnabled(true);
    ui->actionShow_Curves->setEnabled(true);
    ui->menuFile->setEnabled(true);
    if(val)
    {        
        double PI =  3.14159265359;
        f = sender->f;
        lf = sender->lF;
        ff = sender->fF;
        tf = sender->tF;
        vF = sender->vf;
        lvF = sender->lVf;
        tvF = sender->tVf;
        fvF = sender->fVf;
        vG = sender->vg;
        lvG = sender->lVg;
        fvG = sender->fVg;
        tvG = sender->tVg;
        chi = sender->NXIT;
        illRows = sender->illRows;
        avalsVects = sender->eigenValVects;
        lvalsVects = sender->longitudinalEigens;
        fvalsVects = sender->flexuralEigens;
        tvalsVects = sender->torsionalEigens;
        maxX = maxY = 1000;
        sender->close();
        secTime = sender->secTime;
        ui->dispCurvs->clearPlottables();
        gr = new graphics(materialInfo,f.multi_escalar(shearVel/(2*PI)),vF.multi_escalar(shearVel),true,illRows,ui->dispCurvs);
        graphloaded = true;
        connect(gr,SIGNAL(SelectionChanged(int, int, int)),ui->rail3D,SLOT(UpdateRowCol(int,int,int)));
        gr->ShowCurves(false);
        ui->animatePB->setEnabled(true);
        ui->actionShow_Curves->setChecked(false);
    }
}


void MainWindow::on_actionShow_Curves_triggered()
{
    if(ui->actionShow_Curves->isChecked()){
        ui->animatePB->setEnabled(false);
        gr->ShowCurves(true);
        if(animated){//Detener la animacion
            ui->rail3D->Animate(avalsVects,lvalsVects,fvalsVects,tvalsVects);
            ui->animatePB->setText("Animate");
        }
    }
    else{
        gr->ShowCurves(false);
        ui->animatePB->setEnabled(true);
    }
}

void MainWindow::on_animatePB_clicked()
{
    ui->rail3D->Animate(avalsVects,lvalsVects,fvalsVects,tvalsVects);
    if(ui->animatePB->text() == "Animate"){
        ui->animatePB->setText("Stop Animation");
        animated = true;
    }
    else{
        ui->animatePB->setText("Animate");
        animated = false;
    }
    ui->meshRefGB->setEnabled(false);
}



void MainWindow::on_actionAbout_triggered()
{
    Dialog *d = new Dialog(this);
    d->show();
}

void MainWindow::on_actionCreate_Rail_triggered()
{
    DesignRail * d = new DesignRail(this);
    d->show();
    connect(d,SIGNAL(AddThisRail(QList<QPointF>,Mesh*,double,RailDesigned*,bool)),this,SLOT(ShowRail(QList<QPointF>,Mesh*,double,RailDesigned*,bool)));
}

void MainWindow::ShowRail(QList<QPointF> front, Mesh *mesh,double radius,RailDesigned* rail,bool refine){

    ui->spacingSp->setValue(radius/10.5);
    ui->secNumbSp->setValue(15);
    ui->railSec->drawBorder(front);    
    ui->railSec->drawMesh(mesh);        
    railLoaded = true;
    ui->animatePB->setEnabled(false);
    this->radius = radius;
    if(this->mesh.isNull()){
        this->front = front;
        this->mesh = mesh;
        raild = rail;

    }
    if(!refine)
        ui->rail3D->ShowShape(mesh,ui->secNumbSp->value(),ui->spacingSp->value());    
}

void MainWindow::on_pushButton_clicked()
{
    if(railLoaded && !animated){
        raild->target_level++;
        mesh->adapt_refine(0);
        ShowRail(front,mesh,radius,0,true);
    }
}
