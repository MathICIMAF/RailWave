#include "designrail.h"
#include "ui_designrail.h"
#include <QFileDialog>

DesignRail::DesignRail(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::DesignRail)
{
    ui->setupUi(this);
}

DesignRail::~DesignRail()
{
    delete ui;
}

void DesignRail::on_HDlineEdit_textChanged(QString text)
{

}

void DesignRail::on_drawPB_clicked()
{
    double b[12];
    b[0] = toDoubleLength(ui->HDlineEdit->text());
    //b[0] = ui->HDlineEdit->text().toDouble();            //double hd =
    b[1] = toDoubleLength(ui->FDlineEdit->text());           //double fd =
    b[2] = toDoubleLength(ui->BDlineEdit->text());            //double bd =
    b[3] = toDoubleLength(ui->HWlineEdit->text());            //double hw =
    b[4] = toDoubleLength(ui->WlineEdit->text());              //double w = u
    b[5] = toDoubleLength(ui->BWlineEdit->text());            //double bw =
    b[6] = toDoubleLength(ui->R1HlineEdit->text());          //double r1h =
    b[7] = toDoubleLength(ui->R2HlineEdit->text());          //double r2h =
    b[8] = toDoubleLength(ui->RHBlineEdit->text());          //double rhb =
    b[9] = toDoubleAngle(ui->alphalineEdit->text());      //double alpha
    b[10] = toDoubleAngle(ui->betalineEdit->text());        //double beta
    b[11] = toDoubleAngle(ui->thetalineEdit->text());         //double theta

    for(int i = 0; i < 12; i++){
        if(b[i] == 0){
            QMessageBox m;
            m.setText("There are fields equals 0");
            m.exec();
            return;
        }
    }
    pol = new RailPolygon(b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10],b[11],this);
    pol->InitialPolygon();
    pol->SubdividePolygon(5);
    ui->widget->drawpolygon(pol);    
    this->front = pol->Polygon();
    rail = new RailDesigned(pol,28,28);
    mesh = new Mesh(rail);
    //rail->target_level = 1;
}

void DesignRail::on_drawmesh_clicked()
{
    ui->widget->drawmesh(mesh,rail);
    mesh->adapt_refine(0);
}

void DesignRail::on_refineB_clicked()
{
    ui->widget->refine();
}

void DesignRail::on_actionSave_Rail_triggered()
{
    QString file = QFileDialog::getSaveFileName(0,"Save file",QDir::currentPath(),
                                                "Rail files (*.rl)",new QString("Rail files (*.rl)"));
    QFile f(file);
    f.open(QIODevice::WriteOnly);
    QTextStream out(&f);
    out << "HD: " << ui->HDlineEdit->text() << endl;  
    out << "FD: " << ui->FDlineEdit->text() << endl;
    out << "BD: " << ui->BDlineEdit->text() << endl;
    out << "HW: " << ui->HWlineEdit->text() << endl;
    out << "W: " << ui->WlineEdit->text() << endl;
    out << "BW: " << ui->BWlineEdit->text() << endl;
    out << "R1H: " << ui->R1HlineEdit->text() << endl;
    out << "R2H: " << ui->R2HlineEdit->text() << endl;
    out << "RHB: " << ui->RHBlineEdit->text() << endl;
    out << "alpha: " << ui->alphalineEdit->text() << endl;
    out << "beta: " << ui->betalineEdit->text() << endl;
    out << "theta: " << ui->thetalineEdit->text();
   f.close();
}


void DesignRail::on_actionClose_triggered()
{
    this->close();
}

void DesignRail::on_pushButton_clicked()
{
    emit(AddThisRail(front,mesh,pol->GetRadius(),rail,false));
    this->close();
}
