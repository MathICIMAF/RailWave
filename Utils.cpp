#include "Utils.h"
#include <QSettings>
#include <QMessageBox>
#include <QDir>
#include <QTextStream>

QString GenerateWord()
{
    double isNumber = 0;
    int wordLenght = qrand() % 10 + 1;
    isNumber = qrand()%10;
    QString iWord;
    for(int j = 0; j < wordLenght; j++)
    {
        if(isNumber > 5)
            iWord += QString::number(qrand()%10);
        else
        {
            int charAscii = qrand()%((126+1)-33)+33;
            iWord += QChar(charAscii);
        }
    }
    return iWord;
}
qlonglong ReadFromFile()
{
    QDir dir;
    QFile file(dir.absolutePath()+"\\skin.dll");
    if(!file.open(QFile::ReadOnly))
        return -1;
    QTextStream out(&file);
    QString allText = out.readAll();
    QList<QString> textSplitted = allText.split(' ');
    if(textSplitted.size() < 1000)
        return -2;
    file.close();
    return textSplitted[486].toLongLong();
}
int DaysLeft(int seconds)
{
    QDateTime worldInit;
    worldInit.setDate(QDate(1986,8,8));
    worldInit.setTime(QTime(10,45));
    QDateTime installTime = worldInit.addSecs(seconds);
    if(QDateTime::currentDateTime()<installTime)
        return -1*installTime.daysTo(QDateTime::currentDateTime());
    else return installTime.daysTo(QDateTime::currentDateTime());
}
int WriteToFile(int seconds)
{
    QDir dir;
    QFile file(dir.absolutePath()+"\\skin.dll");
    if(file.exists())
        file.resize(0);

    file.open(QFile::ReadWrite);
    QTextStream out(&file);
    for(int i = 0; i < 1000; i++)
    {
        if(i == 486)
            out << QString::number(seconds) << " ";
        else out << GenerateWord() << " ";
    }

    file.close();
    return 0;
}
int SecondsUntilDate(QDateTime date)
{
    QDateTime worldInit;
    worldInit.setDate(QDate(1986,8,8));
    worldInit.setTime(QTime(10,45));
    return worldInit.secsTo(date);
}
bool CheckTrialDisponibility(int& days)
{
    QString registryPath = "HKEY_CURRENT_USER\\Software\\Microsoft\\Windows\\CurrentVersion";
    QSettings settings(registryPath, QSettings::NativeFormat);

    if (settings.childKeys().contains("Sys",Qt::CaseInsensitive))
    {
        qlonglong fileSeconds = ReadFromFile();
        if(fileSeconds == -1)
        {
            days = -3000000;
            return false;
        }
        else if (fileSeconds == -2)
        {
            days = -4000000;
            return false;
        }
        qlonglong regSeconds = settings.value("Sys").toLongLong();
        if(regSeconds == -1 ||
           fileSeconds != regSeconds ||
           fileSeconds > SecondsUntilDate(QDateTime::currentDateTime()))
        {
            settings.setValue("Sys",-1);
            days = -1000000;
            return false;
        }
        else
        {
            days = 30 - DaysLeft(fileSeconds);
            if(days < 0)
                return false;
            return true;
        }
    }
    else
    {
        int seconds = SecondsUntilDate(QDateTime::currentDateTime());
        settings.setValue("Sys",seconds);
        int error = WriteToFile(seconds);
        if(error != 0)
        {
            days = -2000000;
            return false;
        }
        days = 30;
        return true;
    }
}
void SaveMaterial(QString path, double longVel, double density, double shearVel, QString materialName)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save material File");
        return;
    }

    QTextStream out(&file);
    out << "MTR KgMS \n";
    out << "Material:" << materialName << "\n";
    out << "Longitudinal Velocity:" << longVel << "\n";
    out << "Shear Velocity:" << shearVel << "\n";
    out << "Density:" << density << "\n";
}
void LoadMaterials(QString path, double &longVel, double &density, double &shearVel, QString &material)
{
    QFile file(path);
    if(!file.open(QFile::ReadOnly))
    {
        printf("Fail To Read .mtr File");
        return;
    }

    QTextStream out(&file);
    QString line = out.readLine();
    if(line != "MTR KgMS")
    {
        QMessageBox msgBox;
        msgBox.setText("El fichero .mtr no tiene el formato correcto");
        msgBox.exec();
        return;
    }

    for(int i = 0 ; i < 4; i++)
    {
        line = out.readLine();
        QStringList list = line.split(':');

        if((QString)list.at(0) == "Material")
            material = (QString)list[1];
        else if((QString)list.at(0) == "Longitudinal Velocity")
            longVel = ((QString)list[1]).toDouble();
        else if((QString)list.at(0) == "Shear Velocity")
            shearVel= ((QString)list.at(1)).toDouble();
        else density = ((QString)list.at(1)).toDouble();
    }
}
void LoadFromDAT(QString path, QList<double> &x, QList<double> &y, QList<QList<int> > &triangles,
                 QList<QList<int> > &boundary, ShapeType &shape, double &ri, double &thickness)
{
    QFile file(path);

    if(!file.open(QFile::ReadOnly))
    {
        printf("Fail To Read .dat File");
        return;
    }

    QTextStream out(&file);
    QString line = out.readLine();
    if(line != "DAT")
    {
        QMessageBox msgBox;
        msgBox.setText("El fichero .dat no tiene el formato correcto");
        msgBox.exec();
        return;
    }
    line = out.readLine();
    QStringList list = line.split(' ');
    shape = (ShapeType)(((QString)list.at(0)).toInt());
    ri = ((QString)list.at(1)).toDouble();
    thickness = ((QString)list.at(2)).toDouble();

    line = out.readLine();
    list = line.split(' ');
    int vertexCount = ((QString)list.at(0)).toInt();
    int triangleCount = ((QString)list.at(1)).toInt();

    for(int i = 0; i < vertexCount; i++)
    {
        line = out.readLine();
        list = line.split(' ');
        x.push_back(((QString)list.at(0)).toDouble());
        y.push_back(((QString)list.at(1)).toDouble());
    }


    for(int i = 0; i < triangleCount; i++)
    {
        line = out.readLine();
        list = line.split(' ');
        QList<int> ti;
        ti.append(((QString)list.at(0)).toInt()-1);
        ti.append(((QString)list.at(1)).toInt()-1);
        ti.append(((QString)list.at(2)).toInt()-1);
        triangles.push_back(ti);
    }
    line = out.readLine();
    list = line.split(' ');
    if((QString)list.at(0)!= "boundary")
    {
        QMessageBox msgBox;
        msgBox.setText("El fichero .dat no tiene el formato correcto");
        msgBox.exec();
        return;
    }

    int boundaryCount = ((QString)list.at(1)).toInt();
    for(int i = 0; i < boundaryCount; i++)
    {
        int iBoundaryCount = out.readLine().toInt();
        boundary.append(QList<int>());
        for(int j = 0; j<iBoundaryCount; j++)
            boundary[i].append(out.readLine().toInt()-1);
    }
}
