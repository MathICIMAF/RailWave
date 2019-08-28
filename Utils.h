#include <QDateTime>
#include "Constants.h"

QString GenerateWord();
qlonglong ReadFromFile();
int DaysLeft(int seconds);
int WriteToFile(int seconds);
int SecondsUntilDate(QDateTime date);
bool CheckTrialDisponibility(int &days);

void SaveMaterial(QString path, double longVel, double density, double shearVel, QString materialName);
void LoadMaterials(QString path, double &longVel, double &density, double &shearVel, QString &materialName);
void LoadFromDAT(QString path, QList<double> &x, QList<double> &y, QList<QList<int> > &triangles, QList<QList<int> > &boundary, ShapeType &shape, double &ri, double &thickness);
