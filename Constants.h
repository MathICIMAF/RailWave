#ifndef CONSTANTS_H
#define CONSTANTS_H

//Geometry defaults for article
const double defaultRi = 50;//0.001;//0.25;
const double defaultThickness = 50;//0.001;// 0.25;
const int defaultNumberSections = 12;
const int defaultDeltaZ = 0.1;
const double defaultHeight = 1;
const double defaultWidth = 1;

//Discretization defaults for article
const double defaultMaxChi =2;//2;//20;//2000;//
const double defaultStep = 0.05;//0.05;//0.1;//10;//
const int defaultCurves = 6;//6;//20;//6;//

//Material Defaults for article
//char* defaultMtrName = "Steel";
const double defaultDensity = 7.85e3;
const double defaultShearVelocity = 3.26e3;
const double defaultLongitudinalVelocity = 5.996e3;

enum ProblemType{Linear, Quadratic};
enum ShapeType {PipeEnum, BarEnum, RailEnum, PlateEnum};

#endif // CONSTANTS_H
