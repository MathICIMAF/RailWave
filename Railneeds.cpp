#include "RailNeeds.h"

Dictionary::Dictionary(const Dictionary &dic, QObject *parent):QObject(parent){
    values.clear();
    keys.clear();
    for(int i = 0; i < dic.values.count(); i++)
        values.append(dic.values[i]);
    for(int i = 0; i < dic.keys.count(); i++)
        keys.append(dic.keys[i]);
}

void Dictionary::Fill(QString text)
{
    QChar read;
    QString key;
    int length = text.length();

    int index = 0;
    while(index < length)
    {
        try
        {
            read = text[index++];
            if(read == ' ') continue;
            else if(read == ':')
            {
                QString value;
                while(index < length)
                {
                    read = text[index++];
                    if(read == ' ') continue;
                    else if(read != '\n' && read != '\r')
                        value += read;
                    else
                    {
                        Add(key, value.toFloat());
                        key = "";
                        break;
                    }
                    if(index == length-1)
                    {
                        Add(key, value.toFloat());
                        key = "";
                        break;
                    }
                }
            }
            else
            {
                key.append(read);
            }
        }
        catch(std::bad_alloc &ba)
        {

        }
    }
}
float Dictionary::GetValue(QString key)
{
    int pos = -1;
    for(int index = 0; index<keys.count();index++)
    {
        QString aux = keys.at(index);
        if(Equals(aux, key))
        {
            pos = index;
            break;
        }
    }
    return (pos == -1)?NULL:values.at(pos);
}
void Dictionary::Add(QString key, float value)
{
    keys.append(key);
    values.append(value);
}
bool Dictionary::Equals(QString key1, QString key2)
{
    if(key1.length()== key2.length())
    {
        for(int i=0;i<key1.length();i++)
        {
            if(key1.at(i) != key2.at(i))
                return false;
        }
    }
    else return false;
    return true;
}

PointS::PointS(QObject* parent):QObject(parent)
{
    _x = 0;
    _y = 0;
    _tanAsigned = false;
    _tanInf = false;
}

PointS::PointS(const PointS &pt, QObject *parent):QObject(parent){
    _x = pt._x;
    _y = pt._y;
    _tanAsigned = pt._tanAsigned;
    _tan = pt._tan;
    _tanInf = pt._tanInf;
}

PointS& PointS::operator =(const PointS& pt){
    _x = pt._x;
    _y = pt._y;
    _tanAsigned = pt._tanAsigned;
    _tan = pt._tan;
    _tanInf = pt._tanInf;
    return *this;
}

PointS::PointS(double x, double y,QObject* parent):QObject(parent)
{
    _x = x;
    _y = y;
    _tan = NULL;
    _tanAsigned = false;
    _tanInf = false;
}
PointS::PointS(double x, double y, double tan, QObject *parent):QObject(parent)
{
    _x = x;
    _y = y;
    _tan = tan;
    _tanAsigned = true;
    _tanInf = false;
}
void PointS::SetTangent(double value)
{
    _tan = value;
    _tanAsigned = true;
}

Poligon::Poligon(const Poligon &pol,QObject* parent):QObject(parent){
    _excepcion = pol._excepcion;
    _isSimple = pol._isSimple;
    _epsilon = pol._epsilon;
    _vertexCount = pol._vertexCount;
    _endVertex = pol._endVertex;
    _initVertex = pol._initVertex;
    _xCoord = pol._xCoord;
    _yCoord = pol._yCoord;
}

Poligon& Poligon::operator =(const Poligon& pol){
    _excepcion = pol._excepcion;
    _isSimple = pol._isSimple;
    _epsilon = pol._epsilon;
    _vertexCount = pol._vertexCount;
    _endVertex = pol._endVertex;
    _initVertex = pol._initVertex;
    _xCoord = pol._xCoord;
    _yCoord = pol._yCoord;
    return *this;
}

Poligon::Poligon(QList<double> xCoord, QList<double> yCoord,QObject*parent):QObject(parent)
{
    if(xCoord.count() != yCoord.count())
    {
        _excepcion = 1;
    }
    else
    {
        _vertexCount = xCoord.count();
        _xCoord = xCoord;
        _yCoord = yCoord;

        if(xCoord.count() > 2) _isSimple = false;
        else _isSimple = true;

        _initVertex = GVector2(xCoord.at(0), yCoord.at(0));
        _endVertex = GVector2(xCoord.at(xCoord.count()-1), yCoord.at(yCoord.count()-1));

        _epsilon = -1;
        MaximumDistance(xCoord, yCoord);
    }
}
double Poligon::PointToBaseSegmentDistance(GVector2 point)
{
    //Sorting vertex
    GVector2 left;
    GVector2 right;
    if(_initVertex.ax<_endVertex.ax)
    {
        left = _initVertex;
        right = _endVertex;
    }
    else if(_endVertex.ax < _initVertex.ax)
    {
        left = _endVertex;
        right = _initVertex;
    }
    else
    {
        if(_initVertex.ay < _endVertex.ay)
        {
            left = _initVertex;
            right = _endVertex;
        }
        else if(_endVertex.ay < _initVertex.ay)
        {
            left = _endVertex;
            right = _initVertex;
        }
        //Los dos extremos del segmentos son el mismo punto
        else return point.EuclideanDistance(_initVertex);
    }

    GVector2 qMenosP = *right.subs(left, right);
    double normaQMenosP = qMenosP.norm();

    GVector2 p0MenosP = *point.subs(left, point);

    double prodEsc = p0MenosP.dot(p0MenosP, qMenosP);
    double t = prodEsc/pow(normaQMenosP,2);

    GVector2 r = *left.add(left, *(qMenosP.mult(qMenosP,t)));

    GVector2 rMenosQ = *r.subs(right, r);
    double alpha = rMenosQ.norm()/normaQMenosP;

    if(0<=alpha && alpha <= 1)
    {
        GVector2 aux = *point.subs(r, point);
        return aux.norm();
    }
    else
    {
        GVector2 p0MenosQ = *point.subs(right, point);
        return qMin(p0MenosP.norm(), p0MenosQ.norm());
    }
}
void Poligon::MaximumDistance(QList<double> xCoord, QList<double> yCoord)
{
    if(_isSimple)
        _epsilon = 0;
    else{
        for(int i=1;i< xCoord.count()-1;i++)
        {
            GVector2 point(xCoord.at(i), yCoord.at(i));
            double dist = PointToBaseSegmentDistance(point);
            if(dist > _epsilon)
                _epsilon = dist;
        }
    }
}


BinTreePoligon::BinTreePoligon(QObject* parent):QObject(parent)
{
    _left = NULL;
    _right = NULL;
}

BinTreePoligon::BinTreePoligon(const BinTreePoligon& pol,QObject*parent):QObject(parent){
    //if(_left != NULL) delete _left;
    //if(_right != NULL) delete _right;
    _priority = pol._priority;
    _nodePoligon = pol._nodePoligon;
    _left = pol._left;
    _right = pol._right;
    _closestNeighbor = pol._closestNeighbor;
}

BinTreePoligon& BinTreePoligon::operator =(const BinTreePoligon& pol){
    //if(_left != NULL) delete _left;
    //if(_right != NULL) delete _right;//Verificar cdo se elimina un objeto
    _priority = pol._priority;
    _nodePoligon = pol._nodePoligon;
    _left = pol._left;
    _right = pol._right;
    _closestNeighbor = pol._closestNeighbor;
    return *this;
}

BinTreePoligon::BinTreePoligon(Poligon node, QObject *parent):QObject(parent)
{
    _nodePoligon = node;
    if(!_nodePoligon._isSimple)
    {
        QList<double> xCoordLeft;
        QList<double> yCoordLeft;
        QList<double> xCoordRight;
        QList<double> yCoordRight;

        for(int i = 0; i <= _nodePoligon._vertexCount/2;i++)
        {
                xCoordLeft.append(_nodePoligon._xCoord.at(i));
                yCoordLeft.append(_nodePoligon._yCoord.at(i));
        }
        for(int i = _nodePoligon._vertexCount/2; i < _nodePoligon._vertexCount; i++)
        {
                xCoordRight.append(_nodePoligon._xCoord.at(i));
                yCoordRight.append(_nodePoligon._yCoord.at(i));
        }

        Poligon leftPoligon(xCoordLeft, yCoordLeft);
        _left = new BinTreePoligon(leftPoligon);
        Poligon rightPoligon(xCoordRight, yCoordRight);
        _right = new BinTreePoligon(rightPoligon);
    }
    else
    {
        _left = NULL;
        _right = NULL;
    }
}
double BinTreePoligon::ComputePriority(GVector2 point)
{
    _priority =  _nodePoligon.PointToBaseSegmentDistance(point) - _nodePoligon._epsilon;
    return _priority;
}
double BinTreePoligon::PointToPoligonDistance(GVector2 point)
{
        /*
        1: procedure DISTANCIAMÍNIMA(p)
        2: 	distancia  ? infinito
        3: 	nodos ? HEAP( )
        4: 	Añadir raiz a nodos con prioridad DISTANCIACAJA(p, raiz)
        5: 	while nodos no esté vacía do
        6: 		siguiente, cota ? POP(nodos)
        7: 		if cota < distancia then
        8: 			if es una hoja then
        9: 				distancia ? cota
        10: 			else
        11: 				izq, der ?  HIJOS(siguiente)
        12:	 			Añadir izq a nodos, prioridad DISTANCIACAJA(p, izq)
        13:				 Añadir der a nodos, prioridad DISTANCIACAJA(p, der)
                return distancia

        */
    double distance = -1;
    BinHeapPoligon heap(_nodePoligon._vertexCount -1);
    ComputePriority(point);
    heap.add(this);

    while(!heap.is_empty())
    {
        BinTreePoligon next = *heap.remove_min();
        double cote = next.GetPriority();
        if(cote < distance || distance == -1)
        {
            if(next._nodePoligon._isSimple)
            {
                distance = cote;
                _closestNeighbor = next._nodePoligon;
            }
            else
            {
                next.GetLeftChild()->ComputePriority(point);
                next.GetRightChild()->ComputePriority(point);
                heap.add(next.GetLeftChild());
                heap.add(next.GetRightChild());
            }
        }
    }
    return distance;
}


int BinHeapPoligon::add(BinTreePoligon* x)
{
        int i = -1;
        i = BinHeap<BinTreePoligon*>::add(x);
        return i;
}
BinTreePoligon* BinHeapPoligon::remove_min()
{
        BinTreePoligon* x = BinHeap<BinTreePoligon*>::remove_min();
        return x;
}
BinTreePoligon* BinHeapPoligon::remove(int i)
{
        BinTreePoligon *x = BinHeap<BinTreePoligon*>::remove(i);
        return x;
}
void BinHeapPoligon::update(int i)
{
        up_heap(i);
        down_heap(i);
}
void BinHeapPoligon::swap(int i, int j)
{
        /*set_heap(elems[i], j);
        set_heap(elems[j], i);*/
        BinHeap<BinTreePoligon*>::swap(i, j);
}
BinTreePoligon::~BinTreePoligon(){}


QList<float> WReader()
{
    QList<float> ws;
    ws.append(3);
    ws.append(3);
    ws.append(0);
    ws.append(1);
    ws.append(3);
    ws.append(0);
    ws.append(3);
    ws.append(1);
    ws.append(0);
    ws.append(3);
    ws.append(3);
    return ws;
}
QList<PointS> CreatePoints(Dictionary measures)
{
    if(measures.GetValue("B") == NULL || measures.GetValue("H") == NULL || measures.GetValue("S") == NULL ||
       measures.GetValue("D") == NULL || measures.GetValue("F") == NULL || measures.GetValue("K1") == NULL ||
       measures.GetValue("K2") == NULL|| measures.GetValue("C1") == NULL|| measures.GetValue("C2") == NULL)
    {
        return QList<PointS>();
    }

    QList<PointS>  points;
    PointS aux1(-measures.GetValue("B")/2, 0);
    aux1._tanAsigned = true;
    aux1._tanInf = true;
    points.append(aux1);

    PointS aux2(-(2 * measures.GetValue("B")) / 5, measures.GetValue("D"));
    aux2.SetTangent(0.2/1);
    points.append(aux2);

    PointS aux3(-(2 * measures.GetValue("B")) / 5, measures.GetValue("D"));
    aux3.SetTangent(0.2/1);
    points.append(aux3);

    PointS aux4(-measures.GetValue("S") / 2, measures.GetValue("F"));
    aux4._tanAsigned = true;
    aux4._tanInf = true;
    points.append(aux4);

    PointS aux5(-measures.GetValue("S") / 2, measures.GetValue("H") - measures.GetValue("K1"));
    aux5._tanAsigned = true;
    aux5._tanInf = true;
    points.append(aux5);

    PointS aux6(-measures.GetValue("C2") / 2, measures.GetValue("H") - measures.GetValue("K2"));
    aux6.SetTangent(0/-1);
    points.append(aux6);

    PointS aux7(-measures.GetValue("C2") / 2, measures.GetValue("H") - measures.GetValue("K2"));
    aux7._tanAsigned = true;
    aux7._tanInf = true;
    points.append(aux7);

    PointS aux8(-measures.GetValue("C1") / 3, measures.GetValue("H"));
    aux8.SetTangent(0/1);
    points.append(aux8);

    ////////////////////////////////

    PointS aux9(measures.GetValue("C1") / 3, measures.GetValue("H"));
    aux9.SetTangent(0/1);
    points.append(aux9);

    PointS aux10(measures.GetValue("C2") / 2, measures.GetValue("H") - measures.GetValue("K2"));
    aux10._tanAsigned = true;
    aux10._tanInf = true;
    points.append(aux10);

    PointS aux11(measures.GetValue("C2") / 2, measures.GetValue("H") - measures.GetValue("K2"));
    aux11.SetTangent(0/1);
    points.append(aux11);

    PointS aux12(measures.GetValue("S") / 2, measures.GetValue("H") - measures.GetValue("K1"));
    aux12._tanAsigned = true;
    aux12._tanInf = true;
    points.append(aux12);

    PointS aux13(measures.GetValue("S") / 2, measures.GetValue("F"));
    aux13._tanAsigned = true;
    aux13._tanInf = true;
    points.append(aux13);

    PointS aux14((2 * measures.GetValue("B")) / 5, measures.GetValue("D"));
    aux14.SetTangent(-0.2/1);
    points.append(aux14);

    PointS aux15((2 * measures.GetValue("B")) / 5, measures.GetValue("D"));
    aux15.SetTangent(-0.2/1);
    points.append(aux15);

    PointS aux16(measures.GetValue("B") / 2, 0);
    aux16._tanAsigned = true;
    aux16._tanInf = true;
    points.append(aux16);

    return points;
}
QList<PointS> DelimitCurves(int curveNumber, QList<PointS> equis)
{
    QList<PointS> result;
    switch(curveNumber)
    {
        case 1:
        {
            for(int i=0;i<=1;i++)
                result.append(equis.at(i));
            break;
        }
        case 2:
        {
            for(int i=2;i<=5;i++)
                result.append(equis.at(i));
            break;
        }
        case 3:
        {
            for(int i=6;i<=9;i++)
                result.append(equis.at(i));
            break;
        }
        case 4:
        {
            for(int i=10;i<=13;i++)
                result.append(equis.at(i));
            break;
        }
        case 5:
        {
            for(int i=14;i<=15;i++)
                result.append(equis.at(i));
            break;
        }
    }
    return result;
}
QList<PointS> FillSpace(int totalPoints, PointS initPoint, PointS endPoint)
{
    QList<PointS> result;
    double step = (initPoint._x - endPoint._x)/(totalPoints+1);
    for(int index = 0; index <= totalPoints; index++)
    {
        PointS aux(initPoint._x-(index*step),initPoint._y);
        result.append(aux);
    }
    return result;
}
void CreateTriangulation(QList<PointS> result, double *pts, double *pts2, int *fcs)
{
    pts[0] = result[5]._x;     pts[1] = result[5]._y;    pts[2] = 0;
    pts[3] = result[0]._x;     pts[4] = result[0]._y;    pts[5] = 0;
    pts[6] = result[12]._x;    pts[7] = result[12]._y;   pts[8] = 0;
    pts[9] = result[12]._x;    pts[10] = 0;              pts[11] = 0;
    pts[12] = result[32]._x;   pts[13] = result[32]._y;  pts[14] = 0;
    pts[15] = result[32]._x;   pts[16] = 0;              pts[17] = 0;
    pts[18] = result[272]._x;  pts[19] = result[272]._y; pts[20] = 0;
    pts[21] = result[272]._x;  pts[22] = 0;              pts[23] = 0;
    pts[24] = result[292]._x;  pts[25] = result[292]._y; pts[26] = 0;
    pts[27] = result[292]._x;  pts[28] = 0;              pts[29] = 0;
    pts[30] = result[299]._x;  pts[31] = result[299]._y; pts[32] = 0;
    pts[33] = result[304]._x;  pts[34] = result[304]._y; pts[35] = 0;
    pts[36] = result[52]._x;   pts[37] = result[52]._y;  pts[38] = 0;
    pts[39] = result[252]._x;  pts[40] = result[252]._y; pts[41] = 0;
    pts[42] = result[64]._x;   pts[43] = result[64]._y;  pts[44] = 0;
    pts[45] = result[240]._x;  pts[46] = result[240]._y; pts[47] = 0;
    pts[48] = result[80]._x;   pts[49] = result[80]._y;  pts[50] = 0;
    pts[51] = result[224]._x;  pts[52] = result[224]._y; pts[53] = 0;
    pts[54] = result[108]._x;  pts[55] = result[108]._y; pts[56] = 0;
    pts[57] = result[104]._x;  pts[58] = result[104]._y; pts[59] = 0;
    pts[60] = result[80]._x;   pts[61] = result[108]._y; pts[62] = 0;
    pts[63] = result[224]._x;  pts[64] = result[196]._y; pts[65] = 0;
    pts[66] = result[196]._x;  pts[67] = result[196]._y; pts[68] = 0;
    pts[69] = result[200]._x;  pts[70] = result[200]._y; pts[71] = 0;
    pts[72] = result[124]._x;  pts[73] = result[124]._y; pts[74] = 0;
    pts[75] = result[80]._x;   pts[76] = result[136]._y; pts[77] = 0;
    pts[78] = result[224]._x;  pts[79] = result[136]._y; pts[80] = 0;
    pts[81] = result[180]._x;  pts[82] = result[180]._y; pts[83] = 0;

    QFile file("pts.txt");
    file.open(QFile::WriteOnly);
    QTextStream out(&file);
    for(int i = 0; i < 84; i+=3){
        out << pts[i] << " " << pts[i+1] <<" " <<pts[i+2]<< "\n";
    }
    file.close();

    pts2[0] = result[5]._x;    pts2[1] = result[5]._y;
    pts2[2] = result[0]._x;    pts2[3] = result[0]._y;
    pts2[4] = result[12]._x;   pts2[5] = result[12]._y;
    pts2[6] = result[12]._x;   pts2[7] = 0;
    pts2[8] = result[32]._x;   pts2[9] = result[32]._y;
    pts2[10] = result[32]._x;  pts2[11] = 0;
    pts2[12] = result[272]._x; pts2[13] = result[272]._y;
    pts2[14] = result[272]._x; pts2[15] = 0;
    pts2[16] = result[292]._x; pts2[17] = result[292]._y;
    pts2[18] = result[292]._x; pts2[19] = 0;
    pts2[20] = result[299]._x; pts2[21] = result[299]._y;
    pts2[22] = result[304]._x; pts2[23] = result[304]._y;
    pts2[24] = result[52]._x;  pts2[25] = result[52]._y;
    pts2[26] = result[252]._x; pts2[27] = result[252]._y;
    pts2[28] = result[64]._x;  pts2[29] = result[64]._y;
    pts2[30] = result[240]._x; pts2[31] = result[240]._y;
    pts2[32] = result[80]._x;  pts2[33] = result[80]._y;
    pts2[34] = result[224]._x; pts2[35] = result[224]._y;
    pts2[36] = result[108]._x; pts2[37] = result[108]._y;
    pts2[38] = result[104]._x; pts2[39] = result[104]._y;
    pts2[40] = result[80]._x;  pts2[41] = result[108]._y;
    pts2[42] = result[224]._x; pts2[43] = result[196]._y;
    pts2[44] = result[196]._x; pts2[45] = result[196]._y;
    pts2[46] = result[200]._x; pts2[47] = result[200]._y;
    pts2[48] = result[124]._x; pts2[49] = result[124]._y;
    pts2[50] = result[80]._x;  pts2[51] = result[136]._y;
    pts2[52] = result[224]._x; pts2[53] = result[136]._y;
    pts2[54] = result[180]._x; pts2[55] = result[180]._y;



    fcs[0] = 3;   fcs[1] = 0;   fcs[2] = 2;
    fcs[3] = 0;   fcs[4] = 3;   fcs[5] = 1;
    fcs[6] = 3;   fcs[7] = 4;   fcs[8] = 5;
    fcs[9] = 4;   fcs[10] = 3;  fcs[11] = 2;
    fcs[12] = 5;  fcs[13] = 6;  fcs[14] = 7;
    fcs[15] = 6;  fcs[16] = 5;  fcs[17] = 4;
    fcs[18] = 9;  fcs[19] = 6;  fcs[20] = 8;
    fcs[21] = 6;  fcs[22] = 9;  fcs[23] = 7;
    fcs[24] = 9;  fcs[25] = 10; fcs[26] = 11;
    fcs[27] = 10; fcs[28] = 9;  fcs[29] = 8;
    fcs[30] = 4;  fcs[31] = 13; fcs[32] = 6;
    fcs[33] = 13; fcs[34] = 4;  fcs[35] = 12;
    fcs[36] = 12; fcs[37] = 15; fcs[38] = 13;
    fcs[39] = 15; fcs[40] = 12; fcs[41] = 14;
    fcs[42] = 14; fcs[43] = 17; fcs[44] = 15;
    fcs[45] = 17; fcs[46] = 14; fcs[47] = 16;
    fcs[48] = 19; fcs[49] = 20; fcs[50] = 16;
    fcs[51] = 20; fcs[52] = 19; fcs[53] = 18;
    fcs[54] = 16; fcs[55] = 21; fcs[56] = 17;
    fcs[57] = 21; fcs[58] = 16; fcs[59] = 20;
    fcs[60] = 23; fcs[61] = 21; fcs[62] = 22;
    fcs[63] = 21; fcs[64] = 23; fcs[65] = 17;
    fcs[66] = 18; fcs[67] = 25; fcs[68] = 20;
    fcs[69] = 25; fcs[70] = 18; fcs[71] = 24;
    fcs[72] = 20; fcs[73] = 26; fcs[74] = 21;
    fcs[75] = 26; fcs[76] = 20; fcs[77] = 25;
    fcs[78] = 22; fcs[79] = 26; fcs[80] = 27;
    fcs[81] = 26; fcs[82] = 22; fcs[83] = 21;

    /*QFile file2("faces.txt");
    file2.open(QFile::WriteOnly);
    QTextStream out2(&file2);
    for(int i = 0; i < 84; i+=3){
        out2 << fcs[i] << " " << fcs[i+1] <<" " <<fcs[i+2]<< "\n";
    }
    file2.close();

    QFile file3("frontera.txt");
    file3.open(QFile::WriteOnly);
    QTextStream out3(&file3);
    for(int i = 0; i < result.size(); i++){
        out3 << result[i]._x<<" " << result[i]._y << "\n";
    }
    file3.close();*/
}
QList<PointS> AproximateSpline(QList<PointS> &curva1, QList<PointS> &curva2, QList<PointS> &curva3, QList<PointS> &curva4, QList<PointS> &curva5, QList<float> ws)
{
    QList<PointS> result;// = new QList<PointS>();
    int totalPoints = 0;
    result.append(ConicSpline(3, curva1, ws, 0, totalPoints));
    result.append(ConicSpline(5, curva2, ws, 1, totalPoints));
    result.append(ConicSpline(5, curva3, ws, 4, totalPoints));
    result.append(ConicSpline(5, curva4, ws, 7, totalPoints));
    result.append(ConicSpline(3, curva5, ws, 10, totalPoints));

    QList<PointS> toAdd = FillSpace(5*totalPoints,curva5[curva5.count()-1],curva1[0]);
    result.append(toAdd);
    return result;
}


PointS GetQi(PointS p1, PointS p2)
{
    PointS result;
    //y = mx+n
    //Hallando la recta que pasa por los puntos p1 con tangente p1._tan y p2 con tangente p2._tan
    //Entonces es infinito
    if(p1._tanInf)
    {
        //Entonces son dos rectas paralelas y no tienen punto de interseccion
        if(p2._tanInf)
        {
            result._x = (p1._x+p2._x)/2;
            result._y = (p1._y+p2._y)/2;
            result._tanInf = true;
            result._tanAsigned = true;
        }
        else
        {
            double n2 = p2._y - p2._tan * p2._x;
            result._x =p1._x;
            result._y = p2._tan*p1._x+n2;
        }
    }
    //La tangente en p1 no es perpendicular
    else
    {
        //La tangente de p2 es infinito
        if(p2._tanInf)
        {
            double n1 = p1._y - p1._tan * p1._x;
            result._x = p2._x;
            result._y = p1._tan*p2._x+n1;
        }
        else
        {
            if(qAbs(p1._tan - p2._tan)< 0.000000001)
            {
                result._x = (p1._x+p2._x)/2;
                result._y = (p1._y+p2._y)/2;
                result._tanInf = true;
                result._tanAsigned = true;
            }
            else
            {
                double n1 = p1._y - p1._tan * p1._x;
                double n2 = p2._y - p2._tan * p2._x;
                //Hallando el unto donde se intersecan esas rectas, que sera Qi
                double x = (n2 - n1)/(p1._tan - p2._tan);
                double y = p1._tan * x + n1;

                result._x = x;
                result._y = y;
            }
        }
    }
    return result;
}
double GetAutomaticTangent(QList<PointS> coord, int pos)
{
    PointS p1;
    PointS p2;
    PointS p3;
    PointS p4;
    PointS p5;

    if(pos == 0)
    {
        p1 = coord.at(coord.count()-2);
        p2 = coord.at(coord.count()-1);
        p3 = coord.at(pos);
    }
    else if(pos == 1)
    {
        p1 = coord.at(coord.count()-1);
        p2 = coord.at(pos-1);
        p3 = coord.at(pos);
    }
    else
    {
        p1 = coord.at(pos-2);
        p2 = coord.at(pos-1);
        p3 = coord.at(pos);
    }

    if(pos == coord.count()-1)
    {
        p4 = coord.at(0);
        p5 = coord.at(1);
    }
    else if(pos == coord.count()-2)
    {
        p4 = coord.at(coord.count()-1);
        p5 = coord.at(0);
    }
    else
    {
        p4 = coord.at(pos+1);
        p5 = coord.at(pos+2);
    }

    //Si son paralelas p1p2 y p3p4
    if(((p2._y - p1._y)/(p2._x - p1._x) == (p4._y - p3._y)/(p4._x - p3._x)) ||
       (p2._x - p1._x == 0 && p4._x - p3._x == 0))
    {
        PointS b = IntersectionPoint(p2, p3, p4, p5);
        PointS c = IntersectionPoint(b, p2, p1, p5);

        if((c._x == -DBL_MAX && c._y == -DBL_MAX) || c._x - p3._x == 0)
        {
            PointS newP;
            newP._x = coord.at(pos)._x;
            newP._y = coord.at(pos)._y;
            newP._tan = coord.at(pos)._tan;
            newP._tanAsigned = coord.at(pos)._tanAsigned;
            newP._tanInf = true;
            coord.replace(pos, newP);
            return -DBL_MAX;
        }
        else
        {
            PointS newP;
            newP._x = coord.at(pos)._x;
            newP._y = coord.at(pos)._y;
            newP._tan = coord.at(pos)._tan;
            newP._tanAsigned = coord.at(pos)._tanAsigned;
            newP._tanInf = false;
            coord.replace(pos, newP);
            return (c._y - p3._y)/(c._x - p3._y);
        }
    }
    else
    {
        PointS a = IntersectionPoint(p1, p2, p3, p4);
        //son paraleas p2p3 y p4p5
        if(((p3._y - p2._y)/(p3._x - p2._x) == (p5._y - p4._y)/(p5._x - p4._x)) ||
           (p3._x - p2._x == 0 && p5._x - p4._x == 0))
        {
            PointS c = IntersectionPoint(a, p3, p1, p5);
            if((c._x == -DBL_MAX && c._y == -DBL_MAX)|| c._x - p3._x == 0)
            {
                PointS newP;
                newP._x = coord.at(pos)._x;
                newP._y = coord.at(pos)._y;
                newP._tan = coord.at(pos)._tan;
                newP._tanAsigned = coord.at(pos)._tanAsigned;
                newP._tanInf = true;
                coord.replace(pos, newP);
                return -DBL_MAX;
            }
            else
            {
                PointS newP;
                newP._x = coord.at(pos)._x;
                newP._y = coord.at(pos)._y;
                newP._tan = coord.at(pos)._tan;
                newP._tanAsigned = coord.at(pos)._tanAsigned;
                newP._tanInf = false;
                coord.replace(pos, newP);
                return (c._y - p3._y)/(c._x - p3._x);
            }
        }
        else
        {
            PointS b = IntersectionPoint(p2, p3, p4, p5);
            //si son paralelas ab y p1p5
            if(((b._y - a._y)/(b._x - a._x) == (p5._y - p1._y)/(p5._x - p1._x)) ||
               (b._x - a._x == 0 && p5._x - p1._x == 0))
            {
                if(p5._x - p1._x == 0)
                {
                    PointS newP;
                    newP._x = coord.at(pos)._x;
                    newP._y = coord.at(pos)._y;
                    newP._tan = coord.at(pos)._tan;
                    newP._tanAsigned = coord.at(pos)._tanAsigned;
                    newP._tanInf = true;
                    coord.replace(pos, newP);
                    return -DBL_MAX;
                }
                else
                {
                    PointS newP;
                    newP._x = coord.at(pos)._x;
                    newP._y = coord.at(pos)._y;
                    newP._tan = coord.at(pos)._tan;
                    newP._tanAsigned = coord.at(pos)._tanAsigned;
                    newP._tanInf = false;
                    coord.replace(pos, newP);
                    return (p5._y - p1._y)/(p5._x - p1._x);
                }
            }
            else
            {
                PointS c = IntersectionPoint(p1, p5, a, b);
                if((c._x == -DBL_MAX && c._y == -DBL_MAX) || c._x - p3._x == 0)
                {
                    PointS newP;
                    newP._x = coord.at(pos)._x;
                    newP._y = coord.at(pos)._y;
                    newP._tan = coord.at(pos)._tan;
                    newP._tanAsigned = coord.at(pos)._tanAsigned;
                    newP._tanInf = true;
                    coord.replace(pos, newP);
                    return -DBL_MAX;
                }
                else
                {
                    PointS newP;
                    newP._x = coord.at(pos)._x;
                    newP._y = coord.at(pos)._y;
                    newP._tan = coord.at(pos)._tan;
                    newP._tanAsigned = coord.at(pos)._tanAsigned;
                    newP._tanInf = false;
                    coord.replace(pos, newP);
                    return (c._y - p3._y)/(c._x - p3._x);
                }
            }
        }
    }
}
double GetShoulderTan(PointS initPoint, PointS endPoint)
{
    return (endPoint._y-initPoint._y) / (endPoint._x - initPoint._x);
}
int SubdivideEdges(QList<PointS> &result, int iterationNumber, PointS initPoint, PointS endPoint, double w)
{
    if(iterationNumber == 0) return 0;

    PointS qi = GetQi(initPoint, endPoint);
    double shoulderTan = GetShoulderTan(initPoint, endPoint);
    double shoulderX = initPoint._x/(2*(1+w)) + qi._x*w/(1+w) + endPoint._x/(2*(1+w));
    double shoulderY = initPoint._y/(2*(1+w)) + qi._y*w/(1+w) + endPoint._y/(2*(1+w));
    PointS shoulderPoint(shoulderX, shoulderY, shoulderTan);
    if(qAbs(initPoint._x - endPoint._x) < 0.000000001)
        shoulderPoint._tanInf = true;

    double newW = sqrt((1+w)/2);
    int left = 0;
    left = SubdivideEdges(result, iterationNumber-1, initPoint, shoulderPoint, newW);
    result.append(shoulderPoint);
    int right = 0;
    right = SubdivideEdges(result, iterationNumber-1, shoulderPoint, endPoint, newW);
    return left+right+1;
}
PointS IntersectionPoint(PointS p1, PointS p2, PointS p3, PointS p4)
{
    if(p2._x == p1._x)
    {
        //las dos rectas son de la forma x = cte y son paralelas.
        if(p4._x == p3._x)
        {
            return PointS(1/p4._x - p3._x, 1/p4._x - p3._x);
        }
        //p1p2 es de la forma x=cte pero p3p4 no.
        else
        {
            double m2 = (p4._y - p3._y)/(p4._x - p3._x);
            double n2 = p4._y - m2*p4._x;

            return PointS(p2._x, m2*p2._x + n2);
        }
    }
    //p1p2 no es de la forma x = cte.
    else
    {
        double m1 = (p2._y - p1._y)/(p2._x - p1._x);
        double n1 = p2._y - m1*p2._x;

        //p3p4 es de la forma x = cte.
        if(p4._x == p3._x)
        {
            return PointS(p4._x, m1*p4._x + n1);
        }
        //p1p2 no es de la forma x = cte.
        else
        {
            double m2 = (p4._y - p3._y)/(p4._x - p3._x);
            double n2 = p4._y - m2*p4._x;

            if(m2 == m1)
                return PointS();
            double x = (n1 - n2)/(m2 - m1);
            double y = m1*x + n1;

            return PointS(x,y);
        }
    }
}
double GetAutomaticParameter(PointS p1, PointS pq, PointS p2, PointS p)
{
    /* Dadas las coordnadas baricentricas(u, v, 1-u-v) de pos3 se puede calcular
    el w resolviendo la siguiente ecuacion: v^2 - 4w^2*u(1-u-v) = 0 */
    double detA = ((pq._y-p2._y) * (p1._x - p2._x) - (pq._x - p2._x) * (p1._y-p2._y));
    double v = (((p1._x- p2._x) * (p._y - p2._y)) - (p._x - p2._x) * (p1._y - p2._y)) / detA;
    double u = ((p._x - p2._x)*(pq._y - p2._y) - (p._y - p2._y)*(pq._x-p2._x))/detA;

    return qAbs(sqrt(pow(v,2)/(4*(u*(1-u-v)))));
}
//El vector debe tener al menos 3 puntos para llamar al metodo GetAutomaticParameter y 5 para llamar al GetAutomaticTangent
QList<PointS> ConicSpline(int iterationNumber, QList<PointS> &coord, QList<float> ws, int wIndex, int &addedPoints)
{
    QList<PointS> result;
    for(int index = 0; index < coord.count(); index++)
    {
        if(!coord.at(index)._tanAsigned)
        {
            PointS newP;
            newP._x = coord.at(index)._x;
            newP._y = coord.at(index)._y;
            newP._tan = coord.at(index)._tan;
            newP._tanAsigned = coord.at(index)._tanAsigned;
            newP._tanInf = coord.at(index)._tanAsigned;
            newP.SetTangent(GetAutomaticTangent(coord, index));
            coord.replace(index, newP);
        }
    }

    addedPoints = 0;
    for(int index = 0;index < coord.count()-1;index++)
    {
        int pos2 = index+1;
        //PointS qi = GetQi(coord.at(index), coord.at(pos2));
        double w = ws.at(index+wIndex);
        result.append(coord.at(index));
        addedPoints += SubdivideEdges(result, iterationNumber,coord.at(index),coord.at(pos2), w);
    }
    return result;
}

RailPolygon::RailPolygon(double hd, double fd, double bd, double hw, double w,
                         double bw, double r1h, double r2h, double rhb, double alpha,
                         double beta, double theta, QObject *parent):QObject(parent){

    //puntos.assign(3,vector<double>(2));

    this->r1h = r1h;
    this->r2h = r2h;
    this->rhb = rhb;
    a.setX(0.0);
    a.setY(bd+fd+hd);

    b.setY(bd+fd+hd);

    c.setX(hw/2);

    d.setX(w/2);

    e.setX(0);
    e.setY(bd+fd);

    c.setY(e.y() + (hw/2)*tan((alpha/180.0)*3.14));

    b.setX(c.x() - tan((beta/180.0)*3.14)*(a.y() - c.y()));

    d.setY(c.y() + tan((alpha/180.0)*3.14)*(d.x() - c.x()));

    f.setX(w/2); f.setY(bd);
    g.setX(bw/2);
    g.setY(f.y() - tan((theta/180.0)*3.14)*(g.x()-f.x()));

    h.setX(bw/2);h.setY(0.0);

    bp.setX(-b.x());bp.setY(b.y());
    cp.setX(-c.x());cp.setY(c.y());
    dp.setX(-d.x());dp.setY(d.y());
    fp.setX(-f.x());fp.setY(f.y());
    gp.setX(-g.x());gp.setY(g.y());
    hp.setX(-h.x());hp.setY(h.y());
    double maxX = h.x(),minX = hp.x();
    double maxY = b.y(),minY = h.y();

    radius = max((maxX - minX)/2, (maxY - minY)/2)+0.5;
    center.setX( (maxX + minX)/2);
    center.setY( -(maxY + minY)/2);
}


double RailPolygon::norm(QPointF p){
    return sqrt(p.x()*p.x() + p.y()*p.y());
}

QList<QPointF> RailPolygon::Bezcircle(QPointF A, QPointF B, QPointF C, double r, double &w){
    QPointF vBA = A-B;
    vBA/=norm(vBA);
    QPointF vBC = C-B;
    vBC/=norm(vBC);
    QPointF vbis = vBA + vBC;
    vbis/=norm(vbis);

    double tbis = abs(r/vbis.y());
    double tBA = abs(tbis*vbis.x());

    QPointF P1 = B;
    QPointF P0 = B + vBA*tBA;
    QPointF P2 = B + vBC*tBA;

    QPointF r20 = P2-P0;
    r20/=norm(r20);
    QPointF r10 = P1-P0;
    r10/=norm(r10);
    w = r20.x()*r10.x() + r20.y()*r10.y();

    QList<QPointF> res;
    res.append(P0); res.append(P1); res.append(P2);

    return res;
}


void RailPolygon::InitialPolygon(){
    QPointF P0,P1,P2;
    QList<QPointF> aux,Q;
    double w;

    aux = Bezcircle(a,b,c,r1h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0);Q.append(P1);Q.append(P2);
    polig.append(Q[0]);polig.append(P1); polig.append(Q[2]);
    ws.append(w);//w0

    aux = Bezcircle(b,c,d,r2h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0);Q.append(P1);Q.append(P2);
    QPointF C4 = P0;
    QPointF C2 = polig[2];
    QPointF C3 = (C2+C4)/2;
    polig.append(C3);polig.append(C4);ws.append(0);

    polig.append(P1);polig.append(Q[5]);
    ws.append(w);

    aux = Bezcircle(c,d,f,rhb,w); //Arco 3
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0);Q.append(P1);Q.append(P2);
    QPointF C8 = Q[6]; QPointF C6 = Q[5];
    QPointF C7 = (C6+C8)/2;
    polig.append(C7);polig.append(C8);ws.append(0);

    polig.append(P1);polig.append(Q[8]);
    ws.append(w); //w4

    aux = Bezcircle(d,f,g,rhb,w); //Arco 4
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0);Q.append(P1);Q.append(P2);
    QPointF C12 = Q[9]; QPointF C10 = Q[8];
    QPointF C11 = (C10+C12)/2;
    polig.append(C11); polig.append(C12);ws.append(0);//w5

    polig.append(P1); polig.append(Q[11]);
    ws.append(w);//w6

    aux = Bezcircle(f,g,h,r2h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C16 = Q[12]; QPointF C14 = Q[11];
    QPointF C15 = (C14+C16)/2;
    polig.append(C15); polig.append(C16); ws.append(0);

    polig.append(P1); polig.append(Q[14]);
    ws.append(w);

    QPointF orig(0,0);
    aux = Bezcircle(g,h,orig,r2h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C20 = Q[15]; QPointF C18 = Q[14];
    QPointF C19 = (C18+C20)/2;
    polig.append(C19); polig.append(C20); ws.append(0);

    polig.append(P1); polig.append(Q[17]);
    ws.append(w);

    aux = Bezcircle(orig,hp,gp,r2h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C24 = Q[18]; QPointF C22 = Q[17];
    QPointF C23 = (C22+C24)/2;
    polig.append(C23); polig.append(C24); ws.append(0);

    polig.append(P1); polig.append(Q[20]);
    ws.append(w);

    aux = Bezcircle(hp,gp,fp,r2h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C28 = Q[21]; QPointF C26 = Q[20];
    QPointF C27 = (C26+C28)/2;
    polig.append(C27);polig.append(C28); ws.append(0);

    polig.append(P1); polig.append(Q[23]);
    ws.append(w);

    aux = Bezcircle(gp,fp,dp,rhb,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C32 = Q[24]; QPointF C30 = Q[23];
    QPointF C31 = (C30+C32)/2;
    polig.append(C31); polig.append(C32);ws.append(0);

    polig.append(P1);polig.append(Q[26]);
    ws.append(w);

    aux = Bezcircle(fp,dp,cp,rhb,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C36 = Q[27]; QPointF C34 = Q[26];
    QPointF C35 = (C34+C36)/2;
    polig.append(C35);polig.append(C36);ws.append(0);

    polig.append(P1); polig.append(Q[29]);
    ws.append(w);

    aux = Bezcircle(dp,cp,bp,r2h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C40 = Q[30]; QPointF C38 = Q[29];
    QPointF C39 = (C38+C40)/2;
    polig.append(C39); polig.append(C40); ws.append(0);

    polig.append(P1); polig.append(Q[32]);
    ws.append(w);

    aux = Bezcircle(cp,bp,a,r1h,w);
    P0 = aux[0]; P1 = aux[1]; P2 = aux[2];
    Q.append(P0); Q.append(P1);Q.append(P2);
    QPointF C44 = Q[33]; QPointF C42 = Q[32];
    QPointF C43 = (C42+C44)/2;
    polig.append(C43);polig.append(C44);ws.append(0);

    polig.append(P1);polig.append(Q[35]);
    ws.append(w);

    QPointF C48 = Q[0]; QPointF C46 = Q[35];
    QPointF C47 = (C46+C48)/2;
    polig.append(C47);polig.append(C48);ws.append(0);    
    for(int i = 0; i < polig.count(); i++)
        CInitial.append(polig[i]);
    QFile file("prueba.txt");
    file.open(QFile::WriteOnly);
    QTextStream out(&file);
    out << "[";
    for(int i = 0; i < polig.count(); i++){
        out << "["<< polig[i].x() << "," << polig[i].y()<<"] ;" << "\n";
    }
    out << "];";
    file.close();

}



QList<QPointF> RailPolygon::SubdividePoints(QPointF A, QPointF B, QPointF C, int iter, double w, vector<double> &waux){
    waux.push_back(w);
    int exp = (int)pow(2,iter+1);
    QPointF mat[iter+1][exp+1];
    mat[0][0] = A; mat[0][1] = B; mat[0][2] = C;
    for(int j = 1; j <= iter; j++){
        int l = (int)pow(2,j+1);
        waux.push_back( sqrt((1+waux[j-1])/2));

        mat[j][0].setX(mat[j-1][0].x());
        mat[j][0].setY(mat[j-1][0].y());
        mat[j][l].setX(mat[j-1][l/2].x());
        mat[j][l].setY(mat[j-1][l/2].y());

        QPointF p1 = mat[j][0];
        QPointF p2 = mat[j][l];
        QPointF p3;

        for(int k = 1; k <= l; k+=2){
            if(k % 4 == 1){
                mat[j][k].setX((mat[j-1][(k-1)/2].x() + waux[j-1]*mat[j-1][(k+1)/2].x())/(1+waux[j-1]));
                mat[j][k].setY((mat[j-1][(k-1)/2].y() + waux[j-1]*mat[j-1][(k+1)/2].y())/(1+waux[j-1]));

            }
            else{
                mat[j][k].setX((mat[j-1][(k+1)/2].x() + waux[j-1]*mat[j-1][(k-1)/2].x())/(1+waux[j-1]));
                mat[j][k].setY((mat[j-1][(k+1)/2].y() + waux[j-1]*mat[j-1][(k-1)/2].y())/(1+waux[j-1]));
            }
            p3 = mat[j][k];
        }
        for(int k = 2; k < l; k+=2){
            mat[j][k].setX((mat[j][k-1].x() + mat[j][k+1].x())/2);
            mat[j][k].setY((mat[j][k-1].y() + mat[j][k+1].y())/2);
            p3 = mat[j][k];
        }
    }

    int m = (int)pow(2,iter)+1;
    QList<QPointF> result;
    for(int i = 0; i < m; i++){
        QPointF p;       
        p.setX(mat[iter][2*i].x());
        p.setY(mat[iter][2*i].y());
        result.append(p);
    }
    return result;
}

void RailPolygon::SubdividePolygon(int k){
    QList<QPointF> aux;
    int n = polig.size(), c = 0;
    vector<double> waux;    
    for(int i = 0; i < n-2; i+=2){
        QList<QPointF> temp = SubdividePoints(polig[i],polig[i+1],polig[i+2],k,ws[c],waux);
        c++;
        temp.removeAt(0);
        aux.append(temp);
    }
    polig = aux;
}

QList<QPointF> CreateTriangulation(QList<QPointF> CInitial, double *pts, int *fcs,QList<double> ws){
    QList<QPointF> adic;
    /*QFile file1("frontera.txt");
    file1.open(QFile::WriteOnly);
    QTextStream out1(&file1);
    for(int i = 0; i < polig.size(); i++){
        out1 << "c"<<i <<" " << polig[i].x() << " " << polig[i].y() << "\n";
    }*/
    //file1.close();
    pts[3] = CInitial[1].x();pts[4] = CInitial[1].y();//v2
    QPointF temp = (CInitial[1] + CInitial[4])/2;adic.append(temp);
    pts[6] = temp.x(); pts[7] = temp.y();//v3
    pts[9] = CInitial[5].x(); pts[10] = CInitial[5].y();//v4
    QPointF v5 = (CInitial[8] + 2*ws[4]*CInitial[9]+CInitial[10])/(2*(1+ws[4]));adic.append(v5);
    pts[12] = v5.x(); pts[13] = v5.y(); //pts[12] = CInitial[9].x(); pts[13] = CInitial[9].y();\\v5
    temp = (2*CInitial[9] + CInitial[13])/3; adic.append(temp);
    pts[15] = temp.x(); pts[16] = temp.y();//v6
    temp = (CInitial[9] + 2*CInitial[13])/3; adic.append(temp);
    pts[18] = temp.x(); pts[19] = temp.y();//v7
    QPointF v8 = (CInitial[12] + 2*ws[6]*CInitial[13] +CInitial[14])/(2*(1+ws[6]));
    pts[21] = v8.x(); pts[22] = v8.y();adic.append(v8);//pts[21] = CInitial[13].x(); pts[22] = CInitial[13].y();//v8
    temp = (CInitial[13] + CInitial[17])/2; adic.append(temp);
    pts[24] = temp.x(); pts[25] = temp.y();//v9
    pts[27] = CInitial[17].x(); pts[28] = CInitial[17].y();//v10
    pts[30] = CInitial[21].x(); pts[31] = CInitial[21].y();//v11
    QPointF v13(CInitial[9].x(),0),v11 = CInitial[21]; adic.append(v13);
    temp = (v11 + v13)/2; adic.append(temp);
    pts[33] = temp.x(); pts[34] = temp.y();//v12
    pts[36] = v13.x(); pts[37] = v13.y();//v13
    QPointF v14(CInitial[33].x(),0); adic.append(v14);
    pts[39] = v14.x(); pts[40] = v14.y();//v14
    QPointF v16 = CInitial[25];
    QPointF v15 = (v14 + v16)/2; adic.append(v15);
    pts[42] = v15.x(); pts[43] = v15.y();//v15
    pts[45] = CInitial[25].x(); pts[46] = CInitial[25].y();//v16
    pts[48] = CInitial[29].x(); pts[49] = CInitial[29].y();//v17
    temp = (CInitial[29] + CInitial[33])/2; adic.append(temp);
    pts[51] = temp.x(); pts[52] = temp.y();//v18
    pts[54] = -v8.x(); pts[55] = v8.y(); QPointF v19(-v8.x(),v8.y());adic.append(v19); //pts[54] = CInitial[33].x(); pts[55] = CInitial[33].y();//v19
    temp = (2*CInitial[33] + CInitial[37])/3; adic.append(temp);
    pts[57] = temp.x(); pts[58] = temp.y();//v20
    temp = (CInitial[33] + 2*CInitial[37])/3; adic.append(temp);
    pts[60] = temp.x(); pts[61] = temp.y();//v21
    pts[63] = -v5.x(); pts[64] = v5.y();QPointF v22(-v5.x(),v5.y());adic.append(v22);//pts[63] = CInitial[37].x(); pts[64] = CInitial[37].y();//v22
    pts[66] = CInitial[41].x(); pts[67] = CInitial[41].y();//v23
    temp = (CInitial[41] + CInitial[45])/2; adic.append(temp);
    pts[69] = temp.x(); pts[70] = temp.y();//v24
    pts[72] = CInitial[45].x(); pts[73] = CInitial[45].y();//v25
    QPointF v26 = (2*CInitial[45] + CInitial[1])/3;adic.append(v26);
    pts[75] = v26.x(); pts[76] = v26.y();//v26
    QPointF v1 = (CInitial[45] + 2*CInitial[1])/3;adic.append(v1);
    pts[0] = v1.x(); pts[1] = v1.y();//v1
    temp = (v1 + CInitial[9])/2; adic.append(temp);
    pts[78] = temp.x(); pts[79] = temp.y();//v27
    temp = (v26 + CInitial[37])/2; adic.append(temp);//v28
    pts[81] = temp.x(); pts[82] = temp.y();
    for(int i = 2; i < 84; i+=3)
        pts[i] = 0;

    /*QFile file2("vertices.txt");
    file2.open(QFile::WriteOnly);
    QTextStream out2(&file2);
    out2<< "[";
    for(int i = 0; i < 84; i+=3){
        out2   <<"[" << pts[i] << ", " << pts[i+1] << "];" << "\n";
    }
    out2 << "];";
    file2.close();*/

    fcs[0] = 0; fcs[1] = 2; fcs[2] = 1;
    fcs[3] = 0; fcs[4] = 26; fcs[5] = 2;
    fcs[6] = 0; fcs[7] = 25; fcs[8] = 27;
    fcs[9] = 0; fcs[10] = 27; fcs[11] = 26;
    fcs[12] = 25; fcs[13] = 24; fcs[14] = 23;
    fcs[15] = 25; fcs[16] = 23; fcs[17] = 27;
    fcs[18] = 2; fcs[19] = 26; fcs[20] = 3;
    fcs[21] = 26; fcs[22] = 4; fcs[23] = 3;
    fcs[24] = 26; fcs[25] = 27; fcs[26] = 21;
    fcs[27] = 26; fcs[28] = 21; fcs[29] = 4;
    fcs[30] = 27; fcs[31] = 23; fcs[32] = 22;
    fcs[33] = 27; fcs[34] = 22; fcs[35] = 21;
    fcs[36] = 4; fcs[37] = 21; fcs[38] = 20;
    fcs[39] = 4; fcs[40] = 20; fcs[41] = 5;
    fcs[42] = 5; fcs[43] = 20; fcs[44] = 19;
    fcs[45] = 5; fcs[46] = 19; fcs[47] = 6;
    fcs[48] = 6; fcs[49] = 19; fcs[50] = 18;
    fcs[51] = 6; fcs[52] = 18; fcs[53] = 7;
    fcs[54] = 7; fcs[55] = 18; fcs[56] = 13;
    fcs[57] = 7; fcs[58] = 13; fcs[59] = 12;
    fcs[60] = 7; fcs[61] = 12; fcs[62] = 11;
    fcs[63] = 7; fcs[64] = 11; fcs[65] = 8;
    fcs[66] = 9; fcs[67] = 8; fcs[68] = 11;
    fcs[69] = 9; fcs[70] = 11; fcs[71] = 10;
    fcs[72] = 18; fcs[73] = 14; fcs[74] = 13;
    fcs[75] = 18; fcs[76] = 17; fcs[77] = 14;
    fcs[78] = 17; fcs[79] = 16; fcs[80] = 14;
    fcs[81] = 16; fcs[82] = 15; fcs[83] = 14;

    return adic;
}

double toDoubleAngle(QString angle){
    double res = 0;
    int i = 0;QString ang = "";
    while(angle[i].isDigit()){
        ang.append(angle[i++]);
    }
    res += ang.toDouble();
    i++;
    QString minut = "";
    while(angle[i].isDigit()){
        minut.append(angle[i++]);
    }
    double min = minut.toDouble();
    res += (min/60);
    return res;
}

double toDoubleLength(QString length){
    double res = 0;
    QStringList list = length.split('+',QString::SkipEmptyParts);
    res += ((QString)list[0]).toDouble();
    QString right = list[1]; right.remove('('); right.remove(')');
    QStringList frac = right.split('/');
    QString denom = frac[1], num = frac[0];
    double d = denom.toDouble();
    double n = num.toDouble();
    res += (d == 0)?0:n/d;
    return res;
}
