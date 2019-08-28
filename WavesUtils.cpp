#include "WavesUtils.h"
#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QMessageBox>
#include <QtGui>
#include <climits>
#include <ctime>
#include <float.h>

AvalAvect::AvalAvect(const AvalAvect &val, QObject *parent):QObject(parent){
    index = val.GetIndex();
    autoval = val.GetAutoVal();
    autovect = val.GetAutoVect();
    waveNumber = val.GetWaveNumber();
}

AvalAvect::AvalAvect(double pwaveNumber, int pindex, double aVal, QList<complex<double> > aVect,QObject* parent):QObject(parent)
{
    index = pindex;
    autoval = aVal;
    autovect = aVect;
    waveNumber = pwaveNumber;
}

AvalAvect& AvalAvect::operator =(const AvalAvect& val){
    index = val.GetIndex();
    autoval = val.GetAutoVal();
    autovect = val.GetAutoVect();
    waveNumber = val.GetWaveNumber();
    return *this;
}

double AvalAvect::AVectRealMax()
{
    return max(max(autovect[0].real(),autovect[1].real()),autovect[2].real());
}

int AvalAvect::GetIndex()const
{
    return index;
}
double AvalAvect::GetAutoVal()const
{
    return autoval;
}
double AvalAvect::GetWaveNumber()const
{
    return waveNumber;
}
void AvalAvect::Normalize2AVector()
{
    complex<double> *x = ConvertAvectsToPointer();
    double norm = AvectorNorm(x);

    int n = autovect.size();
    for(int i = 0; i < n; i ++)
    {
        autovect[i].real(x[i].real()/norm);
        autovect[i].imag(x[i].imag()/norm);
    }
}
void AvalAvect::Normalize3AVector()
{
    complex<double> *x = ConvertAvectsToPointer();
    complex<double> normX = AvectorNorm(0,3);
    complex<double> normY = AvectorNorm(1,3);
    complex<double> normZ = AvectorNorm(2,3);

    double nx = abs(normX);
    double ny = abs(normY);
    double nz = abs(normZ);
    int n = autovect.size();
    for(int i = 0; i < n; i+=3)
    {
        autovect[i].real(x[i].real()/nx);
        autovect[i].imag(x[i].imag()/nx);
        autovect[i+1].real(x[i+1].real()/ny);
        autovect[i+1].imag(x[i+1].imag()/ny);
        autovect[i+2].real(x[i+2].real()/nz);
        autovect[i+2].imag(x[i+2].imag()/nz);
    }
}
void AvalAvect::SetAutoVal(double aVal)
{
    autoval = aVal;
}
QList<complex<double> > AvalAvect::GetAutoVect()const
{
    return autovect;
}
double AvalAvect::AvectorNorm(complex<double> *x)
{
    int uno = 1;
    int n = autovect.size();
    return dznrm2_(&n, x, &uno);
}
complex<double>* AvalAvect::ConvertAvectsToPointer()
{
    complex<double> *result = new complex<double>[autovect.size()];
    for(int i = 0; i < autovect.size(); i++)
        result[i] = autovect.at(i);
    return result;
}
void AvalAvect::AddValueToAutoVect(complex<double> aVal)
{
    autovect.push_back(aVal);
}
void AvalAvect::SetAutoVect(QList<complex<double> > aVect)
{
    autovect = aVect;
}
complex<double> AvalAvect::AvectorNorm(int initPos, int inc)
{
    complex<double> bigger;
    bigger.real(autovect[initPos].real());
    bigger.imag(autovect[initPos].imag());

    for(int i = initPos; i < autovect.size()/3; i+=inc )
    {
        complex<double> iElem;
        iElem.real(autovect[i].real());
        iElem.imag(autovect[i].imag());
        if(abs(iElem) > abs(bigger))
            bigger = iElem;
    }
    return bigger;
}
void AvalAvect::SetWaveNumber_Index(double pwaveNumber, int pindex)
{
    index = pindex;
    waveNumber = pwaveNumber;
}
void AvalAvect::ChangeValueInAutoVect(complex<double> aVect, int pos)
{
    if(pos >= autovect.count()) return;
    autovect[pos] = aVect;
}
QList<GVector3> AvalAvect::GetDisplacements(double z, double epsilon, int t)
{
    epsilon *= 10;
    QList<GVector3> result;
    for(int k = 0; k < autovect.size()/3; k++)
    {
        double x = epsilon*(autovect[3*k].real()*cos(waveNumber*z-autoval*(t)) - autovect[3*k].imag()*sin(waveNumber*z-autoval*(t)));
        double y = epsilon*(autovect[3*k+1].real()*cos(waveNumber*z-autoval*(t)) - autovect[3*k+1].imag()*sin(waveNumber*z-autoval*(t)));
        double zeta = epsilon*(autovect[3*k+2].real()*cos(waveNumber*z-autoval*(t)) - autovect[3*k+2].imag()*sin(waveNumber*z-autoval*(t)));

        result.push_back(GVector3(x, y, zeta));
    }
    return result;
}


CVertex::CVertex(const CVertex &vert, QObject* parent):QObject(parent){
    sectionId = vert.GetSectionId();
    x = vert.X();
    y = vert.Y();
    z = vert.Z();
    displacement = vert.Displacement();
}

CVertex& CVertex::operator =(const CVertex& vert){
    sectionId = vert.GetSectionId();
    x = vert.X();
    y = vert.Y();
    z = vert.Z();
    displacement = vert.Displacement();
    return *this;
}

CVertex::CVertex(double pX, double pY, double pZ, int pSectionId,QObject* parent):QObject(parent)
{
    x = pX;
    y = pY;
    z = pZ;
    sectionId = pSectionId;
}

double CVertex::X()const{return x;}
double CVertex::Y()const{return y;}
double CVertex::Z()const{return z;}
int CVertex::GetSectionId()const
{
    return sectionId;
}
void CVertex::X(double newX)
{
    x = newX;
}
void CVertex::Y(double newY)
{
    y = newY;
}
void CVertex::Z(double newZ)
{
    z = newZ;
}
GVector3* CVertex::Displacement()const
{
    return displacement;
}
void CVertex::Displacement(GVector3 *value)
{
    displacement = value;
}


Cuad::Cuad(CVertex v0, CVertex v1, CVertex v2, CVertex v3,QObject* parent):QObject(parent)
{
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    vertices[3] = v3;
}

Cuad::Cuad(const Cuad &cuad, QObject *parent):QObject(parent){
    vertices[0] = cuad.Vertice(0);
    vertices[1] = cuad.Vertice(1);
    vertices[2] = cuad.Vertice(2);
    vertices[3] = cuad.Vertice(3);
}

Cuad& Cuad::operator =(const Cuad& cuad){
    vertices[0] = cuad.Vertice(0);
    vertices[1] = cuad.Vertice(1);
    vertices[2] = cuad.Vertice(2);
    vertices[3] = cuad.Vertice(3);
    return *this;
}

CVertex Cuad::Vertice(int index)const
{
    if(index == 0)
        return vertices[0];
    else if(index == 1)
        return vertices[1];
    else if(index == 2)
        return vertices[2];
    else if(index == 3)
        return vertices[3];
    return CVertex();
}
void Cuad::Vertice(int index, CVertex newValue)
{
    if(index == 0)
        vertices[0] = newValue;
    else if(index == 1)
        vertices[1] = newValue;
    else if(index == 2)
        vertices[2] = newValue;
    else if(index == 3)
        vertices[3] = newValue;
}


Section::Section(QList<CVertex*> pVertices, int pId,QObject* parent):QObject(parent)
{
    id = pId;
    vertices = pVertices;
}

int Section::VerticesCount()
{
    return vertices.count();
}
CVertex* Section::Vertice(int index)
{
    return vertices[index];
}


MeshCuad::MeshCuad(QList<QList<GVector2> > bound, QList<QList<int> > index, double pdeltaZ,
                   int pnumberSections, QList<QList<QList<GVector3*> > > pDisplacements,QObject*parent):QObject(parent)
{
    //Crear todos los vertices
    indexes = index;
    deltaZ = pdeltaZ;
    numberSections = pnumberSections;

    CreateMesh(bound, pDisplacements);
}

void MeshCuad::AdjustDeltha()
{
    for(int i = 0 ; i < BoundariesCount(); i++)
    {
        double z = 0;
        for(int j = 0; j < SectionsCount(); j++)
        {
            Section *jSection = BoundaryISectionJ(i,j);
            for(int k = 0; k < jSection->VerticesCount(); k++)
                jSection->Vertice(k)->Z(z);
            z += deltaZ;
        }
    }
}
void MeshCuad::AdjustSections()
{
    if(SectionsCount() == numberSections) return;
    double z = deltaZ*SectionsCount();
    //Sections were added
    if(SectionsCount() < numberSections)
    {
        int sectionsToAdd = numberSections - SectionsCount();
        //Adding new sections
        for(int i = 0; i < BoundariesCount(); i++)
        {
            for(int j = 0; j < sectionsToAdd; j++)
            {
                QList<CVertex*> newSectionVertices;
                for(int k = 0; k < BoundaryISectionJ(i,j)->VerticesCount(); k++)
                {
                    CVertex *cvertexK = new CVertex(BoundaryISectionJ(i,j)->Vertice(k)->X(), BoundaryISectionJ(i,j)->Vertice(k)->Y(), z, SectionsCount());
                    cvertexK->Displacement(new GVector3(0,0,0));
                    newSectionVertices.append(cvertexK);
                }
                sections[i].append(new Section(newSectionVertices,SectionsCount()));
                z += deltaZ;
            }
            z -= deltaZ*sectionsToAdd;
        }
        z = deltaZ*(SectionsCount() - 1);
    }
    //sections were removed
    else
    {
        int sectionsToRemove =  SectionsCount() - numberSections;
        //deleting sections
        for(int i = 0; i < BoundariesCount();i++)
        {
            for(int j = 0 ; j < sectionsToRemove; j++)
                sections[i].removeLast();
        }
    }
}
void MeshCuad::CreateMesh(QList<QList<GVector2> > pbound, QList<QList<QList<GVector3*> > > pdisplacements)
{
    for(int k = 0; k < pbound.count(); k++)
    {
        double z = 0;
        sections.append(QList<Section*>());
        for(int i = 0; i < numberSections; i++)
        {
            QList<CVertex*> vertices;
            for(int j = 0; j < pbound[k].size(); j++)
            {
                CVertex *v = new CVertex(pbound[k][j].ax, pbound[k][j].ay, z, i);
                v->Displacement(pdisplacements[k][i][j]);
                vertices.push_back(v);
            }
            sections[k].append(new Section(vertices,i));
            z += deltaZ;
        }
    }
}

int MeshCuad::SectionsCount()
{
    return sections[0].count();
}
int MeshCuad::BoundariesCount()
{
    return sections.count();
}
void MeshCuad::NoDisplacements()
{
    for(int i = 0; i < sections.count(); i++)
    {
        for(int j = 0; j < sections.at(i).count();j++)
        {
            for(int k = 0; k < sections.at(i).at(j)->VerticesCount(); k++)
            {
                BoundaryISectionJ(i,j)->Vertice(k)->Displacement()->ax =
                BoundaryISectionJ(i,j)->Vertice(k)->Displacement()->ay =
                BoundaryISectionJ(i,j)->Vertice(k)->Displacement()->az = 0;
            }
        }
    }
}
QList<Cuad> MeshCuad::CrossSection()
{
    QList<Cuad> result;
    if(BoundariesCount() == 2)
    {
        Section *innerSection = BoundaryISectionJ(0,0);
        Section *outterSection = BoundaryISectionJ(1,0);
        for(int k = 0; k < innerSection->VerticesCount()-1; k++)
            result.append(Cuad(*innerSection->Vertice(k), *innerSection->Vertice(k+1),
                               *outterSection->Vertice(k+1), *outterSection->Vertice(k)));
        result.append(Cuad(*innerSection->Vertice(innerSection->VerticesCount()-1), *innerSection->Vertice(0),
                           *outterSection->Vertice(0), *outterSection->Vertice(innerSection->VerticesCount()-1)));
    }
    return result;
}
void MeshCuad::SectionsNumber(int val)
{
    numberSections = val;
    AdjustSections();
}
void MeshCuad::DelthaValue(double val)
{
    deltaZ = val;
    AdjustDeltha();
}
Section* MeshCuad::BoundaryISectionJ(int i, int j)
{
    return sections[i][j];
}
QList<Cuad> MeshCuad::BoundaryICuadsSectionJ(int i, int j)
{
    QList<Cuad> result;
    if(j == 0 && BoundariesCount() == 2)
    {
        Section *innerSection = BoundaryISectionJ(0,0);
        Section *outterSection = BoundaryISectionJ(1,0);
        for(int k = 0; k < innerSection->VerticesCount()-1; k++)
            result.append(Cuad(*innerSection->Vertice(k), *innerSection->Vertice(k+1),
                               *outterSection->Vertice(k+1), *outterSection->Vertice(k)));
        result.append(Cuad(*innerSection->Vertice(innerSection->VerticesCount()-1), *innerSection->Vertice(0),
                           *outterSection->Vertice(0), *outterSection->Vertice(innerSection->VerticesCount()-1)));
    }
    Section *jSection = BoundaryISectionJ(i,j);
    Section *jPlusOne = BoundaryISectionJ(i,j+1);
    for(int k = 0; k < jSection->VerticesCount()-1; k++)
        result.append(Cuad(*jSection->Vertice(k), *jPlusOne->Vertice(k), *jPlusOne->Vertice(k+1), *jSection->Vertice(k+1)));
    result.append(Cuad(*jSection->Vertice(jSection->VerticesCount()-1), *jPlusOne->Vertice(jPlusOne->VerticesCount()-1),
                       *jPlusOne->Vertice(0), *jSection->Vertice(0)));
    if(j == SectionsCount()-2 && BoundariesCount() == 2)
    {
        Section *innerSection = BoundaryISectionJ(0,SectionsCount()-1);
        Section *outterSection = BoundaryISectionJ(1,SectionsCount()-1);
        for(int k = 0; k < innerSection->VerticesCount()-1; k++)
            result.append(Cuad(*innerSection->Vertice(k), *innerSection->Vertice(k+1),
                               *outterSection->Vertice(k+1), *outterSection->Vertice(k)));
        result.append(Cuad(*innerSection->Vertice(innerSection->VerticesCount()-1), *innerSection->Vertice(0),
                           *outterSection->Vertice(0), *outterSection->Vertice(innerSection->VerticesCount()-1)));
    }
    return result;
}
void MeshCuad::SetBoundaryISectionJ(int i, int j, Section *value)
{
    sections[i][j] = value;
}
void MeshCuad::ChangeDisplacements(QList<QList<QList<GVector3> > > newDisplacements)
{
    for(int i = 0; i < BoundariesCount(); i++)
    {
        for(int j = 0; j <SectionsCount(); j++)
        {
            for(int k = 0 ; k < BoundaryISectionJ(i,j)->VerticesCount(); k++)
            {
                BoundaryISectionJ(i,j)->Vertice(k)->Displacement()->ax = newDisplacements[i][j][k].ax;
                BoundaryISectionJ(i,j)->Vertice(k)->Displacement()->ay = newDisplacements[i][j][k].ay;
                BoundaryISectionJ(i,j)->Vertice(k)->Displacement()->az = newDisplacements[i][j][k].az;
            }
        }
    }
}


//Versión propuesta por profesor Chistopher Beattie
void WorkerThread::doBeattieWork()
{
    /*
      Last revision September 2015
      This is the main program to solve the wave popagation
      problem using SAFEM with linear elements. In this
      version we apply Christopher Beattie's suggestions to
      tranform the generalized eigenvalue problem into a
      singular value problem....
    */
    double minimumDifference = GetMinimumDifference(triangles, x, y);
    ScaleCoordinates(minimumDifference, x, y);
    SaveVertices(x, y, "C://Users/Portege/Desktop/vertex.txt");
    SaveElems(triangles, "C://Users/Portege/Desktop/elems.txt");

    qDebug()<<"Beattie's Solution";
    clock_t start = clock();
    int NNT = x.size(); //cantidad de vertices
    int NET = triangles.size(); //cantidad de triangulos

    //////////Datos del material elastico//////////
    double c44 = young/2/(1+poisson);
    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
    double c11 = c12+2*c44;
    double vt = sqrt(c44/density);
    double f[3] ;
    f[0] = c11/c44;
    f[1] = c12/c44;
    f[2] = 1;

    FMatrix kk1(9*NET, 9);
    FMatrix kk2(9*NET, 9);
    FMatrix kk3(9*NET, 9);
    FMatrix mg(3*NNT, 3*NNT);
    FMatrix Mup(3*NNT, 3*NNT);
    FMatrix H(9, 3*NNT);

    /*
      Loop over the elements to compute local maytrices K1e,K2e,K3e and Me
      and global mass matrix MG
    */
    FMatrix *k1 = new FMatrix(9, 9);
    FMatrix *k2 = new FMatrix(9, 9);
    FMatrix *k3 = new FMatrix(9, 9);
    FMatrix *m = new FMatrix(9, 9);
    double dj = 0;
    for(int NE = 0; NE < NET; NE++)
    {
        //vertices of the e-th elements are
        QList<int> elemI = triangles.at(NE);
        double coord[6] = {x[elemI.at(0)], y[elemI.at(0)],
                           x[elemI.at(1)], y[elemI.at(1)],
                           x[elemI.at(2)], y[elemI.at(2)]};

        //Compute the local stiffness and mass matrices
        double toCheck = MatrixK1K2K3MT3N(coord, f, k1, k2, k3, m);
        if(NE == 0)
            dj = toCheck;
        if(Sign(dj)!= Sign(toCheck))
        {
            qDebug()<<"ERROR EN EL JACOBIANO";
            emit finished();
        }

        //Fill matrices KK1,KK2 and KK3 where we store local matrices K1e,K2e,K3e
        int i1 = 9*NE;
        for(int i = 0; i < 9; i++)
        {
            for(int j = 0; j < 9; j++)
            {
                kk1.SetElem(i1+i,j,k1->GetElem(i,j));
                kk2.SetElem(i1+i,j,k2->GetElem(i,j));
                kk3.SetElem(i1+i,j,k3->GetElem(i,j));
            }
        }
        k1->WriteMatrix("C://Users/Portege/Desktop/k1.txt");
        k2->WriteMatrix("C://Users/Portege/Desktop/k2.txt");
        k3->WriteMatrix("C://Users/Portege/Desktop/k3.txt");
        /*
          Update global mass matrix MG=MG+Ge'*Me*Ge
          We use the function xtimesGe that multiples a row vector by Ge
          taking into account that Ge'*Me*Ge=((Me*Ge)'*Ge)'=(Me*Ge)'*Ge
          H=Me*Ge
          H1=K1e*Ge, H2=K2e*Ge, H3=K3e*Ge
        */
        for(int i = 0; i < 9; i++)
        {
            QList<double> xtimesge = XTimesGe(m, i, elemI, NNT, false);
            for(int j = 0; j < 3*NNT; j++)
                H.SetElem(i,j,xtimesge[j]);
        }

        for(int i = 0; i < 3*NNT; i++)
        {
            QList<double> xtimesge = XTimesGe(&H, i, elemI, NNT, true);
            for(int j = 0; j < 3*NNT; j++)
                Mup.SetElem(i,j,xtimesge[j]);
        }
        mg = mg + Mup.transpuesta();

        k1->Clear();
        k2->Clear();
        k3->Clear();
        m->Clear();
    }
    // To ensure the symmetry of mass matrix we average the
    // values up and below the diagonal
    mg = (mg + mg.transpuesta()).multi_escalar(0.5);


    kk1.WriteMatrix("C://Users/Portege/Desktop/kk1.txt");
    kk2.WriteMatrix("C://Users/Portege/Desktop/kk2.txt");
    kk3.WriteMatrix("C://Users/Portege/Desktop/kk3.txt");
//    mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");

    // Compute the Cholesky factorization
    // of the global mass matrix MG
    // MG=L^t*L
    // L is upper triangular matrix
    //zpotrf
    //mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");
    char *uplo = new char[1];
    uplo[0] = 'U';
    int info = 0;
    complex<double> *mgValues = mg.ToComplexArray();
    zpotrf_(uplo, &mg.rows, mgValues, &mg.cols, &info);
    if(info != 0)
    {
        qDebug()<<"ERROR EN LA FACTORIZACION DE CHOLESKY(linea 532). info="<<info;
        emit finished();
        return;
    }
    FComplexMatrix L(mgValues, mg.rows, mg.cols);

    NXIh.clear();
    for(double i = step; i <= maxWavesNumber; i+=step)
        NXIh.push_back(i);

    emit maxIteration(NXIh.size());
    for(int i = 0; i < NXIh.size(); i++)
    {
        eigenValsVects.append(QList<AvalAvect>());
        for(int j=0;j<autovalNumber;j++)
            eigenValsVects[i].append(AvalAvect(-1, -1, 0, QList<complex<double> >()));
    }

    // Loop on values of the wave number xi
    // Solving the generalized eigenvalue
    // problem for each value of xi
    int mDim = 9*NET;
    FComplexMatrix *k2giXIh = new FComplexMatrix(9, 9);
    FComplexMatrix *suma = new FComplexMatrix(9, 9);
    FMatrix *k3gXIh2 = new FMatrix(9, 9);
    FComplexMatrix *Ke = new FComplexMatrix(9, 9);
    FComplexMatrix *Us = new FComplexMatrix(mDim, mDim);
    int sCount = min(mDim, 3*NNT);
    double *S = new double[sCount];
    FComplexMatrix *Vs = new FComplexMatrix(3*NNT, 3*NNT);
    int lwork = max(1,2*min(mDim,3*NNT)+max(mDim,3*NNT));
    complex<double> *work = new complex<double>[lwork];
    double *rwork = new double[5*min(mDim,3*NNT)];
    char *jobu = new char[1];
    jobu[0] = 'N';
    char *jobvt = new char[1];
    jobvt[0] = 'S';
    for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
    {
        double XIh = NXIh[INXIh];
        qDebug()<<XIh;

        /* The global stiffness matrix K ( of order 3nv x 3nv)
           satisfies K=R^t*R
           This is NOT the Cholesky factorization of K since R
           is a rectangular matrix (that depends on xi)
        */
        FComplexMatrix R(mDim, 3*NNT);


        // A = K1G+i*xih*K2G+xih^2*K3G;
        // Loop into the elements of the triangulation
        for(int e = 0; e < NET; e++)
        {
            /* Compute the local stiffness matrix Ke
               corresponding to the e-th element
               and the wave numver xih
               Ke=K1e+i*XIh*K2e+XIh^2*K3e
               where i is the imaginary unit
            */
            int i1 = 9*e;
            for(int i = 0; i < 9; i++)
            {
                for(int j = 0; j < 9; j++)
                {
                    k1->SetElem(i,j,kk1.GetElem(i1+i,j));
                    k2->SetElem(i,j,kk2.GetElem(i1+i,j));
                    k3->SetElem(i,j,kk3.GetElem(i1+i,j));
                }
            }

            complex<double> iXIh;
            iXIh.imag(XIh);

            k2->multi_escalar(iXIh, k2giXIh);
            k3->multi_escalar(pow(XIh, 2), k3gXIh2);
            k3gXIh2->Add(k2giXIh, suma);
            k1->Add(suma, Ke);
            Ke->Hermitize();

//            Ke->WriteRealComplexMatrix("C://Users/Portege/Desktop/KeHR.txt");
//            Ke->WriteImagComplexMatrix("C://Users/Portege/Desktop/KeHI.txt");
            // Compute the Cholesky factorization of Ke
            // Ke= Re'*Re
            complex<double> *kevalues = Ke->values;
            zpotrf_(uplo, &Ke->rows, kevalues, &Ke->cols, &info);
            if(info != 0)
            {
                qDebug()<<"ERROR EN LA FACTORIZACIÓN DE CHOLESKY(línea 621). info="<<info <<"chi="<<XIh;
                //emit finished();
                //return;
                continue;
            }
            FComplexMatrix Re(kevalues, Ke->rows, Ke->cols);
//            Re.WriteRealComplexMatrix("C://Users/Portege/Desktop/ReR.txt");
//            Re.WriteImagComplexMatrix("C://Users/Portege/Desktop/ReI.txt");

            // Fill the 9 rows of matrix R corresponding to
            // the e-th element vertices of the e-th
            // elements are
            QList<int> elemI = triangles.at(e);
            for(int i = 0; i < 9; i++)
            {
                QList<complex<double> > xtimesge = XTimesGe(&Re, i, elemI, NNT, false);
                for(int j = 0; j < 3*NNT; j++)
                    R.SetElem(i1+i-1, j, xtimesge[j]);
            }
        }


        /*
          Compute the (9*nt x 3nv ) matrix B = RL^{-1}.
          For k=1,...,9*nt, the k-th row of B, denoted here
          by xk, is computed solving the linear
          system L'*xk' = (R(k,:))' by forwards substitution
        */
//        R.WriteRealComplexMatrix("C://Users/Portege/Desktop/RR.txt");
//        R.WriteImagComplexMatrix("C://Users/Portege/Desktop/RI.txt");
        complex<double> *bValues = new complex<double>[R.rows*R.cols];
        Copy(R.values, bValues, R.rows*R.cols);
        FComplexMatrix B(bValues, mDim, 3*NNT);
        /* Solve the linear system
           L'*xk'=(R(k,:))'
           Assign xk' to the k-th row of matrix B
           of order 9nt x 3nv
        */
        zpotrs_(uplo, &L.rows, &B.cols, L.values, &L.cols, B.values, &B.rows, &info);
        if(info != 0)
        {
            qDebug()<<"ERROR EN LA SOLUCIÓN DEL SISTEMA(línea 644). info="<<info;
            emit finished();
            return;
        }

        /*
           Compute the smallest singular values of B and the
           corresponding right singular vectors
           sqrt(eig(B'*B))=singv(B)
        */

        B.WriteRealComplexMatrix("C://Users/Portege/Desktop/BR.txt");
        B.WriteImagComplexMatrix("C://Users/Portege/Desktop/BI.txt");
        //SVD(B, Us, S, Vs);
        zgesvd_(jobu, jobvt, &B.rows, &B.cols, B.values, &B.rows, S,
                Us->values, &B.rows, Vs->values, &B.cols, work, &lwork,
                rwork, &info);
        if(info != 0)
        {
            qDebug()<<"ERROR EN LA DESCOMPOSICIÓN SVD (línea 664). info="<<info;
            emit finished();
            return;
        }

        QList<double> auxD = ToDoubleList(S, sCount);
        //EIGENVALUES
        for(int i = 0; i < autovalNumber; i++)
        {
            eigenValsVects[INXIh][i].SetAutoVal(S[sCount-1 - i]);
            eigenValsVects[INXIh][i].SetWaveNumber_Index(XIh,i);
            //Compute the eigenvector from the right
            //singular vectors.
            //Compute the eigenvector u_k of the generalized problem
            //K*uk= omega*MG*uk
            //solving the linear system by backward substitution
            //L*uk=Vs(:,3*nv-k+1)
            //zpotrs
            //QList<complex<double> > uk = LinSolver(L, Vs->Diagonal(), 3*NET-i-1);
            complex<double> *ithValues = Vs->GetColumn(sCount-1 - i);
            int one = 1;
            zpotrs_(uplo, &L.rows, &one, L.values, &L.cols, ithValues, &L.rows, &info);
            if(info != 0)
            {
                qDebug()<<"ERROR EN LA SOLUCIÓN DEL SISTEMA(línea 677). info="<<info;
                emit finished();
                return;
            }
            eigenValsVects[INXIh][i].SetAutoVect(ToComplexList(ithValues, L.rows));
        }
        k2giXIh->Clear();
        suma->Clear();
        k3gXIh2->Clear();
        Ke->Clear();
        Us->Clear();
        Vs->Clear();
        emit posChanged(INXIh+1);
    }

    delete(k1);
    delete(k2);
    delete(k3);
    delete(m);
    delete(k2giXIh);
    delete(suma);
    delete(k3gXIh2);
    delete(Ke);
    delete(Us);
    delete(Vs);
    delete(uplo);
    delete(work);
    delete(rwork);
    delete(jobu);
    delete (jobvt);

    QList<int> deleted = PostProcessResults();
    PostProcessEigenValVects(deleted,NXIh.size(),autovalNumber);
    autovalNumber = autovalNumber - deleted.size();
    vf = FMatrix(NXIh.size(),autovalNumber);
    fM = FMatrix(NXIh.size(),autovalNumber);

    for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
    {
        for(int index = 0; index < autovalNumber; index++)
        {
            double val = eigenValsVects[INXIh][index].GetAutoVal();
            double aux1 = val/NXIh[INXIh];
            vf.SetElem(INXIh, index, aux1);     // Velocidad de fase/velocidad transversal (m/s)
            double aux2 = val*vt/(2*3.1415);
            fM.SetElem(INXIh, index, aux2);     // Frecuencia*h (Lado de la Seccion, o diametro  (Hz)
        }
    }
    vg = Diff(NXIh.size(),autovalNumber,eigenValsVects).multi_escalar(vt/step);

    clock_t end = clock();
    double diffticks = end-start;
    secTime =(diffticks)/(CLOCKS_PER_SEC);
    emit finished();
}
//[deprecated]
//M y A FULL Cuadratico
void WorkerThread::doQuadraticWork()
{
    qDebug()<<"Cuadrático [deprecated] (M y A FULL)";
    clock_t start = clock();
    int nnt = x.size(); //cantidad de vertices
    int NET = triangles.size(); //cantidad de triangulos

    //////////Datos del material elastico//////////
    double c44 = young/2/(1+poisson);
    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
    double c11 = c12+2*c44;
    double vt = sqrt(c44/density);
    double f[3];
    f[0] = c11/c44;
    f[1] = c12/c44;
    f[2] = 1;

    int nextIndex = nnt;
    QList<CuadraticVertexData> allMidpoints;
    QList<QList<int> > edges;
    for(int NE = 0; NE < NET; NE++)
    {
        QList<int> elemI = triangles.at(NE);
        int cuadraticElemI[6];
        cuadraticElemI[0] = elemI[0];
        cuadraticElemI[1] = elemI[1];
        cuadraticElemI[2] = elemI[2];

        double coord[12];
        coord[0] = x[elemI.at(0)];
        coord[1] = y[elemI.at(0)];
        coord[2] = x[elemI.at(1)];
        coord[3] = y[elemI.at(1)];
        coord[4] = x[elemI.at(2)];
        coord[5] = y[elemI.at(2)];

        int index = Contains(edges,elemI[0],elemI[1]);
        if(index != -1)
        {
            //Si ya existia esa arista, se guarda en la posicion 3 el indice que referencia al punto medio de ella
            cuadraticElemI[3] = edges[index][2];
        }
        else
        {
            cuadraticElemI[3] = nextIndex++;

            coord[6] = (coord[0]+coord[2])/2;
            coord[7] = (coord[1]+coord[3])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI[3],coord[6], coord[7]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(0));
            edgeI.push_back(elemI.at(1));
            edgeI.push_back(cuadraticElemI[3]);
            edges.push_back(edgeI);
        }
        index = Contains(edges,elemI.at(1),elemI.at(2));
        if(index != -1)
        {
            cuadraticElemI[4] = edges[index][2];
        }
        else
        {
            cuadraticElemI[4] = nextIndex++;
            coord[8] = (coord[2]+coord[4])/2;
            coord[9] = (coord[3]+coord[5])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI[4], coord[8], coord[9]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(1));
            edgeI.push_back(elemI.at(2));
            edgeI.push_back(cuadraticElemI[4]);
            edges.push_back(edgeI);
        }
        index = Contains(edges,elemI.at(2),elemI.at(0));
        if(index != -1)
        {
            cuadraticElemI[5] = edges[index][2];
        }
        else
        {
            cuadraticElemI[5] = nextIndex++;
            coord[10] = (coord[4]+coord[0])/2;
            coord[11] = (coord[5]+coord[1])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI[5], coord[10], coord[11]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(2));
            edgeI.push_back(elemI.at(0));
            edgeI.push_back(cuadraticElemI[5]);
            edges.push_back(edgeI);
        }
    }

    ////////////Matriz Global//////////////
    FMatrix k1g(3*nextIndex,3*nextIndex);
    FMatrix k2g(3*nextIndex,3*nextIndex);
    FMatrix k3g(3*nextIndex,3*nextIndex);
    FMatrix mg(3*nextIndex,3*nextIndex);

    FMatrix *k1 = new FMatrix(18, 18);
    FMatrix *k2 = new FMatrix(18, 18);
    FMatrix *k3 = new FMatrix(18, 18);
    FMatrix *m = new FMatrix(18, 18);
    double dj = 0;
    for(int NE = 0; NE < NET; NE ++)
    {
        QList<int> elemI = triangles.at(NE);
        int cuadraticElemI[6];
        cuadraticElemI[0] = elemI[0];
        cuadraticElemI[1] = elemI[1];
        cuadraticElemI[2] = elemI[2];

        int index3 = Contains(edges,cuadraticElemI[0],cuadraticElemI[1]);
        cuadraticElemI[3] = edges[index3][2];
        int index4 = Contains(edges,cuadraticElemI[1],cuadraticElemI[2]);
        cuadraticElemI[4] = edges[index4][2];
        int index5 = Contains(edges,cuadraticElemI[2],cuadraticElemI[0]);
        cuadraticElemI[5] = edges[index5][2];

        double coord[12];
        coord[0] = x[elemI.at(0)];
        coord[1] = y[elemI.at(0)];
        coord[2] = x[elemI.at(1)];
        coord[3] = y[elemI.at(1)];
        coord[4] = x[elemI.at(2)];
        coord[5] = y[elemI.at(2)];
        coord[6] = allMidpoints[edges[index3][2]- nnt].x;
        coord[7] = allMidpoints[edges[index3][2]- nnt].y;
        coord[8] = allMidpoints[edges[index4][2]- nnt].x;
        coord[9] = allMidpoints[edges[index4][2]- nnt].y;
        coord[10] = allMidpoints[edges[index5][2]- nnt].x;
        coord[11] = allMidpoints[edges[index5][2]- nnt].y;

        /* Formulacion SUBPARAMETRICA */
        double toCheck = MatrixK1K2K3MT6N(coord, f, k1, k2, k3, m) ;
        if(NE == 0)
            dj = toCheck;
        if(Sign(dj)!= Sign(toCheck))
        {
            qDebug()<<"ERROR EN EL JACOBIANO";
            emit finished();
        }

        /*
        % Para triangulo de 6 Nodos
        % ntab contiene las posiciones en el vector de las incognitas y en las
        % matrices globales correspondientes a los nodos del triangulo NE
        */
        double ntab[]={3*cuadraticElemI[0],3*cuadraticElemI[0]+1, 3*cuadraticElemI[0]+2,
                       3*cuadraticElemI[1],3*cuadraticElemI[1]+1, 3*cuadraticElemI[1]+2,
                       3*cuadraticElemI[2],3*cuadraticElemI[2]+1, 3*cuadraticElemI[2]+2,
                       3*cuadraticElemI[3],3*cuadraticElemI[3]+1, 3*cuadraticElemI[3]+2,
                       3*cuadraticElemI[4],3*cuadraticElemI[4]+1, 3*cuadraticElemI[4]+2,
                       3*cuadraticElemI[5],3*cuadraticElemI[5]+1, 3*cuadraticElemI[5]+2};

        int DMK = 18;
        for(int i = 0; i < DMK; i++)
        {
            for(int j=0;j<DMK;j++)
            {
                if(ntab[i] > ntab[j]) continue;
                k1g.SetElem(ntab[i],ntab[j],k1g.GetElem(ntab[i],ntab[j])+k1->GetElem(i,j));
                k2g.SetElem(ntab[i],ntab[j],k2g.GetElem(ntab[i],ntab[j])+k2->GetElem(i,j));
                k3g.SetElem(ntab[i],ntab[j],k3g.GetElem(ntab[i],ntab[j])+k3->GetElem(i,j));
                mg.SetElem(ntab[i],ntab[j],mg.GetElem(ntab[i],ntab[j])+m->GetElem(i,j));
            }
        }

        if(NE < NET-1){k1->Clear(); k2->Clear(); k3->Clear(); m->Clear();}
    }
    delete(k1);
    delete(k2);
    delete(k3);
    delete(m);

    k1g.Symmetrize();
    k2g.AntiSymmetrize();
    k3g.Symmetrize();
    mg.Symmetrize();

    NXIh.clear();
    for(double i = step; i <= maxWavesNumber; i += step)
        NXIh.push_back(i);

    emit maxIteration(NXIh.size());
    for(int i = 0; i < NXIh.size(); i++)
    {
        eigenValsVects.append(QList<AvalAvect>());
        for(int j = 0; j < autovalNumber; j++)
            eigenValsVects[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
    }

    //Variables
    complex<double> *b;
    complex<double> *a;
    int n = k1g.rows;
    int nvc = min(max(2*autovalNumber,20),n);
    complex<double> *w = new complex<double>[autovalNumber+1];
    complex<double> *v = new complex<double>[n*nvc];

    int info = 0;
    b = mg.ToComplexArray();
   // mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");

    FComplexMatrix *k2giXIh = new FComplexMatrix(k2g.rows,k2g.cols);
    FComplexMatrix *suma = new FComplexMatrix(k2g.rows,k2g.cols);
    FComplexMatrix *A = new FComplexMatrix(k2g.rows,k2g.cols);
    FMatrix *k3gXIh2 = new FMatrix(k2g.rows,k2g.cols);
    for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
    {
        double XIh = NXIh[INXIh];
        //  A = K1G+i*XIh*K2G+XIh^2*K3G;
        complex<double> iXIh;
        iXIh.imag(XIh);

        k2g.multi_escalar(iXIh,k2giXIh);
        k3g.multi_escalar(pow(XIh,2),k3gXIh2);
        k3gXIh2->Add(k2giXIh,suma);
        k1g.Add(suma,A);
        a = A->values;
//        A->WriteRealComplexMatrix("C://Users/Portege/Desktop/AR.txt");
//        A->WriteImagComplexMatrix("C://Users/Portege/Desktop/AI.txt");

        eigensolverq_(b,a,&n,&autovalNumber,&nvc,w,v,&info);

        if(info != 0)
        {
           qDebug()<<"info ="<<info<<"chi ="<<XIh;
        }

        SortRealNorm(w, eigenValsVects, INXIh, v, autovalNumber, n, XIh);
        emit posChanged(INXIh+1);

        k2giXIh->Clear();
        suma->Clear();
        A->Clear();
        k3gXIh2->Clear();
    }
    delete a;
    delete b;
    delete w;
    delete v;
    delete A;
    delete suma;
    delete k2giXIh;
    delete k3gXIh2;

    QList<int> deleted;
    deleted = PostProcessResults();
    PostProcessEigenValVects(deleted,NXIh.size(),autovalNumber);
    autovalNumber = autovalNumber - deleted.size();
    vf = FMatrix(NXIh.size(),autovalNumber);
    fM = FMatrix(NXIh.size(),autovalNumber);

    for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
    {
        for(int index = 0; index < autovalNumber; index++)
        {
            double val = sqrt(eigenValsVects[INXIh][index].GetAutoVal());
            double aux1 = val/NXIh[INXIh];
            vf.SetElem(INXIh, index, aux1);
            double aux2 = val;
            fM.SetElem(INXIh, index, aux2);
        }
    }

    vg = Diff(NXIh.size(),autovalNumber, eigenValsVects).multi_escalar(vt/step);

    clock_t end = clock();
    double diffticks = end-start;
    secTime =(diffticks)/(CLOCKS_PER_SEC);
    emit finished();
}
//M y A FULL Simetrico
void WorkerThread::doSymmetricWork()
{
//    qDebug()<<"Simétrico 1/4";
//    clock_t start = clock();
//    int NNT = symmetricalX.size(); //cantidad de vertices
//    int NET = symmetricalTriangles.size(); //cantidad de triangulos

//    //////////Datos del material elastico//////////
//    double c44 = young/2/(1+poisson);
//    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
//    double c11 = c12+2*c44;
//    double vt = sqrt(c44/density);
//    double f[3] ;
//    f[0] = c11/c44;
//    f[1] = c12/c44;
//    f[2] = 1;

//    int max_val = 0;
//    QList<int> neighbors = GetNeighbors(symmetricalTriangles,NNT, max_val);

//    ////////////Matriz Global//////////////
//    CMatrix mgC(3*NNT,3*max_val+3, neighbors);
//    FMatrix k1g(3*NNT,3*NNT);
//    FMatrix k2g(3*NNT,3*NNT);
//    FMatrix k3g(3*NNT,3*NNT);

//    FMatrix *k1 = new FMatrix(9, 9);
//    FMatrix *k2 = new FMatrix(9, 9);
//    FMatrix *k3 = new FMatrix(9, 9);
//    FMatrix *m = new FMatrix(9, 9);
//    for(int NE = 0; NE < NET; NE++)
//    {
//        QList<int> elemI = symmetricalTriangles.at(NE);

//        int ntab [9] = {3*elemI.at(0)-2, 3*elemI.at(0)-1, 3*elemI.at(0),
//                        3*elemI.at(1)-2, 3*elemI.at(1)-1, 3*elemI.at(1),
//                        3*elemI.at(2)-2, 3*elemI.at(2)-1, 3*elemI.at(2)};
//        double coord[6] = {symmetricalX[elemI[0]-1], symmetricalY[elemI[0]-1],
//                           symmetricalX[elemI[1]-1], symmetricalY[elemI[1]-1],
//                           symmetricalX[elemI[2]-1], symmetricalY[elemI[2]-1]};

//        MatrixK1K2K3MT3N(coord, f, k1, k2, k3, m);

//        int DMK = 9;
//        for(int i = 0; i < DMK; i++)
//        {
//            for(int j=0;j<DMK;j++)
//            {
//                int row = ntab[i] -1;
//                int column = ntab[j] -1;
//                k1g.SetElem(row,column,k1g.GetElem(row,column)+k1->GetElem(i,j));
//                k2g.SetElem(row,column,k2g.GetElem(row,column)+k2->GetElem(i,j));
//                k3g.SetElem(row,column,k3g.GetElem(row,column)+k3->GetElem(i,j));
//               // mg.SetElem(row,column,mg.GetElem(row,column)+m->GetElem(i,j));
//            }
//        }
//        mgC.SetElems(*m, elemI[0], elemI[1],elemI[2]);
//        if(NE < NET - 1)
//        {
//            k1->Clear(); k2->Clear(); k3->Clear(); m->Clear();
//        }
//    }
//    delete(k1);
//    delete(k2);
//    delete(k3);
//    delete(m);

//    NXIh.clear();
//    for(double i = step; i <= maxWavesNumber; i +=step)
//        NXIh.push_back(i);
//    emit maxIteration(3*NXIh.size());

//    longitudinalEigens.clear();
//    torsionalEigens.clear();
//    flexionalEigens.clear();
//    for(int i = 0; i < NXIh.size(); i++)
//    {
//        longitudinalEigens.append(QList<AvalAvect>());
//        torsionalEigens.append(QList<AvalAvect>());
//        flexionalEigens.append(QList<AvalAvect>());
//        for(int j = 0; j < autovalNumber; j++)
//        {
//            longitudinalEigens[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
//            torsionalEigens[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
//            flexionalEigens[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
//        }
//    }
//    int counter = 1;
//    for(int i = 0; i < 3; i++)
//    {
//        QList<int> indicesToDelete;

//        if(i == 0)//longitudinal
//        {
//            for(int i = 0 ; i < xZeroBoundaryIndices.count(); i++)
//                indicesToDelete.push_back(3*(xZeroBoundaryIndices[i]+1)-2);
//            for(int i = 0 ; i < yZeroBoundaryIndices.count(); i++)
//                 indicesToDelete.push_back(3*(yZeroBoundaryIndices[i]+1)-3);
//        }
//        else if(i == 1)//torsional
//        {
//            for(int i = 0 ; i < xZeroBoundaryIndices.count(); i++)
//            {
//                indicesToDelete.push_back(3*(xZeroBoundaryIndices[i]+1)-3);
//                indicesToDelete.push_back(3*(xZeroBoundaryIndices[i]+1)-1);
//            }
//            for(int i = 0 ; i < yZeroBoundaryIndices.count(); i++)
//            {
//                 indicesToDelete.push_back(3*(yZeroBoundaryIndices[i]+1)-2);
//                 indicesToDelete.push_back(3*(yZeroBoundaryIndices[i]+1)-1);
//            }
//        }
//        else //flexional
//        {
//            for(int i = 0 ; i < xZeroBoundaryIndices.count(); i++)
//            {
//                indicesToDelete.push_back(3*(xZeroBoundaryIndices[i]+1)-3);
//                indicesToDelete.push_back(3*(xZeroBoundaryIndices[i]+1)-1);
//            }
//            for(int i = 0 ; i < yZeroBoundaryIndices.count(); i++)
//                 indicesToDelete.push_back(3*(yZeroBoundaryIndices[i]+1)-3);
//        }

//        int dimensionX = mgC.rows;// - indicesToDelete.count();
//        int dimensionY = mgC.columns;// - indicesToDelete.count();
//        FMatrix mgSymmetric(dimensionX,dimensionY);
//        FMatrix k1gSymmetric(dimensionX,dimensionY);
//        FMatrix k2gSymmetric(dimensionX,dimensionY);
//        FMatrix k3gSymmetric(dimensionX,dimensionY);
//        mgC.TrimMatrix(indicesToDelete, &mgSymmetric);
//        k1g.TrimMatrix(indicesToDelete, &k1gSymmetric);
//        k2g.TrimMatrix(indicesToDelete, &k2gSymmetric);
//        k3g.TrimMatrix(indicesToDelete, &k3gSymmetric);
//        int nn = mg.rows;

//        //Variables
//        complex<double> *b;
//        complex<double> *a;
//        int n = mgSymmetric.rows;
//        int nvc = min(max(2*autovalNumber,20),n);
//        complex<double> *w = new complex<double>[autovalNumber];
//        complex<double> *v = new complex<double>[n*nvc];
//        int info = 0;
//        b = mgSymmetric.ToComplexArray();

//        FComplexMatrix *k2giXIh = new FComplexMatrix(k2gSymmetric.rows,k2gSymmetric.cols);
//        FComplexMatrix *suma = new FComplexMatrix(k2gSymmetric.rows,k2gSymmetric.cols);
//        FComplexMatrix *A = new FComplexMatrix(k2gSymmetric.rows,k2gSymmetric.cols);
//        FMatrix *k3gXIh2 = new FMatrix(k2gSymmetric.rows,k2gSymmetric.cols);
//        QList<int>symmetryArray;
//        if(i == 0)//longitudinal
//        {
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//        }
//        else if(i == 1)//torsional
//        {
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//        }
//        else //flexional
//        {
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//            symmetryArray.push_back(1);
//            symmetryArray.push_back(-1);
//        }

//        for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
//        {
//            double XIh = NXIh[INXIh];
//            //  A = K1G+i*XIh*K2G+XIh^2*K3G;
//            complex<double> iXIh;
//            iXIh.imag(XIh);

//            k2gSymmetric.multi_escalar(iXIh,k2giXIh);
//            k3gSymmetric.multi_escalar(pow(XIh,2),k3gXIh2);
//            k3gXIh2->Add(k2giXIh,suma);
//            k1gSymmetric.Add(suma,A);
//            a = A->values;

//            eigensolvermxm_(b,a,&n,&cCols,&autovalNumber,&nvc,w,v,neighFortran,&info&numIter,&numOpS,&numOpBX,&numReOR);
//            if(info != 0)
//            {
//                emit finished();
//            }
/////////////////////////AQUI
//            if(i==0)
//                SortRealNorm(w, longitudinalEigens, INXIh, v, autovalNumber, nn, XIh, indicesToDelete, symmetryArray);
//            else if(i==1)
//                SortRealNorm(w, torsionalEigens, INXIh, v, autovalNumber, nn, XIh, indicesToDelete, symmetryArray);
//            else
//                SortRealNorm(w, flexionalEigens, INXIh, v, autovalNumber, nn, XIh, indicesToDelete, symmetryArray);

//            emit posChanged(counter++);
//        }
//        if(i==0)
//        {
//            lVf = FMatrix(NXIh.size(),autovalNumber);
//            lF = FMatrix(NXIh.size(),autovalNumber);
//            for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
//            {
//                for(int index = 0; index < autovalNumber; index++)
//                {
//                    double val = sqrt(longitudinalEigens[INXIh][index].GetAutoVal());
//                    double aux1 = val/NXIh[INXIh];
//                    lVf.SetElem(INXIh, index, aux1);
//                    double aux2 = val;
//                    lF.SetElem(INXIh, index, aux2);
//                }
//            }

//            lVg = Diff(NXIh.size(),autovalNumber,longitudinalEigens).multi_escalar(vt/step);
//        }
//        else if(i==1)
//        {
//            tVf = FMatrix(NXIh.size(),autovalNumber);
//            tF = FMatrix(NXIh.size(),autovalNumber);
//            for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
//            {
//                for(int index = 0; index < autovalNumber; index++)
//                {
//                    double val = sqrt(torsionalEigens[INXIh][index].GetAutoVal());
//                    double aux1 = val/NXIh[INXIh];
//                    tVf.SetElem(INXIh, index, aux1);
//                    double aux2 = val;
//                    tF.SetElem(INXIh, index, aux2);
//                }
//            }

//            tVg = Diff(NXIh.size(),autovalNumber,torsionalEigens).multi_escalar(vt/step);
//        }
//        else
//        {
//            fVf = FMatrix(NXIh.size(),autovalNumber);
//            fF = FMatrix(NXIh.size(),autovalNumber);
//            for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
//            {
//                for(int index = 0; index < autovalNumber; index++)
//                {
//                    double val = sqrt(flexionalEigens[INXIh][index].GetAutoVal());
//                    double aux1 = val/NXIh[INXIh];
//                    fVf.SetElem(INXIh, index, aux1);
//                    double aux2 = val;
//                    fF.SetElem(INXIh, index, aux2);
//                }
//            }

//            fVg = Diff(NXIh.size(),autovalNumber,flexionalEigens).multi_escalar(vt/step);
//        }

//        delete(k2giXIh);
//        delete(suma);
//        delete(A);
//        delete(k3gXIh2);
//        delete a;
//        delete w;
//        delete v;
//    }
    emit finished();
}
//M compacta, A FULL Lineal
void WorkerThread::doLinealCompactedWork()
{
    qDebug()<<"Lineal";
    clock_t start = clock();
    int NNT = x.size(); //cantidad de vertices
    int NET = triangles.size(); //cantidad de triangulos

    //////////Datos del material elastico//////////
    double c44 = young/2/(1+poisson);
    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
    double c11 = c12+2*c44;
    double vt = sqrt(c44/density);
    double f[3] ;
    f[0] = c11/c44;
    f[1] = c12/c44;
    f[2] = 1;

    int max_val = 0;
    QList<int> neighbors = GetNeighbors(triangles,NNT, max_val);

//    clock_t startMatrix = clock();
    ////////////Matriz Global//////////////
    CMatrix mgC(3*NNT,3*max_val+3, neighbors);
//    CMatrix k1gC(3*NNT,3*max_val+3, neighbors);
//    CMatrix k2gC(3*NNT,3*max_val+3, neighbors);
//    CMatrix k3gC(3*NNT,3*max_val+3, neighbors);
    FMatrix k1g(3*NNT,3*NNT);
    FMatrix k2g(3*NNT,3*NNT);
    FMatrix k3g(3*NNT,3*NNT);
//    FMatrix mg(3*NNT,3*NNT);

    FMatrix *k1 = new FMatrix(9,9);
    FMatrix *k2 = new FMatrix(9,9);
    FMatrix *k3 = new FMatrix(9,9);
    FMatrix *m = new FMatrix(9,9);
    double dj = 0;
    for(int NE = 0; NE < NET; NE++)
    {
        QList<int> elemI = triangles.at(NE);

        int ntab [9] = {3*(elemI.at(0)+1)-2, 3*(elemI.at(0)+1)-1, 3*(elemI.at(0)+1),
                        3*(elemI.at(1)+1)-2, 3*(elemI.at(1)+1)-1, 3*(elemI.at(1)+1),
                        3*(elemI.at(2)+1)-2, 3*(elemI.at(2)+1)-1, 3*(elemI.at(2)+1)};
        double coord[6] = {x[elemI.at(0)], y[elemI.at(0)],
                           x[elemI.at(1)], y[elemI.at(1)],
                           x[elemI.at(2)], y[elemI.at(2)]};

        double toCheck = MatrixK1K2K3MT3N(coord, f, k1, k2, k3, m);
        if(NE == 0)
            dj = toCheck;
        if(Sign(dj)!= Sign(toCheck))
        {
            qDebug()<<"ERROR EN EL JACOBIANO";
            emit finished();
        }

        int DMK = 9;
        for(int i = 0; i < DMK; i++)
        {
            for(int j = 0; j < DMK; j++)
            {
                int row = ntab[i] -1;
                int column = ntab[j] -1;

                // Vamos a almacenar solo los valores por encima (y sobre) la diagonal.
                // Al final se llama al método Symmetrize de FMatrix para hacer la matriz simétrica.
                // Este método copia la mitad superior de la matriz en la inferior.
                if(row > column) continue;
                k1g.SetElem(row,column,k1g.GetElem(row,column)+k1->GetElem(i,j));
                k2g.SetElem(row,column,k2g.GetElem(row,column)+k2->GetElem(i,j));
                k3g.SetElem(row,column,k3g.GetElem(row,column)+k3->GetElem(i,j));
//                mg.SetElem(row,column,mg.GetElem(row,column)+m->GetElem(i,j));
            }
        }

        mgC.SetElems(*m, elemI[0], elemI[1],elemI[2]);
//        k1gC.SetElems(*k1, elemI[0], elemI[1],elemI[2]);
//        k2gC.SetElems(*k2, elemI[0], elemI[1],elemI[2]);
//        k3gC.SetElems(*k3, elemI[0], elemI[1],elemI[2]);

        if(NE < NET-1){k1->Clear(); k2->Clear(); k3->Clear(); m->Clear();}
    }

    delete(k1);
    delete(k2);
    delete(k3);
    delete(m);

    k1g.Symmetrize();
    k2g.AntiSymmetrize();
    k3g.Symmetrize();
//    k1gC.Symmetrize();
//    k2gC.AntiSymmetrize();
//    k3gC.Symmetrize();
//    mg.Symmetrize();
    mgC.Symmetrize();
//    mgC.WriteMatrixAndNeighbors("C://Users/Portege/Desktop/mgC.txt", "C://Users/Portege/Desktop/neigh.txt");
//    mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");
//    FMatrix dmg = DecompactMatrix("C://Users/Portege/Desktop/mgC.txt", "C://Users/Portege/Desktop/neigh.txt", "C://Users/Portege/Desktop/mgD.txt");
//    for(int i = 0; i < dmg.rows; i++)
//    {
//        for(int j = i+1; j < dmg.cols; j++)
//        {
//            if(dmg.GetElem(i,j) != dmg.GetElem(j,i) || Abs(dmg.GetElem(i,j) - mg.GetElem(i,j))>1e-6)
//            {
//                qDebug()<<"i ="<<i<<", j ="<<j ;
//                qDebug()<<"D(i,j)= "<< dmg.GetElem(i,j);
//                qDebug()<<"O(i,j)= "<< mg.GetElem(i,j);
//            }

//        }
//    }

//    clock_t endMatrix = clock();
//    double diffticksMatrix = (endMatrix-startMatrix);
//    double secMatrix = (diffticksMatrix)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo en construir las matrices globales" << secMatrix << " segundos.";
    //mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");

    NXIh.clear();
    for(double i = step; i <= maxWavesNumber; i+=step)
        NXIh.push_back(i);

    emit maxIteration(NXIh.size());
    for(int i = 0; i < NXIh.size(); i++)
    {
        eigenValsVects.append(QList<AvalAvect>());
        for(int j=0;j<autovalNumber;j++)
            eigenValsVects[i].append(AvalAvect(-1, -1, 0, QList<complex<double> >()));
    }

    //Variables
    complex<double> *b;
    complex<double> *a;
    int n = mgC.rows;
    int nvc = min(max(2*autovalNumber, 20), n);
    complex<double> *w = new complex<double>[autovalNumber+1];
    complex<double> *v = new complex<double>[n*nvc];
    int cCols = mgC.columns;

    int info = 0;
    b = mgC.ToComplexArray();
    int *neighFortran = ConvertToFortranParam(neighbors, NNT, max_val);
//    mgC.WriteMatrixAndNeighbors("C://Users/Portege/Desktop/mCg.txt","C://Users/Portege/Desktop/neigh.txt");
//    mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");
//    WriteVector("C://Users/Portege/Desktop/neighbors.txt",neighFortran,NNT*max_val);

    FComplexMatrix *k2giXIh = new FComplexMatrix(k2g.rows,k2g.cols);
    FComplexMatrix *suma = new FComplexMatrix(k2g.rows,k2g.cols);
    FComplexMatrix *A = new FComplexMatrix(k2g.rows,k2g.cols);
    FMatrix *k3gXIh2 = new FMatrix(k2g.rows,k2g.cols);

//    CComplexMatrix *k2giXIhC = new CComplexMatrix(k2gC.rows,cCols,neighbors);
//    CComplexMatrix *sumaC = new CComplexMatrix(k2gC.rows,cCols,neighbors);
//    CComplexMatrix *AC = new CComplexMatrix(k2gC.rows,cCols,neighbors);
//    CMatrix *k3gXIh2C = new CMatrix(k2gC.rows,cCols,neighbors);

//    complex<double> *afact = new complex<double>[n*n];
//    complex<double> *work = new complex<double>[2*n];
//    complex<double> *workC = new complex<double>[2*n];
//    double *rwork = new double[2*n];
//    int *ipiv = new int[n];
//    char *norm = new char[1];
//    norm[0] = 'I';
//    double rcond = 0;

//    clock_t endConstructionProcess = clock();
//    double diffticksConstructionProcess = (endConstructionProcess - start);
//    double secConstructionProcess = (diffticksConstructionProcess)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo que demora en preparar todo para comenzar a calcular " << secConstructionProcess << " segundos.";

    int numIter = 0;
    int numOpS = 0;
    int numOpBX = 0;
    int numReOR = 0;
    for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
    {
      //  clock_t startOnePGA = clock();
        double XIh = NXIh[INXIh];
        //  A = K1G+i*XIh*K2G+XIh^2*K3G;
        complex<double> iXIh;
        iXIh.imag(XIh);

        k2g.multi_escalar(iXIh,k2giXIh);
        k3g.multi_escalar(pow(XIh,2),k3gXIh2);
        k3gXIh2->Add(k2giXIh,suma);
        k1g.Add(suma,A);
        a = A->values;

//        for(int i = 0; i < n*n; i++)
//        {
//            double real = a[i].real();
//            double imag = a[i].imag();
//            afact[i].real(real);
//            afact[i].imag(imag);
//        }

//        double anorm = zlange_(norm, &n, &n, afact, &n, work);
//        zgetrf_(&n, &n, afact, &n, ipiv, &info);
//        if(info != 0)
//        {
//            int asas = 0;
//        }
//        zgecon_(norm, &n, afact, &n, &anorm, &rcond, workC, rwork, &info);
//        if(info != 0)
//        {
//            int adaw = 0;
//        }
//        if(1/rcond > 1e2)
//            illConditionatedRowsChis.append(INXIh);

//        k2gC.multi_escalar(iXIh,k2giXIhC);
//        k3gC.multi_escalar(pow(XIh,2),k3gXIh2C);
//        k3gXIh2C->Add(k2giXIhC,sumaC);
//        k1gC.Add(sumaC,AC);

//        A->WriteRealComplexMatrix("C://Users/Portege/Desktop/AR.txt");
//        A->WriteImagComplexMatrix("C://Users/Portege/Desktop/AI.txt");
//        AC->WriteRealComplexMatrix("C://Users/Portege/Desktop/ACR.txt");
//        AC->WriteImagComplexMatrix("C://Users/Portege/Desktop/ACI.txt");
//        FMatrix adR = DecompactMatrix("C://Users/Portege/Desktop/ACR.txt","C://Users/Portege/Desktop/neigh.txt","C://Users/Portege/Desktop/ARD.txt");
//        FMatrix adI = DecompactMatrix("C://Users/Portege/Desktop/ACI.txt","C://Users/Portege/Desktop/neigh.txt","C://Users/Portege/Desktop/AID.txt");

//        for(int i = 0; i < adR.rows; i++)
//        {
//            for(int j = 0; j < adR.cols; j++)
//            {
//                if(Abs(adR.GetElem(i,j) - A->GetElem(i,j).real()) > 1e-6 || Abs(adI.GetElem(i,j) - A->GetElem(i,j).imag()) > 1e-6)
//                {
//                    qDebug()<<"D != O (i,j)=("<<i<<","<<j<<") D = ("<<adR.GetElem(i,j)<<","<<adI.GetElem(i,j)<<")  O = ("<<A->GetElem(i,j).real()<<","<<A->GetElem(i,j).imag()<<")";
//                }
//            }
//        }

        eigensolvermxm_(b,a,&n,&cCols,&autovalNumber,&nvc,w,v,neighFortran,&info,&numIter,&numOpS,&numOpBX,&numReOR);

        if(info != 0)
        {
            qDebug()<<"info ="<<info<<"chi ="<<XIh;
//            qDebug()<<"     Numero Iteraciones:"<<numIter;
//            qDebug()<<"     Numero Operaciones Sistemas:"<<numOpS;
//            qDebug()<<"     Numero Operaciones MX:"<<numOpBX;
//            qDebug()<<"     Numero de pasos de Reortogonalizaciones:"<<numReOR;
        }

        SortRealNorm(w, eigenValsVects, INXIh, v, autovalNumber, n, XIh);
        emit posChanged(INXIh+1);

//        for(int i = 0; i < autovalNumber+1; i++)
//        {
//            w[i].real(0);
//            w[i].imag(0);
//        }
//        for(int i = 0; i < nvc*n; i++)
//        {
//            v[i].real(0);
//            v[i].imag(0);
//        }
        k2giXIh->Clear();
        suma->Clear();
        A->Clear();
        k3gXIh2->Clear();
//        clock_t endOnePGA = clock();
//        double diffticksOnePGA = (endOnePGA - startOnePGA);
//        double secOnePGA = (diffticksOnePGA)/(CLOCKS_PER_SEC);
//        qDebug()<< "Tiempo que demoró en correr el PGA" << secOnePGA << " segundos.";
//        int a = 0;
    }

    delete w;
    delete v;
    delete neighFortran;
    delete a;
    delete b;
    delete k2giXIh;
    delete suma;
    delete A;
    delete k3gXIh2;

    QList<int> deleted = PostProcessResults();
    PostProcessEigenValVects(deleted,NXIh.size(),autovalNumber);
    autovalNumber = autovalNumber - deleted.size();
    vf = FMatrix(NXIh.size(),autovalNumber);
    fM = FMatrix(NXIh.size(),autovalNumber);

    for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
    {
        for(int index = 0; index < autovalNumber; index++)
        {
            double val = sqrt(eigenValsVects[INXIh][index].GetAutoVal());
            double aux1 = val/NXIh[INXIh];
            vf.SetElem(INXIh, index, aux1);     // Velocidad de fase/velocidad transversal (m/s)
            double aux2 = val;
            fM.SetElem(INXIh, index, aux2);     // Frecuencia*h (Lado de la Seccion, o diametro  (Hz)
        }
    }
    vg = Diff(NXIh.size(),autovalNumber,eigenValsVects).multi_escalar(vt/step);

    clock_t end = clock();
    double diffticks = end-start;
    secTime =(diffticks)/(CLOCKS_PER_SEC);
    emit finished();
}
//[allCompacted. Uses CG]
//M y A compactas Lineal
void WorkerThread::doLinealCGCompactedWork()
{
    qDebug()<<"Lineal (M y A compactas) GC";
    clock_t start = clock();
    int NNT = x.size(); //cantidad de vertices
    int NET = triangles.size(); //cantidad de triangulos

    //////////Datos del material elastico//////////
    double c44 = young/2/(1+poisson);
    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
    double c11 = c12+2*c44;
    double vt = sqrt(c44/density);
    double f[3];
    f[0] = c11/c44;
    f[1] = c12/c44;
    f[2] = 1;

    int max_val = 0;
    QList<int> neighbors = GetNeighbors(triangles, NNT, max_val);

    ////////////Matriz Global//////////////
    CMatrix mgC(3*NNT, 3*max_val+3, neighbors);
    CMatrix k1gC(3*NNT, 3*max_val+3, neighbors);
    CMatrix k2gC(3*NNT, 3*max_val+3, neighbors);
    CMatrix k3gC(3*NNT, 3*max_val+3, neighbors);

    FMatrix *k1 = new FMatrix(9, 9);
    FMatrix *k2 = new FMatrix(9, 9);
    FMatrix *k3 = new FMatrix(9, 9);
    FMatrix *m = new FMatrix(9, 9);
    double dj = 0;
    for(int NE = 0; NE < NET; NE++)
    {
        QList<int> elemI = triangles.at(NE);
        double coord[6] = {x[elemI.at(0)], y[elemI.at(0)],
                           x[elemI.at(1)], y[elemI.at(1)],
                           x[elemI.at(2)], y[elemI.at(2)]};

        double toCheck = MatrixK1K2K3MT3N(coord, f, k1, k2, k3, m) ;
        if(NE == 0)
            dj = toCheck;
        if(Sign(dj)!= Sign(toCheck))
        {
            qDebug()<<"ERROR EN EL JACOBIANO";
            emit finished();
        }

        k1gC.SetElems(*k1,elemI[0], elemI[1], elemI[2]);
        k2gC.SetElems(*k2,elemI[0], elemI[1], elemI[2]);
        k3gC.SetElems(*k3,elemI[0], elemI[1], elemI[2]);
        mgC.SetElems(*m,elemI[0], elemI[1], elemI[2]);

        if(NE < NET-1){k1->Clear(); k2->Clear(); k3->Clear(); m->Clear();}
    }
    delete(k1);
    delete(k2);
    delete(k3);
    delete(m);

    k1gC.Symmetrize();
    k2gC.AntiSymmetrize();
    k3gC.Symmetrize();
    mgC.Symmetrize();

    NXIh.clear();
    for(double i = step; i <= maxWavesNumber; i+=step)
        NXIh.push_back(i);

    emit maxIteration(NXIh.size());
    for(int i = 0; i < NXIh.size(); i++)
    {
        eigenValsVects.append(QList<AvalAvect>());
        for(int j = 0; j < autovalNumber; j++)
            eigenValsVects[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
    }

    //Variables
    complex<double> *b;
    complex<double> *a;
    int n = mgC.rows;
    int nvc = min(max(2*autovalNumber,20),n);
    complex<double> *w = new complex<double>[autovalNumber+1];
    complex<double> *v = new complex<double>[n*nvc];
    int cCols = mgC.columns;

    int info = 0;
    b = mgC.ToComplexArray();
    int *neighFortran = ConvertToFortranParam(neighbors, NNT, max_val);
//    mgC.WriteMatrixAndNeighbors("C://Users/Portege/Desktop/MC.txt", "C://Users/Portege/Desktop/neigh.txt");
//    WriteVector("C://Users/Portege/Desktop/neighbors.txt", neighFortran, NNT*max_val);

    CComplexMatrix *k2giXIh = new CComplexMatrix(k2gC.rows, cCols, neighbors);
    CComplexMatrix *suma = new CComplexMatrix(k2gC.rows, cCols, neighbors);
    CComplexMatrix *A = new CComplexMatrix(k2gC.rows, cCols, neighbors);
    CMatrix *k3gXIh2 = new CMatrix(k2gC.rows, cCols, neighbors);
    for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
    {
        double XIh = NXIh[INXIh];
        //  A = K1G+i*XIh*K2G+XIh^2*K3G;
        complex<double> iXIh;
        iXIh.imag(XIh);

        k2gC.multi_escalar(iXIh, k2giXIh);
        k3gC.multi_escalar(pow(XIh,2), k3gXIh2);
        k3gXIh2->Add(k2giXIh, suma);
        k1gC.Add(suma,A);
        a = A->ToArray();
//        A->WriteRealComplexMatrix("C://Users/Portege/Desktop/ACR.txt");
//        A->WriteImagComplexMatrix("C://Users/Portege/Desktop/ACI.txt");

        eigensolvermxmaxb_(b,a,&n,&cCols,&autovalNumber,&nvc,w,v,neighFortran,&info);
        if(info != 0)
        {
            qDebug()<<"info ="<<info<<"chi ="<<XIh;
//            A->WriteRealComplexMatrix("C://Users/Portege/Desktop/ACR.txt");
//            A->WriteImagComplexMatrix("C://Users/Portege/Desktop/ACI.txt");
        }

        SortRealNorm(w, eigenValsVects, INXIh, v, autovalNumber, n, XIh);
        emit posChanged(INXIh+1);

        k2giXIh->Clear();
        suma->Clear();
        A->Clear();
        k3gXIh2->Clear();
    }

    delete w;
    delete v;
    delete neighFortran;
    delete a;
    delete b;
    delete k2giXIh;
    delete suma;
    delete A;
    delete k3gXIh2;

    QList<int> deleted = PostProcessResults();
    PostProcessEigenValVects(deleted,NXIh.size(),autovalNumber);
    autovalNumber = autovalNumber - deleted.size();
    vf = FMatrix(NXIh.size(),autovalNumber);
    fM = FMatrix(NXIh.size(),autovalNumber);

    for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
    {
        for(int index = 0; index < autovalNumber; index++)
        {
            double val = sqrt(eigenValsVects[INXIh][index].GetAutoVal());
            double aux1 = val/NXIh[INXIh];
            vf.SetElem(INXIh, index, aux1); // Velocidad de fase/velocidad transversal (m/s)
            double aux2 = val;
            fM.SetElem(INXIh, index, aux2); // Frecuencia*h (Lado de la Seccion, o diametro  (Hz)
        }
    }
    vg = Diff(NXIh.size(),autovalNumber,eigenValsVects).multi_escalar(vt/step);

    clock_t end = clock();
    double diffticks = end-start;
    secTime =(diffticks)/(CLOCKS_PER_SEC);
    emit finished();
}
//M compacta y A FULL Cuadratico
void WorkerThread::doQuadraticCompactedWork()
{
    qDebug()<< "Cuadrático (MC, A full)";
    clock_t start = clock();
    int nnt = x.size(); //cantidad de vertices
    int NET = triangles.size(); //cantidad de triangulos

    //////////Datos del material elastico//////////
    double c44 = young/2/(1+poisson);
    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
    double c11 = c12+2*c44;
    double vt = sqrt(c44/density);
    double f[3];
    f[0] = c11/c44;
    f[1] = c12/c44;
    f[2] = 1;

//    clock_t startQNodes = clock();
    int nextIndex = nnt;
    QList<CuadraticVertexData> allMidpoints;
    QList<QList<int> > edges;
    QList<QList<int> > quadraticTriangles;
    for(int NE = 0; NE < NET; NE++)
    {
        QList<int> elemI = triangles.at(NE);
        QList<int> cuadraticElemI;
        cuadraticElemI.append(elemI.at(0));
        cuadraticElemI.append(elemI.at(1));
        cuadraticElemI.append(elemI.at(2));
        cuadraticElemI.append(-1);
        cuadraticElemI.append(-1);
        cuadraticElemI.append(-1);

        double coord[12];
        coord[0] = x[elemI.at(0)];
        coord[1] = y[elemI.at(0)];
        coord[2] = x[elemI.at(1)];
        coord[3] = y[elemI.at(1)];
        coord[4] = x[elemI.at(2)];
        coord[5] = y[elemI.at(2)];

        int index = Contains(edges,elemI.at(0),elemI.at(1));
        if(index != -1)
        {
            //Si ya existia esa arista, se guarda en la posicion 3 el indice que referencia al punto medio de ella
            cuadraticElemI.replace(3, edges[index][2]);
        }
        else
        {
            cuadraticElemI.replace(3, nextIndex++);

            coord[6] = (coord[0]+coord[2])/2;
            coord[7] = (coord[1]+coord[3])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI.at(3),coord[6], coord[7]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(0));
            edgeI.push_back(elemI.at(1));
            edgeI.push_back(cuadraticElemI[3]);
            edges.push_back(edgeI);
        }
        index = Contains(edges,elemI.at(1),elemI.at(2));
        if(index != -1)
        {
            cuadraticElemI.replace(4, edges[index][2]);
        }
        else
        {
            cuadraticElemI.replace(4, nextIndex++);
            coord[8] = (coord[2]+coord[4])/2;
            coord[9] = (coord[3]+coord[5])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI.at(4), coord[8], coord[9]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(1));
            edgeI.push_back(elemI.at(2));
            edgeI.push_back(cuadraticElemI.at(4));
            edges.push_back(edgeI);
        }
        index = Contains(edges,elemI.at(2),elemI.at(0));
        if(index != -1)
        {
            cuadraticElemI.replace(5, edges[index][2]);
        }
        else
        {
            cuadraticElemI.replace(5, nextIndex++);
            coord[10] = (coord[4]+coord[0])/2;
            coord[11] = (coord[5]+coord[1])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI.at(5), coord[10], coord[11]));

            QList<int> edgeI;
            edgeI.push_back(elemI[2]);
            edgeI.push_back(elemI[0]);
            edgeI.push_back(cuadraticElemI.at(5));
            edges.push_back(edgeI);
        }
        quadraticTriangles.push_back(cuadraticElemI);
    }

    int max_val  = 0;
    QList<int> neigh = GetQuadraticNeighbors(quadraticTriangles, nextIndex, max_val);

//    clock_t endQNodes = clock();
//    double diffticksQNodes = (endQNodes - startQNodes);
//    double secQNodes = (diffticksQNodes)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo en construir los nodos cuadraticos" << secQNodes << " segundos.";

//    clock_t startQMatrix = clock();
    ////////////Matriz Global//////////////
    FMatrix k1g(3*nextIndex,3*nextIndex);
    FMatrix k2g(3*nextIndex,3*nextIndex);
    FMatrix k3g(3*nextIndex,3*nextIndex);
//    FMatrix mg(3*nextIndex, 3*nextIndex);

    CMatrix mCg(3*nextIndex,3*(max_val+1), neigh);
//    CMatrix k1Cg(3*nextIndex,3*(max_val+1), neigh);
//    CMatrix k2Cg(3*nextIndex,3*(max_val+1), neigh);
//    CMatrix k3Cg(3*nextIndex,3*(max_val+1), neigh);


    FMatrix *k1 = new FMatrix(18, 18);
    FMatrix *k2 = new FMatrix(18, 18);
    FMatrix *k3 = new FMatrix(18, 18);
    FMatrix *m = new FMatrix(18, 18);
    double dj = 0;
    for(int NE = 0; NE < NET; NE ++)
    {
        QList<int> elemI = triangles.at(NE);
        int cuadraticElemI[6];
        cuadraticElemI[0] = elemI.at(0);
        cuadraticElemI[1] = elemI.at(1);
        cuadraticElemI[2] = elemI.at(2);

        int index3 = Contains(edges,cuadraticElemI[0],cuadraticElemI[1]);
        cuadraticElemI[3] = edges[index3][2];
        int index4 = Contains(edges,cuadraticElemI[1],cuadraticElemI[2]);
        cuadraticElemI[4] = edges[index4][2];
        int index5 = Contains(edges,cuadraticElemI[2],cuadraticElemI[0]);
        cuadraticElemI[5] = edges[index5][2];

        double coord[12];
        coord[0] = x[elemI.at(0)];
        coord[1] = y[elemI.at(0)];
        coord[2] = x[elemI.at(1)];
        coord[3] = y[elemI.at(1)];
        coord[4] = x[elemI.at(2)];
        coord[5] = y[elemI.at(2)];
        coord[6] = allMidpoints[edges[index3][2]- nnt].x;
        coord[7] = allMidpoints[edges[index3][2]- nnt].y;
        coord[8] = allMidpoints[edges[index4][2]- nnt].x;
        coord[9] = allMidpoints[edges[index4][2]- nnt].y;
        coord[10] = allMidpoints[edges[index5][2]- nnt].x;
        coord[11] = allMidpoints[edges[index5][2]- nnt].y;

        /* Formulacion SUBPARAMETRICA */
        double toCheck = MatrixK1K2K3MT6N(coord, f, k1, k2, k3, m);
        if(NE == 0)
            dj = toCheck;
        if(Sign(dj)!= Sign(toCheck))
        {
            qDebug()<<"ERROR EN EL JACOBIANO";
            emit finished();
        }

        /*
        % Para triangulo de 6 Nodos
        % ntab contiene las posiciones en el vector de las incognitas y en las
        % matrices globales correspondientes a los nodos del triangulo NE
        */
        double ntab[]={3*cuadraticElemI[0],3*cuadraticElemI[0]+1, 3*cuadraticElemI[0]+2,
                       3*cuadraticElemI[1],3*cuadraticElemI[1]+1, 3*cuadraticElemI[1]+2,
                       3*cuadraticElemI[2],3*cuadraticElemI[2]+1, 3*cuadraticElemI[2]+2,
                       3*cuadraticElemI[3],3*cuadraticElemI[3]+1, 3*cuadraticElemI[3]+2,
                       3*cuadraticElemI[4],3*cuadraticElemI[4]+1, 3*cuadraticElemI[4]+2,
                       3*cuadraticElemI[5],3*cuadraticElemI[5]+1, 3*cuadraticElemI[5]+2};

        int DMK = 18;
        for(int i = 0; i < DMK; i++)
        {
            for(int j=0;j<DMK;j++)
            {
                if(ntab[i] > ntab[j]) continue;
                k1g.SetElem(ntab[i],ntab[j],k1g.GetElem(ntab[i],ntab[j])+k1->GetElem(i,j));
                k2g.SetElem(ntab[i],ntab[j],k2g.GetElem(ntab[i],ntab[j])+k2->GetElem(i,j));
                k3g.SetElem(ntab[i],ntab[j],k3g.GetElem(ntab[i],ntab[j])+k3->GetElem(i,j));
//                mg.SetElem(ntab[i],ntab[j],mg.GetElem(ntab[i],ntab[j])+m->GetElem(i,j));
            }
        }
        mCg.SetElems(*m,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);
//        k1Cg.SetElems(*k1,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
//                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);
//        k2Cg.SetElems(*k2,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
//                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);
//        k3Cg.SetElems(*k3,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
//                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);

        if(NE < NET - 1){k1->Clear(); k2->Clear(); k3->Clear(); m->Clear();}
    }
    delete(k1);
    delete(k2);
    delete(k3);
    delete(m);

//    clock_t endQMatrix = clock();
//    double diffticksQMatrix = (endQMatrix - startQMatrix);
//    double secQMatrix  = (diffticksQMatrix)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo en construir las matrices globales cuadraticas" << secQMatrix << " segundos.";

    k1g.Symmetrize();
    k2g.AntiSymmetrize();
    k3g.Symmetrize();
//    k1Cg.Symmetrize();
//    k2Cg.AntiSymmetrize();
//    k3Cg.Symmetrize();
    mCg.Symmetrize();

    NXIh.clear();
    for(double i = step; i <= maxWavesNumber; i+=step)
        NXIh.push_back(i);
    emit maxIteration(NXIh.size());

    for(int i = 0; i < NXIh.size(); i++)
    {
        eigenValsVects.append(QList<AvalAvect>());
        for(int j=0;j<autovalNumber;j++)
            eigenValsVects[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
    }

    //Variables
    complex<double> *b;
    complex<double> *a;
    int n = mCg.rows;
    int nvc = min(max(2*autovalNumber,20),n);
    complex<double> *w = new complex<double>[autovalNumber+1];
    complex<double> *v = new complex<double>[n*nvc];

    int info = 0;
    b = mCg.ToComplexArray();
    int cCols = mCg.columns;
    int *neighFortran = ConvertToFortranParam(neigh,nextIndex,max_val);
//    mCg.WriteMatrixAndNeighbors("C://Users/Portege/Desktop/mgC.txt", "C://Users/Portege/Desktop/neigh.txt");
//    mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");
//    WriteVector("C://Users/Portege/Desktop/neighbors.txt",neighFortran,nextIndex*max_val);
//    FMatrix dmg = DecompactMatrix("C://Users/Portege/Desktop/mgC.txt", "C://Users/Portege/Desktop/neigh.txt", "C://Users/Portege/Desktop/mgD.txt");
//    for(int i = 0; i < dmg.rows; i++)
//    {
//        for(int j = 0; j < dmg.cols; j++)
//        {
//            if(Abs(dmg.GetElem(i,j) - mg.GetElem(i,j)) > 1e-6)
//                qDebug()<<"i="<<i<<"j="<<j << "D:"<< dmg.GetElem(i,j) << "O:" << mg.GetElem(i,j);
//        }
//    }
    FComplexMatrix *k2giXIh = new FComplexMatrix(k2g.rows,k2g.cols);
    FComplexMatrix *suma = new FComplexMatrix(k2g.rows,k2g.cols);
    FComplexMatrix *A = new FComplexMatrix(k2g.rows,k2g.cols);
    FMatrix *k3gXIh2 = new FMatrix(k2g.rows,k2g.cols);
//    CComplexMatrix *k2giXIhC = new CComplexMatrix(k2Cg.rows,cCols, neigh);
//    CComplexMatrix *sumaC = new CComplexMatrix(k2Cg.rows,cCols, neigh);
//    CComplexMatrix *AC = new CComplexMatrix(k2Cg.rows,cCols, neigh);
//    CMatrix *k3gXIh2C = new CMatrix(k2Cg.rows,cCols, neigh);

//    clock_t endQConstructionProcess = clock();
//    double diffticksQConstructionProcess = (endQConstructionProcess - start);
//    double secQConstrictionProcess = (diffticksQConstructionProcess)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo que demora en preparar todo para comenzar a calcular " << secQConstrictionProcess << " segundos.";
    int numIter = 0;
    int numOpS = 0;
    int numOpBX = 0;
    int numReOR = 0;
    for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
    {
//        clock_t startOnePGA = clock();
        double XIh = NXIh[INXIh];
        //  A = K1G+i*XIh*K2G+XIh^2*K3G;
        complex<double> iXIh;
        iXIh.imag(XIh);

        k2g.multi_escalar(iXIh,k2giXIh);
        k3g.multi_escalar(pow(XIh,2),k3gXIh2);
        k3gXIh2->Add(k2giXIh,suma);
        k1g.Add(suma,A);
        a = A->values;
//        A->WriteRealComplexMatrix("C://Users/Portege/Desktop/AR.txt");
//        A->WriteImagComplexMatrix("C://Users/Portege/Desktop/AI.txt");

//        k2Cg.multi_escalar(iXIh,k2giXIhC);
//        k3Cg.multi_escalar(pow(XIh,2),k3gXIh2C);
//        k3gXIh2C->Add(k2giXIhC,sumaC);
//        k1Cg.Add(sumaC,AC);
//        AC->WriteRealComplexMatrix("C://Users/Portege/Desktop/ACR.txt");
//        AC->WriteImagComplexMatrix("C://Users/Portege/Desktop/ACI.txt");
//        FComplexMatrix da = DecompactCMatrix("C://Users/Portege/Desktop/ACR.txt", "C://Users/Portege/Desktop/ACI.txt", "C://Users/Portege/Desktop/neigh.txt", "C://Users/Portege/Desktop/ADR.txt", "C://Users/Portege/Desktop/ADI.txt");

//        for(int i = 0; i < da.rows; i++)
//        {
//            for(int j = 0; j < da.cols; j++)
//            {
//                if(Abs(da.GetElem(i,j).real() - A->GetElem(i,j).real()) > 1e-6 || Abs(da.GetElem(i,j).imag() - A->GetElem(i,j).imag()) > 1e-6)
//                {
//                    qDebug()<<"REAL i="<<i<<"j="<<j << "D:"<< da.GetElem(i,j).real() << "O:" << A->GetElem(i,j).real();
//                    qDebug()<<"IMAG i="<<i<<"j="<<j << "D:"<< da.GetElem(i,j).imag() << "O:" << A->GetElem(i,j).imag();
//                }
//            }
//        }

//        a = AC->ToArray();
        eigensolvermxm_(b,a,&n,&cCols,&autovalNumber,&nvc,w,v,neighFortran,&info,&numIter,&numOpS,&numOpBX,&numReOR);

//        if(info != 0)
//        {
//            qDebug()<<"info ="<<info<<"chi ="<<XIh;
//            qDebug()<<"     Numero Iteraciones:"<<numIter;
//            qDebug()<<"     Numero Operaciones Sistemas:"<<numOpS;
//            qDebug()<<"     Numero Operaciones MX:"<<numOpBX;
//            qDebug()<<"     Numero de pasos de Reortogonalizaciones:"<<numReOR;
//        }

        SortRealNorm(w, eigenValsVects, INXIh, v, autovalNumber, n, XIh);
        emit posChanged(INXIh+1);
        //qDebug() << qSetRealNumberPrecision( 10 ) << "XIh = " << XIh;
        k2giXIh->Clear();
        suma->Clear();
        A->Clear();
        k3gXIh2->Clear();
//        k2giXIhC->Clear();
//        sumaC->Clear();
//        AC->Clear();
//        k3gXIh2C->Clear();
//        clock_t endOnePGA = clock();
//        double diffticksOnePGA = (endOnePGA - startOnePGA);
//        double secOnePGA = (diffticksOnePGA)/(CLOCKS_PER_SEC);
//        qDebug()<< "Tiempo que demoró en resolver el PGA " << secOnePGA << " segundos.";
//        int a = 0;

    }
    delete k2giXIh;
    delete suma;
    delete A;
    delete k3gXIh2;
//    delete(k2giXIhC);
//    delete(sumaC);
//    delete(AC);
//    delete(k3gXIh2C);
    delete a;
    delete b;
    delete w;
    delete v;
    delete neighFortran;

    QList<int> deleted = PostProcessResults();
    PostProcessEigenValVects(deleted,NXIh.size(),autovalNumber);
    autovalNumber = autovalNumber - deleted.size();
    vf = FMatrix(NXIh.size(),autovalNumber);
    fM = FMatrix(NXIh.size(),autovalNumber);

    for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
    {
        for(int index = 0; index < autovalNumber; index++)
        {
            double val = sqrt(eigenValsVects[INXIh][index].GetAutoVal());
            double aux1 = val/NXIh[INXIh];
            vf.SetElem(INXIh, index, aux1);
            double aux2 = val;
            fM.SetElem(INXIh, index, aux2);
        }
    }

    vg = Diff(NXIh.size(),autovalNumber, eigenValsVects).multi_escalar(vt/step);

    clock_t end = clock();
    double diffticks = end-start;
    secTime =(diffticks)/(CLOCKS_PER_SEC);
    emit finished();
}
//M y A compactas
void WorkerThread::doQuadraticCGCompactedWork()
{
    //FComplexMatrix da = DecompactCMatrix("C://Users/Portege/Desktop/ACR.txt", "C://Users/Portege/Desktop/ACI.txt", "C://Users/Portege/Desktop/neigh.txt", "C://Users/Portege/Desktop/ADR.txt", "C://Users/Portege/Desktop/ADI.txt");
    qDebug()<< "Cuadrático (M y A compactas)";
    clock_t start = clock();
    int nnt = x.size(); //cantidad de vertices
    int NET = triangles.size(); //cantidad de triangulos

    //////////Datos del material elastico//////////
    double c44 = young/2/(1+poisson);
    double c12 = young*poisson/(1+poisson)/(1-2*poisson);
    double c11 = c12+2*c44;
    double vt = sqrt(c44/density);
    double f[3];
    f[0] = c11/c44;
    f[1] = c12/c44;
    f[2] = 1;

//    clock_t startQNodes = clock();
    int nextIndex = nnt;
    QList<CuadraticVertexData> allMidpoints;
    QList<QList<int> > edges;
    QList<QList<int> > quadraticTriangles;
    for(int NE = 0; NE < NET; NE++)
    {
        QList<int> elemI = triangles.at(NE);
        QList<int> cuadraticElemI;
        cuadraticElemI.append(elemI.at(0));
        cuadraticElemI.append(elemI.at(1));
        cuadraticElemI.append(elemI.at(2));
        cuadraticElemI.append(-1);
        cuadraticElemI.append(-1);
        cuadraticElemI.append(-1);

        double coord[12];
        coord[0] = x[elemI.at(0)];
        coord[1] = y[elemI.at(0)];
        coord[2] = x[elemI.at(1)];
        coord[3] = y[elemI.at(1)];
        coord[4] = x[elemI.at(2)];
        coord[5] = y[elemI.at(2)];

        int index = Contains(edges,elemI.at(0),elemI.at(1));
        if(index != -1)
        {
            //Si ya existia esa arista, se guarda en la posicion 3 el indice que referencia al punto medio de ella
            cuadraticElemI.replace(3, edges[index][2]);
        }
        else
        {
            cuadraticElemI.replace(3, nextIndex++);

            coord[6] = (coord[0]+coord[2])/2;
            coord[7] = (coord[1]+coord[3])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI.at(3),coord[6], coord[7]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(0));
            edgeI.push_back(elemI.at(1));
            edgeI.push_back(cuadraticElemI[3]);
            edges.push_back(edgeI);
        }
        index = Contains(edges,elemI.at(1),elemI.at(2));
        if(index != -1)
        {
            cuadraticElemI.replace(4, edges[index][2]);
        }
        else
        {
            cuadraticElemI.replace(4, nextIndex++);
            coord[8] = (coord[2]+coord[4])/2;
            coord[9] = (coord[3]+coord[5])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI.at(4), coord[8], coord[9]));

            QList<int> edgeI;
            edgeI.push_back(elemI.at(1));
            edgeI.push_back(elemI.at(2));
            edgeI.push_back(cuadraticElemI.at(4));
            edges.push_back(edgeI);
        }
        index = Contains(edges,elemI.at(2),elemI.at(0));
        if(index != -1)
        {
            cuadraticElemI.replace(5, edges[index][2]);
        }
        else
        {
            cuadraticElemI.replace(5, nextIndex++);
            coord[10] = (coord[4]+coord[0])/2;
            coord[11] = (coord[5]+coord[1])/2;
            allMidpoints.push_back(CuadraticVertexData(cuadraticElemI.at(5), coord[10], coord[11]));

            QList<int> edgeI;
            edgeI.push_back(elemI[2]);
            edgeI.push_back(elemI[0]);
            edgeI.push_back(cuadraticElemI.at(5));
            edges.push_back(edgeI);
        }
        quadraticTriangles.push_back(cuadraticElemI);
    }

    int max_val  = 0;
    QList<int> neigh = GetQuadraticNeighbors(quadraticTriangles, nextIndex, max_val);

//    clock_t endQNodes = clock();
//    double diffticksQNodes = (endQNodes - startQNodes);
//    double secQNodes = (diffticksQNodes)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo en construir los nodos cuadraticos" << secQNodes << " segundos.";

//    clock_t startQMatrix = clock();
    ////////////Matriz Global//////////////
//    FMatrix k1g(3*nextIndex,3*nextIndex);
//    FMatrix k2g(3*nextIndex,3*nextIndex);
//    FMatrix k3g(3*nextIndex,3*nextIndex);
//    FMatrix mg(3*nextIndex, 3*nextIndex);

    CMatrix mCg(3*nextIndex,3*(max_val+1), neigh);
    CMatrix k1Cg(3*nextIndex,3*(max_val+1), neigh);
    CMatrix k2Cg(3*nextIndex,3*(max_val+1), neigh);
    CMatrix k3Cg(3*nextIndex,3*(max_val+1), neigh);


    FMatrix *k1 = new FMatrix(18, 18);
    FMatrix *k2 = new FMatrix(18, 18);
    FMatrix *k3 = new FMatrix(18, 18);
    FMatrix *m = new FMatrix(18, 18);
    double dj = 0;
    for(int NE = 0; NE < NET; NE ++)
    {
        QList<int> elemI = triangles.at(NE);
        int cuadraticElemI[6];
        cuadraticElemI[0] = elemI.at(0);
        cuadraticElemI[1] = elemI.at(1);
        cuadraticElemI[2] = elemI.at(2);

        int index3 = Contains(edges,cuadraticElemI[0],cuadraticElemI[1]);
        cuadraticElemI[3] = edges[index3][2];
        int index4 = Contains(edges,cuadraticElemI[1],cuadraticElemI[2]);
        cuadraticElemI[4] = edges[index4][2];
        int index5 = Contains(edges,cuadraticElemI[2],cuadraticElemI[0]);
        cuadraticElemI[5] = edges[index5][2];

        double coord[12];
        coord[0] = x[elemI.at(0)];
        coord[1] = y[elemI.at(0)];
        coord[2] = x[elemI.at(1)];
        coord[3] = y[elemI.at(1)];
        coord[4] = x[elemI.at(2)];
        coord[5] = y[elemI.at(2)];
        coord[6] = allMidpoints[edges[index3][2]- nnt].x;
        coord[7] = allMidpoints[edges[index3][2]- nnt].y;
        coord[8] = allMidpoints[edges[index4][2]- nnt].x;
        coord[9] = allMidpoints[edges[index4][2]- nnt].y;
        coord[10] = allMidpoints[edges[index5][2]- nnt].x;
        coord[11] = allMidpoints[edges[index5][2]- nnt].y;

        /* Formulacion SUBPARAMETRICA */
        double toCheck = MatrixK1K2K3MT6N(coord, f, k1, k2, k3, m);
        if(NE == 0)
            dj = toCheck;
        if(Sign(dj)!= Sign(toCheck))
        {
            qDebug()<<"ERROR EN EL JACOBIANO";
            emit finished();
        }

        /*
        % Para triangulo de 6 Nodos
        % ntab contiene las posiciones en el vector de las incognitas y en las
        % matrices globales correspondientes a los nodos del triangulo NE
        */
//        double ntab[]={3*cuadraticElemI[0],3*cuadraticElemI[0]+1, 3*cuadraticElemI[0]+2,
//                       3*cuadraticElemI[1],3*cuadraticElemI[1]+1, 3*cuadraticElemI[1]+2,
//                       3*cuadraticElemI[2],3*cuadraticElemI[2]+1, 3*cuadraticElemI[2]+2,
//                       3*cuadraticElemI[3],3*cuadraticElemI[3]+1, 3*cuadraticElemI[3]+2,
//                       3*cuadraticElemI[4],3*cuadraticElemI[4]+1, 3*cuadraticElemI[4]+2,
//                       3*cuadraticElemI[5],3*cuadraticElemI[5]+1, 3*cuadraticElemI[5]+2};

//        int DMK = 18;
//        for(int i = 0; i < DMK; i++)
//        {
//            for(int j=0;j<DMK;j++)
//            {
//                if(ntab[i] > ntab[j]) continue;
//                k1g.SetElem(ntab[i],ntab[j],k1g.GetElem(ntab[i],ntab[j])+k1->GetElem(i,j));
//                k2g.SetElem(ntab[i],ntab[j],k2g.GetElem(ntab[i],ntab[j])+k2->GetElem(i,j));
//                k3g.SetElem(ntab[i],ntab[j],k3g.GetElem(ntab[i],ntab[j])+k3->GetElem(i,j));
//                mg.SetElem(ntab[i],ntab[j],mg.GetElem(ntab[i],ntab[j])+m->GetElem(i,j));
//            }
//        }
        mCg.SetElems(*m,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);
        k1Cg.SetElems(*k1,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);
        k2Cg.SetElems(*k2,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);
        k3Cg.SetElems(*k3,cuadraticElemI[0], cuadraticElemI[1], cuadraticElemI[2],
                     cuadraticElemI[3], cuadraticElemI[4], cuadraticElemI[5]);

        if(NE < NET - 1){k1->Clear(); k2->Clear(); k3->Clear(); m->Clear();}
    }
    delete(k1);
    delete(k2);
    delete(k3);
    delete(m);

//    clock_t endQMatrix = clock();
//    double diffticksQMatrix = (endQMatrix - startQMatrix);
//    double secQMatrix  = (diffticksQMatrix)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo en construir las matrices globales cuadraticas" << secQMatrix << " segundos.";

    k1Cg.Symmetrize();
    k2Cg.AntiSymmetrize();
    k3Cg.Symmetrize();
    mCg.Symmetrize();

    NXIh.clear();
    for(double i = step; i <= maxWavesNumber; i+=step)
        NXIh.push_back(i);
    emit maxIteration(NXIh.size());

    for(int i = 0; i < NXIh.size(); i++)
    {
        eigenValsVects.append(QList<AvalAvect>());
        for(int j=0;j<autovalNumber;j++)
            eigenValsVects[i].append(AvalAvect(-1,-1,0,QList<complex<double> >()));
    }

    //Variables
    complex<double> *b;
    complex<double> *a;
    int n = mCg.rows;
    int nvc = min(max(2*autovalNumber,20),n);
    complex<double> *w = new complex<double>[autovalNumber+1];
    complex<double> *v = new complex<double>[n*nvc];

    int info = 0;
    b = mCg.ToComplexArray();
    int cCols = mCg.columns;
    int *neighFortran = ConvertToFortranParam(neigh,nextIndex,max_val);
    //mCg.WriteMatrixAndNeighbors("C://Users/Portege/Desktop/mgC.txt", "C://Users/Portege/Desktop/neigh.txt");
//    mg.WriteMatrix("C://Users/Portege/Desktop/mg.txt");
    //WriteVector("C://Users/Portege/Desktop/neighbors.txt",neighFortran,nextIndex*max_val);
//    FMatrix dmg = DecompactMatrix("C://Users/Portege/Desktop/mgC.txt", "C://Users/Portege/Desktop/neigh.txt", "C://Users/Portege/Desktop/mgD.txt");
//    for(int i = 0; i < dmg.rows; i++)
//    {
//        for(int j = 0; j < dmg.cols; j++)
//        {
//            if(Abs(dmg.GetElem(i,j) - mg.GetElem(i,j)) > 1e-6)
//                qDebug()<<"i="<<i<<"j="<<j << "D:"<< dmg.GetElem(i,j) << "O:" << mg.GetElem(i,j);
//        }
//    }
//    FComplexMatrix *k2giXIh = new FComplexMatrix(k2g.rows,k2g.cols);
//    FComplexMatrix *suma = new FComplexMatrix(k2g.rows,k2g.cols);
//    FComplexMatrix *A = new FComplexMatrix(k2g.rows,k2g.cols);
//    FMatrix *k3gXIh2 = new FMatrix(k2g.rows,k2g.cols);
    CComplexMatrix *k2giXIhC = new CComplexMatrix(k2Cg.rows,cCols, neigh);
    CComplexMatrix *sumaC = new CComplexMatrix(k2Cg.rows,cCols, neigh);
    CComplexMatrix *AC = new CComplexMatrix(k2Cg.rows,cCols, neigh);
    CMatrix *k3gXIh2C = new CMatrix(k2Cg.rows,cCols, neigh);

//    clock_t endQConstructionProcess = clock();
//    double diffticksQConstructionProcess = (endQConstructionProcess - start);
//    double secQConstrictionProcess = (diffticksQConstructionProcess)/(CLOCKS_PER_SEC);
//    qDebug()<< "Tiempo que demora en preparar todo para comenzar a calcular " << secQConstrictionProcess << " segundos.";
    for(int INXIh=0; INXIh < NXIh.size(); INXIh++)
    {
//        clock_t startOnePGA = clock();
        double XIh = NXIh[INXIh];
        //  A = K1G+i*XIh*K2G+XIh^2*K3G;
        complex<double> iXIh;
        iXIh.imag(XIh);

//        k2g.multi_escalar(iXIh,k2giXIh);
//        k3g.multi_escalar(pow(XIh,2),k3gXIh2);
//        k3gXIh2->Add(k2giXIh,suma);
//        k1g.Add(suma,A);
//        a = A->values;
//        A->WriteRealComplexMatrix("C://Users/Portege/Desktop/AR.txt");
//        A->WriteImagComplexMatrix("C://Users/Portege/Desktop/AI.txt");

        k2Cg.multi_escalar(iXIh,k2giXIhC);
        k3Cg.multi_escalar(pow(XIh,2),k3gXIh2C);
        k3gXIh2C->Add(k2giXIhC,sumaC);
        k1Cg.Add(sumaC,AC);
        //AC->WriteRealComplexMatrix("C://Users/Portege/Desktop/ACR.txt");
        //AC->WriteImagComplexMatrix("C://Users/Portege/Desktop/ACI.txt");
//        FComplexMatrix da = DecompactCMatrix("C://Users/Portege/Desktop/ACR.txt", "C://Users/Portege/Desktop/ACI.txt", "C://Users/Portege/Desktop/neigh.txt", "C://Users/Portege/Desktop/ADR.txt", "C://Users/Portege/Desktop/ADI.txt");

//        for(int i = 0; i < da.rows; i++)
//        {
//            for(int j = 0; j < da.cols; j++)
//            {
//                if(Abs(da.GetElem(i,j).real() - A->GetElem(i,j).real()) > 1e-6 || Abs(da.GetElem(i,j).imag() - A->GetElem(i,j).imag()) > 1e-6)
//                {
//                    qDebug()<<"REAL i="<<i<<"j="<<j << "D:"<< da.GetElem(i,j).real() << "O:" << A->GetElem(i,j).real();
//                    qDebug()<<"IMAG i="<<i<<"j="<<j << "D:"<< da.GetElem(i,j).imag() << "O:" << A->GetElem(i,j).imag();
//                }
//            }
//        }

        a = AC->ToArray();

        eigensolvermxmaxb_(b,a,&n,&cCols,&autovalNumber,&nvc,w,v,neighFortran,&info);


        if(info != 0)
        {
            qDebug()<<"info ="<<info<<"chi ="<<XIh;
//          AC->WriteRealComplexMatrix("C://Users/Portege/Desktop/ACR.txt");
//          AC->WriteImagComplexMatrix("C://Users/Portege/Desktop/ACI.txt");
        }

        SortRealNorm(w, eigenValsVects, INXIh, v, autovalNumber, n, XIh);
        emit posChanged(INXIh+1);
        //qDebug() << qSetRealNumberPrecision( 10 ) << "XIh = " << XIh;
        k2giXIhC->Clear();
        sumaC->Clear();
        AC->Clear();
        k3gXIh2C->Clear();
//        k2giXIhC->Clear();
//        sumaC->Clear();
//        AC->Clear();
//        k3gXIh2C->Clear();
//        clock_t endOnePGA = clock();
//        double diffticksOnePGA = (endOnePGA - startOnePGA);
//        double secOnePGA = (diffticksOnePGA)/(CLOCKS_PER_SEC);
//        qDebug()<< "Tiempo que demoró en resolver el PGA " << secOnePGA << " segundos.";
//        int a = 0;

    }
    delete k2giXIhC;
    delete sumaC;
    delete AC;
    delete k3gXIh2C;
//    delete(k2giXIhC);
//    delete(sumaC);
//    delete(AC);
//    delete(k3gXIh2C);
    delete a;
    delete b;
    delete w;
    delete v;
    delete neighFortran;

    QList<int> deleted = PostProcessResults();
    PostProcessEigenValVects(deleted,NXIh.size(),autovalNumber);
    autovalNumber = autovalNumber - deleted.size();
    vf = FMatrix(NXIh.size(),autovalNumber);
    fM = FMatrix(NXIh.size(),autovalNumber);

    for(int INXIh = 0 ; INXIh < NXIh.size(); INXIh++)
    {
        for(int index = 0; index < autovalNumber; index++)
        {
            double val = sqrt(eigenValsVects[INXIh][index].GetAutoVal());
            double aux1 = val/NXIh[INXIh];
            vf.SetElem(INXIh, index, aux1);
            double aux2 = val;
            fM.SetElem(INXIh, index, aux2);
        }
    }

    vg = Diff(NXIh.size(),autovalNumber, eigenValsVects).multi_escalar(vt/step);

    clock_t end = clock();
    double diffticks = end-start;
    secTime =(diffticks)/(CLOCKS_PER_SEC);
    emit finished();
}
int WorkerThread::IsReplica(int index)
{
    int times = 0;
    for(int i = 0 ; i < NXIh.size(); i++)
    {
        if(abs(eigenValsVects[i][index-1].GetAutoVal() - eigenValsVects[i][index].GetAutoVal())<0.00001 ||
           abs(eigenValsVects[i][index].GetAutoVal() - eigenValsVects[i][index+1].GetAutoVal())<0.00001)
            times++;
    }
    return times;
}
QList<int> WorkerThread::PostProcessResults()
{
    double percent = NXIh.size()*9.5/10;
    QList<int> toDelete;
    //VERSION 1 quita curvas que no deberia
//    for(int i = 0; i < NXIh.size(); i++)
//    {
//        for(int j = 0; j < autovalNumber-1; j++)
//        {
//            if(toDelete.contains(j))
//                continue;
//            if(abs(eigenValsVects[i][j].GetAutoVal() - eigenValsVects[i][j+1].GetAutoVal())<0.00001)
//            {
//                int index = WhoGoes(j,j+1,eigenValsVects[i]);
//                if(toDelete.contains(index))
//                    continue;
//                toDelete.push_back(index);
//            }
//        }
//    }
    //VERSION 2 quita aquellas curvas tales que la cantidad de
    //sus puntos se solapan con la anterior y la posterior en un
    //porciento mayor que x
    for(int i = 1 ; i < autovalNumber-1; i++)
    {
        int times = IsReplica(i);
        if(times > percent)
            toDelete.push_back(i);
    }
    return toDelete;
}
int WorkerThread::Contains(QList<int> elems, int elem)
{
    for(int i = 0 ; i < elems.size(); i++)
        if(elems[i]==elem)return i;
    return -1;
}
void WorkerThread::ShiftOneLeft(int pos, int steps, int size)
{
    for(int i = 0; i < steps; i++)
        for(int j = 0 ; j < size-pos-1; j++)
            eigenValsVects[i][pos+j] = eigenValsVects[i][pos+j+1];
}
int WorkerThread::Contains(QList<QList<int> > edges, int v0, int v1)
{
    for(int i = 0 ;i < edges.size();i++)
    {
        if((edges[i][0] == v0 && edges[i][1] == v1) || (edges[i][0] == v1 && edges[i][1] == v0))
            return i;
    }
    return -1;
}
int* WorkerThread::ConvertToFortranParam(QList<int> neighbors,int rows, int cols)
{
    int* result = new int[rows*cols];
    for(int i = 0 ; i < rows*cols; i ++)
        result[i] =(neighbors.at(i) == -1)?neighbors.at(i):neighbors.at(i)+1;
    return result;
}
void WorkerThread::PostProcessEigenValVects(QList<int> deleted, int steps, int size)
{
    //En deleted estan ordenadas las curvas de menor a mayor
    for(int i = 0 ; i < deleted.size(); i++)
        ShiftOneLeft(deleted[deleted.size()-1-i],steps,size);
}
QList<int> WorkerThread::GetNeighbors(QList<QList<int> >triangles, int n, int &ncols)
{
    QList<QList<int> >result;
    for(int i = 0 ; i < n; i++)
        result.push_back(QList<int>());

    for(int i = 0; i<triangles.size();i++)
    {
        int i0 = triangles.at(i)[0];
        int i1 = triangles.at(i)[1];
        int i2 = triangles.at(i)[2];
       // Incorporar i2 y i3 a la lista de los vecinos de i1, etc
        if(Contains(result[i0],i1) == -1)
           //% i2 no habia sido incluido en la lista d evecinos de i1, incluirlo!
           result[i0].push_back(i1);
        if(Contains(result[i0],i2)==-1)
           // i3 no habia sido incluido en la lista d evecinos de i1, incluirlo!
           result[i0].push_back(i2);
        if(ncols < result[i0].size())
            ncols = result[i0].size();

       if(Contains(result[i1],i0) == -1)
//           % i1 no habia sido incluido en la lista d evecinos de i2, incluirlo!
           result[i1].push_back(i0);
       if(Contains(result[i1],i2) == -1)
//           % i3 no habia sido incluido en la lista d evecinos de i2, incluirlo!
           result[i1].push_back(i2);
       if(ncols < result[i1].size())
           ncols = result[i2].size();

       if(Contains(result[i2],i0) == -1)
//           % i1 no habia sido incluido en la lista d evecinos de i3, incluirlo!
           result[i2].push_back(i0);
        if(Contains(result[i2],i1) == -1)
//           % i2 no habia sido incluido en la lista d evecinos de i3, incluirlo!
            result[i2].push_back(i1);
        if(ncols < result[i2].size())
            ncols = result[i2].size();
    }

    QList<int> toreturn;
    for(int i = 0 ; i < ncols; i++)
    {
        for(int s = 0 ; s < n; s++)
        {
            if(i >= result[s].size())
                toreturn.append(-1);
            else toreturn.append(result[s][i]);
        }
    }
    return toreturn;
}
FMatrix WorkerThread::Diff(int filas, int columnas, QList<QList<AvalAvect> > eigensValsVects)
{
    FMatrix res(filas-1, columnas);
    for(int i = 0 ; i < filas-1; i ++)
    {
        QList<AvalAvect> iM1Av = eigensValsVects[i+1];
        QList<AvalAvect> iAv = eigensValsVects[i];
        for(int j = 0 ; j < columnas; j++)
            res.SetElem(i,j,sqrt(iM1Av[j].GetAutoVal()) - sqrt(iAv[j].GetAutoVal()));
    }
    return res;
}
QList<double> WorkerThread::XTimesGe(FMatrix *m, int row, QList<int> elems, int size, bool transpose)
{
    /*
      This function computes the product
      y = x*Ge
      where x is a row vector and Ge is a 9xn matrix (associated to the e-th
      triangle) which elements are all zero except
      Ge(1,3*iv1-2)=1, Ge(2,3*iv1-1)=1, Ge(3,3*iv1)=1
      Ge(4,3*iv2-2)=1, Ge(5,3*iv2-1)=1, Ge(6,3*iv2)=1
      Ge(7,3*iv3-2)=1, Ge(8,3*iv3-1)=1, Ge(9,3*iv3)=1
      Input
      m- matrix
      row- row index of the matrix to multiply
      size- number of vertices of the triangulation
      elemI- indices of vertices of the i-th triangle
      transpose- indicates if it shpuld be taken the column of m instead of the row
      Output
      y= x*Ge, a row vector of length 3*n
    */
    QList<double> result;
    for(int i = 0; i < 3*size; i++)
        result.push_back(0);

    if(!transpose)
    {
        result[3*elems[0]] = m->GetElem(row, 0);
        result[3*elems[0]+1] = m->GetElem(row, 1);
        result[3*elems[0]+2] = m->GetElem(row, 2);

        result[3*elems[1]] = m->GetElem(row, 3);
        result[3*elems[1]+1] = m->GetElem(row, 4);
        result[3*elems[1]+2] = m->GetElem(row, 5);

        result[3*elems[2]] = m->GetElem(row, 6);
        result[3*elems[2]+1] = m->GetElem(row, 7);
        result[3*elems[2]+2] = m->GetElem(row, 8);
    }
    else
    {
        result[3*elems[0]] = m->GetElem(0, row);
        result[3*elems[0]+1] = m->GetElem(1, row);
        result[3*elems[0]+2] = m->GetElem(2, row);

        result[3*elems[1]] = m->GetElem(3, row);
        result[3*elems[1]+1] = m->GetElem(4, row);
        result[3*elems[1]+2] = m->GetElem(5, row);

        result[3*elems[2]] = m->GetElem(6, row);
        result[3*elems[2]+1] = m->GetElem(7, row);
        result[3*elems[2]+2] = m->GetElem(8, row);
    }
    return result;
}
void WorkerThread::SortRealNorm(complex<double> *w, QList<QList<AvalAvect> > &eigenvals,
                                int pos, complex<double> *vr, int maxNumber,int n, double waveNumber)
{
    for(int i = 0 ; i < maxNumber;i++)
    {
        eigenvals[pos][i].SetAutoVal(abs(w[i].real()));
        QList<complex<double> > vector;
        for(int j = 0 ; j < n; j ++)
            vector.push_back(vr[i*n+j]);
        eigenvals[pos][i].SetAutoVect(vector);
        eigenvals[pos][i].SetWaveNumber_Index(waveNumber,pos);
        eigenvals[pos][i].Normalize2AVector();
    }

    for(int i = 0 ;i < maxNumber-1; i ++)
    {
        for(int j = i+1; j < maxNumber; j++)
        {
            if(eigenvals[pos][i].GetAutoVal() > eigenvals[pos][j].GetAutoVal())
            {
                AvalAvect aux = eigenvals[pos][i];
                eigenvals[pos][i] = eigenvals[pos][j];
                eigenvals[pos][j] = aux;
            }
        }
    }
}
void WorkerThread::SortRealNorm(complex<double> *w, QList<QList<AvalAvect> > &eigenvals, int pos,
                                complex<double> *vr, int maxNumber, int n,
                                double waveNumber, QList<int> indicesToDelete, QList<int> symmetryElems)
{
    complex<double> zero;
    for(int i = 0 ; i < maxNumber;i++)
    {
        eigenvals[pos][i].SetAutoVal(abs(w[i].real()));

        QList<complex<double> > vector;
        int count = 0;
        //primer cuadrante
        for(int j = 0 ; j < n; j ++)
        {
            if(indicesToDelete.contains(j))
                vector.push_back(zero);
            vector.push_back(vr[count++]);
        }
         eigenvals[pos][i].SetAutoVect(vector);
         eigenvals[pos][i].SetWaveNumber_Index(waveNumber,pos);
    }

    for(int i = 0 ;i < maxNumber-1; i ++)
    {
        for(int j = i+1; j < maxNumber; j++)
        {
            if(eigenvals[pos][i].GetAutoVal() > eigenvals[pos][j].GetAutoVal())
            {
                AvalAvect aux = eigenvals[pos][i];
                eigenvals[pos][i] = eigenvals[pos][j];
                eigenvals[pos][j] = aux;
            }
        }
    }

    for(int i = 0; i < maxNumber; i++)
    {
       QList<complex<double> > vector = eigenvals[pos][i].GetAutoVect();

        complex<double> mult;
        //segundo cuadrante
        for(int j = 0; j < n; j+=3)
        {
            mult.real(vector[j].real()*symmetryElems[0]);
            mult.imag(vector[j].imag()*symmetryElems[0]);
            vector.push_back(mult);


            mult.real(vector[j+1].real()*symmetryElems[1]);
            mult.imag(vector[j+1].imag()*symmetryElems[1]);
            vector.push_back(mult);

            mult.real(vector[j+2].real()*symmetryElems[2]);
            mult.imag(vector[j+2].imag()*symmetryElems[2]);
            vector.push_back(mult);
        }

        //tercer cuadrante
        for(int j = 0; j < n; j+=3)
        {
            mult.real(vector[j].real()*symmetryElems[3]);
            mult.imag(vector[j].imag()*symmetryElems[3]);
            vector.push_back(mult);

            mult.real(vector[j+1].real()*symmetryElems[4]);
            mult.imag(vector[j+1].imag()*symmetryElems[4]);
            vector.push_back(mult);

            mult.real(vector[j+2].real()*symmetryElems[5]);
            mult.imag(vector[j+2].imag()*symmetryElems[5]);
            vector.push_back(mult);
        }

        //cuarto cuadrante
        for(int j = 0; j < n; j+=3)
        {
            mult.real(vector[j].real()*symmetryElems[6]);
            mult.imag(vector[j].imag()*symmetryElems[6]);
            vector.push_back(mult);

            mult.real(vector[j+1].real()*symmetryElems[7]);
            mult.imag(vector[j+1].imag()*symmetryElems[7]);
            vector.push_back(mult);

            mult.real(vector[j+2].real()*symmetryElems[8]);
            mult.imag(vector[j+2].imag()*symmetryElems[8]);
            vector.push_back(mult);
        }

        eigenvals[pos][i].SetAutoVect(vector);
        eigenvals[pos][i].Normalize2AVector();
    }
}
QList<int> WorkerThread::GetQuadraticNeighbors(QList<QList<int> > quadraticTriangles, int n, int &ncols)
{
    QList<QList<int> >result;
    for(int i = 0 ; i < n; i++)
        result.push_back(QList<int>());

    for(int i = 0; i < quadraticTriangles.size();i++)
    {
        int i0 = quadraticTriangles.at(i)[0];
        int i1 = quadraticTriangles.at(i)[1];
        int i2 = quadraticTriangles.at(i)[2];
        int i3 = quadraticTriangles.at(i)[3];
        int i4 = quadraticTriangles.at(i)[4];
        int i5 = quadraticTriangles.at(i)[5];

        //Incorporar i1, i2, i3, i4 y i5 a la lista de los vecinos de i0
        if(Contains(result[i0],i1) == -1)
           result[i0].push_back(i1);
        if(Contains(result[i0],i2)==-1)
           result[i0].push_back(i2);
        if(Contains(result[i0],i3)==-1)
           result[i0].push_back(i3);
        if(Contains(result[i0],i4)==-1)
           result[i0].push_back(i4);
        if(Contains(result[i0],i5)==-1)
           result[i0].push_back(i5);
        if(ncols < result[i0].size())
            ncols = result[i0].size();

        // Incorporar i0, i2, i3, i4 y i5 a la lista de los vecinos de i1
        if(Contains(result[i1],i0) == -1)
           result[i1].push_back(i0);
       if(Contains(result[i1],i2) == -1)
           result[i1].push_back(i2);
       if(Contains(result[i1],i3) == -1)
           result[i1].push_back(i3);
       if(Contains(result[i1],i4) == -1)
           result[i1].push_back(i4);
       if(Contains(result[i1],i5) == -1)
           result[i1].push_back(i5);
       if(ncols < result[i1].size())
           ncols = result[i1].size();

       // Incorporar i0, i1, i3, i4 y i5 a la lista de los vecinos de i2
       if(Contains(result[i2],i0) == -1)
           result[i2].push_back(i0);
       if(Contains(result[i2],i1) == -1)
           result[i2].push_back(i1);
       if(Contains(result[i2],i3) == -1)
           result[i2].push_back(i3);
       if(Contains(result[i2],i4) == -1)
           result[i2].push_back(i4);
       if(Contains(result[i2],i5) == -1)
           result[i2].push_back(i5);
       if(ncols < result[i2].size())
           ncols = result[i2].size();

       // Incorporar i0, i1, i2, i4 y i5 a la lista de los vecinos de i3
       if(Contains(result[i3],i0) == -1)
           result[i3].push_back(i0);
       if(Contains(result[i3],i1) == -1)
           result[i3].push_back(i1);
       if(Contains(result[i3],i2) == -1)
           result[i3].push_back(i2);
       if(Contains(result[i3],i4) == -1)
           result[i3].push_back(i4);
       if(Contains(result[i3],i5) == -1)
           result[i3].push_back(i5);
       if(ncols < result[i3].size())
           ncols = result[i3].size();

       // Incorporar i0, i1, i2, i3 y i5 a la lista de los vecinos de i4
       if(Contains(result[i4],i0) == -1)
           result[i4].push_back(i0);
       if(Contains(result[i4],i1) == -1)
           result[i4].push_back(i1);
       if(Contains(result[i4],i2) == -1)
           result[i4].push_back(i2);
       if(Contains(result[i4],i3) == -1)
           result[i4].push_back(i3);
       if(Contains(result[i4],i5) == -1)
           result[i4].push_back(i5);
       if(ncols < result[i4].size())
           ncols = result[i4].size();

       // Incorporar i0, i1, i2, i3 y i4 a la lista de los vecinos de i5
       if(Contains(result[i5],i0) == -1)
           result[i5].push_back(i0);
       if(Contains(result[i5],i1) == -1)
           result[i5].push_back(i1);
       if(Contains(result[i5],i2) == -1)
           result[i5].push_back(i2);
       if(Contains(result[i5],i3) == -1)
           result[i5].push_back(i3);
       if(Contains(result[i5],i4) == -1)
           result[i5].push_back(i4);
       if(ncols < result[i5].size())
           ncols = result[i5].size();
    }

    QList<int> toreturn;
    for(int i = 0 ; i < ncols; i++)
    {
        for(int s = 0 ; s < n; s++)
        {
            if(i >= result[s].count())
                toreturn.append(-1);
            else toreturn.append(result[s][i]);
        }
    }
    return toreturn;
}
double WorkerThread::MatrixK1K2K3MT3N(double coord[], double f[], FMatrix *k1, FMatrix *k2, FMatrix *k3, FMatrix *m)
{
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Geometria                                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    double x1 = coord[0];     double y1 = coord[1];
    double x2 = coord[2];     double y2 = coord[3];
    double x3 = coord[4];     double y3 = coord[5];

/*
% La base de funciones lineales que vamos a utilizar sobre cada
% triangulo esta formada por 3 funciones Ni(xi,eta), i=1,...,3 tales que
% Ni(xi_j,eta_j)=0 para j diferente de i, j=1,..,3
% donde (xi_j,eta_j) son las coordenadas locales (xi,eta) del j-esimo nodo.
% Mas precisamente
% N1=xi
% N2=eta
% N3=(1-xi-eta)
% Estas funciones se organizan en una matriz N(xi,eta)(que depende de xi
% y eta) de orden 3x 9 (3 columnas por cada funcion de la base)
% Un punto con coordenadas (x,y) en el interior del triangulo se puede
% escribir como
% x=N1(xi,eta)*x1+N2(xi,eta)*x2+N3(xi,eta)*x3 (1)
% y=N1(xi,eta)*y1+N2(xi,eta)*y2+N3(xi,eta)*y3 (2)

% Esta expresion nos da el paso de las coordenadas locales (xi,eta) a las
% globales (x,y). A traves de su inversa podriamos calcular la expresion de
% (xi,eta) en terminos de (x,y).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Integracion Numerica                                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La formula de integracion numerica que vamos a utilizar sobre el triangulo
% canonico es
% int int f(xi,eta) d xi d eta = 1/3(f(1/2,1/2)+f(0,1/2)+f(1/2,0))
% Esta formula tiene error O(h^3) y es exacta para polinomios cuadraticos
% como los que tenemos que integrar para calcular las matrices K3 y M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    double gxi[] = {1.0/2, 0, 1.0/2 };
    double geta[] = { 1.0/2, 1.0/2, 0 };

    double peso = 1.0/3;
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determinate del Jacobiano de la transformacion de coordenadas (1)-(2)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    double djToCheck = x1*y2 - x1*y3 - x3*y2 - y1*x2 + y1*x3 + y3*x2;
    double dj = abs(djToCheck);
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Matrix Elastica                                                %%%
%%%     C = [F(1) F(2) F(2)  0   0   0;                             %%%
%%%          F(2) F(1) F(2)  0   0   0;                             %%%
%%%          F(2) F(2) F(1)  0   0   0;                             %%%
%%%            0    0    0 F(3)  0   0;                             %%%
%%%            0    0    0   0 F(3)  0;                             %%%
%%%            0    0    0   0   0 F(3)];                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */
    FMatrix c(6, 6);
    c.SetElem(0, 0, f[0]);
    c.SetElem(0, 1, f[1]);
    c.SetElem(0, 2, f[1]);
    c.SetElem(1, 0, f[1]);
    c.SetElem(1, 1, f[0]);
    c.SetElem(1, 2, f[1]);
    c.SetElem(2, 0, f[1]);
    c.SetElem(2, 1, f[1]);
    c.SetElem(2, 2, f[0]);
    c.SetElem(3, 3, f[2]);
    c.SetElem(4, 4, f[2]);
    c.SetElem(5, 5, f[2]);
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Matriz B1                                 %
% Observar que esta matriz NO depende de xi y eta.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    FMatrix b1(6, 9);

   //////////////////////////////////////////////////////////
   ////                   B1                             ////
   //////////////////////////////////////////////////////////
   //(y2 - y3) / DJ 0 0 -(y1 - y3) / DJ 0 0 -(y2 - y3) / DJ + (y1 - y3) / DJ 0 0
   b1.SetElem(0, 0, (y2 - y3) / dj);
   b1.SetElem(0, 3, -(y1 - y3) / dj);
   b1.SetElem(0, 6, -(y2 - y3) / dj + (y1 - y3) / dj);

   //0 -(x2 - x3) / DJ 0 0 (x1 - x3) / DJ 0 0 (x2 - x3) / DJ - (x1 - x3) / DJ 0;
   b1.SetElem(1, 1, -(x2 - x3) / dj);
   b1.SetElem(1, 4, (x1 - x3) / dj);
   b1.SetElem(1, 7, (x2 - x3) / dj - (x1 - x3) / dj);

   //0 0 0 0 0 0 0 0 0

   //0 0 -(x2 - x3) / DJ 0 0 (x1 - x3) / DJ 0 0 (x2 - x3) / DJ - (x1 - x3) / DJ
   b1.SetElem(3, 2, -(x2 - x3) / dj);
   b1.SetElem(3, 5, (x1 - x3) / dj);
   b1.SetElem(3, 8, (x2 - x3) / dj - (x1 - x3) / dj);

   //0 0 (y2 - y3) / DJ 0 0 -(y1 - y3) / DJ 0 0 -(y2 - y3) / DJ + (y1 - y3) / DJ;
   b1.SetElem(4, 2, (y2 - y3) / dj);
   b1.SetElem(4, 5, -(y1 - y3) / dj);
   b1.SetElem(4, 8, -(y2 - y3) / dj + (y1 - y3) / dj);

   // -(x2 - x3) / DJ (y2 - y3) / DJ 0 (x1 - x3) / DJ -(y1 - y3) / DJ 0 (x2 - x3) / DJ - (x1 - x3) / DJ -(y2 - y3) / DJ + (y1 - y3) / DJ 0
   b1.SetElem(5, 0, -(x2 - x3) / dj);
   b1.SetElem(5, 1, (y2 - y3) / dj);
   b1.SetElem(5, 3, (x1 - x3) / dj);
   b1.SetElem(5, 4, -(y1 - y3) / dj);
   b1.SetElem(5, 6, (x2 - x3) / dj - (x1 - x3) / dj);
   b1.SetElem(5, 7, -(y2 - y3) / dj + (y1 - y3) / dj);
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Matriz B1t*C*B1  orden  9x9                                  %%%
% Observar que esta matriz NO depende de xi y eta. Por tanto                    %
%  K1=Integral_T (B1t*C*B1 ) dxdy = Integral_Tcanonico (B1t*C*B1 )*DJ dxideta
%                                 = DJ*(B1t*C*B1)=Area(T)*(B1t*C*B1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
   //determinacion de b1'*c*b1
   FMatrix b1Tc = b1.transpuesta() * c;
   FMatrix b1Tcb1 = b1Tc* b1;
   b1Tcb1.multi_escalar(dj,k1);

   FMatrix n(3, 9);
   FMatrix b2(6, 9);
   FMatrix cb1 = c*b1;

    for (int i = 0; i < 3; i++)
    {
        double xi = gxi[i];
        double eta = geta[i];

        n.SetElem(0, 0, xi);
        n.SetElem(1, 1, xi);
        n.SetElem(2, 2, xi);
        n.SetElem(0, 3, eta);
        n.SetElem(1, 4, eta);
        n.SetElem(2, 5, eta);
        n.SetElem(0, 6, 1 - xi - eta);
        n.SetElem(1, 7, 1 - xi - eta);
        n.SetElem(2, 8, 1 - xi - eta);

        //////////////////////////////////////////////////////////
        ////                   B2                             ////
        //////////////////////////////////////////////////////////
        b2.SetElem(2, 2, xi);
        b2.SetElem(2, 5, eta);
        b2.SetElem(2, 8, 1 - xi - eta);
        b2.SetElem(3, 1, xi);
        b2.SetElem(3, 4, eta);
        b2.SetElem(3, 7, 1 - xi - eta);
        b2.SetElem(4, 0, xi);
        b2.SetElem(4, 3, eta);
        b2.SetElem(4, 6, 1 - xi - eta);

        //determinacion de b1'*c*b2
        FMatrix b12 = b1Tc * b2;
        //determinacion de b1'*c*b2
        FMatrix b21 = b2.transpuesta() * cb1;
        //determinacion de b2'*c*b2
        FMatrix b22 = b2.transpuesta() * c * b2;
        //determinacion de n'*n
        FMatrix n2 = n.transpuesta() * n;

        k2->AddThis((b12-b21).multi_escalar(peso));
        k3->AddThis(b22.multi_escalar(peso));
        m->AddThis(n2.multi_escalar(peso));
    }
/*
% OJO Hay que multiplicar K2, K3 y M por el determinante de la Jacobiana de la
% transformacion (1)-(2), que es constante como funcion de xi y eta

% OJO- tenemos que multiplicar cada sumando por el determinante de la
% Jacobiana de la transformacion, para cambiar la integral sobre un
% triangulo arbitrario T por la integral sobre el triangulo canonico.
% Como el determinante de la  Jacobiana DJ NO DEPENDE de xi y eta,
% podemos multiplicar la integral final por DJ.
*/
    k2->multi_escalar_this(dj);
    k3->multi_escalar_this(dj);
    m->multi_escalar_this(dj);

    return djToCheck;
}
double WorkerThread::MatrixK1K2K3MT6N(double coord[], double f[], FMatrix *k1, FMatrix *k2, FMatrix *k3, FMatrix *m)
{
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Julio 2012 Elementos finitos cuadraticos (6 nodos, 6 funciones basicas:
% N1, N2, ...N6)
% Formulacion SUBPARAMETRICA
% o sea las funciones de la base (o de forma) son cuadraticas y la
% transformacion de las coordenadas afines a las baricentricas es LINEAL.
% Esto ultimo tiene la ventaja de que la Jacobiana de la transformacion de
% coordenadas es una funcion constante (no depende de xi y eta) y su
% determinante tambien. Por lo tanto, las integrales que debemos calcular
% para hallar las matrices K1,K2,K3 y M  son integrales de polinomios de
% grado menor o igual que 4 que se integran exactamente con la cuadratura
% de 7 puntos que estamos utilizando. Esto no ocurre en la version
% isoparametrica, donde las integrales pueden involucrar funciones racionales con
% un polinomio de grado 4 en el numerador y uno de grado 2 en el
% denominador (el determinante de la Jacobiana, que es una funcion
% cuadratica de xi y eta).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Inicializacion de las Matrices                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
    FMatrix N(3, 18);
    FMatrix dNdx(3,18);
    FMatrix dNdy(3,18);
    FMatrix L1(6,3);
    FMatrix L2(6,3);
    FMatrix L3(6,3);

    L1.SetElem(0,0,1);
    L1.SetElem(4,2,1);
    L1.SetElem(5,1,1);
    L2.SetElem(1,1,1);
    L2.SetElem(3,2,1);
    L2.SetElem(5,0,1);
    L3.SetElem(2,2,1);
    L3.SetElem(3,1,1);
    L3.SetElem(4,0,1);

    /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Geometria                                                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     (xi,yi) son las coordenadas de los nodos del triangulo que se esta
     analizando. Como estamos utilizando elementos finitos cuadraticos los
     nodos son los vertices del triangulo y un punto intermedios por cada arista
    */
    double x1 = coord[0];     double y1 = coord[1];
    double x2 = coord[2];     double y2 = coord[3];
    double x3 = coord[4];     double y3 = coord[5];

    /*
    % La base de funciones cuadraticas que vamos a utilizar sobre cada
    % triangulo esta formada por 6 funciones Ni(xi,eta), i=1,...,6 tales que
    % Ni(xi_j,eta_j)=0 para j diferente de i, j=1,..,6
    % donde (xi_j,eta_j) son las coordenadas locales (xi,eta) del j-esimo nodo.
    % Mas precisamente
    % N1=xi*(2*xi-1)
    % N2=eta*(2*eta-1)
    % N3=(1-xi-eta)*(2*(1-xi-eta)-1)
    % N4=4*xi*eta
    % N5=4*(1-xi-eta)*eta
    % N6=4*xi*(1-xi-eta)
    % Estas funciones se organizan en una matriz N(xi,eta)(que depende de xi
    % y eta) de orden 3x 18 (3 columnas por cada funcion de la base)

    % Un punto con coordenadas (x,y) en el interior del triangulo se puede
    % escribir como
    % x=N1(xi,eta)*x1+N2(xi,eta)*x2+N3(xi,eta)*x3      (1)
    % y=N1(xi,eta)*y1+N2(xi,eta)*y2+N3(xi,eta)*y3      (2)

    % Esta expresion nos da el paso de las coordenadas locales (xi,eta) a las
    % globales (x,y). A traves de su inversa podriamos calcular la expresion de
    % (xi,eta) en terminos de (x,y). Estas expresiones se pueden sustituir en
    % N(xi,eta) para obtener una matriz N(x,y). En la practica necesitamos las
    % derivadas parciales de N(x,y) respecto a (x,y). Una forma mas facil de
    % hacerlo es utilizar la regla de la cadena, es decir calcular las derivadas
    % parciales de N(xi,eta) respecto a xi y eta y multiplicar el resultado por
    % las derivadas parciales de xi y eta respecto a x y y. Estas ultimas se
    % pueden obtener calculando la matriz inversa de la Jacobiana de la
    % transformacion dada por (1) y (2).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determinate del Jacobiano de la transformacion de coordenadas (1)-(2)
    % (es constante como funcion de xi y eta y es igual a 2 veces el area del
    % triangulo)
    % Transpuesta del Jacobiano
    */
    double j11 = x1 - x3;
    double j12 = y1 - y3;
    double j21 = x2 - x3;
    double j22 = y2 - y3;

    /*
    % Determinante del Jacobiano, es una funcion que NO depende de xi y eta,
    % solo depende del triangulo
    */
    double djToCheck = j11*j22-j12*j21;
    double dj = abs(djToCheck);

    /*
    % Inversa de la transpuesta del Jacobiano
    */
    double invj11 = j22/dj;
    double invj12 = -j12/dj;
    double invj21 = -j21/dj;
    double invj22 = j11/dj;

    /*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Integracion Numerica                                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % La formula de integracion numerica que vamos a utilizar sobre el triangulo
    % canonico es una formula basada en 7 puntos
    % Esta formula tiene error O(h^6) y es exacta para polinomios de grado 5.

    % int int f(xi,eta) d xi d eta = Area(T)* sum_{i=1,..,7} wi*f(xi_i,eta_i)
    % donde
    %(xi, eta) --> peso
    %(1/3,1/3)  --> 0.225
    %(  0.0597158717 ,  0.4701420641 ),( 0.4701420641 ,  0.0597158717 ), (0.4701420641  , 0.4701420641 ) --> 0.1323941527
    %(  0.7974269853  , 0.1012865073), ( 0.1012865073 ,  0.7974269853 ), ( 0.1012865073 , 0.1012865073) --> 0.1259391805

    */
    double gxi[] = {1.0/3, 0.0597158717, 0.4701420641, 0.4701420641, 0.7974269853, 0.1012865073, 0.1012865073};
    double geta[] = {1.0/3, 0.4701420641, 0.0597158717, 0.4701420641, 0.1012865073, 0.7974269853, 0.1012865073};
    double peso[] = {0.225,0.1323941527,0.1323941527,0.1323941527,0.1259391805,0.1259391805,0.1259391805};

    /*
    %%%  Matrix Elastica                                                %%%
        C = [F(1) F(2) F(2)  0   0   0;
             F(2) F(1) F(2)  0   0   0;
             F(2) F(2) F(1)  0   0   0;
               0    0    0 F(3)  0   0;
               0    0    0   0 F(3)  0;
               0    0    0   0   0 F(3)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    */
    FMatrix C(6,6);
    C.SetElem(0,0,f[0]);
    C.SetElem(0,1,f[1]);
    C.SetElem(0,2,f[1]);
    C.SetElem(1,0,f[1]);
    C.SetElem(1,1,f[0]);
    C.SetElem(1,2,f[1]);
    C.SetElem(2,0,f[1]);
    C.SetElem(2,1,f[1]);
    C.SetElem(2,2,f[0]);
    C.SetElem(3,3,f[2]);
    C.SetElem(4,4,f[2]);
    C.SetElem(5,5,f[2]);

    for(int i = 0; i < 7; i++)
    {
        double xi = gxi[i];
        double eta = geta[i];

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Matriz N                                                    %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %N(1,1:18) = [xi*(2*xi-1) 0 0 eta*(2*eta-1) 0 0 (1-xi-eta)*(1-2*xi-2*eta) 0 0 4*xi*eta 0 0 4*(1-xi-eta)*eta 0 0 4*xi*(1-xi-eta) 0 0];
        %N(2,1:18) = [ 0 xi*(2*xi-1) 0 0 eta*(2*eta-1) 0 0 (1-xi- eta)*(1-2*xi-2*eta) 0 0 4*xi*eta 0 0 4*(1-xi-eta)*eta 0 0 4*xi*(1-xi-eta) 0];
        %N(3,1:18) = [ 0 0 xi*(2*xi-1) 0 0 eta*(2*eta-1) 0 0 (1-xi-eta)*(1-2*xi-2*eta) 0 0 4*xi*eta 0 0 4*(1-xi-eta)*eta 0 0 4*xi*(1-xi-eta)];

        % Calculo de las derivadas parciales de N(x,y) respecto a (x,y).
        % Utilizamos la regla de la cadena, es decir calculamos las derivadas
        % parciales de cada funcion Nj(xi,eta) j=1,...,6 respecto a xi y eta y
        % multiplicamos el resultado por las derivadas parciales de xi y eta
        % respecto a x y y. Estas ultimas estan en la inversa de la matriz
        % de la Jacobiana de la transformacion dada por (1) y (2).

        % N1=xi*(2*xi-1)
        */
        double n1 = xi*(2*xi-1);
        double dn1dxi = 4*xi-1;
        double dn1dx = dn1dxi*invj11;
        double dn1dy = dn1dxi*invj21;

        double n2 = eta*(2*eta-1);
        double dn2deta = 4*eta-1;
        double dn2dx = dn2deta*invj12;
        double dn2dy = dn2deta*invj22;

        double n3 = (1-xi-eta)*(2*(1-xi-eta)-1);
        double dn3dxi = -3+4*xi+4*eta;
        double dn3deta = -3+4*xi+4*eta;
        double dn3dx = dn3dxi*invj11+dn3deta*invj12;
        double dn3dy = dn3dxi*invj21+dn3deta*invj22;

        double n4 = 4*xi*eta;
        double dn4dxi = 4*eta;
        double dn4deta = 4*xi;
        double dn4dx = dn4dxi*invj11+dn4deta*invj12;
        double dn4dy = dn4dxi*invj21+dn4deta*invj22;

        double n5 = 4*(1-xi-eta)*eta;
        double dn5dxi = -4*eta;
        double dn5deta = -8*eta+4-4*xi;
        double dn5dx = dn5dxi*invj11+dn5deta*invj12;
        double dn5dy = dn5dxi*invj21+dn5deta*invj22;

        double n6 = 4*xi*(1-xi-eta);
        double dn6dxi = 4-8*xi-4*eta;
        double dn6deta = -4*xi;
        double dn6dx = dn6dxi*invj11+dn6deta*invj12;
        double dn6dy = dn6dxi*invj21+dn6deta*invj22;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % N es una matriz de orden 3x18 con las estructura
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % [ N1 0  0  N2 0  0 .... N6 0  0 ]
        % [ 0  N1 0  0  N2 0 .... 0  N6 0 ]
        % [ 0  0  N1 0  0  N2.... 0  0  N6]
        */
        N.SetElem(0,0,n1); N.SetElem(0,3,n2); N.SetElem(0,6,n3); N.SetElem(0,9,n4);  N.SetElem(0,12,n5); N.SetElem(0,15,n6);
        N.SetElem(1,1,n1); N.SetElem(1,4,n2); N.SetElem(1,7,n3); N.SetElem(1,10,n4); N.SetElem(1,13,n5); N.SetElem(1,16,n6);
        N.SetElem(2,2,n1); N.SetElem(2,5,n2); N.SetElem(2,8,n3); N.SetElem(2,11,n4); N.SetElem(2,14,n5); N.SetElem(2,17,n6);

        /*
        % dNdx es la matriz q
        % [ dN1dx  0      0       dNdx2 0      0      .... dN6dx  0        0 ]
        % [ 0      dN1dx  0       0     dN2dx  0       .... 0     dN6dx    0 ]
        % [ 0      0      dN1dx   0     0      dN2dx   .... 0     0     dN6dx]
        */
        dNdx.SetElem(0,0,dn1dx); dNdx.SetElem(0,3,dn2dx); dNdx.SetElem(0,6,dn3dx); dNdx.SetElem(0,9,dn4dx);  dNdx.SetElem(0,12,dn5dx); dNdx.SetElem(0,15,dn6dx);
        dNdx.SetElem(1,1,dn1dx); dNdx.SetElem(1,4,dn2dx); dNdx.SetElem(1,7,dn3dx); dNdx.SetElem(1,10,dn4dx); dNdx.SetElem(1,13,dn5dx); dNdx.SetElem(1,16,dn6dx);
        dNdx.SetElem(2,2,dn1dx); dNdx.SetElem(2,5,dn2dx); dNdx.SetElem(2,8,dn3dx); dNdx.SetElem(2,11,dn4dx); dNdx.SetElem(2,14,dn5dx); dNdx.SetElem(2,17,dn6dx);

        /*
        % dNdy es similar
        */
        dNdy.SetElem(0,0,dn1dy); dNdy.SetElem(0,3,dn2dy); dNdy.SetElem(0,6,dn3dy); dNdy.SetElem(0,9,dn4dy);  dNdy.SetElem(0,12,dn5dy); dNdy.SetElem(0,15,dn6dy);
        dNdy.SetElem(1,1,dn1dy); dNdy.SetElem(1,4,dn2dy); dNdy.SetElem(1,7,dn3dy); dNdy.SetElem(1,10,dn4dy); dNdy.SetElem(1,13,dn5dy); dNdy.SetElem(1,16,dn6dy);
        dNdy.SetElem(2,2,dn1dy); dNdy.SetElem(2,5,dn2dy); dNdy.SetElem(2,8,dn3dy); dNdy.SetElem(2,11,dn4dy); dNdy.SetElem(2,14,dn5dy); dNdy.SetElem(2,17,dn6dy);

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Matriz B1 de orden 6x18                                    %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % B1=L1*dN/dx +L2* dN/dy
        */
        FMatrix B1=L1*dNdx + L2*dNdy;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Matriz B2 de orden 6x18                                      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % B2=L3*N
        */
        FMatrix B2=L3*N;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     DETERMINACION DE K1=B1T*C*B1 de orden 18x18                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        */
        FMatrix B11 = B1.transpuesta()*C*B1;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     DETERMINACION DE B1T*C*B2   orden 18x18                     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        */
        FMatrix B12 = B1.transpuesta()*C*B2;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     DETERMINACION DE B2T*C*B1   orden 18x18                     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        */
        FMatrix B21 = B2.transpuesta()*C*B1;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %****** DETERMINACION DE K3=B2T*C*B2   orden 18x18                  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        */
        FMatrix B22 = B2.transpuesta()*C*B2;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     DETERMINACION DE M=NT*N  orden 18x18                        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        */
        FMatrix N2 = N.transpuesta()*N;

        /*
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % evaluacion del i-esimo sumando para la integracion numerica
        % OJO- tenemos que multiplicar cada sumando por el determinante de la
        % Jacobiana de la transformacion, para cambiar la integral sobre un
        % triangulo arbitrario T por la integral sobre el triangulo canonico.
        % Como el determinante de la  Jacobiana DJ NO DEPENDE de xi y eta,
        % podemos multiplicar la integral final por DJ.
        %
        */
        k1->AddThis(B11.multi_escalar(peso[i]));
        k2->AddThis((B12-B21).multi_escalar(peso[i]));
        k3->AddThis(B22.multi_escalar(peso[i]));
        m->AddThis(N2.multi_escalar(peso[i]));
    }
    k1->multi_escalar_this(dj);
    k2->multi_escalar_this(dj);
    k3->multi_escalar_this(dj);
    m->multi_escalar_this(dj);

    return djToCheck;
}
QList<complex<double> > WorkerThread::XTimesGe(FComplexMatrix *m, int row, QList<int> elems, int size, bool transpose)
{
    /*
      This function computes the product
      y = x*Ge
      where x is a row vector and Ge is a 9xn matrix (associated to the e-th
      triangle) which elements are all zero except
      Ge(1,3*iv1-2)=1, Ge(2,3*iv1-1)=1, Ge(3,3*iv1)=1
      Ge(4,3*iv2-2)=1, Ge(5,3*iv2-1)=1, Ge(6,3*iv2)=1
      Ge(7,3*iv3-2)=1, Ge(8,3*iv3-1)=1, Ge(9,3*iv3)=1
      Input
      m- matrix
      row- row index of the matrix to multiply
      size- number of vertices of the triangulation
      elemI- indices of vertices of the i-th triangle
      transpose- indicates if it shpuld be taken the column of m instead of the row
      Output
      y= x*Ge, a row vector of length 3*n
    */
    QList<complex<double> > result;
    for(int i = 0; i < 3*size; i++)
    {
        complex<double> zero;
        zero.real(0);
        zero.imag(0);
        result.push_back(zero);
    }

    if(!transpose)
    {
        result[3*elems[0]] = m->GetElem(row, 0);
        result[3*elems[0]+1] = m->GetElem(row, 1);
        result[3*elems[0]+2] = m->GetElem(row, 2);

        result[3*elems[1]] = m->GetElem(row, 3);
        result[3*elems[1]+1] = m->GetElem(row, 4);
        result[3*elems[1]+2] = m->GetElem(row, 5);

        result[3*elems[2]] = m->GetElem(row, 6);
        result[3*elems[2]+1] = m->GetElem(row, 7);
        result[3*elems[2]+2] = m->GetElem(row, 8);
    }
    else
    {
        result[3*elems[0]] = m->GetElem(0, row);
        result[3*elems[0]+1] = m->GetElem(1, row);
        result[3*elems[0]+2] = m->GetElem(2, row);

        result[3*elems[1]] = m->GetElem(3, row);
        result[3*elems[1]+1] = m->GetElem(4, row);
        result[3*elems[1]+2] = m->GetElem(5, row);

        result[3*elems[2]] = m->GetElem(6, row);
        result[3*elems[2]+1] = m->GetElem(7, row);
        result[3*elems[2]+2] = m->GetElem(8, row);
    }
    return result;
}

void WorkerThread::run()
{
    poisson = 0.5*(pow(longVel,2)-2*pow(shearVel,2))/(pow(longVel,2)-pow(shearVel,2));
    density = density;
    young = 2*density*(1+poisson)*(pow(shearVel,2));
//    shearVel = sqrt(young/(2*density*(1+poisson)));
//    longVel = sqrt((2*(young/(2*density*(1+poisson)))*(1-poisson))/(1-2*poisson));

//    doBeattieWork();
//    if(symmetry)
//        doSymmetricWork();
//    else
//    {
        if(problemType == Linear)
            doLinealCompactedWork();//CG
            //doLinealCGCompactedWork();
        else
            doQuadraticCGCompactedWork();//CG
//    }
}


IOThread::IOThread()
{
    save = false;
    error = 0;
}
IOThread::IOThread(QString ppath, QString pmaterialName,
                   double pshearVel, double plongVel, double pdensity,
                   FMatrix pf, FMatrix pfF, FMatrix plF, FMatrix ptF,
                   FMatrix pvF, FMatrix pfvF, FMatrix plvF, FMatrix ptvF,
                   FMatrix pvg, FMatrix pfvg, FMatrix plvg, FMatrix ptvg,
                   QList<QList<AvalAvect> > pdisplacements,
                   QList<QList<AvalAvect> > pfdisplacements,
                   QList<QList<AvalAvect> > pldisplacements,
                   QList<QList<AvalAvect> > ptdisplacements,
                   int pmaxCoordXFvsVF, int pmaxCoordYFvsVF,
                   int pmaxCoordXFvsVG, int pmaxCoordYFvsVG,
                   int pmaxCoordXKvsF, int pmaxCoordYKvsF,
                   double pdelthaZ, int psectionsCount,
                   QList<QList<GVector2> > pboundaries,
                   QList<QList<int> > pbIndices,
                   ShapeType pstructureType, double pir,
                   double pthickness, double pratio)
{
    save = true;
    path = ppath;
    materialName = pmaterialName;
    shearVel = pshearVel;
    longVel = plongVel;
    density = pdensity;
    f = pf;
    fF = pfF;
    lF = plF;
    tF = ptF;
    vF = pvF;
    fvF = pfvF;
    lvF = plvF;
    tvF = ptvF;
    vg = pvg;
    fvg = pfvg;
    lvg = plvg;
    tvg = ptvg;
    displacements = pdisplacements;
    fdisplacements = pfdisplacements;
    ldisplacements = pldisplacements;
    tdisplacements = ptdisplacements;
    maxCoordXFvsVF = pmaxCoordXFvsVF;
    maxCoordYFvsVF = pmaxCoordYFvsVF;
    maxCoordXFvsVG = pmaxCoordXFvsVG;
    maxCoordYFvsVG = pmaxCoordYFvsVG;
    maxCoordXKvsF = pmaxCoordXKvsF;
    maxCoordYKvsF = pmaxCoordYKvsF;
    delthaZ = pdelthaZ;
    sectionsCount = psectionsCount;
    boundaries = pboundaries;
    bIndices = pbIndices;
    structureType = pstructureType;
    ir = pir;
    thickness = pthickness;
    ratio = pratio;
    error = 0;
}

void IOThread::SaveExample()
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save .xpl File");
        return;
    }
    QTextStream out(&file);
    out << "XPL\n";

    SaveToFile(out, f, vF, vg, displacements);
    SaveToFile(out, fF, fvF, fvg, fdisplacements);
    SaveToFile(out, lF, lvF, lvg, ldisplacements);
    SaveToFile(out, tF, tvF, tvg, tdisplacements);

    out << materialName << "~" << shearVel << "~" << longVel << "~" << density << "~\n";
    out << maxCoordXFvsVF << "~" << maxCoordYFvsVF << "~";
    out << maxCoordXFvsVG << "~" << maxCoordYFvsVG << "~";
    out << maxCoordXKvsF << "~" << maxCoordYKvsF << "~\n";
    out << structureType << "~" << ir << "~" << thickness << "~" << ratio << "~\n";

    if(boundaries.count() > 0)
    {
        out << delthaZ << "~";
        out << sectionsCount << "~";

        out << boundaries.count() << "~\n";
        for(int i = 0; i < boundaries.count(); i++)
        {
            out << boundaries[i].count() << "~";
            for(int j = 0 ; j < boundaries[i].count(); j++)
            {
                out << boundaries[i][j].ax << "~";
                out << boundaries[i][j].ay << "~";
            }
            out << "\n";
        }

        for(int i = 0; i < boundaries.count(); i++)
        {
            out << bIndices[i].count() << "~";
            for(int j = 0; j < bIndices[i].count(); j++)
                out << bIndices[i][j] << "~";
            out << "\n";
        }
    }
    else
        out << "!~\n";
    out << "!;";

    file.close();
}
void IOThread::LoadExample()
{
    QFile file(path);
    if(!file.open(QFile::ReadOnly))
    {
        printf("Fail To Read .xpl File");
        error = -1;
        return;
    }

    QTextStream out(&file);
    QString line = out.readLine();
    if(line != "XPL")
    {
        QMessageBox msgBox;
        msgBox.setText("File .xpl is corrupt");
        msgBox.exec();
        error = -1;
        return;
    }
    error = ReadFromFile(out, f, vF, vg, displacements, chi);
    if(error != 0) return;
    error = ReadFromFile(out, fF, fvF, fvg, fdisplacements, chi);
    if(error != 0)return;
    error = ReadFromFile(out, lF, lvF, lvg, ldisplacements, chi);
    if(error != 0)return;
    error = ReadFromFile(out, tF, tvF, tvg, tdisplacements, chi);
    if(error != 0)return;

    line = out.readLine();
    QStringList allData = line.split('~');
    if(allData.count()!= 5)
    {
        error = -1;
        return;
    }
    int index = 0;
    materialName = allData[index++];
    shearVel = ((QString)allData[index++]).toDouble();
    longVel = ((QString)allData[index++]).toDouble();
    density = ((QString)allData[index++]).toDouble();

    line = out.readLine();
    allData = line.split('~');
    if(allData.count()!= 7)
    {
        error = -1;
        return;
    }
    index = 0;
    maxCoordXFvsVF = ((QString)allData[index++]).toDouble();
    maxCoordYFvsVF = ((QString)allData[index++]).toDouble();
    maxCoordXFvsVG = ((QString)allData[index++]).toDouble();
    maxCoordYFvsVG = ((QString)allData[index++]).toDouble();
    maxCoordXKvsF = ((QString)allData[index++]).toDouble();
    maxCoordYKvsF = ((QString)allData[index++]).toDouble();

    line = out.readLine();
    allData = line.split('~');
    if(allData.count()!= 5)
    {
        error = -1;
        return;
    }
    index = 0;
    structureType = (ShapeType)((QString)allData[index++]).toInt();
    ir = ((QString)allData[index++]).toDouble();
    thickness = ((QString)allData[index++]).toDouble();
    ratio = ((QString)allData[index++]).toDouble();

    line = out.readLine();
    allData = line.split('~');
    index = 0;
    if((QString)allData[index] != "!")
    {
        if(allData.count()!= 4)
        {
            error = -1;
            return;
        }
        delthaZ = ((QString)allData[index++]).toDouble();
        sectionsCount = ((QString)allData[index++]).toInt();
        int boundariesCount = ((QString)allData[index++]).toInt();

        int verticesCount = -1;
        for(int i = 0; i < boundariesCount; i++)
        {
            line = out.readLine();
            allData = line.split('~');
            index = 0;
            verticesCount = ((QString)allData[index++]).toInt();
            if(allData.count()!= verticesCount*2+2)
            {
                error = -1;
                return;
            }
            boundaries.append(QList<GVector2>());
            for(int j = 0 ; j < verticesCount; j++)
            {
                double x = ((QString)allData[index++]).toDouble();
                double y = ((QString)allData[index++]).toDouble();
                boundaries[i].append(GVector2(x,y));
            }
        }

        for(int i = 0; i < boundariesCount; i++)
        {
            bIndices.append(QList<int>());
            line = out.readLine();
            allData = line.split('~');
            index = 0;
            verticesCount = ((QString)allData[index++]).toInt();
            if(allData.count()!= verticesCount+2)
            {
                error = -1;
                return;
            }

            for(int j = 0; j < verticesCount; j++)
                bIndices[i].append(((QString)allData[index++]).toInt());
        }
    }

    line = out.readLine();
    allData = line.split('~');
    if(allData.count() != 1)
    {
        error = -1;
    }
}
void IOThread::SaveToFile(QTextStream &out , FMatrix f, FMatrix vF, FMatrix vg, QList<QList<AvalAvect> > displacements)
{
    if(f.rows > 0)
    {
        out << f.rows << "~" << f.cols << "~\n";
        for(int i = 0; i < f.rows; i++)
        {
            for(int j = 0; j < f.cols; j++)
                out << f.GetElem(i,j) << "~";
            out << "\n";
        }
        for(int i = 0; i < vF.rows; i++)
        {
            for(int j = 0; j < vF.cols; j++)
                out << vF.GetElem(i,j) << "~";
            out << "\n";
        }

        for(int i = 0; i < vg.rows; i++)
        {
            for(int j = 0; j < vg.cols; j++)
                out << vg.GetElem(i,j) << "~";
            out << "\n";
        }

        out << displacements.count() << "~\n" ;
        for(int i = 0; i < displacements.count(); i++)
        {
            out << displacements[i].count() << "~\n";
            for(int j = 0; j < displacements[i].count(); j++)
            {
                out << displacements[i][j].GetWaveNumber() << "~";
                out << displacements[i][j].GetIndex() << "~";
                out << displacements[i][j].GetAutoVal() << "~";
                out << displacements[i][j].GetAutoVect().count() << "~";
                for(int k = 0; k < displacements[i][j].GetAutoVect().count(); k++)
                {
                    out <<  displacements[i][j].GetAutoVect()[k].real() << "~";
                    out <<  displacements[i][j].GetAutoVect()[k].imag() << "~";
                }
                out << "\n";
            }
        }
    }
    else
        out << "!~\n";
}
int IOThread::ReadFromFile(QTextStream &out, FMatrix &f, FMatrix &vF, FMatrix &vg,
                  QList<QList<AvalAvect> > &displacements, QList<double> &chi)
{
    QString line = out.readLine();
    QStringList allData = line.split('~');
    if(allData[0] == "!") return 0;
    else
    {
        if(allData.count()!= 3) return -1;
        int row = ((QString)allData[0]).toInt();
        int col = ((QString)allData[1]).toInt();
        f.ChangeDimensions(row,col);
        line = out.readLine();
        allData = line.split('~');
        int count = 0;
        for(int i = 0; i < row; i++)
        {
            if(allData.count()!= col+1) return -1;
            for(int j = 0; j < col; j++)
                f.SetElem(i,j,((QString)allData[count++]).toDouble());
            line = out.readLine();
            allData = line.split('~');
            count = 0;
        }

        vF.ChangeDimensions(row,col);
        for(int i = 0; i < row; i++)
        {
            if(allData.count()!= col+1) return -1;
            for(int j = 0; j < col; j++)
                vF.SetElem(i,j,((QString)allData[count++]).toDouble());
            count = 0;
            line = out.readLine();
            allData = line.split('~');
        }

        vg.ChangeDimensions(row-1,col);
        for(int i = 0; i < row-1; i++)
        {
            if(allData.count()!= col+1) return -1;
            for(int j = 0; j < col; j++)
                vg.SetElem(i,j,((QString)allData[count++]).toDouble());
            count = 0;
            line = out.readLine();
            allData = line.split('~');
        }

        int dispRow = ((QString)allData[0]).toInt();
        for(int i = 0; i < dispRow; i++)
        {
            displacements.append(QList<AvalAvect>());
            count = 0;
            line = out.readLine();
            allData = line.split('~');
            int dispCol = ((QString)allData[0]).toInt();
            for(int j = 0; j < dispCol; j++)
            {
                line = out.readLine();
                allData = line.split('~');
                count = 0;
                if(j == 0)
                    chi.append(((QString)allData[count]).toDouble());
                double waveNumber = ((QString)allData[count++]).toDouble();
                int avalIndex = ((QString)allData[count++]).toInt();
                double aval = ((QString)allData[count++]).toDouble();

                QList<complex<double> > avects;
                int avectcount = ((QString)allData[count++]).toInt();
                if(allData.count()!= avectcount*2+5) return -1;
                for(int k = 0; k < avectcount; k++)
                {
                    complex<double> iAvect;
                    iAvect.real(((QString)allData[count++]).toDouble());
                    iAvect.imag(((QString)allData[count++]).toDouble());
                    avects.append(iAvect);
                }
                displacements[i].append(AvalAvect(waveNumber, avalIndex, aval, avects));
            }
        }
    }
    return 0;
}

void IOThread::run()
{
    if(save)
        SaveExample();
    else LoadExample();
    emit(finished());
}


double Abs(double val)
{
    return (val<0)?-1*val:val;
}
QList<QList<int> > ReadNeighbors(QString path)
{
    QList<QList<int> > neigh;
    QFile file(path);
    if(!file.open(QFile::ReadOnly))
    {
        return neigh;
    }

    QTextStream out(&file);
    QString line = out.readLine();
    QStringList list = line.split(' ');
    QString rw = list[0];
    QString cl = list[1];
    int rows = rw.toInt();
    int cols = cl.toInt();

    line = out.readLine();
    list = line.split(' ');
    for(int r = 0; r < rows; r++)
    {
        int c = 0;
        neigh.append(QList<int>());
        while(c < cols)
        {
            for(int i = 0; i < list.count(); i++)
            {
                if(list.at(i) == " " || list.at(i) == "") continue;
                neigh[r].append(((QString)list.at(i)).toInt());
                c++;
            }
            line = out.readLine();
            list = line.split(' ');
        }
    }

    return neigh;
}
void WriteVector(QString path, int *x, int size)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Vector");
        return;
    }

    QTextStream out(&file);
    out << size << " (rows*cols)"<< "\n";
    for(int i = 0; i < size; i++)
    {
            out << x[i]<<"\n";
    }
    file.close();
}
void WriteVector(QString path, QList<int>elems)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Vector");
        return;
    }

    QTextStream out(&file);
    for(int i = 0; i < elems.size(); i++)
    {
        if(elems.at(i) != -1)
            out << elems.at(i)+1<<"\n";
        else out << elems.at(i)<<"\n";
    }
    file.close();
}
void WriteVector(QString path, QList<double>elems)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Vector");
        return;
    }

    QTextStream out(&file);
    for(int i = 0; i < elems.size(); i++)
    {
        out << elems.at(i)<<"\n";
    }
    file.close();
}
void SaveVertices(QList<double> x, QList<double> y);
FMatrix DecompactMatrix(CMatrix A, QString pathToSave)
{
    FMatrix result(A.rows, A.rows);
    for(int i = 0; i < A.nrows; i++)
    {
        result.SetElem(3*i,3*i,A.elems[3*i][0]);
        result.SetElem(3*i,3*i+1,A.elems[3*i][1]);
        result.SetElem(3*i,3*i+2,A.elems[3*i][2]);

        result.SetElem(3*i+1,3*i,A.elems[3*i+1][0]);
        result.SetElem(3*i+1,3*i+1,A.elems[3*i+1][1]);
        result.SetElem(3*i+1,3*i+2,A.elems[3*i+1][2]);

        result.SetElem(3*i+2,3*i,A.elems[3*i+2][0]);
        result.SetElem(3*i+2,3*i+1,A.elems[3*i+2][1]);
        result.SetElem(3*i+2,3*i+2,A.elems[3*i+2][2]);


        QList<int> neigh = A.GetElems(i);
        for(int j = 0; j < neigh.count(); j++)
        {
            int ind = neigh[j];
            if(ind == -1) continue;

            result.SetElem(3*i, 3*ind, A.elems[3*i][3*(j+1)]);
            result.SetElem(3*i, 3*ind+1, A.elems[3*i][3*(j+1)+1]);
            result.SetElem(3*i, 3*ind+2, A.elems[3*i][3*(j+1)+2]);

            result.SetElem(3*i+1, 3*ind, A.elems[3*i+1][3*(j+1)]);
            result.SetElem(3*i+1, 3*ind+1, A.elems[3*i+1][3*(j+1)+1]);
            result.SetElem(3*i+1, 3*ind+2, A.elems[3*i+1][3*(j+1)+2]);

            result.SetElem(3*i+2, 3*ind, A.elems[3*i+2][3*(j+1)]);
            result.SetElem(3*i+2, 3*ind+1, A.elems[3*i+2][3*(j+1)+1]);
            result.SetElem(3*i+2, 3*ind+2, A.elems[3*i+2][3*(j+1)+2]);
        }
    }
    result.WriteMatrix(pathToSave);
    return result;
}
FMatrix* ReadCMatrix(int *rows, int *cols, QString path)
{
    QFile file(path);
    if(!file.open(QFile::ReadOnly))
        return new FMatrix(*rows, *cols);
    QTextStream out(&file);
    QString line = out.readLine();
    QStringList list = line.split(' ');
    QString rw = list[0];
    QString cl = list[1];
    rows[0] = rw.toInt();
    cols[0] = cl.toInt();
    FMatrix* result = new FMatrix(*rows, *cols);

    line = out.readLine();
    list = line.split(' ');
    for(int r = 0; r < *rows; r++)
    {
        int c = 0;
        for(int i = 0; i < list.count(); i++)
        {
            if(list.at(i) == " " || list.at(i) == "") continue;
            double d = ((QString)list.at(i)).toDouble();
            result->SetElem(r,c++,d);
        }
        line = out.readLine();
        list = line.split(' ');
    }

    file.close();
    return result;
}
void SaveElems(QList<QList<int> > triangles, QString path)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Elems");
        return;
    }

    QTextStream out(&file);
    for(int i = 0; i < triangles.count(); i++)
    {
        for(int j = 0; j < triangles[i].count(); j++)
            out << triangles[i][j]+1<<" ";
        out << "\n";
    }
    file.close();
}
double Distance(double x1, double y1, double x2, double y2)
{
    return sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
}
void WriteVector(QString path, complex<double> *x, int size)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Vector");
        return;
    }

    QTextStream out(&file);
    for(int i = 0; i < size; i++)
    {
        out << x[i].real();
        if(x[i].imag() < 0)
            out <<" -";
        else out <<" +";
        out << abs(x[i].imag()) <<"*i\n";
    }
    file.close();
}
void WriteVector(QString path, int *x, int vertexNu, int cols)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Vector");
        return;
    }

    QTextStream out(&file);
    for(int i = 0; i < vertexNu/3; i++)
    {
        for(int j =0; j < cols; j++)
            out << x[j*vertexNu/3+i] ;
        out<< "\n";
    }
    file.close();
}
QList<double> ToDoubleList(double *S, int sCount)
{
    QList<double> result;
    for(int i = 0; i < sCount; i++)
        result.push_back(S[i]);
    return result;
}
void Copy(complex<double> *from, complex<double> *to, int size)
{
    for(int i = 0; i < size; i++)
    {
        to[i].real(from[i].real());
        to[i].imag(from[i].imag());
    }
}
void SaveVertices(QList<double> x, QList<double> y, QString path)
{
    QFile file(path);
    if(!file.open(QFile::WriteOnly))
    {
        printf("Fail To Save Vertices");
        return;
    }

    QTextStream out(&file);
    for(int i = 0; i < x.count(); i++)
        out<< x[i] << " " << y[i] << "\n";
    file.close();
}
QList<complex<double> > ToComplexList(complex<double> *array, int size)
{
    QList<complex<double> > result;
    for(int i = 0; i < size; i++)
        result.push_back(array[i]);

    return result;
}
FComplexMatrix* ReadCMatrix(int *rows, int *cols, QString pathR, QString pathI)
{
    QFile fileR(pathR);
    QFile fileI(pathI);
    if(!fileR.open(QFile::ReadOnly) || !fileI.open(QFile::ReadOnly))
        return new FComplexMatrix(*rows, *cols);
    QTextStream outR(&fileR);
    QString lineR = outR.readLine();
    QStringList listR = lineR.split(' ');
    QTextStream outI(&fileI);
    QString lineI = outI.readLine();
    QStringList listI = lineI.split(' ');
    QString rw = listR[0];
    QString cl = listR[1];
    *rows = rw.toInt();
    *cols = cl.toInt();
    FComplexMatrix* result = new FComplexMatrix(*rows, *cols);

    lineR = outR.readLine();
    listR = lineR.split(' ');
    lineI = outI.readLine();
    listI = lineI.split(' ');
    for(int r = 0; r < *rows; r++)
    {
        int c = 0;
        int indR = 0;
        int indI = 0;
        while(indR < listR.count() && indI < listI.count())
        {
            while(indR < listR.count() && (listR.at(indR) == " " || listR.at(indR) == ""))
            {indR++;}
            while(indI < listI.count() && (listI.at(indI) == " " || listI.at(indI) == ""))
            {indI++;}
            if(indR >= listR.count() || indI >= listI.count())
                break;

            double dR = ((QString)listR.at(indR++)).toDouble();
            double dI = ((QString)listI.at(indI++)).toDouble();
            complex<double> d;
            d.real(dR);
            d.imag(dI);
            result->SetElem(r,c++,d);
        }
        lineR = outR.readLine();
        listR = lineR.split(' ');

        lineI = outI.readLine();
        listI = lineI.split(' ');
    }

    fileR.close();
    fileI.close();
    return result;
}
void ScaleCoordinates(int minimumDifference, QList<double> &x, QList<double> &y)
{
    for(int i = 0; i < x.count(); i++)
    {
        double ix = x.at(i);
        double iy = y.at(i);
        x[i] = ix/minimumDifference;
        y[i] = iy/minimumDifference;
    }
}
FMatrix DecompactMatrix(QString pathToDecompact, QString neighPath, QString pathToSave)
{
    int *rows = new int[1];
    int *cols = new int[1];
    FMatrix *A = ReadCMatrix(rows, cols, pathToDecompact);
    QList<QList<int> > neigh = ReadNeighbors(neighPath);
    FMatrix result(A->rows, A->rows);

    for(int i = 0; i < neigh.count(); i++)
    {
        result.SetElem(3*i,3*i,A->GetElem(3*i,0));
        result.SetElem(3*i,3*i+1,A->GetElem(3*i,1));
        result.SetElem(3*i,3*i+2,A->GetElem(3*i,2));

        result.SetElem(3*i+1,3*i,A->GetElem(3*i+1,0));
        result.SetElem(3*i+1,3*i+1,A->GetElem(3*i+1,1));
        result.SetElem(3*i+1,3*i+2,A->GetElem(3*i+1,2));

        result.SetElem(3*i+2,3*i,A->GetElem(3*i+2,0));
        result.SetElem(3*i+2,3*i+1,A->GetElem(3*i+2,1));
        result.SetElem(3*i+2,3*i+2,A->GetElem(3*i+2,2));
        for(int j = 0; j < neigh[0].count(); j++)
        {
            int ind = neigh[i][j]-1;
            //Se le resta 1 porque el fichero guarda
            //los indices para fortran, es decir cpmenzando en 1 y no en 0.
            //Cuando ind = -2 indica que no tiene mas vecinos (-1-1=-2)
            if(ind == -2) continue;

            result.SetElem(3*i, 3*ind, A->GetElem(3*i, 3*(j+1)));
            result.SetElem(3*i, 3*ind+1, A->GetElem(3*i, 3*(j+1)+1));
            result.SetElem(3*i, 3*ind+2, A->GetElem(3*i, 3*(j+1)+2));

            result.SetElem(3*i+1, 3*ind, A->GetElem(3*i+1, 3*(j+1)));
            result.SetElem(3*i+1, 3*ind+1, A->GetElem(3*i+1, 3*(j+1)+1));
            result.SetElem(3*i+1, 3*ind+2, A->GetElem(3*i+1, 3*(j+1)+2));

            result.SetElem(3*i+2, 3*ind, A->GetElem(3*i+2, 3*(j+1)));
            result.SetElem(3*i+2, 3*ind+1, A->GetElem(3*i+2, 3*(j+1)+1));
            result.SetElem(3*i+2, 3*ind+2, A->GetElem(3*i+2, 3*(j+1)+2));
        }
    }
    result.WriteMatrix(pathToSave);
    return result;
}
double GetMinimumDifference(QList<QList<int> > triangles, QList<double> x, QList<double> y)
{
    //double.maxValue
    double minVal = DBL_MAX;
    for(int i = 0; i < triangles.count(); i++)
    {
        minVal = min(minVal, min(
                                 Distance(x[triangles[i][2]], x[triangles[i][2]], x[triangles[i][0]], y[triangles[i][0]]), min(
                                                                                                                                 Distance(x[triangles[i][0]], y[triangles[i][0]], x[triangles[i][1]], y[triangles[i][1]]), Distance(x[triangles[i][1]], y[triangles[i][1]], x[triangles[i][2]], y[triangles[i][2]])
                                                                                                                               )
                                )
                     );
    }
    return minVal;
}
FComplexMatrix DecompactCMatrix(QString pathToDecompactR, QString pathToDecompactI, QString neighPath, QString pathToSaveR, QString pathToSaveI)
{
    int *rows = new int[1];
    int *cols = new int[1];
    FComplexMatrix *A = ReadCMatrix(rows,cols, pathToDecompactR, pathToDecompactI);
    QList<QList<int> > neigh = ReadNeighbors(neighPath);
    FComplexMatrix result(A->rows, A->rows);

    for(int i = 0; i < neigh.count(); i++)
    {
        result.SetElem(3*i,3*i,A->GetElem(3*i,0));
        result.SetElem(3*i,3*i+1,A->GetElem(3*i,1));
        result.SetElem(3*i,3*i+2,A->GetElem(3*i,2));

        result.SetElem(3*i+1,3*i,A->GetElem(3*i+1,0));
        result.SetElem(3*i+1,3*i+1,A->GetElem(3*i+1,1));
        result.SetElem(3*i+1,3*i+2,A->GetElem(3*i+1,2));

        result.SetElem(3*i+2,3*i,A->GetElem(3*i+2,0));
        result.SetElem(3*i+2,3*i+1,A->GetElem(3*i+2,1));
        result.SetElem(3*i+2,3*i+2,A->GetElem(3*i+2,2));

        for(int j = 0; j < neigh[0].count(); j++)
        {
            int ind = neigh[i][j]-1;
            if(ind == -2) continue;

            result.SetElem(3*i, 3*ind, A->GetElem(3*i, 3*(j+1)));
            result.SetElem(3*i, 3*ind+1, A->GetElem(3*i, 3*(j+1)+1));
            result.SetElem(3*i, 3*ind+2, A->GetElem(3*i, 3*(j+1)+2));

            result.SetElem(3*i+1, 3*ind, A->GetElem(3*i+1, 3*(j+1)));
            result.SetElem(3*i+1, 3*ind+1, A->GetElem(3*i+1, 3*(j+1)+1));
            result.SetElem(3*i+1, 3*ind+2, A->GetElem(3*i+1, 3*(j+1)+2));

            result.SetElem(3*i+2, 3*ind, A->GetElem(3*i+2, 3*(j+1)));
            result.SetElem(3*i+2, 3*ind+1, A->GetElem(3*i+2, 3*(j+1)+1));
            result.SetElem(3*i+2, 3*ind+2, A->GetElem(3*i+2, 3*(j+1)+2));
        }
    }
    result.WriteRealComplexMatrix(pathToSaveR);
    result.WriteImagComplexMatrix(pathToSaveI);

    return result;
}
