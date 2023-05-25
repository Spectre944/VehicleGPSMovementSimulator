#ifndef COORD_WORK_H
#define COORD_WORK_H


#include <QObject>
#include <QString>
#include <QList>
#include <QPair>

#include <qmath.h>

#include <algorithm>
#include <vector>
#include <QVariantList>

#include <math.h>



#define ID_CAR          "1"
#define NUCLEAR_RADIUS  3000

#define SEC_IN_DAY      86400
#define MAX_TIME        1e12

#define R  6378137   // ~ Радиус Земли

#define DEGREE   "° "
#define MINUTE   "' "
#define SECOND   "'' "


class Point {

public:

    Point() {
        x = 0.0;
        y = 0.0;
    }

    Point(double _x, double _y) {
        x = _x;
        y = _y;
    }

    bool operator < (const Point &p) const {
        return x < p.x || (x == p.x && y < p.y);
    }

    bool operator == (const Point &p) const {

        return (x == p.x && y == p.y) ;
    }

    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    void setXY(double _x, double _y) {
        x = _x;
        y = _y;
    }

    const  double getX() const {
        return static_cast<const double>(x);
    }

    const double  getY() const {
        return static_cast<const double>(y);
    }


private:

    double x, y; // x - широта , y - долгота выбранной точки
};


class Work_Point {

public:

    Work_Point() {
        this->clearAll();
    }

    Work_Point(double latid, double longtid, int _move) {
        set_params(latid, longtid, _move);
    }

    void set_move(int value) {  // 1 – машина движется, 0 – машина стоит
        move = value;
    }

    int get_move() {
        return move;
    }

    double get_latid() {
        return point.getX();
    }

    double get_longtid() {
        return point.getY();
    }

    void set_params(double latid, double longtid, int _move) {
        point.setXY(latid, longtid);
        move = _move;
    }

    void save_old_coords() {
        old_point.setXY(point.getX(), point.getY());
    }

    void set_save_point() {
        save_point.setXY(point.getX(), point.getY());
    }

    Point get_point() {
        return point;
    }

    Point get_old_point() {
        return old_point;
    }

    Point get_save_point() {
        return save_point;
    }

    void set_distance(int _distance) {
        save_distance = _distance;
    }

    int get_distance() {
        return save_distance;
    }


    void set_time(qint64 _time) {
        save_time = _time;
    }

    qint64 get_time() {
        return save_time;
    }

    void set_interval_time(qint64 _time) {
        interval_time = _time;
    }

    qint64 get_interval_time() {
        return interval_time;
    }

    void clearAll() {
        point.setXY(0.0, 0.0);
        old_point.setXY(0.0, 0.0);
        save_point.setXY(0.0, 0.0);
        save_distance = 0;
        save_time = 0;
        interval_time = 0;
        move = 0;
        localContourRoute.clear();
        maxTermReconn = 0;
        start_reconn_time = 0;
    }

    void set_localContourRoute(std::vector<Point> vector) {
        localContourRoute.clear();

        foreach(Point p, vector)
            localContourRoute.push_back(p);
    }

    std::vector<Point> get_localContourRoute() {
        return localContourRoute;
    }

    int erase_elem(Point p) {

        if(localContourRoute.size()) {

            for (auto it = localContourRoute.begin(); it != localContourRoute.end(); ++it) {
                if (*it == p) {
                    localContourRoute.erase(it);
                    return 0;
                }
            }
        }

        if(!localContourRoute.size())
            return 1;

        return 1;

    }

    void set_maxTermReconn(qint64 term){
        maxTermReconn = term;
    }

    qint64 get_maxTermReconn() {
        return maxTermReconn;
    }

    void set_start_reconn_time(qint64 _time) {
        start_reconn_time = _time;
    }

    qint64 get_start_reconn_time() {
        return start_reconn_time;
    }





private:

    Point point;

    Point old_point;

    Point save_point;

    int save_distance;

    qint64 interval_time;

    qint64 save_time;

    int move;

    qint64 maxTermReconn;

    qint64 start_reconn_time;

    std::vector<Point> localContourRoute;
};




////////////////////////////////////////////////////////////////////////////////////////////


class Coord_work {

public:

    /*
                 ФУНКЦИЯ ОПРЕДЕЛЕНИЯ ВЫПУКЛОГО КОНТУРА ОБЛАСТИ
    Аргумент ф-ции : vector<Point>   (Point - структура, описывающая координаты точки)
    Результат работы - vector<Point> набор точек, описывающий выпуклую область

    */
    std::vector<Point> static convex_hull(std::vector<Point>);  // Coord_work::convex_hull();


    /*
                 ФУНКЦИЯ ОПРЕДЕЛЕНИЯ ПЛОЩАДИ ОБЛАСТИ
    Аргумент ф-ции : vector<Point>   (Point - структура, описывающая координаты точки)
    Результат работы - площадь области, описанной точками вектора Р (в кв. км)
    */
    double static calcPolygonArea (std::vector<Point>);


    /*
            ФУНКЦИЯ ПРЕОБРАЗОВАНИЯ НАБОРА ПАР КООРДИНАТ В "УПАКОВАННУЮ" СТРОКУ
    Аргумент ф-ции : vector<Point>   (Point - структура, описывающая координаты точки)
    Результат работы - "упакованная" строка для записи в БД
    */
    QString static transPairCoordString(std::vector<Point>);

    /*
        ФУНКЦИЯ ПРЕОБРАЗОВАНИЯ "УПАКОВАННОЙ" СТРОКИ В НАБОР ПАР КООРДИНАТ
    Аргумент ф-ции : "упакованная" строка из БД
    Результат работы - vector<Point>   (Point - структура, описывающая координаты точки)
    */
    std::vector<Point> static transStringPairCoord(QString);

           // Coord_work::transStringPairCoord(str);
    /*
        ФУНКЦИЯ ОПРЕДЕЛЕНИЯ ВИРТУАЛЬНОГО ЦЕНТРА ОБЛАСТИ
    Аргумент ф-ции : vector<Point>   (Point - структура, описывающая координаты точки)
    Результат работы - QPair<latid, longtid>
    */
    QPair<double, double> static calcVirtualCentreArea (std::vector<Point>);


    /*
        ФУНКЦИЯ ОПРЕДЕЛЕНИЯ ВХОЖДЕНИЯ ТОЧКИ В ОБЛАСТЬ
    Аргумент ф-ции :  Point точка, vector<Point> P - ветор Р описывает область
    Результат работы - true  ,  false
    */
    bool static pointIntoArea(Point , std::vector<Point> );

    /*
        ФУНКЦИЯ РАСЧЕТА РАССТОЯНИЯ МЕЖДУ 2 ТОЧКАМИ
    Аргументы ф-ции : координаты 2 точек
    Результат работы - расстояние между точками в м
    */
    double static dist_2_Points(Point , Point);


    /*
       ФУНКЦИЯ ОПРЕДЕЛЕНИЯ КООРД. ТОЧКИ ПО КООРД, АЗИМУТУ И РАССТОЯНИЮ
    Аргумент ф-ции :  Point точка, double азимут , double расстояние
    Результат работы - double координата
    */
    Point static pointByCoordAzimDist(Point , double , double );


    /*
        ФУНКЦИЯ ФОРМИРОВАНИЯ КРУГА
    Аргумент ф-ции :  Point точка - центр, радиус круга в км
    Результат работы - вектор с точками круга
    */
    std::vector<Point> static circleFormation(Point , double); // Coord_work::circleFormation(Point , double);

    /*
        ФУНКЦИЯ ФОРМИРОВАНИЯ КРУГА
    Аргумент ф-ции :  Point точка - центр, радиус круга в км, нач. угол, конечный угол (в градусах)
    Результат работы - вектор с точками круга
    */
    std::vector<Point> static partCircleFormation(Point , double, double, double);


    /*
        ФУНКЦИЯ ФОРМИРОВАНИЯ ТРЕУГОЛЬНИКА
    Аргумент ф-ции :  Point точка - центр, R1, R2, wind_dir, wind_speed
    Результат работы - вектор с точками круга
    */
    std::vector<Point> static triangleFormation(Point , double R1, double R2, double wind_dir);


    /*
        ФУНКЦИЯ РАСЧЕТА УПРОЩЕННОЙ РАДИАЦИОННОЙ ОБЛАСТИ  QList<QVariantList>
    Аргумент ф-ции :  центральная точка (Point), радиус R1 , радиус R2,   угол ветра
    Результат работы - контур и площадь зоны long-term
    */
    QPair<QString, double> static nuclearSimpleAreaCalc(Point, double, double,  double);


    /*
        ФУНКЦИЯ РАСЧЕТА УПРОЩЕННОЙ ХИМИЧЕСКОЙ ОБЛАСТИ  QList<QVariantList> nuclearlist
    Аргумент ф-ции :  список точек QList<QVariantList> с химической угрозой
    Результат работы - QPair<QString, double> - пара контур, площадь зоны
    */
    QPair<QString, double> static chemSimpleAreaCalc(QList<QVariantList> &);

    /*
        ФУНКЦИЯ РАСЧЕТА УПРОЩЕННОЙ БИОЛОГИЧЕСКОЙ ОБЛАСТИ  QList<QVariantList> nuclearlist
    Аргумент ф-ции :  список точек QList<QVariantList> с химической угрозой
    Результат работы - QPair<QString, double> - пара контур, площадь
    */
    QPair<QString, double> static biolSimpleAreaCalc(QList<QVariantList> &);


    double static Angle2D(double , double , double , double );

private :

    double static cross(const Point &, const Point &, const Point &);

};


inline std::vector<Point> Coord_work::convex_hull(std::vector<Point> P)
{
    size_t n = P.size(),
           k = 0;

    std::vector<Point> Empty;
    Empty.clear();

    if (n <= 3)
        return Empty;

    std::vector<Point> H(2 * n);

    sort(P.begin(), P.end());  // Сортировать точки лексикографически

    for (size_t i = 0; i < n; ++i) {    // Build lower hull
        while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
        H[k++] = P[i];
    }

    for (size_t i = n-1, t = k+1; i > 0; --i) {   // Build upper hull
        while (k >= t && cross(H[k-2], H[k-1], P[i-1]) <= 0) k--;
        H[k++] = P[i-1];
    }

    H.resize(k-1);

    return H;
}

inline double Coord_work::calcPolygonArea(std::vector<Point> P)
{
    double area = 0;

    Point p = P.at(0); // Первую точку вставляем в конец вектора
    P.push_back(p);

    for(size_t i = 0; i < P.size() - 1; ++i) {
        Point p1 = P.at(i);
        Point p2 = P.at(i + 1);
        area += qDegreesToRadians(p2.getY() - p1.getY()) * (2 + qSin(qDegreesToRadians(p1.getX()))+ qSin(qDegreesToRadians(p2.getX())));
    }

    return (area * R * R / 2 / 1000000);
}

inline QString Coord_work::transPairCoordString(std::vector<Point> vector)
{
    QString result;
    result.clear();

    foreach(Point p, vector) {
        result += QString::number(p.getX()) + "-" + QString::number(p.getY());
        result += "+";
    }

    return result.remove(result.size() - 1, 1);
}

inline std::vector<Point> Coord_work::transStringPairCoord(QString str)
{
    std::vector<Point> vector;
    vector.clear();

    QStringList outerList = str.split("+"),
                innerList;

    foreach(QString help, outerList) {
        innerList = help.split("-");
        vector.push_back(Point(innerList.at(0).trimmed().toDouble(), innerList.at(1).trimmed().toDouble()));
    }

    return vector;
}

inline QPair<double, double> Coord_work::calcVirtualCentreArea(std::vector<Point> vector)
{
    double latid = 0.0, longtid = 0.0;
    double size = static_cast<double>(vector.size());

    foreach(Point point, vector) {
        latid += point.getX();
        longtid += point.getY();
    }

    latid /= size;
    longtid /= size;

    return qMakePair(latid, longtid);
}

inline bool Coord_work::pointIntoArea(Point p, std::vector<Point> P)
{
    double angle = 0,
           point1_lat,
           point1_long,
           point2_lat,
           point2_long,
           latitude = p.getX(),
            longitude = p.getY();

    size_t n = P.size();

    for (int i = 0; i < n; ++i) {
        point1_lat  = P.at(i).getX() - latitude;
        point1_long = P.at(i).getY() - longitude;
        point2_lat  = P.at((i +1) % n).getX() - latitude;
        point2_long = P.at((i+1) % n).getY() - longitude;
        angle += Angle2D(point1_lat, point1_long, point2_lat, point2_long);
    }

    return (qFabs(angle) > M_PI);

}

inline double Coord_work::dist_2_Points(Point point1, Point point2)
{
    double F1, F2, dF, dL, a, c, res;

    F1 = (point1.getX()) * M_PI/180;
    F2 = (point2.getX()) * M_PI/180;

    dF = (point2.getX() - point1.getX()) * M_PI/180;
    dL = (point2.getY()  - point1.getY()) * M_PI/180;

    a =  qSin(dF/2) * qSin(dF/2) + qCos(F1)   * qCos(F2) * qSin(dL/2) * qSin(dL/2);
    c =  2 * qAtan2(qSqrt(a), qSqrt(1-a));
    res =  R * c;

    return fabs(res);
}

inline Point Coord_work::pointByCoordAzimDist(Point point, double azimut, double distance) // azimut   - в градусах   distance - в км
{
    double lat1 = (point.getX()) * M_PI / 180.0;
    double lon1 = (point.getY()) * M_PI / 180.0;
    double brng = azimut  * M_PI / 180.0;
    double dist = distance   * 1e3;


    double lat2 = qAsin( qSin(lat1) * qCos(dist/R) + qCos(lat1) * qSin(dist/R) * qCos(brng) );

    double lon2 = lon1 + qAtan2( qSin(brng) * qSin(dist/R) * qCos(lat1),
                                 qCos(dist/R) - qSin(lat1) * qSin(lat2));

    lat2 *= 180.0 / M_PI;
    lon2 *= 180.0 / M_PI;

    Point p(lat2, lon2);  // Было

    return p;
}

inline std::vector<Point> Coord_work::circleFormation(Point centralPoint, double radius)
{
    std::vector<Point> circleVec;
    circleVec.clear();
    Point p;

    for(int angle = 0; angle < 360; angle += 10) {
        p = pointByCoordAzimDist(centralPoint, static_cast<double>(angle), radius);
        circleVec.push_back(p);
    }

    return circleVec;
}

inline std::vector<Point> Coord_work::partCircleFormation(Point centralPoint, double radius, double startAngle, double endAngle)
{
    //  endAngle > startAngle - заполнение по часовой
    //  endAngle < startAngle - заполнение против часовой

    std::vector<Point> circleVec;
    circleVec.clear();
    Point p;

    int start  = static_cast<int>(startAngle) ,
        finish = static_cast<int>(endAngle);

    if(endAngle == startAngle)
        return circleVec;

    else if(endAngle > startAngle) {
        for(int angle = start; angle <= finish; angle += 10) {
            p = pointByCoordAzimDist(centralPoint, static_cast<double>(angle), radius);
            circleVec.push_back(p);
        }
    }

    else if(endAngle < startAngle) {
        for(int angle = start; angle >= finish; angle -= 10) {
            p = pointByCoordAzimDist(centralPoint, static_cast<double>(angle), radius);
            circleVec.push_back(p);
        }
    }

    return circleVec;
}

inline std::vector<Point> Coord_work::triangleFormation(Point basic, double R1, double R2, double wind_dir)
{
    std::vector<Point> circleVec;
    circleVec.clear();
    Point A, B, C, D;

    C = pointByCoordAzimDist(basic, wind_dir, R1);
    A = pointByCoordAzimDist(basic, wind_dir + 180, 2 * R2);
    B = pointByCoordAzimDist(A, wind_dir - 30, (R1  + R2 * 2)/0.866 );
    D = pointByCoordAzimDist(A, wind_dir + 30, (R1  + R2 * 2)/0.866 );

    circleVec.push_back(A);
    circleVec.push_back(B);
    circleVec.push_back(C);
    circleVec.push_back(D);

    return circleVec;

}

inline QPair<QString, double> Coord_work::nuclearSimpleAreaCalc(Point centre, double R1 ,double R2, double wind_direct)
{
    /*
        1. Заполнить полукруг
        2. Найти отдаленную по напр ветра точку
        3. Найти точки +- 90 град
        4. Заполнить вектор контура
        5. Найти зону - строку
        6. Найти зону площадь
        7. Вернуть пару
    */

    std::vector<Point> zoneVector;
    zoneVector.clear();
//    std::vector<Point> resultVector;
//    resultVector.clear();

    QString str_zone;
    double square_zone;

    Point apexTriangle  = pointByCoordAzimDist(centre , wind_direct + 180 , R2 * 2);
    Point leftTriangle  = pointByCoordAzimDist(apexTriangle , wind_direct - 30 , (R1  + R2 * 2) / 0.866 );
    Point baseTriangle  = pointByCoordAzimDist(centre , wind_direct , R1);
    Point rightTriangle = pointByCoordAzimDist(apexTriangle , wind_direct + 30 , (R1  + R2 * 2) / 0.866 );

    zoneVector.push_back(apexTriangle);
    zoneVector.push_back(leftTriangle);
    zoneVector.push_back(baseTriangle);
    zoneVector.push_back(rightTriangle);

//    resultVector = convex_hull(zoneVector);
    str_zone     = transPairCoordString(zoneVector);
    square_zone  = calcPolygonArea(zoneVector);

    return qMakePair(str_zone, square_zone);

}



inline QPair<QString, double> Coord_work::chemSimpleAreaCalc(QList<QVariantList> &mainList)
{
    Q_UNUSED(mainList);

    QPair<QString, double> pair;


    return pair;
}

inline QPair<QString, double> Coord_work::biolSimpleAreaCalc(QList<QVariantList> &mainList)
{
     Q_UNUSED(mainList);

    QPair<QString, double> pair;


    return pair;
}

inline double Coord_work::Angle2D(double y1, double x1, double y2, double x2)
{
    double dtheta, theta1, theta2;
    theta1 = qAtan2(y1,x1);
    theta2 = qAtan2(y2,x2);
    dtheta = theta2 - theta1;
    while (dtheta > M_PI)
        dtheta -= 2 * M_PI;
    while (dtheta < -M_PI)
        dtheta += 2 * M_PI;

    return  dtheta;
}


inline double Coord_work::cross(const Point &O, const Point &A, const Point &B)
{
    return ((A.getX() - O.getX()) * (B.getY() - O.getY())) - ((A.getY() - O.getY()) * (B.getX() - O.getX()));
}










































#endif // COORD_WORK_H
