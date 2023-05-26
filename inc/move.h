#ifndef MOVE_H
#define MOVE_H

#include <QMainWindow>
#include <QTimer>
#include <QPainter>

#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>
#include <QJsonArray>

#include <QDir>
#include <QFile>

#include <inc/coord_work.h>


#define TIM_INTER 100

QT_BEGIN_NAMESPACE
namespace Ui { class Move; }
QT_END_NAMESPACE

class Move : public QMainWindow
{
    Q_OBJECT

public:
    Move(QWidget *parent = nullptr);
    ~Move();

    void displayData();
    void movementProc();
    void directrionProc();
    void roatateArrow();
    void savePath();
    void setFirstCoord();
    QString savePathJSON(const QVector<Point>& coordinates);
    Point calculateNewPos(Point, int, int);

    void startTimer() { timer->start(TIM_INTER); }

    //vehicle params
    int speed;
    int directionWheel;

    int distancePassed;
    int pointsPassed;

    double speedMultiplier = 7.5;

    //gps params
    int directionCompas;
    Point currCoordinate;

    QVector<Point> coordList;


private slots:

    void on_pushButtonMoveForward_clicked();
    void on_pushButtonMoveBackward_clicked();
    void on_pushButtonMoveLeft_clicked();
    void on_pushButtonMoveRight_clicked();
    void on_pushButtonResetDir_clicked();
    void on_pushButtonHandBrake_clicked();

private:
    Ui::Move *ui;

    QTimer *timer;
};
#endif // MOVE_H
