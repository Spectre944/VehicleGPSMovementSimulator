#include "inc/move.h"
#include "ui_move.h"

Move::Move(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::Move)
{
    ui->setupUi(this);

    timer = new QTimer();
    timer->setInterval(TIM_INTER);

    connect(timer, &QTimer::timeout, this, QOverload<>::of(&Move::movementProc));
    timer->start();

    //Connect QMenu actions
    connect(ui->actionSavePath, &QAction::triggered, this, QOverload<>::of(&Move::savePath));
    connect(ui->actionWriteFirstCoord, &QAction::triggered, this, QOverload<>::of(&Move::setFirstCoord));

    //stop start timer in QMenu
    connect(ui->actionStopTimer, &QAction::triggered, timer, &QTimer::stop);
    connect(ui->actionStartTimer, &QAction::triggered, this, QOverload<>::of(&Move::startTimer));

    currCoordinate = Point(49.8081, 24.0189);
    directionCompas = 0;
    coordList.clear();

    speed = 0;
    directionWheel = 0;

    distancePassed = 0;
    pointsPassed = 0;

}

Move::~Move()
{
    delete ui;
}

void Move::displayData()
{

    //clear all
    ui->labelSpeedNeutral->setStyleSheet("");

    ui->labelSpeed1->setStyleSheet("");
    ui->labelSpeed2->setStyleSheet("");
    ui->labelSpeed3->setStyleSheet("");
    ui->labelSpeed4->setStyleSheet("");

    ui->labelSpeedRev1->setStyleSheet("");
    ui->labelSpeedRev2->setStyleSheet("");

    ui->lineEditCoord->setText(QString::number(currCoordinate.getX()) + " " + QString::number(currCoordinate.getY()));
    ui->labelSpeedKPH->setText(QString::number(speed*speedMultiplier*3.6));
    ui->labelMoveDir->setText(QString::number(directionCompas));

    roatateArrow();

    ui->labelDistancePassed->setText(QString::number(distancePassed) + " м");
    ui->labelPointsPassed->setText(QString::number(pointsPassed) + " шт");


    switch (speed) {

    case 0:

        ui->labelSpeedNeutral->setStyleSheet("background-color: rgba(255, 247, 24, 150);\n	border: 1px solid black;\n	border-radius: 3px;\n");

        break;

    case 1:

        ui->labelSpeed1->setStyleSheet("background-color: rgba(0, 255, 127, 140);\nborder: 1px solid black;\nborder-radius: 3px;\n");

        break;

    case 2:

        ui->labelSpeed2->setStyleSheet("background-color: rgba(0, 255, 127, 140);\nborder: 1px solid black;\nborder-radius: 3px;\n");

        break;

    case 3:

        ui->labelSpeed3->setStyleSheet("background-color: rgba(0, 255, 127, 140);\nborder: 1px solid black;\nborder-radius: 3px;\n");

        break;

    case 4:

        ui->labelSpeed4->setStyleSheet("background-color: rgba(0, 170, 255, 200);\nborder: 1px solid black;\nborder-radius: 3px;\n");

        break;

    case -1:

        ui->labelSpeedRev1->setStyleSheet("background-color: rgba(255, 141, 26, 170);\nborder: 1px solid black;\nborder-radius: 3px;\n");

        break;

    case -2:

        ui->labelSpeedRev2->setStyleSheet("background-color: rgba(255, 141, 26, 170);\nborder: 1px solid black;\nborder-radius: 3px;\n");

        break;
    default:
        break;
    }

}

void Move::movementProc()
{

    //update position on front of the car comparing to North
    //also imitating wheel rotation if going backward invere changing of compass
    directrionProc();

    //Calculate new position depending of speed and direction
    currCoordinate = calculateNewPos(currCoordinate, speed, directionCompas);

    //update labels
    displayData();

    //dont append points if we stop
    if(speed == 0) return;

    //save array of waypoints
    coordList.append(currCoordinate);

    //write ammout of points
    pointsPassed++;

}

void Move::directrionProc()
{
    //update position on front of the car comparing to North
    //also imitating wheel rotation if going backward invere changing of compass
    if (speed == 0) return;

    if(speed < 0)
        directionCompas -= directionWheel;
    else
        directionCompas += directionWheel;

    directionCompas > 360 ? directionCompas -= 360 : 0;

    directionCompas < 0 ? directionCompas += 360 : 0;
}

void Move::roatateArrow()
{

    QPixmap arrow(":/ico/images/navigation-arrow.png");

    QPixmap rotatedPixmap(arrow.size());
    rotatedPixmap.fill(Qt::transparent);

    QPainter painter(&rotatedPixmap);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::SmoothPixmapTransform);

    // Translate the painter to the center of the pixmap
    painter.translate(arrow.width() / 2, arrow.height() / 2);
    // Rotate the painter by the specified angle
    painter.rotate(directionCompas);
    // Translate the painter back to the original position
    painter.translate(-arrow.width() / 2, -arrow.height() / 2);

    // Draw the rotated pixmap onto the rotatedPixmap
    painter.drawPixmap(0, 0, arrow);

    ui->labelMoveArrow->setPixmap(rotatedPixmap);

}

void Move::savePath()
{

    QString json = savePathJSON(coordList);

    QDir dir = QDir::current();
    QString filePath = dir.absoluteFilePath("ImitationPath.json");

    QFile file(filePath);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QTextStream out(&file);
        out << json;
        file.close();
        qDebug() << "JSON saved to file: " << filePath;
    }
    else
    {
        qDebug() << "Failed to open file for writing: " << filePath;
    }

    coordList.clear();
    distancePassed = 0;
    pointsPassed = 0;

}

void Move::setFirstCoord()
{
    //split arrayt by space or comma
    QStringList coord (ui->lineEditCoord->text().split(QRegExp("[\\s,]+")));

    if (coord.count() >= 2){
         currCoordinate.setXY(coord.at(0).toDouble(), coord.at(1).toDouble());
    }
}

QString Move::savePathJSON(const QVector<Point>& coordinates)
{
    QJsonObject featureObject;
    featureObject["type"] = "Feature";
    featureObject["id"] = "sm98711698";

    QJsonObject geometryObject;
    geometryObject["type"] = "LineString";

    QJsonArray coordinateArray;
    for (const auto& coordinate : coordinates)
    {
        QJsonArray pointArray;
        pointArray.append(coordinate.getY());
        pointArray.append(coordinate.getX());
        coordinateArray.append(pointArray);
    }
    geometryObject["coordinates"] = coordinateArray;

    featureObject["geometry"] = geometryObject;
    featureObject["properties"] = QJsonObject();

    QJsonDocument document(featureObject);

    return document.toJson(QJsonDocument::Indented);;
}

Point Move::calculateNewPos(Point pos, int speed, int direction)
{

    if(speed == 0) return pos;

    Point newPos(0, 0);

    // Convert direction from degrees to radians
    double direction_rad = M_PI / 180.0 * direction;

    // Speed in meters per second
    double speed_mps = speed * speedMultiplier;

    // Calculate the distance traveled in 1 second
    distancePassed += abs(speed_mps);

    double x = speed_mps * sin(direction_rad);
    double y = speed_mps * cos(direction_rad);

    // Calculate the new latitude
    double lat = pos.getX() + 180 / M_PI * y / R;
    double lon = pos.getY() + 180 / M_PI / sin(pos.getX()*M_PI/180) * x / R;

    newPos.setXY(lat, lon);

    return newPos;
}


void Move::on_pushButtonMoveForward_clicked()
{
    (speed < 4) ?  speed++ : 0;
    displayData();
}


void Move::on_pushButtonMoveBackward_clicked()
{
    (speed > -2) ?  speed-- : 0;
    displayData();
}


void Move::on_pushButtonMoveLeft_clicked()
{
    directionWheel -= 5;
    ui->horizontalSliderDir->setValue(directionWheel);
}


void Move::on_pushButtonMoveRight_clicked()
{
    directionWheel += 5;
    ui->horizontalSliderDir->setValue(directionWheel);
}


void Move::on_pushButtonResetDir_clicked()
{
    directionWheel = 0;
    ui->horizontalSliderDir->setValue(directionWheel);
}


void Move::on_pushButtonHandBrake_clicked()
{
    speed = 0;
    displayData();
}

