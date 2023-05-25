#include "move.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Move w;
    w.show();
    return a.exec();
}
