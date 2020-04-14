#include "geometry.h"

void warn(std::string message) {
    std::cout << "Не прошли тест: " << message << std::endl;
}

void PointTest() {
    std::cout << "Тест структуры точка" << std::endl;
    Point x(3.0, 4.0);
    Point y(3.0, 4.0);
    Point t(3.0, -4.0);
    Point z(-1.0, 0);

    if (!(x == y) || (x == z)) {
        warn("Не прошел тест ==");
    }

    if (x != y || !(x != z)) {
        warn("Не прошел тест !=");
    }

    if (!(std::abs(x.getLength() - 5.0) < eps &&
        std::abs(t.getLength() - 5.0) < eps))
    {
        warn("Не прошел тест getLength()");
    }

    z.rotate(10);
    z.rotate(-10);
    if (!(std::abs(z.x - (-1.0)) < eps && (std::abs(z.y) < eps))) {
        warn("Не прошел тест rotate() туда сюда");
    }

    z.rotate(3.14159265358979);
    if(z != Point(1,0)) {
        warn("Не прошел тест rotate() на угол пи");
    }

    Point a(1,0);
    Point b(0, 1);
    if(a + b != Point(1,1)) {
        warn("Не прошел тест +");
    }

    if(a - a != Point(0,0)) {
        warn("Не прошел тест -");
    }

    if(0.5*a != Point(0.5, 0)) {
        warn("Не прошел тест double*Point");
    }

    if(a*0.5 != Point(0.5, 0)) {
        warn("Не прошел тест Point*double");
    }

    if(std::abs(GetAngle(Point(1,0), Point(0,1)) - pi*0.5) > eps) {
        warn("Не прошел тест на угол");
    }

    if(a*0.5*2 != a) {
        warn("Не прошел тест умножения");
    }

    std::cout << "===============" << std::endl;
}


void LineTest() {
    std::cout << "Тест структуры прямая" << std::endl;


    std::cout << "===============" << std::endl;
}


int main() {
    PointTest();
    Polygon poly(Point(0,0), Point(1,1), Point(1,0), Point(1,1));
    return 0;
}