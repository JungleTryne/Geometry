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

    if(std::abs(Point(1, 1).getLength() - sqrt(2)) >= eps) {
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

    Point f(-1, 1);
    Point h(1,1);
    if(GetScaledPoint(h, f, 2) != Point(3,1)) {
        warn("Не прошел тест GetScaledPoint");
    }

    f = Point(0, 0);
    h = Point(1,1);
    if(GetScaledPoint(h, f, 2) != Point(2,2)) {
        warn("Не прошел тест GetScaledPoint");
    }

    std::cout << "===============" << std::endl;
}


void LineTest() {
    std::cout << "Тест структуры прямая" << std::endl;
    Line one(Point(1,1), Point(2,2));
    Line two(Point(0,0), Point(3,3));
    Line three(Point(1,0), Point(0,0));
    if(! ((one == two) && (one != three))) {
        warn("Прошел тест == и !=");
    }

    Line four(Point(0,0), Point(4, 3));
    Point normalVector = four.getNormalVector();
    double angle = GetAngle(normalVector, Point(-3, 4));
    if(!((std::abs(angle) < eps) || (std::abs(angle - pi) < eps))) {
        warn("Не прошел тест getNormalVector");
        std::cout << normalVector.x << ' ' << normalVector.y << std::endl;
        std::cout << angle << std::endl;
    }

    Point a(2,2);
    Point c(-2, 2);
    Line axis(Point(0,0), Point(0,2));
    Point b = GetReflectedPoint(a, axis);
    Point d = GetReflectedPoint(c, axis);
    if(b != Point(-2, 2) || d != Point(2,2)) {
        warn("Не прошел тест GetReflectedPoint");
        std::cout << b.x << ' ' << b.y << std::endl;
    }

    Point z(2, 0);
    Line axis2(Point(0,0), Point(1,1));
    Point x = GetReflectedPoint(z, axis2);
    if(x != Point(0, 2)) {
        warn("Не прошел тест GetReflectedPoint");
    }

    std::cout << "===============" << std::endl;
}

void PolygonTest() {
    std::cout << "Тест структуры многоугольник" << std::endl;

    Polygon poly1(Point(-2,2), Point(1,2), Point(6,1) ,Point(3,-1), Point(-1,-1));
    Polygon poly2(Point(1,2), Point(6,1) ,Point(3,-1), Point(-1,-1), Point(-2,2));
    poly1.rotate(Point(0,0), 360);
    std::cout << "==" << std::endl;
    std::cout << (poly1 == poly2);
    std::cout << (poly2 == poly1) << std::endl;
    std::cout << "congruent" << std::endl;
    //std::cout << poly1.isCongruentTo(poly2);
    //std::cout << poly2.isCongruentTo(poly1) << std::endl;
    std::cout << "similar" << std::endl;
    std::cout << poly1.isSimilarTo(poly2);
    std::cout << poly2.isSimilarTo(poly1) << std::endl;
    std::cout << "contains point and convex" << std::endl;

    std::cout << poly1.containsPoint(Point(0,0));
    std::cout << poly1.isConvex();
    std::cout << std::endl;

    Polygon poly3(Point(0,0), Point(1,1), Point(2, 0));
    Triangle poly4(Point(0,0), Point(2,0), Point(1, -1));
    //std::cout << (poly3 == poly4);
    //std::cout << (poly4 == poly3);
    //std::cout << poly3.isCongruentTo(poly4);
    //std::cout << poly4.isCongruentTo(poly3);
    std::cout << poly3.isSimilarTo(poly4);
    std::cout << poly4.isSimilarTo(poly3);
    std::cout << std::endl;

    Triangle poly5(Point(0,0), Point(0.5,0.5), Point(1,0));
    std::cout << (poly3 == poly5);
    std::cout << (poly5 == poly3);
    //std::cout << poly3.isCongruentTo(poly5);
    //std::cout << poly5.isCongruentTo(poly3);
    std::cout << poly3.isSimilarTo(poly5);
    std::cout << poly5.isSimilarTo(poly3);
    std::cout << std::endl;

    Polygon poly6(Point(0,0), Point(2,0), Point(2,2), Point(0,2));
    Polygon poly7(Point(0,0), Point(-1,0), Point(-1,1), Point(0,1));
    std::cout << poly6.isSimilarTo(poly7) << std::endl;

    Polygon poly8(Point(0,0), Point(2,0), Point(2,2), Point(1,1), Point(0,2));
    Polygon poly9(Point(0,0), Point(2,0), Point(2,2), Point(1,3), Point(0,2));
    //std::cout << poly8.isCongruentTo(poly9) << std::endl;

    Triangle poly10(Point(1,1), Point(3,1), Point(1,3));
    Triangle poly11(Point(2,1), Point(1,1), Point(1,2));
    std::cout << poly10.isSimilarTo(poly11) << std::endl;

    Polygon poly12(Point(0,2), Point(1,3), Point(2,2), Point(2,0), Point(0,0));
    std::cout << poly9.isSimilarTo(poly12) << std::endl;

    Polygon poly13(Point(1,0), Point(0,2), Point(-1,0), Point(0,-2));
    Polygon poly14(Point(2,0), Point(0, -1), Point(-2,0), Point(0,1));
    std::cout << poly13.isSimilarTo(poly14) << std::endl;

    Polygon poly15(Point(0,0), Point(2,0), Point(0,4));
    Polygon poly16(Point(0,0), Point(0,4), Point(-2,0));
    Polygon poly17(Point(0,4), Point(0,0), Point(2,0));

    //std::cout << poly16.isCongruentTo(poly15) << std::endl;
    //std::cout << poly15.isCongruentTo(poly17);

    Polygon poly18(Point(2,2), Point(3,4), Point(5,5),Point(7,2), Point(5,3));
    Polygon poly19(Point(5,-5), Point(7,-2), Point(5,-3), Point(2,-2),Point(3,-4));
    //std::cout << poly18.isCongruentTo(poly19) << std::endl;

    Polygon poly20(Point(-3,1), Point(-5,1), Point(-5,3 ), Point(-4,2), Point(-3,3));
    Polygon poly21(Point(-4, -4), Point(-5,-3), Point(-5,-1), Point(-3,-1), Point(-3,-3));
    std::cout << (poly20 == poly19) << std::endl;

    std::cout << "similarity tests" << std::endl;
    std::cout << poly18.isSimilarTo(poly19) << std::endl;

    std::cout << "===============" << std::endl;
}

void EllipseTest() {
    std::cout << "Тест структуры эллипс" << std::endl;
    Ellipse elps1(Point(1,0), Point(-1,0), 2);
    Ellipse elps2(Point(-1, 0), Point(1, 0), 2);
    std::cout << (elps1 == elps2) << std::endl;

    std::cout << "===============" << std::endl;
}

void DifferentShapesTest() {
    std::cout << "Тест взаимодействия разных структур" << std::endl;
    Ellipse el(Point(0,0), Point(0,0), 1);
    Circle cl(Point(0,0), 1);
    std::cout << (el == cl) << std::endl;
    std::cout << (cl == el) << std::endl;
    std::cout << (el.isSimilarTo(cl)) << std::endl;
    std::cout << (cl.isSimilarTo(el)) << std::endl;
    std::cout << (el.isCongruentTo(cl)) << std::endl;
    std::cout << (cl.isCongruentTo(el)) << std::endl;

    Ellipse el2(Point(0,0), Point(0,0), 2);
    Circle cl2(Point(1,1), 1);
    std::cout << (el2 == cl2) << std::endl;
    std::cout << (cl2 == el2) << std::endl;
    std::cout << (el2.isSimilarTo(cl2)) << std::endl;
    std::cout << (cl2.isSimilarTo(el2)) << std::endl;
    std::cout << (el2.isCongruentTo(cl2)) << std::endl;
    std::cout << (cl2.isCongruentTo(el2)) << std::endl;
    std::cout << "===============" << std::endl;
}

void PolygonTest2() {
    std::cout << "Тест структуры Полигон" << std::endl;
    Polygon poly(Point(1,4), Point(3, 6), Point(4,4), Point(5,6), Point(6,3), Point(3,2));
    std::cout << poly.containsPoint(Point(3,3));
    std::cout << poly.containsPoint(Point(3,2));
    std::cout << poly.containsPoint(Point(2,3));
    std::cout << poly.containsPoint(Point(3.5,5));
    std::cout << std::endl;

    Circle c(Point(0,0), 1);
    std::cout << c.containsPoint(Point(1/sqrt(2), 1/sqrt(2)));
    std::cout << "===============" << std::endl;
}

void RectangleTest() {
    Rectangle r(Point(0,0), Point(4,3), 1.3333333333333);
    Rectangle r1(Point(1,1), Point(3,6), 2.5);
    Rectangle r2(Point(1,1), Point(3,6), 0.4);
    Square sq(Point(0,0), Point(1,1));

}

void TriangleTest() {

    Ellipse el(Point(0,0), Point(0,0), 10);
    Ellipse el1(Point(3,3), Point(3,3), 10);
    Circle c(Point(3,3), 5);
    std::cout << (el == el1) << std::endl;
    std::cout << el1.isCongruentTo(el) << std::endl;
    std::cout << el1.isSimilarTo(el) << std::endl;
    std::cout << (el1 == c) << std::endl;

}

int main() {
    //PointTest();
    //LineTest();
    //PolygonTest();
    //EllipseTest();
    //DifferentShapesTest();
    //PolygonTest2();
    //RectangleTest();
    TriangleTest();
    return 0;
}