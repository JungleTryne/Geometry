#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

const double eps = 1e-6;
const double pi = 3.1415926535;

struct Point {
    double x;
    double y;

    explicit Point(double x=0, double y=0);
    Point(const Point& other);

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;

    double getLength() const;
    void rotate(double angle);
};

Point::Point(double x, double y) {
    this->x = x;
    this->y = y;
}

Point::Point(const Point &other) {
    this->x = other.x;
    this->y = other.y;
}

bool Point::operator==(const Point &other) const {
    return std::abs(this->x - other.x) <= eps && std::abs(this->y - other.y) <= eps;
}

bool Point::operator!=(const Point &other) const {
    return !(*this==other);
}

double Point::getLength() const {
    return sqrt(this->x * this->x + this->y * this->y);
}

void Point::rotate(double angle) {
    //angle в радианах
    double newX = this->x*cos(angle) - this->y*sin(angle);
    double newY = this->x*sin(angle) + this->y*cos(angle);
    this->x = newX;
    this->y = newY;
}

Point operator-(const Point& one, const Point& two) {
    Point newPoint(one.x - two.x, one.y - two.y);
    return newPoint;
}

Point operator+ (const Point& one, const Point& two) {
    Point newPoint(one.x + two.x, one.y + two.y);
    return newPoint;
}

Point operator*(const Point& one, double coefficient) {
    Point newPoint(one.x * coefficient, one.y*coefficient);
    return newPoint;
}

Point operator*(double coefficient, Point& one) {
    return one*coefficient;
}

double GetVectorsAngle(const Point& one, const Point& two) {
    return acos((one.x*two.x + one.y*two.y)/(one.getLength()*two.getLength()));
}

class Line {
private:
    Point one;
    Point two;
public:
    Line(const Point& first, const Point& second);
    Line(double angle, double height);
    Line(const Point& point, double angle);

    bool operator==(const Line& other) const;
    bool operator!=(const Line& other) const;

    Point getNormalVector() const;
};

Line::Line(const Point &first, const Point &second) {
    this->one = first;
    this->two = second;
}

Line::Line(double angle, double height) {
    this->one = Point(0, height);
    this->two = Point(-(height/angle), 0);
}

Line::Line(const Point &point, double angle) {
    this->one = point;
    this->two = Point(point.x + 1, point.y + angle);
}

bool Line::operator==(const Line& other) const {
    Point pointOne = other.one;
    Point pointTwo = other.two;

    bool condOne = (this->one.x - this->two.x)*(pointOne.y - this->two.y)
                    - (pointOne.x - this->two.x)*(this->one.y - this->two.y) < eps;

    bool condTwo = (this->one.x - this->two.x)*(pointTwo.y - this->two.y)
                   - (pointTwo.x - this->two.x)*(this->one.y - this->two.y) < eps;

    return condOne && condTwo;
}

bool Line::operator!=(const Line& other) const {
    return !(*this == other);
}

Point Line::getNormalVector() const {
    return Point(this->one.x - this->two.x, this->two.y - this->one.y);
}

Point GetReflectedPoint(const Point& point, const Line& line) {
    Point normalVector = line.getNormalVector();
    double normalVectorLength = normalVector.getLength();
    normalVector.x /= normalVectorLength;
    normalVector.y /= normalVectorLength;
    double ortolineLength = 2*(normalVector.x*point.x + normalVector.y*point.y);

    Point reflected(point.x - ortolineLength*normalVector.x, point.y - ortolineLength*normalVector.y);
    return reflected;
}

Point GetScaledPoint(const Point& point, const Point& center, double coefficient) {
    Point vector = point - center;
    vector = vector*coefficient;
    Point newPoint = vector + center;
    return newPoint;
}

double GetDeterminant(double x11, double x12, double x21, double x22) {
    return x11*x22 - x12*x21;
}

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    //virtual bool operator==(const Shape& other) const {return false;}
    //virtual bool isCongruentTo(const Shape& other) const {return false;}
    //virtual bool isSimilarTo(const Shape& other) const {return false;}
    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflex(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
};

class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
public:
    //TODO: конструктор с переменным количество параметров
    Polygon(const Polygon& other);
    Polygon(const std::vector<Point>& points);

    double perimeter() const override;
    double area() const override;
    bool operator==(const Polygon& other) const;
    bool isCongruentTo(const Polygon& other) const;
    bool isSimilarTo(const Polygon& other) const;
    bool containsPoint(const Point& point) const override;

    void rotate(const Point& center, double angle) override;
    void reflex(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;
};

double Polygon::perimeter() const {
    double sum = 0;
    for(size_t i = 0; i < this->vertices.size()-1; ++i) {
        sum += (this->vertices[i+1] - this->vertices[i]).getLength();
    }
    return sum;
}

double Polygon::area() const {
    size_t pointer = 0;
    double area = 0;
    for(size_t i = 0; i < this->vertices.size(); ++i) {
        area += GetDeterminant(this->vertices[pointer].x,
                this->vertices[pointer].y,
                this->vertices[(pointer + 1) % this->vertices.size()].x,
                this->vertices[(pointer + 1) % this->vertices.size()].y
                );
    }
    return std::abs(area)*0.5;
}

bool Polygon::operator==(const Polygon &other) const {
    if(this->vertices.size() != other.vertices.size()) {
        return false;
    }
    for(size_t i = 0; i < this->vertices.size(); ++i) {
        if(this->vertices[i] != other.vertices[i]) {
            return false;
        }
    }
    return true;
}

bool Polygon::isSimilarTo(const Polygon &other) const {
    std::vector<double> anglesOne;
    std::vector<double> anglesTwo;

    for(size_t i = 0; i < this->vertices.size()-1; ++i) {
        anglesOne.push_back(GetVectorsAngle(this->vertices[i], this->vertices[i+1]));
        anglesTwo.push_back(GetVectorsAngle(other.vertices[i], other.vertices[i+1]));
    }

    for(size_t i = 0; i < anglesOne.size() + 1; ++i) {
        if(anglesOne == anglesTwo) {
            return true;
        }
        std::rotate(anglesTwo.begin(), anglesTwo.begin() + 1, anglesTwo.end());
    }

    return false;
}

bool Polygon::isCongruentTo(const Polygon &other) const {
    if(!this->isSimilarTo(other)) {
        return false;
    }
    std::vector<double> lengthsOne;
    std::vector<double> lengthsTwo;

    for(size_t i = 0; i < this->vertices.size()-1; ++i) {
        lengthsOne.push_back((this->vertices[i+1] - this->vertices[i]).getLength());
        lengthsTwo.push_back((other.vertices[i+1] - other.vertices[i]).getLength());
    }

    for(size_t i = 0; i < lengthsOne.size() + 1; ++i) {
        if(lengthsOne == lengthsTwo) {
            return true;
        }
        std::rotate(lengthsTwo.begin(), lengthsTwo.begin() + 1, lengthsTwo.end());
    }

    return false;

}

bool Polygon::containsPoint(const Point& point) const {
    bool inside = false;
    for(size_t i = 0, j = this->vertices.size()-1; i < this->vertices.size(); j = i++) {
        if ((this->vertices[i].y > point.y) != (this->vertices[j].y > point.y) &&
        (point.x <
        (this->vertices[j].x-this->vertices[i].x) * (point.y-this->vertices[i].y) / (this->vertices[j].y-this->vertices[i].y) +this->vertices[i].x)
        )
        {
            inside = !inside;
        }
    }
    return inside;
}

void Polygon::rotate(const Point &center, double angle) {
    for(Point& vertex : this->vertices) {
        Point vector = vertex - center;
        vector.rotate((angle/360)*pi);
        vertex = vector + center;
    }
}

void Polygon::reflex(const Line &axis) {
    for(Point& vertex : this->vertices) {
        vertex = GetReflectedPoint(vertex, axis);
    }
}

void Polygon::scale(const Point& center, double coefficient) {
    for(Point& vertex : this->vertices) {
        vertex = GetScaledPoint(vertex, center, coefficient);
    }
}

Polygon::Polygon(const Polygon& other) {
    this->vertices = other.vertices;
}

Polygon::Polygon(const std::vector<Point>& points) {
    this->vertices = points;
}

class Rectangle : public Polygon {
public:
    //TODO: реализовать конструктор по двум точкам
    explicit Rectangle(const std::vector<Point>& vertices);
    Rectangle(const Rectangle& other);

    Point center() const;
    std::pair<Line, Line> diagonals() const;
};

Rectangle::Rectangle(const std::vector<Point> &vertices) : Polygon(vertices) {}

Point Rectangle::center() const {
    return (this->vertices[2] - this->vertices[0])*0.5 + this->vertices[0];
}

std::pair<Line, Line> Rectangle::diagonals() const {
    return std::make_pair(
            Line(this->vertices[0], this->vertices[2]),
            Line(this->vertices[1], this->vertices[3])
            );
}

Rectangle::Rectangle(const Rectangle& other) = default;


class Ellipse : public Shape {
protected:
    std::pair<Point, Point> _focuses;
    double constSum;
    std::pair<double, double> getAxis() const;
public:
    Ellipse(const std::pair<Point, Point>& focuses, double constSum);
    double perimeter() const override;
    double area() const override;
    bool operator==(const Ellipse& other) const;
    bool isCongruent(const Ellipse& other) const;
    bool isSimilarTo(const Ellipse& other) const;
    bool containsPoint(const Point& point) const override;

    void rotate(const Point& center, double angle) override;
    void reflex(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;

    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;
};

Ellipse::Ellipse(const std::pair<Point, Point>& focuses, double constSum) {
    this->_focuses = focuses;
    this->constSum = constSum;
}

std::pair<double, double> Ellipse::getAxis() const {
    //Функция получания полуосей эллепса
    double cathet = (_focuses.first - _focuses.second).getLength() / 2;
    double hypotenuse = constSum/2;
    double smallAxis = sqrt(hypotenuse*hypotenuse - cathet*cathet);
    double bigAxis = constSum/2;
    return std::make_pair(smallAxis, bigAxis);
}

double Ellipse::perimeter() const {
    std::pair<double, double> axis = this->getAxis();
    double smallAxis = axis.first;
    double bigAxis = axis.second;
    double finalPerimeter = 2*pi*sqrt((smallAxis*smallAxis + bigAxis*bigAxis)/2);
    return finalPerimeter;
}

double Ellipse::area() const {
    std::pair<double, double> axis = this->getAxis();
    double smallAxis = axis.first;
    double bigAxis = axis.second;
    double finalArea = pi*smallAxis*bigAxis;
    return finalArea;
}

bool Ellipse::operator==(const Ellipse &other) const {
    return this->_focuses == other._focuses && this->constSum == other.constSum;
}

bool Ellipse::isCongruent(const Ellipse &other) const {
    return (this->_focuses.first - this->_focuses.second).getLength() == (other._focuses.first - other._focuses.second).getLength() &&
            this->constSum == other.constSum;
}

bool Ellipse::isSimilarTo(const Ellipse &other) const {
    double lengthOne = (this->_focuses.first - this->_focuses.second).getLength();
    double lengthTwo = (other._focuses.first - other._focuses.second).getLength();
    double coefficient = lengthOne / lengthTwo;
    double newConstSum = other.constSum * coefficient;
    return std::abs(newConstSum - this->constSum) < eps;
}

bool Ellipse::containsPoint(const Point& point) const {
    return (this->_focuses.first - point).getLength() + (this->_focuses.second - point).getLength() < this->constSum;
}

void Ellipse::rotate(const Point &center, double angle) {
    angle = (angle / 360) * 2 * pi; //переводим в радианы
    Point vectorOne = this->_focuses.first - center;
    Point vectorTwo = this->_focuses.second - center;
    vectorOne.rotate(angle);
    vectorTwo.rotate(angle);
    this->_focuses.first = vectorOne + center;
    this->_focuses.second = vectorTwo + center;
}

void Ellipse::reflex(const Line &axis) {
    Point newFocusOne = GetReflectedPoint(this->_focuses.first, axis);
    Point newFocusTwo = GetReflectedPoint(this->_focuses.second, axis);
    this->_focuses = std::make_pair(newFocusOne, newFocusTwo);
}

void Ellipse::scale(const Point &center, double coefficient) {
    Point newFocusOne = GetScaledPoint(this->_focuses.first, center, coefficient);
    Point newFocusTwo = GetScaledPoint(this->_focuses.second, center, coefficient);
    this->_focuses = std::make_pair(newFocusOne, newFocusTwo);
}

std::pair<Point, Point> Ellipse::focuses() const {
    return this->_focuses;
}

std::pair<Line, Line> Ellipse::directrices() const {
    std::pair<double, double> axis = this->getAxis();
    double smallAxis = axis.first;
    double bigAxis = axis.second;
    double radiusDistance = bigAxis/this->eccentricity();
    Point centerVector = this->center();

    Point unitVector = this->_focuses.first - this->_focuses.second;
    unitVector.x /= unitVector.getLength();
    unitVector.y /= unitVector.getLength();

    Point unitNormalVector(-unitVector.y, unitVector.x);

    Point firstLineFirstPoint = unitVector * this->eccentricity() + centerVector;
    Point firstLineSecondPoint = unitVector*this->eccentricity() + centerVector + unitNormalVector;
    Line firstLine = Line(firstLineFirstPoint, firstLineSecondPoint);

    Point secondLineFirstPoint = unitVector * (-this->eccentricity()) + centerVector;
    Point secondLineSecondPoint = unitVector * (-this->eccentricity()) + centerVector;
    Line secondLine = Line(secondLineFirstPoint, secondLineSecondPoint);

    return std::make_pair(firstLine, secondLine);
}

double Ellipse::eccentricity() const {
    std::pair<double, double> axis = this->getAxis();
    double smallAxis = axis.first;
    double bigAxis = axis.second;
    return sqrt(1 - (smallAxis*smallAxis)/(bigAxis*bigAxis));
}

Point Ellipse::center() const {
    return (this->_focuses.first - this->_focuses.second) * 0.5 + this->_focuses.second;
}


class Circle : public Ellipse {
private:
public:
    double radius() const;
    Circle(const Point& center, double radius);
};

double Circle::radius() const {
    return this->constSum/2;
}

Circle::Circle(const Point& _center, double radius) : Ellipse(std::make_pair(_center, _center), 2*radius) {}

class Square : public Rectangle{
public:
    explicit Square(const std::pair<Point, Point>& points);
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

Circle Square::circumscribedCircle() const {
    return Circle(this->center(), ((this->vertices[2] - this->vertices[0]).getLength())/2);
}

Circle Square::inscribedCircle() const {
    return Circle(this->center(), (((this->vertices[2] - this->vertices[1]).getLength())/2));
}

std::vector<Point> GetSquare(const std::pair<Point, Point>& points) {
    Point vector = points.second - points.first;
    vector.rotate(pi/4);
    Point vector2 = points.second - points.first;
    vector2.rotate(-pi/4);
    std::vector<Point> result = {points.first, vector, points.second, vector2};
    return result;
}

Square::Square(const std::pair<Point, Point>& points) : Rectangle(GetSquare(points)){}

class Triangle : public Polygon {
public:
    explicit Triangle(const std::vector<Point>& points);
};

Triangle::Triangle(const std::vector<Point>& points) : Polygon(points) {

}

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
