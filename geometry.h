#include <algorithm>
#include <cmath>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <vector>

const double eps = 1e-6;
const double pi = 3.1415926535;

class NoLineException : public std::exception {
    const char* what() const noexcept override;
};

const char *NoLineException::what() const noexcept {
    return "Не существует прямой эйлера для данного треугольника";
}

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
    return Point(this->one.y - this->two.y, this->two.x - this->one.x);
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

double GetAngle(const Point& firstVector, const Point& secondVector) {
    double angle = acos((firstVector.x*secondVector.x + firstVector.y*secondVector.y)/(firstVector.getLength()* secondVector.getLength()));
    return angle;
}

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool containsPoint(const Point& point) const = 0;
    virtual ~Shape() = default;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflex(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
};

class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
public:
    Polygon(const Polygon& other);
    explicit Polygon(const std::vector<Point>& points);
    explicit Polygon(const Point& first);

    template<typename... T>
    explicit Polygon(const Point& first, const T&... other);

    double perimeter() const override;
    double area() const override;
    bool operator==(const Polygon& other) const;
    bool operator!=(const Polygon& other) const;
    bool isCongruentTo(const Polygon& other) const;
    bool isSimilarTo(const Polygon& other) const;
    bool containsPoint(const Point& point) const override;
    bool isConvex() const;
    std::vector<Point> getVertices() const;

    void rotate(const Point& center, double angle) override;
    void reflex(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;
};

double Polygon::perimeter() const {
    double sum = 0.0;
    for(size_t i = 0; i < this->vertices.size(); ++i) {
        sum += (this->vertices[(i+1) % this->vertices.size()] - this->vertices[i]).getLength();
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
        pointer++;
    }
    return std::abs(area)*0.5;
}

bool Polygon::operator==(const Polygon &other) const {
    //TODO: Можно оптимизировать!
    if(this->vertices.size() != other.vertices.size()) {
        return false;
    }
    std::vector<Point> thisCopy = this->vertices;
    for(size_t j = 0; j < thisCopy.size(); ++j) {
        bool equal = true;
        for(size_t i = 0; i < thisCopy.size(); ++i) {
            if(thisCopy[i] != other.vertices[i]) {
                equal = false;
                break;
            }
        }
        std::rotate(thisCopy.begin(), thisCopy.begin() + 1, thisCopy.end());
        if(equal) {
            return true;
        }
    }
    return false;
}

bool Polygon::isSimilarTo(const Polygon &other) const {
    std::vector<double> anglesOne;
    std::vector<double> anglesTwo;

    for(size_t i = 0; i < this->vertices.size(); ++i) {
        anglesOne.push_back(GetVectorsAngle(this->vertices[i], this->vertices[(i+1) % this->vertices.size()]));
        anglesTwo.push_back(GetVectorsAngle(other.vertices[i], other.vertices[(i+1) % this->vertices.size()]));
    }

    for(size_t i = 0; i < anglesOne.size() + 1; ++i) {
        bool equal = true;
        for(size_t j = 0; j < anglesOne.size(); ++j) {
            if(std::abs(anglesOne[j] - anglesTwo[j]) >= eps) {
                equal = false;
                break;
            }
        }
        if(equal) {
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

    for(size_t i = 0; i < this->vertices.size(); ++i) {
        lengthsOne.push_back((this->vertices[(i+1) % this->vertices.size()] - this->vertices[i]).getLength());
        lengthsTwo.push_back((other.vertices[(i+1) % this->vertices.size()] - other.vertices[i]).getLength());
    }

    for(size_t i = 0; i < lengthsOne.size() + 1; ++i) {
        bool equal = true;
        for(size_t j = 0; j < lengthsOne.size(); ++j) {
            if(std::abs(lengthsOne[j] - lengthsTwo[j]) >= eps) {
                equal = false;
                break;
            }
        }
        if(equal) {
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
        vector.rotate((angle/360)*2*pi);
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

bool Polygon::operator!=(const Polygon& other) const {
    return !(*this == other);
}

bool Polygon::isConvex() const {
    bool positive = false;
    bool negative = false;
    for(size_t i = 0; i < this->vertices.size(); ++i) {
        size_t firstIndex = i;
        size_t secondIndex = (i + 1) % this->vertices.size();
        size_t thirdIndex = (i + 2) % this->vertices.size();

        double deltaxOne = this->vertices[secondIndex].x - this->vertices[firstIndex].x;
        double deltaxTwo = this->vertices[thirdIndex].x - this->vertices[secondIndex].x;
        double deltayOne = this->vertices[secondIndex].y - this->vertices[firstIndex].y;
        double deltayTwo = this->vertices[thirdIndex].y - this->vertices[secondIndex].y;
        double zComponentProduct = deltaxOne*deltayTwo - deltayOne*deltaxTwo;
        if(zComponentProduct > 0) {
            positive = true;
        }
        if(zComponentProduct < 0) {
            negative = true;
        }
    }
    return positive ^ negative;
}

std::vector<Point> Polygon::getVertices() const {
    return this->vertices;
}

Polygon::Polygon(const Point& first) {
    this->vertices.insert(this->vertices.begin(), first);
}

template<typename... T>
Polygon::Polygon(const Point &first, const T &... other) : Polygon(other...) {
    this->vertices.insert(this->vertices.begin(), first);
}

class Rectangle : public Polygon {
public:
    explicit Rectangle(const std::vector<Point>& points);
    Rectangle(const Rectangle& other);
    Rectangle(const Point& one, const Point& three, double coefficient);

    Point center() const;
    std::pair<Line, Line> diagonals() const;
};

Point Rectangle::center() const {
    return (this->vertices[2] - this->vertices[0])*0.5 + this->vertices[0];
}

std::pair<Line, Line> Rectangle::diagonals() const {
    return std::make_pair(
            Line(this->vertices[0], this->vertices[2]),
            Line(this->vertices[1], this->vertices[3])
            );
}

std::vector<Point> getRectangleFromTwoPoints(const Point& one, const Point& three, double coefficient) {
    Point vector = three - one;
    double angle = atan(coefficient);
    Point side = vector;
    side.rotate(angle);
    side = side * (1/side.getLength()) * sqrt(std::pow(vector.getLength(),2)/(1+coefficient*coefficient));
    Point two = one + side;
    Point four = two - side;
    return std::vector<Point>{one, two, three, four};
}

Rectangle::Rectangle(const Point& one, const Point& three, double coefficient) : Polygon(getRectangleFromTwoPoints(one, three, coefficient)){}

Rectangle::Rectangle(const std::vector<Point> &points) : Polygon(points) {}

Rectangle::Rectangle(const Rectangle& other) = default;


class Ellipse : public Shape {
protected:
    std::pair<Point, Point> _focuses;
    double constSum;
    std::pair<double, double> getAxis() const;
public:
    Ellipse(const Point& one, const Point& two, double constSum);
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

Ellipse::Ellipse(const Point& one, const Point& two, double constSum) {
    this->_focuses = std::make_pair(one, two);
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
    return (this->_focuses == other._focuses ||
            this->_focuses == std::make_pair(other._focuses.second, other._focuses.first)) &&
    (this->constSum - other.constSum) < eps;
}

bool Ellipse::isCongruent(const Ellipse &other) const {
    return std::abs((this->_focuses.first - this->_focuses.second).getLength() - (other._focuses.first - other._focuses.second).getLength()) < eps &&
            std::abs(this->constSum - other.constSum) < eps;
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

Circle::Circle(const Point& _center, double radius) : Ellipse(_center, _center, 2*radius) {}

class Square : public Rectangle{
public:
    explicit Square(const Point& one, const Point& two);
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

Circle Square::circumscribedCircle() const {
    Circle circle(this->center(), ((this->vertices[2] - this->vertices[0]).getLength())/2);
    return circle;
}

Circle Square::inscribedCircle() const {
    Circle circle(this->center(), (((this->vertices[2] - this->vertices[1]).getLength())/2));
    return circle;
}

std::vector<Point> GetSquare(const Point& one, const Point& two) {
    Point vector = two - one;
    vector.rotate(pi/4);
    Point vector2 = two - one;
    vector2.rotate(-pi/4);
    std::vector<Point> result = {one, vector, two, vector2};
    return result;
}

Square::Square(const Point& one, const Point& two) : Rectangle(GetSquare(one, two)){}

class Triangle : public Polygon {
public:
    explicit Triangle(const Point& one, const Point& two, const Point& three);
    Point orthocenter() const;
    Point centroid() const;
    Line EulerLine() const;
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Circle ninePointsCircle() const;
};

Triangle::Triangle(const Point& one, const Point& two, const Point& three) : Polygon(std::vector<Point>{one, two, three}) {}

Point Triangle::orthocenter() const {
    //Формула взята с сайта
    // https://math.stackexchange.com/questions/478069/how-to-calculate-the-coordinates-of-orthocentre
    Point vectorOne = this->vertices[1] - this->vertices[0];
    Point vectorTwo = this->vertices[2] - this->vertices[0];
    Point vectorThree = this->vertices[2] - this->vertices[1];
    double tanA = tan(GetAngle(vectorOne, vectorTwo));
    double tanB = tan(GetAngle(-1*vectorOne, vectorThree));
    double tanC = tan(GetAngle(-1*vectorTwo, -1*vectorThree));
    double orthoX = (vertices[0].x*tanA + vertices[1].x*tanB + vertices[2].x*tanC)/(tanA + tanB + tanC);
    double orthoY = (vertices[0].y*tanA + vertices[1].y*tanB + vertices[2].y*tanC)/(tanA + tanB + tanC);
    return Point(orthoX, orthoY);
}

Point Triangle::centroid() const {
    return (this->vertices[0] + this->vertices[1] + this->vertices[2])*(0.3333333333333);
}

Line Triangle::EulerLine() const {
    Point _orthocenter = this->orthocenter();
    Point _centroid = this->centroid();
    if(_orthocenter == _centroid) {
        throw NoLineException();
    }
    Line line(this->orthocenter(), this->centroid());
    return line;
}

Circle Triangle::circumscribedCircle() const {
    Point vectorOne = this->vertices[1] - this->vertices[0];
    Point vectorTwo = this->vertices[2] - this->vertices[0];
    Point vectorThree = this->vertices[2] - this->vertices[1];
    double radius = vectorOne.getLength()*vectorTwo.getLength()*vectorThree.getLength()/(4*this->area());
    Circle circle(this->orthocenter(), radius);
    return circle;
}

Circle Triangle::inscribedCircle() const {
    //Формула взята с сайта
    // https://mathworld.wolfram.com/Incenter.html
    Point vectorOne = this->vertices[1] - this->vertices[0];
    Point vectorTwo = this->vertices[2] - this->vertices[0];
    Point vectorThree = this->vertices[2] - this->vertices[1];

    double oneLength = vectorOne.getLength();
    double twoLength = vectorTwo.getLength();
    double threeLength = vectorThree.getLength();

    double halfPerimeter = 0.5*(oneLength + twoLength + threeLength);
    double radius = this->area()/halfPerimeter;

    Point center(
            (oneLength*this->vertices[2].x + twoLength*this->vertices[0].x + threeLength*this->vertices[1].x)/
                    (oneLength + twoLength + threeLength),
            ((oneLength*this->vertices[2].y + twoLength*this->vertices[0].y + threeLength*this->vertices[1].y)/
             (oneLength + twoLength + threeLength))
             );

    Circle circle(center, radius);
    return circle;
}

Circle Triangle::ninePointsCircle() const {
    Point vectorOne = this->vertices[1] - this->vertices[0];
    Point vectorTwo = this->vertices[2] - this->vertices[0];
    Point vectorThree = this->vertices[2] - this->vertices[1];

    double oneLength = vectorOne.getLength();
    double twoLength = vectorTwo.getLength();
    double threeLength = vectorThree.getLength();

    Point inscribedCircleCenter(
            (oneLength*this->vertices[2].x + twoLength*this->vertices[0].x + threeLength*this->vertices[1].x)/
            (oneLength + twoLength + threeLength),
            ((oneLength*this->vertices[2].y + twoLength*this->vertices[0].y + threeLength*this->vertices[1].y)/
             (oneLength + twoLength + threeLength))
    );

    Point _orthocenter = this->orthocenter();
    Point center = (_orthocenter + inscribedCircleCenter)*0.5;
    double radius = vectorOne.getLength()*vectorTwo.getLength()*vectorThree.getLength()/(4*this->area());
    radius /= 2;
    Circle circle(center, radius);
    return circle;
}
