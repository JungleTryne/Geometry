#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

const double eps = 1e-6;

struct Point {
    double x;
    double y;

    explicit Point(double x=0, double y=0);
    Point(const Point& other);

    bool operator==(const Point& other) const;
    bool operator!=(const Point& other) const;

    double getLength() const;
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

Point operator-(const Point& one, const Point& two) {
    Point newPoint(one.x - two.x, one.y - two.y);
    return newPoint;
}

double GetTriangleArea(const Point& one, const Point& two, const Point& three) {
    double lengthOne = (two - one).getLength();
    double lengthTwo = (three - two).getLength();
    double lengthThree = (one - three).getLength();
    double halfPerimeter = (lengthOne + lengthTwo + lengthThree) / 2;
    return sqrt(halfPerimeter *
            (halfPerimeter - lengthOne) *
            (halfPerimeter - lengthTwo) *
            (halfPerimeter - lengthThree));
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

class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& other) const = 0;
    virtual bool isCongruentTo(const Shape& other) const = 0;
    virtual bool isSimilarTo(const Shape& other) const = 0;
    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflex(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
};

class Polygon : public Shape {
private:
    std::vector<Point> vertices;
public:
    double perimeter() const override;
    double area() const override;
    bool operator==(const Polygon& other) const;
    bool isCongruentTo(const Polygon& other) const;
    bool isSimilarTo(const Polygon& other) const;
    bool containsPoint(const Point& point) const override;

};

double Polygon::perimeter() const {
    double sum = 0;
    for(size_t i = 0; i < this->vertices.size()-1; ++i) {
        sum += (this->vertices[i+1] - this->vertices[i]).getLength();
    }
    return sum;
}

double Polygon::area() const {
    return -1.0;
    //TODO: implement
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

bool Polygon::containsPoint(const Point &point) const {
    return false;
}


int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
