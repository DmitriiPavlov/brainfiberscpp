#pragma once
#include <cmath>

class Vec3 {
public:
    Vec3(double x = 0, double y = 0, double z = 0) : x_(x), y_(y), z_(z) {}
    double getX() const { return x_; }
    double getY() const { return y_; }
    double getZ() const { return z_; }

    Vec3& operator+=(const Vec3& o) { x_ += o.x_; y_ += o.y_; z_ += o.z_; return *this; }
    Vec3& operator-=(const Vec3& o) { x_ -= o.x_; y_ -= o.y_; z_ -= o.z_; return *this; }
    Vec3& operator*=(double s)      { x_ *= s; y_ *= s; z_ *= s; return *this; }
    Vec3& operator/=(double s)      { x_ /= s; y_ /= s; z_ /= s; return *this; }

    friend Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
    friend Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
    friend Vec3 operator*(Vec3 a, double s)      { return a *= s; }
    friend Vec3 operator*(double s, Vec3 a)      { return a *= s; }
    friend Vec3 operator/(Vec3 a, double s)      { return a /= s; }
    friend Vec3 operator-(const Vec3& a)         { return Vec3(-a.x_, -a.y_, -a.z_); }
    friend bool operator==(const Vec3&, const Vec3&) = default;

    double dot(const Vec3& o) const { return x_ * o.x_ + y_ * o.y_ + z_ * o.z_; }
    double norm2() const { return dot(*this); }
    double norm() const { return std::sqrt(norm2()); }

private:
    double x_, y_, z_;
};

class Vec2 {
public:
    Vec2(double x = 0, double y = 0) : x_(x), y_(y) {}
    double getX() const { return x_; }
    double getY() const { return y_; }

    Vec2& operator+=(const Vec2& o) { x_ += o.x_; y_ += o.y_; return *this; }
    Vec2& operator-=(const Vec2& o) { x_ -= o.x_; y_ -= o.y_; return *this; }
    Vec2& operator*=(double s)      { x_ *= s; y_ *= s; return *this; }
    Vec2& operator/=(double s)      { x_ /= s; y_ /= s; return *this; }

    friend Vec2 operator+(Vec2 a, const Vec2& b) { return a += b; }
    friend Vec2 operator-(Vec2 a, const Vec2& b) { return a -= b; }
    friend Vec2 operator*(Vec2 a, double s)      { return a *= s; }
    friend Vec2 operator*(double s, Vec2 a)      { return a *= s; }
    friend Vec2 operator/(Vec2 a, double s)      { return a /= s; }
    friend Vec2 operator-(const Vec2& a)         { return Vec2(-a.x_, -a.y_); }
    friend bool operator==(const Vec2&, const Vec2&) = default;

    double dot(const Vec2& o) const { return x_ * o.x_ + y_ * o.y_; }
    double norm2() const { return dot(*this); }
    double norm() const { return std::sqrt(norm2()); }

private:
    double x_, y_;
};

