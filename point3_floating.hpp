#include <cmath>

using namespace std;

typedef long long itype;
typedef long double ftype;

const ftype pi = acosl(-1); // 3.1415926535897932384626433832795l
const ftype radian = 180 / pi; // 57.295779513082320876798154814105l
const ftype eps = 1e-12l;
const int inf = 0x7f7f7f7f;
const long long infll = 0x7f7f7f7f7f7f7f7fll;
const ftype infl = 1e20l;

template <class T>
inline int sgn(const T & x)
{
	return (x > eps) - (x < -eps);
}

template <class T>
inline T non_negative(const T & x)
{
	if (x < 0)
		return 0;
	else
		return x;
}

template <class T>
inline T trigonometric(const T & x)
{
	if (x < -1)
		return -1;
	else if (x < 1)
		return x;
	else
		return 1;
}

template <class T>
inline T sqr(const T & x)
{
	return x * x;
}

#define Vector3 Point3

template <class T>
class Point3
{
public:
	T x, y, z;

	Point3(void) : x(0), y(0), z(0)
	{

	}

	Point3(const T & x, const T & y, const T & z) : x(x), y(y), z(z)
	{

	}

	template <class S>
	Point3(const Point3<S> & src) : x(src.x), y(src.y), z(src.z)
	{

	}

	Vector3 operator + (const Vector3 & rhs) const
	{
		return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
	}

	Vector3 operator - (const Vector3 & rhs) const
	{
		return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
	}

	T length(void) const
	{
		return sqrtl(x * x + y * y + z * z);
	}

	T length2(void) const
	{
		return x * x + y * y + z * z;
	}

	T distance(const Point3 & rhs) const
	{
		return (rhs - *this).length();
	}

	T distance2(const Point3 & rhs) const
	{
		return (rhs - *this).length2();
	}

	Vector3 operator * (const T & rhs) const
	{
		return Vector3(x * rhs, y * rhs, z * rhs);
	}

	Vector3 operator / (const T & rhs) const
	{
		return Vector3(x / rhs, y / rhs, z / rhs);
	}

	Vector3 operator * (const Vector3 & rhs) const
	{
		return Vector3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x);
	}

	T operator & (const Vector3 & rhs) const
	{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}

	Vector3 cross(const Point3 & lhs, const Point3 & rhs) const
	{
		return (lhs - *this) * (rhs - *this);
	}

	T dot(const Point3 & lhs, const Point3 & rhs) const
	{
		return (lhs - *this) & (rhs - *this);
	}

	T operator ^ (const Point3 & rhs) const
	{
		T l1 = length(), l2 = rhs.length();
		if (sgn(l1) == 0 || sgn(l2) == 0)
			return 0;
		else
			return acosl(trigonometric((*this & rhs) / l1 / l2));
	}

	T angle(const Point3 & lhs, const Point3 & rhs) const
	{
		return (lhs - *this) ^ (rhs - *this);
	}

	bool operator == (const Point3 & rhs) const
	{
		return sgn(x - rhs.x) == 0 && sgn(y - rhs.y) == 0 && sgn(z - rhs.z) == 0;
	}

	bool operator != (const Point3 & rhs) const
	{
		return sgn(x - rhs.x) != 0 || sgn(y - rhs.y) != 0 || sgn(z - rhs.z) != 0;
	}

	bool operator < (const Point3 & rhs) const
	{
		if (sgn(z - rhs.z) != 0)
			return z < rhs.z;
		else if (sgn(y - rhs.y) != 0)
			return y < rhs.y;
		else
			return x < rhs.x;
	}

	Vector3 unit(const T & n = 1) const
	{
		T l = length();
		return sgn(l) == 0 ? Vector3(n, 0, 0) : Vector3(x * n / l, y * n / l, z * n / l);
	}

	Vector3 normal(void) const
	{
		return sgn(x) == 0 ? Vector3(1, 0, 0) : Vector3(-y, x, 0).unit();
	}

	Vector3 moderate(void) const
	{
		int xs = sgn(x), ys = sgn(y), zs = sgn(z);
		return zs > 0 || zs == 0 && (ys > 0 || ys == 0 && xs >= 0) ? Vector3(x, y, z) : Vector3(-x, -y, -z);
	}

	Point3 transform(const Vector3 & X, const Vector3 & Y, const Vector3 & Z, const Point3 & O = Point3(0, 0, 0)) const
	{
		Vector3 YZ = Y * Z, ZX = Z * X, XY = X * Y;
		Vector3 P = *this - O;
		return Point3((P & YZ) / (X & YZ), (P & ZX) / (Y & ZX), (P & XY) / (Z & XY));
	}

	Point3 project(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return p;
		else
			return p + (q - p) * p.dot(q, *this) / p.distance2(q);
	}

	Point3 mirror(const Point3 & p, const Point3 & q) const
	{
		return project(p, q) * 2 - *this;
	}

	Point3 drop(const Point3 & o, const Vector3 & n) const
	{
		if (n == Vector3(0, 0, 0))
			return o;
		else
			return *this - n * ((*this - o) & n) / n.length2();
	}

	Vector3 rotate(const T & a) const
	{
		T c = cosl(a), s = sinl(a);
		return Point3(x * c - y * s, x * s + y * c, z);
	}

	Point3 rotate(const Point3 & p, const Point3 & q, const T & a) const
	{
		if (p == q)
			return *this;
		else
		{
			Point3 t = project(p, q), x = *this - t, y = (q - p).unit() * x;
			return t + x * cosl(a) + y * sinl(a);
		}
	}

	T distance_to_line(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return distance(p);
		else
			return cross(p, q).length() / p.distance(q);
	}

	T distance_to_halfline(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return distance(p);
		else if (p.dot(q, *this) < 0)
			return distance(p);
		else
			return distance_to_line(p, q);
	}

	T distance_to_segment(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return distance(p);
		else if (p.dot(q, *this) < 0)
			return distance(p);
		else if (q.dot(p, *this) < 0)
			return distance(q);
		else
			return distance_to_line(p, q);
	}

	T distance_to_plane(const Point3 & o, const Vector3 & n) const
	{
		if (n == Vector3(0, 0, 0))
			return distance(o);
		else
			return ((*this - o) & n) / n.length();
	}

	bool on_line(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return *this == p;
		else
			return sgn(cross(p, q).length()) == 0;
	}

	bool on_halfline(const Point3 & p, const Point3 & q, const bool & inclusive = true) const
	{
		if (*this == p)
			return inclusive;
		else
			return on_line(p, q) && sgn(p.dot(q, *this)) >= 0;
	}

	bool on_segment(const Point3 & p, const Point3 & q, const bool & inclusive = true) const
	{
		if (*this == p || *this == q)
			return inclusive;
		else
			return on_line(p, q) && sgn(dot(p, q)) <= 0;
	}

	bool in_triangle(const Point3 & u, const Point3 & v, const Point3 & w, const bool & inclusive = true) const
	{
		Point3 p[3] = {u, v, w};
		Point3 pn = p[0].cross(p[1], p[2]);
		if (sgn(pn.length()) == 0)
			return inclusive && (on_segment(v, w) || on_segment(w, u) || on_segment(u, v));
		else
			pn = pn.unit();
		if (sgn((*this - p[0]) & pn) != 0)
			return false;
		for (int i = 0; i < 3; i++)
		{
			if (on_segment(p[i], p[(i + 1) % 3]))
				return inclusive;
			else if (sgn(cross(p[i], p[(i + 1) % 3]) & pn) < 0)
				return false;
		}
		return true;
	}
};

template <class T>
inline bool xyzcmp(const Point3<T> & lhs, const Point3<T> & rhs)
{
	if (sgn(lhs.x - rhs.x) != 0)
		return lhs.x < rhs.x;
	else if (sgn(lhs.y - rhs.y) != 0)
		return lhs.y < rhs.y;
	else
		return lhs.z < rhs.z;
}

template <class T>
inline bool zyxcmp(const Point3<T> & lhs, const Point3<T> & rhs)
{
	if (sgn(lhs.z - rhs.z) != 0)
		return lhs.z < rhs.z;
	else if (sgn(lhs.y - rhs.y) != 0)
		return lhs.y < rhs.y;
	else
		return lhs.x < rhs.x;
}

