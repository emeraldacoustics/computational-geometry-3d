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

template <class T, class U = ftype>
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

	U length(void) const
	{
		return sqrtl(x * x + y * y + z * z);
	}

	T length2(void) const
	{
		return x * x + y * y + z * z;
	}

	U distance(const Point3 & rhs) const
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

	U operator ^ (const Point3 & rhs) const
	{
		U l1 = length(), l2 = rhs.length();
		if (sgn(l1) == 0 || sgn(l2) == 0)
			return 0;
		else
			return acosl(trigonometric(U(*this & rhs) / l1 / l2));
	}

	U angle(const Point3 & lhs, const Point3 & rhs) const
	{
		return (lhs - *this) ^ (rhs - *this);
	}

	bool operator == (const Point3 & rhs) const
	{
		return x == rhs.x && y == rhs.y && z == rhs.z;
	}

	bool operator != (const Point3 & rhs) const
	{
		return x != rhs.x || y != rhs.y || z != rhs.z;
	}

	bool operator < (const Point3 & rhs) const
	{
		if (z != rhs.z)
			return z < rhs.z;
		else if (y != rhs.y)
			return y < rhs.y;
		else
			return x < rhs.x;
	}

	Vector3 unit(const T & n = 1) const
	{
		T l = abs(__gcd(x, __gcd(y, z)));
		return l == 0 ? Vector3(n, 0) : Vector3(x / l * n, y / l * n, z / l * n);
	}

	Vector3 normal(void) const
	{
		return sgn(x) == 0 ? Vector3(1, 0, 0) : Vector3(-y, x, 0).unit();
	}

	Vector3 moderate(void) const
	{
		return z > 0 || z == 0 && (y > 0 || y == 0 && x >= 0) ? Vector3(x, y, z) : Vector3(-x, -y, -z);
	}

	U distance_to_line(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return distance(p);
		else
			return cross(p, q).length() / p.distance(q);
	}

	U distance_to_halfline(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return distance(p);
		else if (p.dot(q, *this) < 0)
			return distance(p);
		else
			return distance_to_line(p, q);
	}

	U distance_to_segment(const Point3 & p, const Point3 & q) const
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

	U distance_to_plane(const Point3 & o, const Point3 & n) const
	{
		if (n == Vector3(0, 0, 0))
			return distance(o);
		else
			return U((*this - o) & n) / n.length();
	}

	bool on_line(const Point3 & p, const Point3 & q) const
	{
		if (p == q)
			return *this == p;
		else
			return cross(p, q) == Vector3(0, 0, 0);
	}

	bool on_halfline(const Point3 & p, const Point3 & q, const bool & inclusive = true) const
	{
		if (*this == p)
			return inclusive;
		else
			return on_line(p, q) && p.dot(q, *this) >= 0;
	}

	bool on_segment(const Point3 & p, const Point3 & q, const bool & inclusive = true) const
	{
		if (*this == p || *this == q)
			return inclusive;
		else
			return on_line(p, q) && dot(p, q) <= 0;
	}
};

template <class T>
inline bool xyzcmp(const Point3<T> & lhs, const Point3<T> & rhs)
{
	if (lhs.x != rhs.x)
		return lhs.x < rhs.x;
	else if (lhs.y != rhs.y)
		return lhs.y < rhs.y;
	else
		return lhs.z < rhs.z;
}

template <class T>
inline bool zyxcmp(const Point3<T> & lhs, const Point3<T> & rhs)
{
	if (lhs.z != rhs.z)
		return lhs.z < rhs.z;
	else if (lhs.y != rhs.y)
		return lhs.y < rhs.y;
	else
		return lhs.x < rhs.x;
}

