template <class T>
T triple(const Point<T> & a, const Point<T> & b, const Point<T> & c)
{
	return a & b * c;
}

template <class T>
Point<T> line_plane_intersection(const Point<T> & p, const Point<T> & q, const Point<T> & o, const Point<T> & n)
{
	T sp = (p - o) & n, sq = (q - o) & n;
	return (p * -sq + q * sp) / (sp - sq);
}

template <class T>
Point<T> plane_intersection(const Point<T> & a1, const Point<T> & n1, const Point<T> & a2, const Point<T> & n2, const Point<T> & a3, const Point<T> & n3)
{
	Point<T> x(n1.x, n2.x, n3.x);
	Point<T> y(n1.y, n2.y, n3.y);
	Point<T> z(n1.z, n2.z, n3.z);
	Point<T> d(a1 & n1, a2 & n2, a3 & n3);
	return Point<T>(triple(d, y, z), triple(x, d, z), triple(x, y, d)) / triple(n1, n2, n3);
}

