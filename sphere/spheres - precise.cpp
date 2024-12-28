template <class T>
T interior_angle(const T & a, const T & b, const T & c)
{
	return acosl(trigonometric((b * b + c * c - a * a) / b / c / 2));
}

template <class T>
vector<pair<Point<T>, Point<T> > > sphere_sphere_intersection(const Point<T> & o1, const T & r1, const Point<T> & o2, const T & r2)
{
	vector<pair<Point<T>, Point<T> > > ans;
	const T d = o1.distance(o2);
	if (sgn(d) != 0)
	{
		Point<T> v = o2 - o1;
		T a = ((r1 * r1 - r2 * r2) / d + d) / 2;
		int s = sgn(fabsl(a) - r1);
		if (s <= 0)
		{
			long double h = sqrtl(r1 * r1 - a * a);
			ans.emplace_back(o1 + v.unit(a), v.unit(h));
		}
	}
	return ans;
}

template <class T>
vector<pair<Point<T>, Point<T> > > tangents(const Point<T> & o, const T & r, const Point<T> & p)
{
	const T d = p.distance(o);
	return sphere_sphere_intersection(o, r, (o + p) / 2, d / 2);
}

template <class T>
vector<Point<T> > sphere_line_intersection(const Point<T> & o, const T & r, const Point<T> & p, const Point<T> & q, const bool & multiple = false)
{
	vector<Point<T> > ans;
	const Point<T> u = o.project(p, q);
	const T d = o.distance(u);
	const int s = sgn(d - r);
	if (s == 0)
	{
		ans.push_back(u);
		if (multiple)
			ans.push_back(ans[0]);
	}
	else if (s < 0)
	{
		const T l = sqrtl(r * r - d * d);
		const Point<T> v = p == q ? (u - o).normal().unit(l) : (q - p).unit(l);
		ans.push_back(u - v);
		ans.push_back(u + v);
	}
	return ans;
}

template <class T>
T sphere_sphere_intersection_volume(const Point<T> & o1, const T & r1, const Point<T> & o2, const T & r2)
{
	T d = o1.distance(o2);
	if (sgn(r1 + r2 - d) <= 0)
		return 0;
	else if (sgn(fabsl(r2 - r1) - d) >= 0)
	{
		T r = min(r1, r2);
		return r * r * r * pi * 4 / 3;
	}
	else
	{
		T ans = 0, d2 = o1.distance2(o2), a;
		a = interior_angle(r2, r1, d); //acosl((d2 + sqr(r1) - sqr(r2)) / d / r1 / 2);
		ans += r1 * r1 * r1 * (pi * (1 - cosl(a)) / 2 - cosl(a) * sqr(sinl(a))) / 3;
		a = interior_angle(r1, r2, d); //acosl((d2 + sqr(r2) - sqr(r1)) / d / r2 / 2);
		ans += r2 * r2 * r2 * (pi * (1 - cosl(a)) / 2 - cosl(a) * sqr(sinl(a))) / 3;
		return ans;
	}
}

