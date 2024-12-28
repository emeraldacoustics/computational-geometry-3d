template <class T>
bool seg_seg_intersect(const Point<T> & p1, const Point<T> & p2, const Point<T> & q1, const Point<T> & q2, const bool & inclusive = true)
{
	Point<T> n = ((p2 - p1) * (q2 - q1)).unit();
	if (sgn(n & (q1 - p1)) != 0)
		return false;
	Point<T> qv1 = p1.cross(p2, q1), qv2 = p1.cross(p2, q2);
	Point<T> pv1 = q1.cross(q2, p1), pv2 = q1.cross(q2, p2);
	int sp = sgn(pv1 & pv2), sq = sgn(qv1 & qv2);
	if (sp > 0 || sq > 0)
		return false;
	else if (sp == 0 || sq == 0)
		return inclusive;
	else
		return true;
}

template <class T>
T line_line_distance(const Point<T> & p1, const Point<T> & p2, const Point<T> & q1, const Point<T> & q2)
{
	Point<T> pv = p2 - p1, qv = q2 - q1;
	if (sgn((pv * qv).length()) == 0)
		return p1.distance_to_line(q1, q2);
	else
		return fabsl((pv * qv).unit() & (q1 - p1));
}

template <class T>
T seg_seg_distance(const Point<T> & p1, const Point<T> & p2, const Point<T> & q1, const Point<T> & q2)
{
	Point<T> pv = p2 - p1, qv = q2 - q1;
	if (sgn((pv * qv).length()) != 0)
	{
		Point<T> n = (pv * qv).unit();
		T d = (q1 - p1) & n;
		Point<T> u1 = q1 - n * d, u2 = q2 - n * d;
		if (seg_seg_intersect(p1, p2, u1, u2))
			return fabsl(d);
	}
	return min(min(p1.distance_to_segment(q1, q2), p2.distance_to_segment(q1, q2)),
			   min(q1.distance_to_segment(p1, p2), q2.distance_to_segment(p1, p2)));
}

