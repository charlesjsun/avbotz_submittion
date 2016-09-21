// Copyright (c) 2016, Charles JL Sun. All rights reserved.

#include <iostream>
#include <vector>

#define DEBUG 1

// Returns whether a == b with a given precision (epsilon)
bool double_compare(double a, double b, double epsilon)
{
	return std::fabs(a - b) <= epsilon;
}

// Returns whether a and b has the same sign
bool same_sign(double a, double b)
{
	return (a > 0.0 && b > 0.0) || (a < 0.0 && b < 0.0);
}

// Returns the y value given an x value for a hyperbola with foci at <x1, y1> and <x2, y2>, and the constant distant in between them, t.
void hyp_func(const double x_in, double &y1_out, double &y2_out, const double x1, const double y1, const double x2, const double y2, const double t)
{

	// Offset everything by <x1, y1>
	const double nx2 = x2 - x1;
	const double ny2 = y2 - y1;
	const double x = x_in - x1;

	const double u = -nx2;
	const double v = -ny2;
	const double w = nx2 * nx2 + ny2 * ny2 - t * t;

	const double a = 4.0 * (v - t) * (v + t);
	const double b = 4.0 * v * (2.0 * u * x + w);
	//const double c = (4.0 * (u - t) * (u + t)) * x * x + (4.0 * u * w) * x + (w * w);
	const double c = 4.0 * x * ((u - t) * (u + t) * x + u * w) + (w * w);

	// Discriminant
	const double d = b * b - 4.0 * a * c;

	// Negative discriminant is complex numbers, no real solutions
	if (d < 0.0)
		return;

	const double sqrt_d = std::sqrt(d);

	double y_pos = (sqrt_d - b) / (2.0 * a);
	double y_neg = (-sqrt_d - b) / (2.0 * a);

	// Test in original equation
	const double origin_pos = std::sqrt(x * x + y_pos * y_pos);
	const double origin_neg = std::sqrt(x * x + y_neg * y_neg);
	const double t_pos = std::sqrt((x_in - x2) * (x_in - x2) + (y_pos - ny2) * (y_pos - ny2)) - origin_pos;
	const double t_neg = std::sqrt((x_in - x2) * (x_in - x2) + (y_neg - ny2) * (y_neg - ny2)) - origin_neg;

	// Re-add the offset so it's now in world space
	y_pos += y1;
	y_neg += y1;

	y1_out = double_compare(t_pos, t, 1e-6) ? y_pos : INFINITY;
	y2_out = double_compare(t_neg, t, 1e-6) ? y_neg : INFINITY;

}

// Check the result to see if the x_in value correspond with the actual intersection
bool check_result(const double x_in, double &y_out, const double x1, const double y1, const double x2, const double y2, const double x3, const double y3, const double t1, const double t2)
{
	double out1 = INFINITY, out2 = INFINITY, out3 = INFINITY, out4 = INFINITY;

	hyp_func(x_in, out1, out2, x1, y1, x2, y2, t1);
	hyp_func(x_in, out3, out4, x1, y1, x3, y3, t2);

	double closest_to_zero = INFINITY;
	double y_value = INFINITY;
	double d1, d2, d3, d4;

#define SET_D_CLOSEST_TO_ZERO(d_num, a, b) \
		{ \
			if (out##a != INFINITY && out##b != INFINITY) \
			{ \
				d##d_num = out##a - out##b; \
				if (std::fabs(d##d_num) < std::fabs(closest_to_zero)) { \
					closest_to_zero = d##d_num; \
					y_out = (out##a + out##b) / 2.0; \
				}\
			} \
		}

	SET_D_CLOSEST_TO_ZERO(1, 1, 3);
	SET_D_CLOSEST_TO_ZERO(2, 1, 4);
	SET_D_CLOSEST_TO_ZERO(3, 2, 3);
	SET_D_CLOSEST_TO_ZERO(4, 2, 4);

	if (double_compare(closest_to_zero, 0.0, 1e-5))
	{
		return true;
	}

	return false;

}

// The 4th degree polynomial with the roots depicting the intersection of the two hyperbolas
double intersect_func(const double x_in, const double x1, const double y1, const double x2, const double y2, const double x3, const double y3, const double t1, const double t2)
{
	const double x = x_in - x1;

	#define MAKE_VARS(n, coord) \
		const double nx##coord = x##coord - x1; \
		const double ny##coord = y##coord - y1; \
		const double u##n = -nx##coord; \
		const double v##n = -ny##coord; \
		const double w##n = nx##coord * nx##coord + ny##coord * ny##coord - t##n * t##n; \
		const double a##n = 4.0 * (v##n - t##n) * (v##n + t##n); \
		const double b##n = 4.0 * v##n * (2.0 * u##n * x + w##n); \
		const double c##n = 4.0 * x * ((u##n - t##n) * (u##n + t##n) * x + u##n * w##n) + (w##n * w##n); 

	MAKE_VARS(1, 2);
	MAKE_VARS(2, 3);

	return - (a2*b1*b1*c2) - (a1*b2*b2*c1) - (a2*a2*c1*c1) - (a1*a1*c2*c2) + (a2*b1*b2*c1) + (a1*b1*b2*c2) + (2.0*a1*a2*c1*c2);
}

#if DEBUG
int num_finder = 0;
double real_x1;
double real_y1;
#endif

// A root finder that uses midpoint method (binary search) to find the root of a polynomial given an upper and lower bound.
class RootFinder
{

private:

	double low_bound;
	double high_bound;
	double low_bound_y;
	double high_bound_y;

	double x1, y1, x2, y2, x3, y3, t1, t2;

#if DEBUG
	int id;
#endif

public:

	RootFinder(double low_bound, double high_bound, double low_bound_y, double high_bound_y, double x1, double y1, double x2, double y2, double x3, double y3, double t1, double t2)
		: low_bound(low_bound), high_bound(high_bound), low_bound_y(low_bound_y), high_bound_y(high_bound_y), x1(x1), y1(y1), x2(x2), y2(y2), x3(x3), y3(y3), t1(t1), t2(t2)
	{
#if DEBUG
		id = ++num_finder;
#endif
	}

	// Called each iteration of the main loop, returns true if a root is found and sets root_out to the found root.
	bool update(double &root_out)
	{
		const double midpoint = (high_bound + low_bound) / 2.0; // get the mid point x value
		const double y_value = intersect_func(midpoint, x1, y1, x2, y2, x3, y3, t1, t2); // find the y from the x value

#if DEBUG
		std::cout << id <<  " - low: " << low_bound + real_x1 << " high: " << high_bound + real_x1 << ", y: " << y_value << std::endl;
#endif

		// If the y value is 0, then exit the loop because a root is found.
		if (double_compare(0.0, y_value, 1e-7))
		{
			//system("PAUSE");
			root_out = midpoint;
			return true;
		}

		// if the sign of the lower bound y value is the same as the mid point y value, then set midpoint as the lower bound
		// else set the midpoint as the higher bound
		if (same_sign(low_bound_y, y_value))
		{
			low_bound = midpoint;
			low_bound_y = y_value;
		}
		else
		{
			high_bound = midpoint;
			high_bound_y = y_value;
		}

		return false;
	}

};

void tdoa(const double x1, const double y1, const double x2, const double y2, const double x3, const double y3, const double t1, const double t2)
{

	// One loop going in the positive direction and another going in the negative direction, each one with an up and down curve
	const double i = std::fmin(t1 / std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)), t2 / std::sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3)));

#if DEBUG
	real_x1 = x1;
	real_y1 = y1;
#endif

	// offset everything by <x1, y1> so the <x1, y1> is the origin
	const double x1_o = 0;
	const double y1_o = 0;
	const double x2_o = x2 - x1;
	const double y2_o = y2 - y1;
	const double x3_o = x3 - x1;
	const double y3_o = y3 - y1;

	// Two i's looping positive and negative values
	double x_pos = 0;
	double x_neg = 0;
	double y_pos = intersect_func(0, x1_o, y1_o, x2_o, y2_o, x3_o, y3_o, t1, t2);
	double y_neg = y_pos;

	double prev_x_pos = INFINITY;
	double prev_x_neg = INFINITY;
	double prev_y_pos = INFINITY;
	double prev_y_neg = INFINITY;

	// Final values
	double final_x = 0;
	double final_y = 0;

	// The array of RootFinders to update every loop
	std::vector<RootFinder*> root_finders;

	bool done = false;

	while (!done)
	{
		// Start at 0, then for each of the iteration, increase/decrease the current x value by i.
		// If the y value of the previous x is of a different sign (signifying the polynomial crossing the x axis),
		// then create a RootFinder to find that root.
		{
			prev_x_pos = x_pos;
			prev_y_pos = y_pos;
			x_pos += i;
			y_pos = intersect_func(x_pos, x1_o, y1_o, x2_o, y2_o, x3_o, y3_o, t1, t2);
			if (!same_sign(prev_y_pos, y_pos))
			{
#if DEBUG
				std::cout << x_pos + x1 << " low: " << prev_x_pos + x1 << " high: " << x_pos + x1 << ", low y: " << prev_y_pos + y1 << ", high y: " << y_pos + y1 << std::endl;
#endif
				root_finders.push_back(new RootFinder(prev_x_pos, x_pos, prev_y_pos, y_pos, x1_o, y1_o, x2_o, y2_o, x3_o, y3_o, t1, t2));
			}

			prev_x_neg = x_neg;
			prev_y_neg = y_neg;
			x_neg -= i;
			y_neg = intersect_func(x_neg, x1_o, y1_o, x2_o, y2_o, x3_o, y3_o, t1, t2);
			if (!same_sign(prev_y_neg, y_neg))
			{
#if DEBUG
				std::cout << x_neg + x1 << " low: " << prev_x_neg + x1 << " high: " << x_neg + x1 << ", low y: " << prev_y_neg + y1 << ", high y: " << y_neg + y1 << std::endl;
#endif
				root_finders.push_back(new RootFinder(x_neg, prev_x_neg, y_neg, prev_y_neg, x1_o, y1_o, x2_o, y2_o, x3_o, y3_o, t1, t2));
			}
		}

		// For each of the root finders, update, and if the root finder finds a root of the polynomial,
		// then check to see if the root is correct, since only one root is the actual x intercept value of the 2 half hyperbolas
		for (int i = 0; i < root_finders.size(); ++i)
		{
			double x;
			if (root_finders[i]->update(x))
			{
				double y_out = INFINITY;
				if (check_result(x, y_out, x1_o, y1_o, x2_o, y2_o, x3_o, y3_o, t1, t2))
				{
					// Exit out of the loop. Remember to re-offset the coordinates
					final_x = x + x1;
					final_y = y_out + y1;
					done = true;
					break;
				}
			}

		}

	}

	std::cout << std::endl << "Final: " << final_x << ", " << final_y << std::endl;

	// Remove for Linux
	system("PAUSE");

	// Delete all the RootFinders from the heap.
	for (int i = 0; i < root_finders.size(); ++i)
	{
		delete root_finders[i];
	}
	root_finders.clear();
}

int main()
{
	/// -0.0802, -29.428
	const double x1 = 0.4;
	const double y1 = 0.03;
	const double x2 = -1.28;
	const double y2 = 0.4;
	const double x3 = -0.9;
	const double y3 = 1.07;

	/// 1.178, -0.5055
	//const double x1 = 1.615;
	//const double y1 = 0.1;
	//const double x2 = 0.16;
	//const double y2 = 0.0;
	//const double x3 = 0.144;
	//const double y3 = 0.96;

	//const double x1 = 1.74;
	//const double y1 = -0.02;
	//const double x2 = 0.16;
	//const double y2 = 0.0;
	//const double x3 = 0.144;
	//const double y3 = 0.96;

	/// 2.7821, 0.1368
	//const double x1 = 0.005;
	//const double y1 = 0.008;
	//const double x2 = 1.54;
	//const double y2 = -2.78;
	//const double x3 = -0.02;
	//const double y3 = -2.47;

	/// 3.50348, 4.00356
	//const double x1 = 0.0;
	//const double y1 = 1.0;
	//const double x2 = 0.5;
	//const double y2 = 0.0;
	//const double x3 = -0.5;
	//const double y3 = 0.0;

	//const double x1 = -0.32;
	//const double y1 = 2.9;
	//const double x2 = 0.67;
	//const double y2 = 2.41;
	//const double x3 = -0.46;
	//const double y3 = 1.76;

	//const double x1 = -0.73;
	//const double y1 = 0.15;
	//const double x2 = -1.69;
	//const double y2 = 2.07;
	//const double x3 = -0.46;
	//const double y3 = 1.76;

	//const double x1 = 0.0;
	//const double y1 = 0.0;
	//const double x2 = -1.07;
	//const double y2 = -2.26;
	//const double x3 = -1.32;
	//const double y3 = -4.03;

	const double s = 1481.0;
	const double t2 = 2.63474 * 0.0001 * s;
	const double t3 = 7.07023 * 0.0001 * s;

	tdoa(x1, y1, x2, y2, x3, y3, t2, t3);

	return 0;

}
