// Copyright (c) 2016, Charles JL Sun. All rights reserved.

#include <iostream>

#define DEBUG 1

bool double_compare(double a, double b, double epsilon = 0.0001)
{
	return std::fabs(a - b) < epsilon;
}

// x1, y1 = starting focus, x2, y2 = ending focus, t = time interval from p1 to p2
void hyp_func_old(const double x_in, double &y1_out, double &y2_out, const double x1, const double y1, const double x2, const double y2, const double t)
{

	const double u = x1 - x2;
	const double v = y1 - y2;
	const double w = x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 - t * t;

	const double a = 4.0 * (v - t) * (v + t);
	const double b = (8.0 * u * v) * x_in + (4.0 * v * w + 8.0 * t * t * y1);
	const double c = (4.0 * (u - t) * (u + t)) * x_in * x_in + (4.0 * u * w + 8 * t * t * x1) * x_in + (w * w - 4 * t * t * x1 * x1 - 4 * t * t * y1 * y1);

	// Discriminant
	const double d = b * b - 4 * a * c;

	// Negative discriminant is complex numbers, no real solutions
	if (d < 0.0)
		return;

	const double sqrt_d = std::sqrt(d);

	double y_pos = (sqrt_d - b) / (2.0 * a);
	double y_neg = (-sqrt_d - b) / (2.0 * a);

	// Test in original equation
	const double origin_pos = std::sqrt((x_in - x1) * (x_in - x1) + (y_pos - y1) * (y_pos - y1));
	const double origin_neg = std::sqrt((x_in - x1) * (x_in - x1) + (y_neg - y1) * (y_neg - y1));
	const double t_pos = std::sqrt((x_in - x2) * (x_in - x2) + (y_pos - y2) * (y_pos - y2)) - origin_pos;
	const double t_neg = std::sqrt((x_in - x2) * (x_in - x2) + (y_neg - y2) * (y_neg - y2)) - origin_neg;

	y1_out = double_compare(t_pos, t) ? y_pos : INFINITY;
	y2_out = double_compare(t_neg, t) ? y_neg : INFINITY;

}

// NOTE: And then I realised that I could just set <x1, y1> as the origin <0, 0> and simplify a lot of the equations.
//       Afterwards I just add the offset back to the resulting value.
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

	y1_out = double_compare(t_pos, t) ? y_pos : INFINITY;
	y2_out = double_compare(t_neg, t) ? y_neg : INFINITY;

}

bool find_intersection(bool top, double &x_in, double &final_x, double &final_y, double &i, double &prev_d, double x1, double y1, double x2, double y2, double x3, double y3, double t2, double t3)
{
	double out1 = INFINITY;
	double out2 = INFINITY;
	double out3 = INFINITY;
	double out4 = INFINITY;

	hyp_func(x_in, out1, out2, x1, y1, x2, y2, t2);
	hyp_func(x_in, out3, out4, x1, y1, x3, y3, t3);

	double closest_to_zero = INFINITY;
	double y_value = INFINITY;
	double d1, d2, d3, d4;

	#define SET_D_CLOSEST_TO_ZERO(d_num, a, b) \
		{ \
			if (out##a != INFINITY && out##b != INFINITY) \
			{ \
				d##d_num = out##a - out##b; \
				if (std::fabs(d##d_num) < std::fabs(closest_to_zero)) { \
					y_value = (out##a + out##b) / 2.0; \
					closest_to_zero = d##d_num; \
				}\
			} \
			else \
			{ \
				d##d_num = INFINITY; \
			} \
		}

	if (top)
	{
		SET_D_CLOSEST_TO_ZERO(1, 1, 3);
		SET_D_CLOSEST_TO_ZERO(2, 1, 4);
		d3 = INFINITY;
		d4 = INFINITY;
	}
	else
	{
		SET_D_CLOSEST_TO_ZERO(3, 2, 3);
		SET_D_CLOSEST_TO_ZERO(4, 2, 4);
		d1 = INFINITY;
		d2 = INFINITY;
	}

	if (prev_d == -INFINITY)
	{
		prev_d = closest_to_zero;
	}

#if DEBUG
	std::cout << "x: " << x_in << " y: " << out1 << ", " << out2 << ", " << out3 << ", " << out4 << std::endl;
	std::cout << "i: " << i << " d: " << d1 << ", " << d2 << ", " << d3 << ", " << d4 << " prev: " << prev_d << " ctz: " << closest_to_zero << std::endl;
	std::cout << std::endl;
	system("PAUSE");
#endif

	if ((prev_d < 0 && closest_to_zero > 0) || (prev_d > 0 && closest_to_zero < 0) 
		|| (prev_d != INFINITY && closest_to_zero == INFINITY) || (prev_d == INFINITY && closest_to_zero != INFINITY))
	{
		i /= -2.0;
	}

	if (double_compare(closest_to_zero, 0.0, 0.0000001))
	{
		final_x = x_in;
		final_y = y_value;
		return true;
	}
	
	prev_d = closest_to_zero;
	x_in += i;

	return false;
}

int main()
{
	/// -0.0802, -29.428
	//const double x1 = 0.4;
	//const double y1 = 0.03;
	//const double x2 = -1.28;
	//const double y2 = 0.4;
	//const double x3 = -0.9;
	//const double y3 = 1.07;

	/// 1.178, -0.5055
	//const double x1 = 1.615;
	//const double y1 = 0.1;
	//const double x2 = 0.16;
	//const double y2 = 0.0;
	//const double x3 = 0.144;
	//const double y3 = 0.96;

	/// Doesn't work
	//const double x1 = 1.74;
	//const double y1 = -0.02;
	//const double x2 = 0.16;
	//const double y2 = 0.0;
	//const double x3 = 0.144;
	//const double y3 = 0.96;

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

	const double s = 1481.0;
	const double t2 = 2.63474 * 0.0001 * s;
	const double t3 = 7.07023 * 0.0001 * s;

	// One loop going in the positive direction and another going in the negative direction, each one with an up and down curve
	double i_pos1 = 100.0;
	double i_pos2 = 100.0;
	double i_neg1 = -100.0;
	double i_neg2 = -100.0;

	double x_pos1 = 0;
	double x_pos2 = 0;
	double x_neg1 = 0;
	double x_neg2 = 0;

	double prev_d_pos1 = -INFINITY;
	double prev_d_pos2 = -INFINITY;
	double prev_d_neg1 = -INFINITY;
	double prev_d_neg2 = -INFINITY;

	double final_x = 0;
	double final_y = 0;

	while (true)
	{
		if (find_intersection(true,  x_pos1, final_x, final_y, i_pos1, prev_d_pos1, x1, y1, x2, y2, x3, y3, t2, t3)) break;
		if (find_intersection(false, x_pos2, final_x, final_y, i_pos2, prev_d_pos2, x1, y1, x2, y2, x3, y3, t2, t3)) break;
		if (find_intersection(true,  x_neg1, final_x, final_y, i_neg1, prev_d_neg1, x1, y1, x2, y2, x3, y3, t2, t3)) break;
		if (find_intersection(false, x_neg2, final_x, final_y, i_neg2, prev_d_neg2, x1, y1, x2, y2, x3, y3, t2, t3)) break;
	}

	std::cout << "Final: " << final_x << ", " << final_y << std::endl;

	system("PAUSE");

	return 0;

}
