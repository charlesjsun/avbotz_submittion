// Copyright (c) 2016, Charles JL Sun. All rights reserved.

#include <fstream>
#include <iostream>

int main()
{
	// Change the file name if it's different
	std::ifstream in("test.in");
	std::ofstream out("test.out");

	if (!in.is_open() || !out.is_open())
	{
		std::cout << "COULD NOT OPEN FILE" << std::endl;
		return 0;
	}

	int num_larger_than_200 = 0;

	int x, y, z;
	while (in >> x >> y >> z)
	{
		int dist = (int)std::floor(sqrt(x * x + y * y + z * z));
		out << dist << std::endl;

		if (dist > 200)
		{
			++num_larger_than_200;
		}
	}

	out << num_larger_than_200 << std::endl;

	in.close();
	out.close();

	return 0;
}
