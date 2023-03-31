#ifndef LOGO_H
#define LOGO_H

#include "segments.h"

struct Logo
{
	int row_min;
	int row_max;
	int col_min;
	int col_max;
	Segment yellow_segment;
	std::vector<Segment> blue_segments;
	std::vector<Segment> red_segments;
};

#endif