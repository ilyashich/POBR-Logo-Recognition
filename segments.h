#ifndef SEGMENTS_H
#define SEGMENTS_H

#include <vector>
#include <deque>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

enum SegmentType
{
	Undefined,
	Letter_L,
	Letter_D,
	Letter_I,
	Red_Dot,
	Letter_I_With_Dot
};


struct Segment
{
	int row_min;
	int row_max;
	int col_min;
	int col_max;
	std::vector<std::pair<int, int>> pixel_coordinates;
	std::vector<std::pair<int, int>> border_pixel_coordinates;
	SegmentType type;

	int get_width()
	{
		return col_max - col_min + 1;
	}

	int get_height()
	{
		return row_max - row_min + 1;
	}

	bool contains(Segment other)
	{
		return (other.row_min > row_min && other.row_max < row_max&& other.col_min > col_min && other.col_max < col_max);
	}
};

bool compare_segments_by_x(const Segment& a, const Segment& b)
{
	return a.col_max < b.col_min;
}

bool compare_segments_by_y(const Segment& a, const Segment& b)
{
	return a.row_max < b.row_min;
}

bool is_border_color(cv::Mat& mask, int y, int x)
{
	if (x == 0 || x == mask.cols - 1 || y == 0 || y == mask.rows - 1)
	{
		return true;
	}
	return mask.at<cv::Vec3b>(y - 1, x) == cv::Vec3b(0, 0, 0) ||
		mask.at<cv::Vec3b>(y + 1, x) == cv::Vec3b(0, 0, 0) ||
		mask.at<cv::Vec3b>(y, x - 1) == cv::Vec3b(0, 0, 0) ||
		mask.at<cv::Vec3b>(y, x + 1) == cv::Vec3b(0, 0, 0);
}

Segment flood_fill_segment(cv::Mat mask, std::pair<int, int> seed)
{
	int row_min = mask.rows;
	int row_max = 0;
	int col_min = mask.cols;
	int col_max = 0;
	std::vector<std::pair<int, int>> pixel_coordinates;
	std::vector<std::pair<int, int>> border_pixel_coordinates;
	std::deque<std::pair<int, int>> stack;

	stack.push_back(seed);
	while (!stack.empty())
	{
		std::pair<int, int> current_pixel_coords = stack.front();
		stack.pop_front();
		if (current_pixel_coords.first >= 0 && current_pixel_coords.second >= 0 && current_pixel_coords.first < mask.rows && current_pixel_coords.second < mask.cols)
		{

			cv::Vec3b pixel = mask.at<cv::Vec3b>(current_pixel_coords.first, current_pixel_coords.second);
			if (pixel != cv::Vec3b(0, 0, 0))
			{
				if (current_pixel_coords.first < row_min)
				{
					row_min = current_pixel_coords.first;
				}
				if (current_pixel_coords.first > row_max)
				{
					row_max = current_pixel_coords.first;
				}
				if (current_pixel_coords.second < col_min)
				{
					col_min = current_pixel_coords.second;
				}
				if (current_pixel_coords.second > col_max)
				{
					col_max = current_pixel_coords.second;
				}
				pixel_coordinates.push_back(current_pixel_coords);
				mask.at<cv::Vec3b>(current_pixel_coords.first, current_pixel_coords.second) = cv::Vec3b(0, 0, 0);

				if (is_border_color(mask, current_pixel_coords.first, current_pixel_coords.second))
				{
					border_pixel_coordinates.push_back(current_pixel_coords);
				}

				stack.push_back(std::make_pair(current_pixel_coords.first - 1, current_pixel_coords.second));
				stack.push_back(std::make_pair(current_pixel_coords.first + 1, current_pixel_coords.second));
				stack.push_back(std::make_pair(current_pixel_coords.first, current_pixel_coords.second - 1));
				stack.push_back(std::make_pair(current_pixel_coords.first, current_pixel_coords.second + 1));

			}
		}
	}

	return Segment{ row_min, row_max, col_min, col_max, pixel_coordinates, border_pixel_coordinates, Undefined };
}

std::vector<Segment> segment_mask(cv::Mat image)
{
	cv::Mat mask = image.clone();
	std::vector<Segment> segments;
	for (int i = 0; i < mask.rows; i++)
	{
		for (int j = 0; j < mask.cols; j++)
		{
			if (mask.at<cv::Vec3b>(i, j) != cv::Vec3b(0, 0, 0))
			{
				Segment segment = flood_fill_segment(mask, std::make_pair(i, j));
				segments.push_back(segment);
			}
		}
	}
	return segments;
}

std::vector<Segment> filter_out_segments(std::vector<Segment> segments, int min_height, int min_width, int max_height, int max_width)
{
	std::vector<Segment> filtered_segments;
	for (Segment segment : segments)
	{
		if (segment.get_height() >= min_height && segment.get_height() <= max_height && segment.get_width() >= min_width && segment.get_width() <= max_width)
		{
			filtered_segments.push_back(segment);
		}
	}
	return filtered_segments;
}

#endif