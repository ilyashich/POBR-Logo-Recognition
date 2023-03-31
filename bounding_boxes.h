#ifndef BOUNDING_BOXES_H
#define BOUNDING_BOXES_H


#include <cassert>
#include <vector>
#include <tuple>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "segments.h"
#include "logo.h"


#define BOX_COLOR cv::Vec3b(4, 255, 16);

void horizontal_line(cv::Mat& image, int y, int x_start, int x_end)
{
	assert(x_start <= x_end);
	for (int x = x_start; x <= x_end; x++) 
	{
		image.at<cv::Vec3b>(y, x) = BOX_COLOR;
		image.at<cv::Vec3b>(y + 1, x) = BOX_COLOR;
		image.at<cv::Vec3b>(y - 1, x) = BOX_COLOR;
	}
}

void vertical_line(cv::Mat& image, int x, int y_start, int y_end)
{
	assert(y_start <= y_end);
	for (int y = y_start; y <= y_end; y++) 
	{
		image.at<cv::Vec3b>(y, x) = BOX_COLOR;
		image.at<cv::Vec3b>(y, x - 1) = BOX_COLOR;
		image.at<cv::Vec3b>(y, x + 1) = BOX_COLOR;
	}
}

cv::Mat draw_bounding_boxes_for_segments(cv::Mat& image, std::vector<Segment> segments)
{
	cv::Mat result = image.clone();
	for (const auto& segment : segments)
	{
		int x_start = segment.col_min;
		int y_start = segment.row_min;
		int x_end = segment.col_max;
		int y_end = segment.row_max;
		horizontal_line(result, y_start, x_start, x_end);
		horizontal_line(result, y_end, x_start, x_end);
		vertical_line(result, x_start, y_start, y_end);
		vertical_line(result, x_end, y_start, y_end);
	}
	return result;
}

cv::Mat draw_bounding_boxes_for_logos(cv::Mat& image, std::vector<Logo> logos)
{
	cv::Mat result = image.clone();
	for (const auto& logo : logos)
	{
		double a = 0.04 * (logo.col_max - logo.col_min);
		double b = 0.04 * (logo.row_max - logo.row_min);
		
		int x_start = logo.col_min - a;
		int y_start = logo.row_min - b;
		int x_end = logo.col_max + a;
		int y_end = logo.row_max + b;
		
		horizontal_line(result, y_start, x_start, x_end);
		horizontal_line(result, y_end, x_start, x_end);
		vertical_line(result, x_start, y_start, y_end);
		vertical_line(result, x_end, y_start, y_end);
	}
	return result;
}

#endif
