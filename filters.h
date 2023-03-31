#ifndef FILTERS_H
#define FILTERS_H

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

enum FilterType
{
	Erosion,
	Median,
	Dilation
};

int calculate_operation(std::vector<int> row, FilterType type)
{
	switch (type)
	{
	case Erosion:
		return *std::min_element(row.begin(), row.end());
	case Median:
		std::sort(row.begin(), row.end());
		return row[row.size() / 2];
	case Dilation:
		return *std::max_element(row.begin(), row.end());
	default:
		return 0;
	}
}

cv::Mat rank_filter(cv::Mat& src, int filter_size, FilterType type)
{
	cv::Mat dst = cv::Mat::zeros(src.rows, src.cols, src.type());
	int offset = filter_size / 2;
	std::vector<std::vector<int>> acc(3, std::vector<int>(filter_size * filter_size, 0));
	for (int x = offset; x < src.rows - offset; ++x)
	{
		for (int y = offset; y < src.cols - offset; ++y)
		{
			for (int a = 0; a < filter_size; ++a)
			{
				for (int b = 0; b < filter_size; ++b)
				{
					int xn = x + a - offset;
					int yn = y + b - offset;
					cv::Vec3b p = src.at<cv::Vec3b>(xn, yn);
					acc[0][a * filter_size + b] = p[0];
					acc[1][a * filter_size + b] = p[1];
					acc[2][a * filter_size + b] = p[2];
				}

			}
			dst.at<cv::Vec3b>(x, y)[0] = calculate_operation(acc[0], type);
			dst.at<cv::Vec3b>(x, y)[1] = calculate_operation(acc[1], type);
			dst.at<cv::Vec3b>(x, y)[2] = calculate_operation(acc[2], type);
		}
	}
	return dst;
}

cv::Mat erosion_filter(cv::Mat& src, int filter_size, int num_iter)
{
	cv::Mat result = src.clone();
	for (int i = 0; i < num_iter; ++i)
	{
		result = rank_filter(result, filter_size, Erosion);
	}
	return result;
}

cv::Mat dilation_filter(cv::Mat& src, int filter_size, int num_iter)
{
	cv::Mat result = src.clone();
	for (int i = 0; i < num_iter; ++i)
	{
		result = rank_filter(result, filter_size, Dilation);
	}
	return result;
}

#endif
