#include <iostream>
#include <deque>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "colors.h"
#include "segments.h"
#include "filters.h"
#include "shape_matching.h"
#include "logo.h"
#include "bounding_boxes.h"


int main()
{
	std::vector<std::string> files{
			"Resources/1.jpg",
			"Resources/2.jpg",
			"Resources/3.jpg",
			"Resources/4.jpg",
			"Resources/5.jpg",
			"Resources/6.jpg",
			"Resources/7.jpg"
	};

	for (std::string filename : files)
	{
		cv::Mat image = cv::imread(filename);

		cv::Mat hsv_image = bgr2hsv(image);


		cv::Mat blue_mask = inRange(hsv_image, cv::Vec3b(80, 40, 30), cv::Vec3b(130, 255, 225));
		std::vector<Segment> blue_segments = segment_mask(blue_mask);
		blue_segments = filter_out_segments(blue_segments, 7, 5, 150, 150);
		std::sort(blue_segments.begin(), blue_segments.end(), compare_segments_by_x);

		cv::Mat red_mask_1 = inRange(hsv_image, cv::Vec3b(0, 50, 100), cv::Vec3b(15, 255, 255));
		cv::Mat red_mask_2 = inRange(hsv_image, cv::Vec3b(160, 50, 50), cv::Vec3b(179, 255, 255));
		cv::Mat red_mask = mask_or(red_mask_1, red_mask_2);
		std::vector<Segment> red_segments = segment_mask(red_mask);
		red_segments = filter_out_segments(red_segments, 5, 5, 150, 150);
		std::sort(red_segments.begin(), red_segments.end(), compare_segments_by_y);


		cv::Mat yellow_mask = inRange(hsv_image, cv::Vec3b(20, 100, 100), cv::Vec3b(30, 255, 255));
		cv::Mat yellow_mask_filtered = dilation_filter(yellow_mask, 3, 1);
		std::vector<Segment> yellow_segments = segment_mask(yellow_mask_filtered);
		yellow_segments = filter_out_segments(yellow_segments, 15, 30, 500, 500);


		std::vector<Logo> found_logos = build_logos(yellow_segments, blue_segments, red_segments);

		cv::Mat result = draw_bounding_boxes_for_logos(image, found_logos);
		cv::imwrite("out/" + filename.substr(10, filename.length() - 10), result);
	}
	
	return 0;
}