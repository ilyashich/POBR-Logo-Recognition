#ifndef COLORS_H
#define COLORS_H

#include <cassert>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

cv::Vec3b pixel_bgr2hsv(cv::Vec3b& bgr_pixel)
{
    double r = bgr_pixel[2];
    double g = bgr_pixel[1];
    double b = bgr_pixel[0];
    double Cmax = std::max(std::max(r, g), b);
    double Cmin = std::min(std::min(r, g), b);
    double delta = Cmax - Cmin;
    double h = 0.0;
    double s = 0.0;
    double v = 0.0;

    if (Cmax == Cmin)
    {
        h = 0.0;
    }
    else if (Cmax == r && g >= b)
    {
        h = fmod((60.0 * ((g - b) / delta)), 360.0);
    }		
    else if (Cmax == r && g < b)
    {
        h = fmod((60.0 * ((g - b) / delta)) + 360.0, 360.0);
    }
    else if (Cmax == g)
    {
        h = fmod((60.0 * ((b - r) / delta)) + 120.0, 360.0);
    }
    else if (Cmax == b)
    {
        h = fmod((60.0 * ((r - g) / delta)) + 240.0, 360.0);
    }

    if (Cmax == 0)
    {
        s = 0;
    }
    else
    {
        s = (delta * 255) / Cmax;
    }

    v = Cmax;
    h = h / 2.0;
    return cv::Vec3b(h, s, v);
}

cv::Mat bgr2hsv(cv::Mat& image)
{
    cv::Mat hsv = cv::Mat::zeros(image.rows, image.cols, image.type());
    for (int i = 0; i < image.rows; i++)
    {
        for (int j = 0; j < image.cols; j++)
        {
            cv::Vec3b bgr_pixel = image.at<cv::Vec3b>(i, j);
            hsv.at<cv::Vec3b>(i, j) = pixel_bgr2hsv(bgr_pixel);
        }
    }
    return hsv;
}

bool inRangeInner(cv::Vec3b pixel, cv::Vec3b lower, cv::Vec3b upper)
{
    return pixel[0] >= lower[0] && pixel[0] <= upper[0] &&
        pixel[1] >= lower[1] && pixel[1] <= upper[1] &&
        pixel[2] >= lower[2] && pixel[2] <= upper[2];
}

cv::Mat inRange(cv::Mat& image, cv::Vec3b lower, cv::Vec3b upper)
{
    cv::Mat result = cv::Mat::zeros(image.rows, image.cols, image.type());
    for (int i = 0; i < image.rows; i++)
    {
        for (int j = 0; j < image.cols; j++)
        {
            if (inRangeInner(image.at<cv::Vec3b>(i, j), lower, upper))
                result.at<cv::Vec3b>(i, j) = cv::Vec3b(255, 255, 255);
        }
    }
    return result;
}

cv::Mat mask_or(cv::Mat& mask1, cv::Mat& mask2)
{
    assert(mask1.rows == mask2.rows && mask1.cols == mask2.cols);
    cv::Mat result = cv::Mat::zeros(mask1.rows, mask1.cols, mask1.type());
    cv::Vec3b zeroes(0, 0, 0);
    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.cols; j++)
        {
            if (mask1.at<cv::Vec3b>(i, j) != zeroes || mask2.at<cv::Vec3b>(i, j) != zeroes)
                result.at<cv::Vec3b>(i, j) = cv::Vec3b(255, 255, 255);
        }
    }
    return result;
}

cv::Mat mask_and(cv::Mat& mask1, cv::Mat& mask2)
{
    assert(mask1.rows == mask2.rows && mask1.cols == mask2.cols);
    cv::Mat result = cv::Mat::zeros(mask1.rows, mask1.cols, mask1.type());
    cv::Vec3b zeroes(0, 0, 0);
    for (int i = 0; i < result.rows; i++)
    {
        for (int j = 0; j < result.cols; j++)
        {
            if (mask1.at<cv::Vec3b>(i, j) != zeroes && mask2.at<cv::Vec3b>(i, j) != zeroes)
                result.at<cv::Vec3b>(i, j) = cv::Vec3b(255, 255, 255);
        }
    }
    return result;
}

#endif