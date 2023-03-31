#ifndef SHAPE_MATCHING_H
#define SHAPE_MATCHING_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <optional>

#include "logo.h"

struct CentralMoments
{
	long double mu00;
	long double mu01;
	long double mu10;
	long double mu11;
	long double mu20;
	long double mu02;
	long double mu21;
	long double mu12;
	long double mu30;
	long double mu03;
};

struct ScaleInvariants
{
	long double eta11;
	long double eta20;
	long double eta02;
	long double eta21;
	long double eta12;
	long double eta30;
	long double eta03;
};

struct RotationInvariants
{
	long double M1;
	long double M2;
	long double M3;
	long double M4;
	long double M5;
	long double M6;
	long double M7;
};

long double m(int i, int j, std::vector<std::pair<int, int>> pixels)
{
	long double result = 0;
	for (int k = 0; k < pixels.size(); k++)
	{
		long double x = (long double)pixels[k].second;
		long double y = (long double)pixels[k].first;
		result = result + (long double)(std::pow(x, i) * std::pow(y, j));

	}
	return result;
}

CentralMoments mu_table(std::vector<std::pair<int, int>> pixels)
{
	CentralMoments moments;
	long double m_00 = m(0, 0, pixels);
	long double m_10 = m(1, 0, pixels);
	long double m_01 = m(0, 1, pixels);
	long double m_11 = m(1, 1, pixels);
	long double m_20 = m(2, 0, pixels);
	long double m_02 = m(0, 2, pixels);
	long double m_21 = m(2, 1, pixels);
	long double m_12 = m(1, 2, pixels);
	long double m_30 = m(3, 0, pixels);
	long double m_03 = m(0, 3, pixels);
	long double x_ = m_10 / m_00;
	long double y_ = m_01 / m_00;
	moments.mu00 = m_00;
	moments.mu01 = 0.0;
	moments.mu10 = 0.0;
	moments.mu11 = m_11 - x_ * m_01;
	moments.mu20 = m_20 - x_ * m_10;
	moments.mu02 = m_02 - y_ * m_01;
	moments.mu21 = m_21 - 2.0 * x_ * m_11 - y_ * m_20 + 2.0 * x_ * x_ * m_01;
	moments.mu12 = m_12 - 2.0 * y_ * m_11 - x_ * m_02 + 2.0 * y_ * y_ * m_10;
	moments.mu30 = m_30 - 3.0 * x_ * m_20 + 2.0 * x_ * x_ * m_10;
	moments.mu03 = m_03 - 3.0 * y_ * m_02 + 2.0 * y_ * y_ * m_01;
	return moments;
}

long double calculate_eta(int i, int j, long double mu_ij, long double mu_00)
{
	long double power = (((long double)i + (long double)j + 2.0) / 2.0);
	long double eta = mu_ij / (std::pow(mu_00, power));
	return eta;
}

ScaleInvariants eta_table(std::vector<std::pair<int, int>> pixels)
{
	CentralMoments moments = mu_table(pixels);
	ScaleInvariants eta_table;
	eta_table.eta11 = calculate_eta(1, 1, moments.mu11, moments.mu00);
	eta_table.eta20 = calculate_eta(2, 0, moments.mu20, moments.mu00);
	eta_table.eta02 = calculate_eta(0, 2, moments.mu02, moments.mu00);
	eta_table.eta21 = calculate_eta(2, 1, moments.mu21, moments.mu00);
	eta_table.eta12 = calculate_eta(1, 2, moments.mu12, moments.mu00);
	eta_table.eta30 = calculate_eta(3, 0, moments.mu30, moments.mu00);
	eta_table.eta03 = calculate_eta(0, 3, moments.mu03, moments.mu00);
	return eta_table;
}

RotationInvariants hu_moments(std::vector<std::pair<int, int>> pixels)
{
	std::vector<long double> hu(7, 0.0);
	ScaleInvariants e = eta_table(pixels);
	RotationInvariants i;

	i.M1 = e.eta20 + e.eta02;
	i.M2 = std::pow((e.eta20 - e.eta02), 2) + 4.0 * std::pow(e.eta11, 2);
	i.M3 = std::pow((e.eta30 - 3.0 * e.eta12), 2) + std::pow((3.0 * e.eta21 - e.eta03), 2);
	i.M4 = std::pow((e.eta30 + e.eta12), 2) + std::pow((e.eta21 + e.eta03), 2);
	i.M5 = (e.eta30 - 3.0 * e.eta12) * (e.eta30 + e.eta12) * (std::pow((e.eta30 + e.eta12), 2) - 3.0 * std::pow((e.eta21 + e.eta03), 2)) + (3.0 * e.eta21 - e.eta03) * (e.eta21 + e.eta03) * (3.0 * std::pow((e.eta30 + e.eta12), 2) - std::pow((e.eta21 + e.eta03), 2));
	i.M6 = (e.eta20 - e.eta02) * (std::pow((e.eta30 + e.eta12), 2) - std::pow((e.eta21 + e.eta03), 2)) + 4.0 * e.eta11 * (e.eta30 + e.eta12) * (e.eta21 + e.eta03);
	i.M7 = e.eta20 * e.eta02 - e.eta11 * e.eta11;

	return i;
}

bool is_letter_l(RotationInvariants r)
{
	if (r.M1 < 0.232904 * 0.7 || r.M1 > 0.452047 * 1.2)
		return false;
	if (r.M2 < 0.001229 * 0.7 || r.M2 > 0.23478 * 1.2)
		return false;
	if (r.M3 < 0.002978 * 0.7 || r.M3 > 0.01445 * 1.2)
		return false;
	if (r.M4 < 0.000188 * 0.7 || r.M4 > 0.006359 * 1.2)
		return false;
	if (r.M5 < -0.000001 * 1.3 || r.M5 > 0.00006 * 1.2)
		return false;
	if (r.M6 < -0.000068 * 1.3 || r.M6 > 0.002069 * 1.2)
		return false;
	if (r.M7 < 0.011456 * 0.7 || r.M7 > 0.0165816 * 1.2)
		return false;
	return true;
}

bool is_letter_d(RotationInvariants r)
{
	if (r.M1 < 0.197702 * 0.8 || r.M1 > 0.237674 * 1.2)
		return false;
	if (r.M2 < 0.00011 * 0.7 || r.M2 > 0.02854 * 1.2)
		return false;
	if (r.M3 < 0.00028 * 0.8 || r.M3 > 0.0011 * 1.2)
		return false;
	if (r.M4 < -0.000001 * 1.2 || r.M4 > 0.00035 * 1.2)
		return false;
	if (r.M5 < -0.000001 * 1.2 || r.M5 > 0.003797 * 1.2)
		return false;
	if (r.M6 < -0.000001 * 1.2 || r.M6 > 0.000064 * 1.2)
		return false;
	if (r.M7 < 0.0089269 * 0.7 || r.M7 > 0.0109863 * 1.2)
		return false;
	return true;
}

bool is_letter_i(RotationInvariants r)
{
	if (r.M1 < 0.166381 * 0.8 || r.M1 > 0.19249 * 1.2)
		return false;
	if (r.M2 < 0.000004 * 0.8 || r.M2 > 0.0085 * 1.2)
		return false;
	if (r.M3 < 0.00001 * 0.8 || r.M3 > 0.000852 * 1.2)
		return false;
	if (r.M4 < -0.000001 * 1.2 || r.M4 >  0.000075 * 1.2)
		return false;
	if (r.M5 < -0.000001 * 1.2 || r.M5 > 0.000001 * 1.2)
		return false;
	if (r.M6 < -0.000007 * 1.2 || r.M6 >  0.0000013 * 1.2)
		return false;
	if (r.M7 < 0.00805614 * 0.8 || r.M7 >  0.00962518 * 1.2)
		return false;
	return true;
}

bool is_yellow_circle(RotationInvariants r)
{
	if (r.M1 < 0.210226 * 0.8 || r.M1 > 0.252874 * 1.2)
		return false;
	if (r.M2 < 0.00012 * 0.8 || r.M2 >  0.021667 * 1.2)
		return false;
	if (r.M3 < -0.000001 * 1.2 || r.M3 > 0.000058 * 1.2)
		return false;
	if (r.M4 < -0.000001 * 1.2 || r.M4 > 0.000018 * 1.2)
		return false;
	if (r.M5 < -0.000001 * 1.2 || r.M5 > 0.000001 * 1.2)
		return false;
	if (r.M6 < -0.000001 * 1.2 || r.M6 > 0.000001 * 1.2)
		return false;
	if (r.M7 < 0.00891144 * 0.8 || r.M7 > 0.0102189 * 1.2)
		return false;
	return true;
}

bool is_red_dot(RotationInvariants r)
{
	if (r.M1 < 0.159161 * 0.8 || r.M1 > 0.195348 * 1.2)
		return false;
	if (r.M2 < -0.000001 * 1.2 || r.M2 > 0.013819 * 1.2)
		return false;
	if (r.M3 < -0.000001 * 1.2 || r.M3 > 0.000105 * 1.2)
		return false;
	if (r.M4 < -0.000001 * 1.2 || r.M4 >  0.00001 * 1.2)
		return false;
	if (r.M5 < -0.000001 * 1.2 || r.M5 > 0.000001 * 1.2)
		return false;
	if (r.M6 < -0.000001 * 1.2 || r.M6 > 0.000001 * 1.2)
		return false;
	if (r.M7 < 0.0062891 * 0.8 || r.M7 > 0.00634452 * 1.2)
		return false;
	return true;
}

bool is_i_with_dot(RotationInvariants r)
{
	if (r.M1 < 0.2564 * 0.8 || r.M1 > 0.46 * 1.2)
		return false;
	if (r.M2 < 0.0214 * 0.8 || r.M2 > 0.1525 * 1.2)
		return false;
	if (r.M3 < 0.00379 * 0.8 || r.M3 > 0.02201 * 1.2)
		return false;
	if (r.M4 < 0.00081339 * 0.8 || r.M4 >  0.01332 * 1.2)
		return false;
	if (r.M5 < 0.00000134 * 0.8 || r.M5 > 0.0002276 * 1.2)
		return false;
	if (r.M6 < 0.0001176 * 0.8 || r.M6 > 0.005164 * 1.2)
		return false;
	if (r.M7 < 0.011 * 0.8 || r.M7 > 0.0153 * 1.2)
		return false;
	return true;
}
bool is_correct_logo(std::vector<Segment> blue_segments, std::vector<Segment> red_segments) 
{
	int min_y = std::min(std::min(blue_segments[0].row_min, blue_segments[1].row_min), blue_segments[2].row_min);
	int max_y = std::max(std::max(blue_segments[0].row_min, blue_segments[1].row_min), blue_segments[2].row_min);
	if (max_y - min_y > 30)
		return false;
	return (
		blue_segments.size() == 3 && red_segments.size() == 2 &&
		blue_segments[0].type == Letter_L &&
		blue_segments[1].type == Letter_D &&
		blue_segments[2].type == Letter_L &&
		red_segments[0].type == Red_Dot &&
		red_segments[1].type == Letter_I &&
		blue_segments[0].col_max <= blue_segments[1].col_min &&
		blue_segments[1].col_max <= blue_segments[2].col_min &&
		red_segments[0].col_min >= red_segments[1].col_min &&
		red_segments[0].col_max <= red_segments[1].col_max &&
		red_segments[1].col_min >= blue_segments[0].col_max &&
		red_segments[1].col_min <= blue_segments[1].col_min) ||
		(
			blue_segments.size() == 3 && red_segments.size() == 1 &&
			blue_segments[0].type == Letter_L &&
			blue_segments[1].type == Letter_D &&
			blue_segments[2].type == Letter_L &&
			red_segments[0].type == Letter_I_With_Dot &&
			blue_segments[0].col_max <= blue_segments[1].col_min &&
			blue_segments[1].col_max <= blue_segments[2].col_min &&
			red_segments[0].col_min >= blue_segments[0].col_max &&
			red_segments[0].col_min <= blue_segments[1].col_min
		);
}

std::optional<Logo> build_logo(Segment yellow_segment, std::vector<Segment> blue_segments, std::vector<Segment> red_segments)
{
	std::vector<Segment> matched_blue_segments;
	std::vector<Segment> matched_red_segments;
	bool yellow_circle = false;
	if (is_yellow_circle(hu_moments(yellow_segment.pixel_coordinates)))
	{
		yellow_circle = true;
		for (auto& blue_segment : blue_segments)
		{
			if (yellow_segment.contains(blue_segment))
			{
				RotationInvariants invariants = hu_moments(blue_segment.pixel_coordinates);
				if (is_letter_l(invariants))
				{
					blue_segment.type = Letter_L;
					matched_blue_segments.push_back(blue_segment);
				}
				if (is_letter_d(invariants))
				{
					blue_segment.type = Letter_D;
					matched_blue_segments.push_back(blue_segment);
				}

			}
		}

		for (auto& red_segment : red_segments)
		{
			if (yellow_segment.contains(red_segment))
			{
				RotationInvariants invariants = hu_moments(red_segment.pixel_coordinates);
				if (is_red_dot(invariants))
				{
					red_segment.type = Red_Dot;
					matched_red_segments.push_back(red_segment);
				}
				if (is_letter_i(invariants))
				{
					red_segment.type = Letter_I;
					matched_red_segments.push_back(red_segment);
				}
				if (is_i_with_dot(invariants))
				{
					red_segment.type = Letter_I_With_Dot;
					matched_red_segments.push_back(red_segment);
				}

			}
		}
	}

	if (yellow_circle && is_correct_logo(matched_blue_segments, matched_red_segments))
	{
		return Logo{
			yellow_segment.row_min,
			yellow_segment.row_max,
			yellow_segment.col_min,
			yellow_segment.col_max,
			yellow_segment,
			matched_blue_segments,
			matched_red_segments
		};
	}
	return std::nullopt;
}

std::vector<Logo> build_logos(std::vector<Segment> yellow_segments, std::vector<Segment> blue_segments, std::vector<Segment> red_segments)
{
	std::vector<Logo> logos;
	for (const auto& yellow_segment : yellow_segments)
	{
		auto logo = build_logo(yellow_segment, blue_segments, red_segments);
		if (logo)
		{
			logos.push_back(*logo);
		}
	}
	return logos;
}

#endif
