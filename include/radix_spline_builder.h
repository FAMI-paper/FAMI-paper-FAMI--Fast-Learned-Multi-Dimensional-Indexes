#pragma once

#include <cassert>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include "common.h"
#include "radix_spline.h"

namespace rs
{
	double error_cnt = 0;
	template <class KeyType>
	class Builder
	{
	public:
		Builder(KeyType min_key, KeyType max_key, size_t group_size, size_t cell_size, size_t num_radix_bits = 18)
			: min_key_(min_key),
			  max_key_(max_key),
			  group_size_(group_size / cell_size),
			  cell_size_(cell_size),
			  num_radix_bits_(num_radix_bits),
			  num_shift_bits_(GetNumShiftBits(max_key - min_key, num_radix_bits_)),
			  curr_num_keys_(0),
			  prev_prefix_(0)
		{
			// Initialize radix table, needs to contain all prefixes up to the largest key + 1.
			const uint32_t max_prefix = (max_key - min_key) >> num_shift_bits_;
			radix_table_.resize(max_prefix + 2, 0);
		}

		void AddKey(KeyType key)
		{
			group_points_.push_back({key, curr_num_keys_});
			if (group_points_.size() == group_size_)
			{
				DoFitting();
				group_points_.clear();
			}
			curr_num_keys_++;
		}

		// Finalizes the construction and returns a read-only `GridFile`.
		RadixSpline<KeyType> Finalize()
		{
			if (group_points_.size() != 0)
			{
				DoFitting();
				group_points_.clear();
			}

			// Maybe even size the radix based on max key right from the start
			FinalizeRadixTable();
			return RadixSpline<KeyType>(min_key_,
										max_key_,
										curr_num_keys_,
										num_radix_bits_,
										num_shift_bits_,
										group_size_,
										std::move(radix_table_),
										std::move(spline_points_),
										std::move(fit_params_));
		}

	private:
		// Returns the number of shift bits based on the `diff` between the largest and the smallest key.
		// KeyType == uint32_t.
		static size_t GetNumShiftBits(uint32_t diff, size_t num_radix_bits)
		{
			const uint32_t clz = __builtin_clz(diff);
			if ((32 - clz) < num_radix_bits)
				return 0;
			return 32 - num_radix_bits - clz;
		}
		// KeyType == uint64_t.
		static size_t GetNumShiftBits(uint64_t diff, size_t num_radix_bits)
		{
			const uint32_t clzl = __builtin_clzl(diff);
			if ((64 - clzl) < num_radix_bits)
				return 0;
			return 64 - num_radix_bits - clzl;
		}

		void AddKeyToSpline(KeyType key)
		{
			spline_points_.push_back(key);
			PossiblyAddKeyToRadixTable(key);
		}

		void PossiblyAddKeyToRadixTable(KeyType key)
		{
			const KeyType curr_prefix = (key - min_key_) >> num_shift_bits_;
			if (curr_prefix != prev_prefix_)
			{
				const uint32_t curr_index = spline_points_.size() - 1;
				for (KeyType prefix = prev_prefix_ + 1; prefix <= curr_prefix; ++prefix)
					radix_table_[prefix] = curr_index;
				prev_prefix_ = curr_prefix;
			}
		}

		void FinalizeRadixTable()
		{
			++prev_prefix_;
			const uint32_t num_spline_points = spline_points_.size();
			for (; prev_prefix_ < radix_table_.size(); ++prev_prefix_)
				radix_table_[prev_prefix_] = num_spline_points;
		}

		double ComputeAverageRes(double slope, double x_mean, double y_mean)
		{
			double residual = 0;
			double y_pred, key, pos;
			int n = group_points_.size();
			for (size_t i = 0; i < n; i++)
			{
				key = group_points_[i].x;
				pos = group_points_[i].y;
				y_pred = std::fma(key - x_mean, slope, y_mean);
				residual += abs(int(pos / cell_size_) - int(y_pred / cell_size_));
				if (int(pos / cell_size_) != int(y_pred / cell_size_))
				{
					error_cnt++;
				}
			}
			return residual / n;
		}

		void DoFitting()
		{
			int n = group_points_.size();
			double sum_x = 0, sum_xx = 0, sum_y = 0, sum_xy = 0;
			for (auto point : group_points_)
			{
				sum_x += (double)point.x;
				sum_y += (double)point.y;
				sum_xy += (double)point.x * point.y;
				sum_xx += (double)point.x * point.x;
			}
			const double slope = (sum_x * sum_y - n * sum_xy) / (sum_x * sum_x - n * sum_xx);
			const double intercept = sum_y / n - slope * sum_x / n;
			ComputeAverageRes(slope, sum_x / n, sum_y / n);
			fit_params_.push_back({slope, intercept});
			AddKeyToSpline(group_points_.back().x);
		}

		const KeyType min_key_;
		const KeyType max_key_;
		const size_t num_radix_bits_;
		const size_t num_shift_bits_;
		std::vector<uint32_t> radix_table_;
		std::vector<KeyType> spline_points_;

		// Vector of model parameters
		std::vector<FitParam> fit_params_;

		size_t curr_num_keys_;
		KeyType prev_prefix_;

		// Number of CDF points in one group
		size_t group_size_;

		size_t cell_size_;

		// vector of a CDF points group
		std::vector<Coord<KeyType>> group_points_;
	};
} // namespace rs