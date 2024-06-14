#pragma once
#include "libmorton/morton.h"
#include "radix_spline.h"
#include "common.h"
#include <vector>
#include <iostream>
namespace gf
{
  template <class KeyType>
  class GridFile
  {
  public:
    GridFile() = default;
    // GridFile(std::vector<KeyType> min_keys, std::vector<KeyType> max_keys, int dim, double grid_size)
    //     : min_keys_(min_keys), max_keys_(max_keys), dim_(dim), grid_size_(grid_size) {}
    GridFile(KeyType min_key, KeyType max_key, size_t num_radix_bits, size_t num_shift_bits, size_t group_size, size_t cell_size,
             std::vector<uint32_t> radix_table,
             std::vector<KeyType> spline_points,
             std::vector<rs::RadixSpline<KeyType>> rss)
            //  std::vector<MLP> mlps)
        : min_key_(min_key),
          max_key_(max_key),
          num_radix_bits_(num_radix_bits),
          num_shift_bits_(num_shift_bits),
          group_size_(group_size),
          cell_size_(cell_size),
          dim_(2),
          radix_table_(std::move(radix_table)),
          spline_points_(std::move(spline_points)),
          rss_(std::move(rss)){}
          // mlps_(std::move(mlps)) {}

    uint_fast64_t GetGridIndex(std::vector<KeyType> key)
    {
      std::vector<uint_fast32_t> coord = GetCoord(key);
      // std::cout << coord[0] << std::endl;
      // std::cout << coord[1] << std::endl;
      if (dim_ == 2)
      {
        return libmorton::morton2D_64_encode(coord[0], coord[1]);
      }
      return 0;
    }

    std::vector<uint_fast64_t> GetGridIndex(std::vector<KeyType> begin, std::vector<KeyType> end)
    {
      std::vector<uint_fast32_t> begin_coord = GetCoord(begin);
      std::vector<uint_fast32_t> end_coord = GetCoord(end);
      std::vector<uint_fast64_t> grid_indexes;
      if (dim_ == 2)
      {
        uint_fast64_t begin_zaddr = libmorton::morton2D_64_encode(begin_coord[0], begin_coord[1]);
        uint_fast64_t end_zaddr = libmorton::morton2D_64_encode(end_coord[0], end_coord[1]);
        for (uint_fast64_t zaddr = begin_zaddr; zaddr <= end_zaddr; zaddr++)
        {
          if (!IsRelevant2D(begin_coord, end_coord, zaddr))
          {
            zaddr = NextJumpIn2D(begin_zaddr, end_zaddr, zaddr);
          }
          grid_indexes.push_back(zaddr);
        }
        return grid_indexes;
      }
      return std::vector<uint_fast64_t>();
    }

    // Returns the size in bytes.
    size_t GetSizeInByte() const
    {
      size_t rs_size = 0;
      int cnt=0;
      for (rs::RadixSpline<KeyType> rs : rss_)
      {
        rs_size += rs.GetSize();
        cnt++;
      }
      return sizeof(*this) + radix_table_.size() * sizeof(uint32_t) + spline_points_.size() * sizeof(KeyType) + rs_size;
    }

    std::vector<uint_fast32_t> GetCoord(std::vector<KeyType> key)
    {
      if (key[0] < min_key_)
        return {0, 0};
      if (key[0] > max_key_)
        return {spline_points_.size() - 1, 0};
      const size_t index = GetSplineSegment(key[0]);
      
      // torch::Tensor input = torch::tensor({float(key[1])});
      // input = torch::unsqueeze(input, 1);
      // // std::cout << input << std::endl;
      // double ret = mlps_[index].forward(input).item().toDouble();
      // if (ret < 0)
      //   return {index, 0};
      // if (ret >= (group_size_ - 1) / cell_size_)
      //   return {index, (group_size_ - 1) / cell_size_};
      // return {index, ret};

      // double ret = rss_[index].GetEstimatedPosition(key[1]);
      // if (ret < 0)
      //   return {index, 0};
      // if (ret >= group_size_)
      //   return {index, (group_size_ - 1) / cell_size_};
      // return {index, ret / cell_size_};

      double ret = rss_[index].GetEstimatedPosition(key[1]);
      return {index, ret};
    }

    uint_fast32_t GetCoordX(KeyType keyx)
    {
      if (keyx < min_key_)
        return 0;
      if (keyx > max_key_)
        return spline_points_.size() - 1;
      return GetSplineSegment(keyx);
    }

    uint_fast32_t GetCoordY(uint_fast32_t x,KeyType keyy)
    {
      double ret = rss_[x].GetEstimatedPosition(keyy);
      return ret;
    }

  private:
    // std::vector<uint_fast32_t> GetCoord(std::vector<KeyType> key)
    // {
    //     std::vector<uint_fast32_t> coord(dim_);
    //     for (int i = 0; i < dim_; i++)
    //     {
    //         coord[i] = (key[i] - min_keys_[i]) / grid_size_;
    //     }
    //     return coord;
    // }
    
    size_t GetSplineSegment(const KeyType key) const
    {
      // Narrow search range using radix table.
      const KeyType prefix = (key - min_key_) >> num_shift_bits_;
      assert(prefix + 1 < radix_table_.size());
      const uint32_t begin = radix_table_[prefix];
      const uint32_t end = radix_table_[prefix + 1];

      if (end - begin < 32)
      {
        // Do linear search over narrowed range.
        uint32_t current = begin;
        while (spline_points_[current] < key)
          ++current;
        return current;
      }

      // Do binary search over narrowed range.
      const auto lb = std::lower_bound(spline_points_.begin() + begin,
                                       spline_points_.begin() + end,
                                       key);
      return std::distance(spline_points_.begin(), lb);
    }
    bool IsRelevant2D(std::vector<uint_fast32_t> begin, std::vector<uint_fast32_t> end, uint64_t testad)
    {
      uint64_t testx, testy;
      libmorton::morton2D_64_decode(testad, testx, testy);
      return testx >= begin[0] && testx <= end[0] && testy >= begin[1] && testy <= end[1];
    }
    uint_fast64_t NextJumpIn2DSimple(std::vector<uint_fast32_t> begin, std::vector<uint_fast32_t> end, uint_fast64_t from_zaddr)
    {
      uint_fast64_t to_zaddr;
      uint64_t end_zaddr = libmorton::morton2D_64_encode(end[0], end[1]);
      for (to_zaddr = from_zaddr + 1; to_zaddr <= end_zaddr; to_zaddr++)
      {
        if (IsRelevant2D(begin, end, to_zaddr))
        {
          break;
        }
      }
      return to_zaddr;
    }
    uint_fast64_t NextJumpIn2D(uint_fast64_t begin_zaddr, uint_fast64_t end_zaddr, uint_fast64_t from_zaddr)
    {
      uint_fast64_t min_v, max_v, z, bit_position, load_mask, load_ones;
      uint_fast64_t bigmin;
      min_v = begin_zaddr;
      max_v = end_zaddr;
      z = from_zaddr;
      bit_position = bit_position_init_;
      load_mask = load_mask_init_;
      load_ones = load_ones_init_;
      while (bit_position)
      {
        uint_fast64_t z_bit, min_bit, max_bit;
        z_bit = z & bit_position;
        min_bit = min_v & bit_position;
        max_bit = max_v & bit_position;
        int op = (z_bit != 0) << 2 | (min_bit != 0) << 1 | (max_bit != 0);
        switch (op)
        {
        case 0:
          break;
        case 1:
          bigmin = min_v & load_mask | bit_position;
          max_v = max_v & load_mask | load_ones;
          break;
        case 3:
          return min_v;
          break;
        case 4:
          return bigmin;
          break;
        case 5:
          min_v = min_v & load_mask | bit_position;
          break;
        case 7:
          break;
        default:
          std::cerr << "Z-order index search failed. Something wrong..." << std::endl;
          break;
        }
        bit_position >>= 1;
        load_ones >>= 1;
        load_mask >>= 1;
        load_mask |= bit_position_init_;
      }
      return bigmin;
    }
    // std::vector<KeyType> min_keys_;
    // std::vector<KeyType> max_keys_;
    // double grid_size_;

    int dim_;

    KeyType min_key_;
    KeyType max_key_;
    size_t num_radix_bits_;
    size_t num_shift_bits_;
    size_t group_size_;
    size_t cell_size_;

    std::vector<uint32_t> radix_table_;
    std::vector<KeyType> spline_points_;
    std::vector<rs::RadixSpline<KeyType>> rss_;

    // std::vector<MLP> mlps_;

    static const uint_fast64_t bit_position_init_ = 0x8000000000000000;
    static const uint_fast64_t load_mask_init_ = 0x5555555555555555;
    static const uint_fast64_t load_ones_init_ = 0x2aaaaaaaaaaaaaaa;
  };
}