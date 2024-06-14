#pragma once

#include <cassert>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
// #include "MLP.h"
#include "common.h"
#include "grid_file.h"
#include "radix_spline.h"
#include "radix_spline_builder.h"

namespace gf
{
  int error_cnt = 0;
  template <class KeyType>
  class Builder
  {
  public:
    Builder(KeyType min_key, KeyType max_key, size_t group_size,size_t rs_group_size, size_t num_radix_bits = 18, size_t cell_size = 100)
        : min_key_(min_key),
          max_key_(max_key),
          group_size_(group_size),
          rs_group_size_(rs_group_size),
          num_radix_bits_(num_radix_bits - log(group_size / cell_size + 1) / log(2)),
          num_shift_bits_(GetNumShiftBits(max_key - min_key, num_radix_bits_)),
          cell_size_(cell_size),
          prev_prefix_(0)
    {
      // Initialize radix table, needs to contain all prefixes up to the largest key + 1.
      const uint32_t max_prefix = (max_key - min_key) >> num_shift_bits_;
      radix_table_.resize(max_prefix + 2, 0);
    }

    void AddKey(KeyType key)
    {
      group_points_.push_back(key);
      if (group_points_.size() == group_size_)
      {
        DoFitting();
        group_points_.clear();
      }
    }

    void AddKeyToSpline(KeyType key)
    {
      spline_points_.push_back(key);
      PossiblyAddKeyToRadixTable(key);
    }

    // Finalizes the construction and returns a read-only `GridFile`.
    GridFile<KeyType> Finalize()
    {
      if (group_points_.size() != 0)
      {
        DoFitting();
        group_points_.clear();
      }

      // Maybe even size the radix based on max key right from the start
      FinalizeRadixTable();
      return GridFile<KeyType>(min_key_,
                               max_key_,
                               num_radix_bits_,
                               num_shift_bits_,
                               group_size_,
                               cell_size_,
                               std::move(radix_table_),
                               std::move(spline_points_),
                               std::move(rss_));
                              //  std::move(mlps_));
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

    void DoFitting()
    {
      std::sort(group_points_.begin(), group_points_.end());
      int n = group_points_.size();

      // int batch_size = 10000;
      // std::vector<float> index;
      // std::vector<float> keys;
      // int max_index = (n - 1) / cell_size_;
      // int diff = group_points_.back() - group_points_.front();
      // for (size_t i = 0; i < n; i++)
      // {
      //   // keys.push_back(group_points_[i]);
      //   keys.push_back((group_points_[i] - group_points_.front()) / 1.0 / diff);
      //   index.push_back(i / cell_size_);
      // }
      // // torch::Tensor input = torch::tensor(keys);
      // // input = torch::unsqueeze(input, 1);
      // // torch::Tensor target = torch::tensor(index);
      // // target = torch::unsqueeze(target, 1);
      // auto dataset = myDataset(keys, index).map(torch::data::transforms::Stack<>());
      // auto dataloader = torch::data::make_data_loader<torch::data::samplers::RandomSampler>(std::move(dataset), batch_size);
      // auto mlp = MLP(1, 1);
      // bool suc = false;
      // float learning_rate = 0.01;
      // float pre_loss_test = 0;
      // for (int epoch = 0; epoch < 5000; epoch++)
      // {
      //   int batch_index = 0;
      //   float acc_train = 0;
      //   float loss_train = 0;
      //   float acc_test = 0;
      //   float loss_test = 0;
      //   std::cout<<learning_rate<<std::endl;
      //   torch::optim::Adam optimizer_mlp(mlp.parameters(), learning_rate);
      //   mlp.train();
      //   for (auto batch : *dataloader)
      //   {
      //     optimizer_mlp.zero_grad();
      //     auto data = batch.data;
      //     auto target = batch.target;
      //     // std::cout << data << std::endl;
      //     // std::cout << target<< std::endl;
      //     auto prediction = mlp.forward(data);
      //     auto loss = torch::mse_loss(prediction, target);
      //     loss_train += loss.item<float>();
      //     auto acc = ((prediction + 0.5).toType(torch::ScalarType::Int) == target.toType(torch::ScalarType::Int)).sum();
      //     // if (epoch == 32)
      //     // {
      //     //   std::cout << prediction << std::endl;
      //     //   std::cout << target << std::endl;
      //     //   std::cout << (prediction + 0.5).toType(torch::ScalarType::Int) << std::endl;
      //     //   std::cout << target.toType(torch::ScalarType::Int) << std::endl;
      //     // }
      //     acc_train += acc.template item<float>() / data.size(0);
      //     loss.backward();
      //     optimizer_mlp.step();
      //     batch_index++;
      //   }
      //   // std::cout << "Epoch: " << epoch << " |Loss_train: " << loss_train / batch_index << " |Acc_train:" << acc_train / batch_index << std::endl;
      //   batch_index = 0;
      //   mlp.eval();
      //   for (auto batch : *dataloader)
      //   {
      //     auto data = batch.data;
      //     auto target = batch.target;
      //     auto prediction = mlp.forward(data);
      //     auto loss = torch::mse_loss(prediction, target);
      //     loss_test += loss.item<float>();
      //     auto acc = ((prediction + 0.5).toType(torch::ScalarType::Int) == target.toType(torch::ScalarType::Int)).sum();
      //     acc_test += acc.template item<float>() / data.size(0);
      //     batch_index++;
      //   }
      //   std::cout << "Epoch: " << epoch << " |Loss_test: " << loss_test / batch_index << " |Acc_test:" << acc_test / batch_index << std::endl;
      //   if (acc_test / batch_index >= 0.7)
      //   {
      //     suc = true;
      //     std::cout << std::endl;
      //     break;
      //   }
      //   if (epoch != 0 && loss_test > pre_loss_test)
      //   {
      //     if(rand()%10==0)
      //     {
      //       learning_rate /= 5;
      //     }
      //   }
      //   pre_loss_test = loss_test;
      // }
      // if (!suc)
      // {
      //   std::cerr << "train failed" << std::endl;
      // }
      // mlp.eval();
      // mlps_.push_back(mlp);
      // // auto out = mlp.forward(input);
      // // std::cout << out << std::endl;
      // // auto out2 = mlp.forward(torch::unsqueeze(input[20513], 1));
      // // std::cout << out2 << std::endl;

      rs::Builder<KeyType> rsb(group_points_.front(), group_points_.back(), rs_group_size_, cell_size_, num_radix_bits_);
      
      // for (KeyType key : group_points_)
      // {
      //   rsb.AddKey(key);
      // }

      for (size_t i = 0; i < n; i++)
      {
        if (i % cell_size_ == 0)
        {
          rsb.AddKey(group_points_[i]);
        }
      }
      rs::RadixSpline<KeyType> rs = rsb.Finalize();
      rss_.push_back(rs);
      double residual = 0;
      double y_pred, key, pos;
      for (size_t i = 0; i < n; i++)
      {
        key = group_points_[i];
        pos = i;
        y_pred = rs.GetEstimatedPosition(key);
        if (int(pos / cell_size_) != int(y_pred))
        {
          error_cnt++;
        }
      }
    }

    const KeyType min_key_;
    const KeyType max_key_;
    const size_t num_radix_bits_;
    const size_t num_shift_bits_;
    std::vector<uint32_t> radix_table_;
    std::vector<KeyType> spline_points_;
    std::vector<rs::RadixSpline<KeyType>> rss_;

    // std::vector<MLP> mlps_;

    KeyType prev_prefix_;

    // Number of CDF points in one group
    size_t group_size_;

    size_t rs_group_size_;

    // Number of CDF points in one cell
    size_t cell_size_;

    // vector of a CDF points group
    std::vector<KeyType> group_points_;
  };
} // namespace gf