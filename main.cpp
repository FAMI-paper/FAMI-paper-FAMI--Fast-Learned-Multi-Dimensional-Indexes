#include <vector>
#include <chrono>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include "libmorton/morton.h"
#include "include/multiD_map.h"
#include <stx/btree_multimap.h>
using namespace std;

namespace util
{
  class MyVisitor : public SpatialIndex::IVisitor
  {
  public:
    void visitNode(const SpatialIndex::INode & /* n */) override {}

    void visitData(const SpatialIndex::IData &d) override
    {
      rstar_res.push_back(d.getIdentifier());
      // the ID of this data entry is an answer to the query. I will just print it to stdout.
    }

    void visitData(std::vector<const SpatialIndex::IData *> &v) override
    {
    }
  };

  template <class KeyType>
  vector<vector<KeyType>> LoadData(const string &filename, int dim, int data_size)
  {
    vector<vector<KeyType>> data;
    ifstream in(filename);
    if (!in.is_open())
    {
      cerr << "unable to open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    vector<KeyType> item(dim);
    float key;
    while (data_size--)
    {
      if (!(in >> key))
      {
        cerr << "Not enough data entries" << endl;
        throw;
      }
      item[0] = fabs(key);
      for (int i = 1; i < dim; i++)
      {
        in >> key;
        item[i] = fabs(key);
      }
      data.push_back(item);
    }
    return data;
  }

  template <class KeyType>
  vector<vector<KeyType>> LoadData(const string &filename, int dim)
  {
    vector<vector<KeyType>> data;
    ifstream in(filename);
    if (!in.is_open())
    {
      cerr << "unable to open " << filename << endl;
      exit(EXIT_FAILURE);
    }
    vector<KeyType> item(dim);
    float key;
    while (in >> key)
    {
      item[0] = fabs(key);
      for (int i = 1; i < dim; i++)
      {
        in >> key;
        item[i] = fabs(key);
      }
      data.push_back(item);
    }
    return data;
  }

  template <class KeyType>
  vector<KeyType> GetMinKey(const vector<vector<KeyType>> &data, int dim)
  {
    vector<KeyType> min_keys(dim);
    int size = data.size();
    for (int i = 0; i < dim; i++)
    {
      min_keys[i] = data[0][i];
      for (int j = 1; j < size; j++)
      {
        min_keys[i] = min(min_keys[i], data[j][i]);
      }
    }
    return min_keys;
  }
  template <class KeyType>
  vector<KeyType> GetMaxKey(const vector<vector<KeyType>> &data, int dim)
  {
    vector<KeyType> max_keys(dim);
    int size = data.size();
    for (int i = 0; i < dim; i++)
    {
      max_keys[i] = data[0][i];
      for (int j = 1; j < size; j++)
      {
        max_keys[i] = max(max_keys[i], data[j][i]);
      }
    }
    return max_keys;
  }

  float gridsize[2];
  template <class KeyType>
  vector<vector<KeyType>> GetUpperData(const vector<vector<KeyType>> &data, int dim, vector<float> max_keys)
  {
    std::default_random_engine e;
    std::normal_distribution<float> n(0, 1);
    vector<vector<float>> upper_data;
    vector<KeyType> item(dim);
    int size = data.size();
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < dim; j++)
      {
        // item[j] = data[i][j] + max(0.0, gridsize[j] + n(e));
        item[j] = min(data[i][j] + gridsize[j], max_keys[j]);
      }
      upper_data.push_back(item);
    }
    return upper_data;
  }
}

std::vector<Polygon<uint64_t>> pols_;
int main(int argc, char **argv)
{
  if (argc != 13)
  {
    cerr << "usage: " << argv[0] << "<data_file> <data_num> <query_num> <dimension> <exp> <index_type> <num_radix_bits> <cell_size> <rs_group_size> <query_type> <range_size> <update_ratio>" << endl;
    throw;
  }
  srand(921224);
  const string data_file = argv[1];
  const int data_size = atoi(argv[2]);
  const int query_size = atoi(argv[3]);
  const int dim = atoi(argv[4]);
  const int exp = atoi(argv[5]);
  string index_type = argv[6];
  const int num_radix_bits = atoi(argv[7]);
  const int cell_size = atoi(argv[8]);
  const int rs_group_size = atoi(argv[9]);
  const int query_type = atoi(argv[10]);
  int range_size;
  double update_ratio;
  if (query_type)
  {
    range_size = atoi(argv[11]);
  }
  if (query_type == 2)
  {
    update_ratio = atof(argv[12]);
  }
  size_t ratio = pow(10, exp);
  vector<vector<float>> data = util::LoadData<float>(data_file, dim, data_size);
  // vector<vector<float>> data = util::LoadData<float>(data_file, dim);
  // std::random_shuffle(data.begin(), data.end());
  // data.resize(data_size);
  int n = data.size();
  const int t_size = std::sqrt(n / cell_size) + 1;
  vector<float> min_keys = util::GetMinKey<float>(data, dim);
  vector<float> max_keys = util::GetMaxKey<float>(data, dim);
  util::gridsize[0] = (max_keys[0] - min_keys[0]) / range_size;
  util::gridsize[1] = (max_keys[1] - min_keys[1]) / range_size;
  uint64_t x_cell_size = uint64_t((max_keys[0] - min_keys[0]) * ratio) / (t_size - 1);
  uint64_t y_cell_size = uint64_t((max_keys[1] - min_keys[1]) * ratio) / (t_size - 1);
  uint64_t minx = round(min_keys[0] * ratio);
  uint64_t miny = round(min_keys[1] * ratio);

  vector<vector<float>> insert_data;
  if (query_type == 2)
  {
    insert_data.insert(insert_data.end(), data.begin() + int(update_ratio * data_size), data.end());
    data.resize(int(update_ratio * data_size));
    cout << data.size() << endl;
  }
  std::sort(data.begin(), data.end(), [](const std::vector<float> &a, const std::vector<float> &b)
            { return a[0] < b[0]; });
  n = data.size();

  // grid+r-tree
  mdm::MultiDMap<float> *map;

  // avggrid+r-tree
  RTree<uint64_t> ***rtree = (RTree<uint64_t> ***)malloc(sizeof(RTree<uint64_t> **) * t_size);
  for (int i = 0; i < t_size; i++)
  {
    rtree[i] = (RTree<uint64_t> **)malloc(sizeof(RTree<uint64_t> *) * t_size);
  }

  // r*-tree
  SpatialIndex::IStorageManager *memoryfile = SpatialIndex::StorageManager::createNewMemoryStorageManager();
  SpatialIndex::StorageManager::IBuffer *file = SpatialIndex::StorageManager::createNewRandomEvictionsBuffer(*memoryfile, 10, false);
  int64_t rtree_id;
  SpatialIndex::ISpatialIndex *tree = SpatialIndex::RTree::createNewRTree(*file, 0.7, 20, 20, 2, SpatialIndex::RTree::RV_RSTAR, rtree_id);

  // avggrid+b-tree
  stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                      stx::btree_default_map_traits<uint64_t, uint64_t>>
      ***bt = (stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                   stx::btree_default_map_traits<uint64_t, uint64_t>> ***)malloc(sizeof(stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                                                                                                            stx::btree_default_map_traits<uint64_t, uint64_t>> **) *
                                                                                                 t_size);
  for (int i = 0; i < t_size; i++)
  {
    bt[i] = (stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                 stx::btree_default_map_traits<uint64_t, uint64_t>> **)malloc(sizeof(stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                                                                                                         stx::btree_default_map_traits<uint64_t, uint64_t>> *) *
                                                                                              t_size);
  }

  auto build_begin = chrono::high_resolution_clock::now();

  if (index_type == "grid+r-tree")
  {
    // grid+r-tree
    map = new mdm::MultiDMap<float>(data, num_radix_bits, cell_size, rs_group_size, exp);
  }

  if (index_type == "avggrid+r-tree")
  {
    // avggrid+r-tree
    for (int i = 0; i < t_size; i++)
    {
      for (int j = 0; j < t_size; j++)
      {
        rtree[i][j] = new RTree<uint64_t>(100);
      }
    }
    for (int i = 0; i < n; i++)
    {
      uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
      pols_.push_back(Polygon<uint64_t>(Point<uint64_t>(x, y)));
    }
    for (int i = 0; i < n; i++)
    {
      uint64_t x = pols_[i].get_Pmin().get_X(), y = pols_[i].get_Pmin().get_Y();
      rtree[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert_polygon(&pols_[i], &pols_[i]);
    }
  }

  if (index_type == "r*-tree")
  {
    // r*-tree
    for (int i = 0; i < n; i++)
    {
      uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
      std::string value = std::to_string(x);
      double p[2];
      p[0] = x;
      p[1] = y;
      SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
      tree->insertData((uint32_t)(value.size() + 1), reinterpret_cast<const uint8_t *>(value.c_str()), r, x);
    }
  }

  if (index_type == "avggrid+b-tree")
  {
    // avggrid+b-tree
    for (int i = 0; i < t_size; i++)
    {
      for (int j = 0; j < t_size; j++)
      {
        bt[i][j] = new stx::btree_multimap<uint64_t, uint64_t, std::less<uint64_t>,
                                           stx::btree_default_map_traits<uint64_t, uint64_t>>();
      }
    }
    for (int i = 0; i < n; i++)
    {
      uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
      uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
      bt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert2(id, x);
      // cnt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]++;
    }
  }

  auto build_end = chrono::high_resolution_clock::now();
  double build_s = chrono::duration_cast<chrono::nanoseconds>(build_end - build_begin).count() / 1.0 / 1e9;
  cout << "build/s: " << build_s << endl;

  std::random_shuffle(data.begin(), data.end());
  data.resize(query_size);
  vector<vector<float>> upper_data = util::GetUpperData<float>(data, dim, max_keys);
  n = data.size();

  auto lookup_begin = chrono::high_resolution_clock::now();
  uint_fast64_t k = 0;
  if (query_type)
  {
    // range query

    if (index_type == "grid+r-tree")
    {
      // grid+r-tree
      for (int i = 0; i < n; i++)
      {
        map->range_search(data[i], upper_data[i]);
      }
    }

    if (index_type == "avggrid+r-tree")
    {
      // avggrid+r-tree
      for (int i = 0; i < n; i++)
      {
        std::vector<uint64_t> plow(2);
        plow[0] = round(data[i][0] * ratio);
        plow[1] = round(data[i][1] * ratio);
        std::vector<uint64_t> phigh(2);
        phigh[0] = round(upper_data[i][0] * ratio);
        phigh[1] = round(upper_data[i][1] * ratio);
        std::vector<Polygon<uint64_t> *> total_answ;

        int bx = (plow[0] - minx) / x_cell_size;
        int ex = (phigh[0] - minx) / x_cell_size;
        int by = (plow[1] - miny) / y_cell_size;
        int ey = (phigh[1] - miny) / y_cell_size;
        for (int i = bx; i <= ex; i++)
        {
          for (int j = by; j <= ey; j++)
          {
            std::vector<Polygon<uint64_t> *> answ;
            rtree[i][j]->range_search(Polygon<uint64_t>(Point<uint64_t>(plow[0], plow[1]), Point<uint64_t>(phigh[0], phigh[1])), answ);
            total_answ.insert(total_answ.end(), answ.begin(), answ.end());
          }
        }
        total_cnt += total_answ.size();
        // for (int i = 0; i < total_answ.size(); i++)
        // {
        //   uint64_t x = total_answ[i]->get_Pmin().get_X();
        //   uint64_t y = total_answ[i]->get_Pmin().get_Y();
        //   if (x < plow[0] || y < plow[1] || x > phigh[0] || y > phigh[1])
        //   {
        //     printf("error\n");
        //     exit(-1);
        //   }
        // }
        k |= total_answ[0]->get_Pmin().get_X();
      }
    }

    if (index_type == "r*-tree")
    {
      // r*-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t xlow = round(data[i][0] * ratio), ylow = round(data[i][1] * ratio);
        uint64_t xhigh = round(upper_data[i][0] * ratio), yhigh = round(upper_data[i][1] * ratio);
        double plow[2], phigh[2];
        plow[0] = xlow;
        plow[1] = ylow;
        phigh[0] = xhigh;
        phigh[1] = yhigh;
        util::MyVisitor vis;
        SpatialIndex::Region r = SpatialIndex::Region(plow, phigh, 2);
        tree->intersectsWithQuery(r, vis);
        total_cnt += rstar_res.size();
        // for (uint64_t res : rstar_res)
        // {
        //   if (res < xlow || res > xhigh)
        //   {
        //     printf("error\n");
        //     // exit(-1);
        //   }
        // }
        k |= rstar_res[0];
        rstar_res.clear();
      }
    }

    if (index_type == "avggrid+b-tree")
    {
      // avggrid+b-tree
      for (int i = 0; i < n; i++)
      {
        std::vector<uint64_t> plow(2);
        plow[0] = round(data[i][0] * ratio);
        plow[1] = round(data[i][1] * ratio);
        std::vector<uint64_t> phigh(2);
        phigh[0] = round(upper_data[i][0] * ratio);
        phigh[1] = round(upper_data[i][1] * ratio);
        std::vector<Polygon<uint64_t> *> total_answ;

        uint64_t bx = (plow[0] - minx) / x_cell_size;
        uint64_t ex = (phigh[0] - minx) / x_cell_size;
        uint64_t by = (plow[1] - miny) / y_cell_size;
        uint64_t ey = (phigh[1] - miny) / y_cell_size;
        for (uint64_t i = bx; i <= ex; i++)
        {
          for (uint64_t j = by; j <= ey; j++)
          {
            std::vector<Polygon<uint64_t> *> answ;
            uint64_t x = max(i * x_cell_size + minx, plow[0]);
            uint64_t y = max(j * y_cell_size + miny, plow[1]);
            uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
            auto iter = bt[i][j]->lower_bound(id);
            uint64_t resx, resy;
            while (iter != bt[i][j]->end())
            {
              id = iter->first;
              libmorton::morton2D_64_decode(id, resx, resy);
              resx += i * x_cell_size + minx;
              resy += j * y_cell_size + miny;
              if (resx >= plow[0] && resx <= phigh[0] && resy >= plow[1] && resy <= phigh[1])
                total_cnt++;
              iter++;
            }
          }
        }
      }
    }
  }
  else
  {
    // point query

    if (index_type == "grid+r-tree")
    {
      // grid+r-tree
      for (int i = 0; i < n; i++)
      {
        k |= map->find(data[i]);
      }
    }

    if (index_type == "avggrid+r-tree")
    {
      // avggrid+r-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
        std::vector<Polygon<uint64_t> *> answ;
        rtree[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->range_search(Polygon<uint64_t>(Point<uint64_t>(x, y)), answ);
        total_cnt += answ.size();
        if (answ[0]->get_Pmin().get_X() != x || answ[0]->get_Pmin().get_Y() != y)
        {
          printf("error result\n");
          exit(-1);
        }
        k |= answ[0]->get_Pmin().get_X();
      }
    }

    if (index_type == "r*-tree")
    {
      // r*-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
        double p[2];
        p[0] = x;
        p[1] = y;
        util::MyVisitor vis;
        SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
        tree->intersectsWithQuery(r, vis);
        total_cnt += rstar_res.size();
        if (rstar_res[0] != x)
        {
          printf("error\n");
          // exit(-1);
        }
        k |= rstar_res[0];
        rstar_res.clear();
      }
    }

    if (index_type == "avggrid+b-tree")
    {
      // avggrid+b-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(data[i][0] * ratio), y = round(data[i][1] * ratio);
        uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
        auto iter = bt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->lower_bound(id);
        uint64_t res = iter->second;
        if (res != x)
        {
          printf("error\n");
          // exit(-1);
        }
        while (iter->first == id)
        {
          total_cnt++;
          iter++;
        }
        k |= res;
      }
    }
  }
  auto lookup_end = chrono::high_resolution_clock::now();
  uint64_t lookup_ns = chrono::duration_cast<chrono::nanoseconds>(lookup_end - lookup_begin).count();
  cout << "lookup/ns: " << lookup_ns / n << endl;
  cout << "total_cnt: " << total_cnt << endl;
  cout << k << endl;

  if (query_type == 2)
  {
    auto insert_begin = chrono::high_resolution_clock::now();
    n = insert_data.size();
    if (index_type == "grid+r-tree")
    {
      // grid+r-tree
      for (int i = 0; i < n; i++)
      {
        map->insert(insert_data[i]);
      }
    }
    if (index_type == "avggrid+r-tree")
    {
      // avggrid+r-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(insert_data[i][0] * ratio), y = round(insert_data[i][1] * ratio);
        pols_.push_back(Polygon<uint64_t>(Point<uint64_t>(x, y)));
      }
      for (int i = 0; i < n; i++)
      {
        uint64_t x = pols_[i].get_Pmin().get_X(), y = pols_[i].get_Pmin().get_Y();
        rtree[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert_polygon(&pols_[i], &pols_[i]);
      }
    }
    if (index_type == "r*-tree")
    {
      // r*-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(insert_data[i][0] * ratio), y = round(insert_data[i][1] * ratio);
        std::string value = std::to_string(x);
        double p[2];
        p[0] = x;
        p[1] = y;
        SpatialIndex::Region r = SpatialIndex::Region(p, p, 2);
        tree->insertData((uint32_t)(value.size() + 1), reinterpret_cast<const uint8_t *>(value.c_str()), r, x);
      }
    }

    if (index_type == "avggrid+b-tree")
    {
      // avggrid+b-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t x = round(insert_data[i][0] * ratio), y = round(insert_data[i][1] * ratio);
        uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
        bt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]->insert2(id, x);
        // cnt[(x - minx) / x_cell_size][(y - miny) / y_cell_size]++;
      }
    }

    auto insert_end = chrono::high_resolution_clock::now();
    uint64_t insert_ns = chrono::duration_cast<chrono::nanoseconds>(insert_end - insert_begin).count();
    cout << "insert/ns: " << insert_ns / n << endl;

    n = data.size();
    total_cnt = 0;
    auto afterlookup_begin = chrono::high_resolution_clock::now();
    k = 0;

    if (index_type == "grid+r-tree")
    {
      // grid+r-tree
      for (int i = 0; i < n; i++)
      {
        map->range_search(data[i], upper_data[i]);
      }
    }

    if (index_type == "avggrid+r-tree")
    {
      // avggrid+r-tree
      for (int i = 0; i < n; i++)
      {
        std::vector<uint64_t> plow(2);
        plow[0] = round(data[i][0] * ratio);
        plow[1] = round(data[i][1] * ratio);
        std::vector<uint64_t> phigh(2);
        phigh[0] = round(upper_data[i][0] * ratio);
        phigh[1] = round(upper_data[i][1] * ratio);
        std::vector<Polygon<uint64_t> *> total_answ;

        int bx = (plow[0] - minx) / x_cell_size;
        int ex = (phigh[0] - minx) / x_cell_size;
        int by = (plow[1] - miny) / y_cell_size;
        int ey = (phigh[1] - miny) / y_cell_size;
        for (int i = bx; i <= ex; i++)
        {
          for (int j = by; j <= ey; j++)
          {
            std::vector<Polygon<uint64_t> *> answ;
            rtree[i][j]->range_search(Polygon<uint64_t>(Point<uint64_t>(plow[0], plow[1]), Point<uint64_t>(phigh[0], phigh[1])), answ);
            total_answ.insert(total_answ.end(), answ.begin(), answ.end());
          }
        }
        total_cnt += total_answ.size();
        // for (int i = 0; i < total_answ.size(); i++)
        // {
        //   uint64_t x = total_answ[i]->get_Pmin().get_X();
        //   uint64_t y = total_answ[i]->get_Pmin().get_Y();
        //   if (x < plow[0] || y < plow[1] || x > phigh[0] || y > phigh[1])
        //   {
        //     printf("error\n");
        //     exit(-1);
        //   }
        // }
        k |= total_answ[0]->get_Pmin().get_X();
      }
    }

    if (index_type == "r*-tree")
    {
      // r*-tree
      for (int i = 0; i < n; i++)
      {
        uint64_t xlow = round(data[i][0] * ratio), ylow = round(data[i][1] * ratio);
        uint64_t xhigh = round(upper_data[i][0] * ratio), yhigh = round(upper_data[i][1] * ratio);
        double plow[2], phigh[2];
        plow[0] = xlow;
        plow[1] = ylow;
        phigh[0] = xhigh;
        phigh[1] = yhigh;
        util::MyVisitor vis;
        SpatialIndex::Region r = SpatialIndex::Region(plow, phigh, 2);
        tree->intersectsWithQuery(r, vis);
        total_cnt += rstar_res.size();
        // for (uint64_t res : rstar_res)
        // {
        //   if (res < xlow || res > xhigh)
        //   {
        //     printf("error\n");
        //     // exit(-1);
        //   }
        // }
        k |= rstar_res[0];
        rstar_res.clear();
      }
    }

    if (index_type == "avggrid+b-tree")
    {
      // avggrid+b-tree
      for (int i = 0; i < n; i++)
      {
        std::vector<uint64_t> plow(2);
        plow[0] = round(data[i][0] * ratio);
        plow[1] = round(data[i][1] * ratio);
        std::vector<uint64_t> phigh(2);
        phigh[0] = round(upper_data[i][0] * ratio);
        phigh[1] = round(upper_data[i][1] * ratio);
        std::vector<Polygon<uint64_t> *> total_answ;

        uint64_t bx = (plow[0] - minx) / x_cell_size;
        uint64_t ex = (phigh[0] - minx) / x_cell_size;
        uint64_t by = (plow[1] - miny) / y_cell_size;
        uint64_t ey = (phigh[1] - miny) / y_cell_size;
        for (uint64_t i = bx; i <= ex; i++)
        {
          for (uint64_t j = by; j <= ey; j++)
          {
            std::vector<Polygon<uint64_t> *> answ;
            uint64_t x = max(i * x_cell_size + minx, plow[0]);
            uint64_t y = max(j * y_cell_size + miny, plow[1]);
            uint64_t id = libmorton::morton2D_64_encode((x - minx) % x_cell_size, (y - miny) % y_cell_size);
            auto iter = bt[i][j]->lower_bound(id);
            uint64_t resx, resy;
            while (iter != bt[i][j]->end())
            {
              id = iter->first;
              libmorton::morton2D_64_decode(id, resx, resy);
              resx += i * x_cell_size + minx;
              resy += j * y_cell_size + miny;
              if (resx >= plow[0] && resx <= phigh[0] && resy >= plow[1] && resy <= phigh[1])
                total_cnt++;
              iter++;
            }
          }
        }
      }
    }
    auto afterlookup_end = chrono::high_resolution_clock::now();
    uint64_t afterlookup_ns = chrono::duration_cast<chrono::nanoseconds>(afterlookup_end - afterlookup_begin).count();
    cout << "afterlookup/ns: " << afterlookup_ns / n << endl;
    cout << "total_cnt: " << total_cnt << endl;
    cout << k << endl;
  }

  return 1;
}