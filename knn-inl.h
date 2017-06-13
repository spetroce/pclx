#ifndef __KNN_INL_H__
#define __KNN_INL_H__

#include <flann/flann.hpp>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include "mio/altro/opencv.hpp"
#include "mio/altro/casting.hpp"
#include "mio/altro/rand.hpp"
#include "pclx/core.h"


namespace pclx{

namespace knn{
  static const int ALL = -1;
  static int BAD_ALLOC = -2;
  static int SUCCESS = 0;
}

//IMPORTANT NOTE: the values of the vector rnn_table[i] are cloud.pnts indices
template <typename CloudType>
int FlannRNN(const CloudType &cloud, pclx::table32S_t &rnn_table, pclx::table32F_t &rnn_table_dist,
             const float radius, const uint32_t max_neighbors, const bool print_info = false){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_E(cloud.indices.size() > 3)

  mio::CSysTimer sysTime("FlannRNN()");
  mio::CRUTimer ruTime("buildIndex()");
  sysTime.start();   

  const uint32_t num_pnt = cloud.indices.size();
  const float radius_sqrd = radius * radius;

  std::vector<typename CloudType::PointType> temp;
  cloud.CreateContinuousPntArray(temp);
  flann::Matrix<float> flann_pnt_mat( mio::StaticCastPtr<float>( temp.data() ), num_pnt, 3 );
  
  ruTime.start();
  flann::Index< flann::L2_3D<float> > flann_index( flann_pnt_mat, flann::KDTreeSingleIndexParams(10) );
  flann_index.buildIndex();
  ruTime.stop(print_info);
  
  //clear index and distance tables
  rnn_table.clear();
  rnn_table_dist.clear();

  flann::SearchParams search_params(32);
  search_params.max_neighbors = max_neighbors;

  if(true){
    flann_index.radiusSearch(flann_pnt_mat, rnn_table, rnn_table_dist, radius_sqrd, search_params);
    for(auto &vec : rnn_table)
      for(auto &i : vec)
        i = cloud.indices[i];
  }
  else{
    try{
      rnn_table.resize(num_pnt);
      rnn_table_dist.resize(num_pnt);
    }
    catch(std::bad_alloc const &e){
      FRIENDLY_RETHROW(e)
    }

    auto tbb_func = [&](const tbb::blocked_range<uint32_t> &br){
      //flann::Matrix(T* data_, uint32_t rows_, uint32_t cols_, uint32_t stride_ = 0)
      float query_pnt[3];
      flann::Matrix<float> query(query_pnt, 1, 3);
      flann::Matrix<int> indices(new int[num_pnt], 1, num_pnt);
      flann::Matrix<float> dists(new float[num_pnt], 1, num_pnt);

      const uint32_t br_end = br.end();
      for(uint32_t i = br.begin(); i != br_end; ++i){
        query_pnt[0] = temp[i].x;
        query_pnt[1] = temp[i].y;
        query_pnt[2] = temp[i].z;
        //the closest point to a point is itself; indices[0][0] = query indice
        uint32_t num_indices = flann_index.radiusSearch(query, indices, dists, radius_sqrd, search_params);

        try{
          rnn_table[i].resize(num_indices);
          rnn_table_dist[i].resize(num_indices);
        }
        catch(std::bad_alloc const &e){
          FRIENDLY_RETHROW(e)
        }
        for(uint32_t j = 0; j < num_indices; ++j){
          rnn_table[i][j] = cloud.indices[ indices[0][j] ];
          rnn_table_dist[i][j] = dists[0][j]; //distances are arranged smallest to largest
        }
      }

      delete[] indices.ptr();
      delete[] dists.ptr();
    };

    tbb::parallel_for(tbb::blocked_range<uint32_t>(0, num_pnt), tbb_func);
  }

  sysTime.stop(print_info);
  return 0;
}


//IMPORTANT NOTE: the values of the vector knn_table[i] are cloud.pnts indices
template <typename CloudType>
int FlannKNN(const CloudType &cloud, pclx::table32S_t &knn_table, pclx::table32F_t &knn_table_dist,
             const int k_neighbors, const bool print_info = false){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_E(cloud.indices.size() > 3)

  mio::CRUTimer ruTime("FlannKNN()"), ruTime_2("FlannKNN buildIndex()");
  ruTime.start();   

  const uint32_t num_pnt = cloud.indices.size();
  const uint32_t k = k_neighbors > num_pnt ? num_pnt : k_neighbors;
  typename CloudType::PointType *cont_pnt_array = nullptr;
  cloud.CreateContinuousPntArray(cont_pnt_array);
  float *cont_pnt_array_flt = mio::StaticCastPtr<float>(cont_pnt_array);
  flann::Matrix<float> flann_pnt_mat(cont_pnt_array_flt, num_pnt, 3);
  
  //clear and reallocate index and distance tables
  knn_table.clear();
  knn_table_dist.clear();

  ruTime_2.start();
  flann::Index< flann::L2_3D<float> > flann_index( flann_pnt_mat, flann::KDTreeSingleIndexParams(10) );
  flann_index.buildIndex();
  ruTime_2.stop(print_info);

  //the closest point to a point is itself; flann_indices[0][0] = query indice
  flann_index.knnSearch( flann_pnt_mat, knn_table, knn_table_dist, k, flann::SearchParams(32) );
  for(auto &vec : knn_table)
    for(auto &i : vec)
      i = cloud.indices[i];
  
  delete[] cont_pnt_array;
  ruTime.stop(print_info);
  return 0;
}


template <typename CloudType>
inline float AverageResolution(const CloudType &cloud, uint32_t num_samples, const pclx::table32F_t *nn_dist = nullptr){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_E(cloud.indices.size() > 0)

  float dist_total = 0;

  if(nn_dist == nullptr){
    pclx::table32S_t &knn_table;
    pclx::table32F_t &knn_table_dist;

    if(cloud.indices.size() < num_samples)
      num_samples = cloud.indices.size();

    std::vector<uint32_t> rand_indices(num_samples);
    mio::RandomIntVec<uint32_t>(rand_indices, 0, cloud.indices.size()-1);

    pclx::FlannKNN(cloud, knn_table, knn_table_dist, 1);

    for(auto &i : rand_indices)
      dist_total += knn_table_dist[i][0];
  }
  else{
    STD_INVALID_ARG_E(nn_dist->size() == 0)
    const auto &nn_dist_ = *nn_dist;
    if(nn_dist_.size() < num_samples)
      num_samples = nn_dist_.size();

    std::vector<uint32_t> rand_indices(num_samples);
    mio::RandomIntVec<uint32_t>(rand_indices, 0, nn_dist_.size()-1);

    for(auto &i : rand_indices)
      dist_total += nn_dist_[i][0];
  }

  return (dist_total/num_samples);
}

}

#endif //__KNN_INL_H__
