#ifndef __PCLX_CLUSTER_H__
#define __PCLX_CLUSTER_H__

#include <iostream>
#include <queue>
#include "mio/altro/error.h"

namespace pclx{

template <typename CloudType>
void EuclideanCluster(CloudType &src_cloud, const table32S_t &nn_table,
                      std::vector<CloudType> &cloud_vec, const bool print_flag){
  STD_INVALID_ARG_E(src_cloud.indices.size() > 0)
  STD_INVALID_ARG_E(src_cloud.pnts != nullptr)
  STD_INVALID_ARG_E( nn_table.size() == src_cloud.indices.size() )
  mio::CRUTimer ruTime("pclx::Cluster()");
  ruTime.start();

  cloud_vec.clear();
  const uint32_t num_src_pnt = src_cloud.indices.size();
  uint32_t label = 0, seed_idx, lowest_labeled = 0;
  bool exit_flag = false;
  std::vector<bool> is_labeled(num_src_pnt, false);
  std::queue<uint32_t> index_queue; //holds unlabeled surface point indices that need to be radius searched
  if(src_cloud.local_idx_lu.size() == 0)
    src_cloud.BuildLocalIdxLU();

  for(;;){ //run until all points are labeled (added to a cluster)
    seed_idx = lowest_labeled; //remember where the last loop started
    for(;;){ //find an unlabeled point
      if(seed_idx == num_src_pnt){
        exit_flag = true;
        break;
      }
      if(!is_labeled[seed_idx])
        break;
      seed_idx++;
    }
    if(exit_flag)
      break;

    //found an unlabeled point. 
    //mark it as labeled, add it to a new typename CloudType, and push it onto the queue
    is_labeled[seed_idx] = true;
    cloud_vec.push_back( CloudType() );
    CloudType &cur_cloud = cloud_vec.back();
    cur_cloud.SetData(src_cloud.pnts);
    cur_cloud.indices.push_back( src_cloud.indices[seed_idx] );
    index_queue.push(seed_idx);
    lowest_labeled = seed_idx;

    //pop all neighbors onto queue, add them to the current cluster, and then check all their neighbors.
    //ie. exhaustively check and add neighbors onto the queue until it empties.
    while( !index_queue.empty() ){
      const uint32_t seed_idx_ = index_queue.front();
      index_queue.pop();
      const uint32_t num_nn = nn_table[seed_idx_].size();
      for(uint32_t i = 1; i < num_nn; i++){
        const uint32_t nn_idx = src_cloud.local_idx_lu[ nn_table[seed_idx_][i] ];
        if(!is_labeled[nn_idx]){
          index_queue.push(nn_idx);
          is_labeled[nn_idx] = true;
          cur_cloud.indices.push_back( nn_table[seed_idx_][i] );
        }
      }
    }

    label++;
  }

  ruTime.stop(print_flag);
  if(print_flag)
    printf("pclx::EuclideanCluster(): Found %d clusters\n", label);
}


template <typename CloudType>
void CombineClusters(const std::vector<CloudType> &cloud_vec, CloudType &dst_cloud,
                     const uint32_t size_threshold = 0, const bool sort_indices = false){
  STD_INVALID_ARG_E(cloud_vec.size() > 0)
  //check that all the clouds have the same 'pnts' vector
  const auto pnts_ = cloud_vec.front().pnts;
  for(const auto &cloud : cloud_vec){
    STD_INVALID_ARG_E(cloud.pnts != nullptr && cloud.pnts == pnts_)
  }
  dst_cloud.Clear();
  dst_cloud.SetData(cloud_vec.front().pnts);
  std::vector<bool> include_cloud;

  //sum all cloud_vec[i].indices.size() to get size for dst_cloud.indices
  uint32_t num_pnts = 0, i = 0, j;
  if(size_threshold > 0){
    include_cloud.resize( cloud_vec.size() );
    for(const auto &cloud : cloud_vec)
      if(cloud.indices.size() >= size_threshold){
        num_pnts += cloud.indices.size();
        include_cloud[i++] = true;
      }
      else
        include_cloud[i++] = false;
  }
  else
    for(auto &cloud : cloud_vec)
      num_pnts += cloud.indices.size();
  STD_INVALID_ARG_E(num_pnts > 0);
  dst_cloud.indices.resize(num_pnts);

  //copy indices
  i = j = 0;
  if(size_threshold > 0){
    for(const auto &cloud : cloud_vec)
      if(include_cloud[i++])
        for(const auto &idx : cloud.indices)
          dst_cloud.indices[j++] = idx;
  }
  else
    for(const auto &cloud : cloud_vec)
      for(const auto &idx : cloud.indices)
        dst_cloud.indices[j++] = idx;

  if(sort_indices)
    std::sort( dst_cloud.indices.begin(), dst_cloud.indices.end() );
}

} //namespace pclx

#endif //__PCLX_CLUSTER_H__
