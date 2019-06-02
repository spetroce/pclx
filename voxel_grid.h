#ifndef __PCLX_VOXEL_GRID_H__
#define __PCLX_VOXEL_GRID_H__

#include <iostream>
#include <unordered_map>
#include <vector>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include "mio/math/math.h"
#include "pclx/core.h"


namespace pclx{

struct Leaf{
  std::vector<uint32_t> indices;
  vertex3f_t centroid; //computed/actual centroid
  uint32_t centroid_index; //the index of a cloud.pnts point that is closest to the centroid
};

template <typename CloudType>
class CVoxelGrid{
  private:
    std::unordered_map<int, Leaf> m_leaf_map;
    
  public:
    void ComputeVoxelGrid(CloudType &src_cloud, CloudType &dst_cloud, 
                          vertex3f_t leaf_size, const bool print_info = false);
};


template <typename CloudType>
void CVoxelGrid<CloudType>::ComputeVoxelGrid(CloudType &src_cloud, CloudType &dst_cloud, 
                                            vertex3f_t leaf_size, const bool print_info){
  STD_INVALID_ARG_E(src_cloud.pnts != nullptr)
  if(src_cloud.indices.size() == 0)
    return;

  mio::CRUTimer ruTime("CVoxelGrid::ComputeVoxelGrid");
  ruTime.start();
  m_leaf_map.clear(); //clear out any old data in unordered_map m_leaf_map
  const auto &pnts = *src_cloud.pnts;

  vertex3f_t inv_leaf_size(1.0f / leaf_size.x, 1.0f / leaf_size.y, 1.0f / leaf_size.z);
  src_cloud.ComputeDimensions();
  //get voxel coordinates for the min and max points
  vertex3n_t min_vox( static_cast<int>( std::floor(src_cloud.min.x * inv_leaf_size.x) ),
                      static_cast<int>( std::floor(src_cloud.min.y * inv_leaf_size.y) ),
                      static_cast<int>( std::floor(src_cloud.min.z * inv_leaf_size.z) ) );
  vertex3n_t max_vox( static_cast<int>( std::floor(src_cloud.max.x * inv_leaf_size.x) ),
                      static_cast<int>( std::floor(src_cloud.max.y * inv_leaf_size.y) ),
                      static_cast<int>( std::floor(src_cloud.max.z * inv_leaf_size.z) ) );

  //number of divisions needed for each axis
  const vertex3n_t vox_div(max_vox.x - min_vox.x + 1, max_vox.y - min_vox.y + 1, max_vox.z - min_vox.z + 1);
  //voxel grid division multiplier
  const vertex3n_t vox_div_mul(1.0f, vox_div.x, vox_div.x * vox_div.y);

  vertex3n_t vox_coord;
  for(auto &i : src_cloud.indices){
    vox_coord.x = static_cast<int>(std::floor(pnts[i].x * inv_leaf_size.x) - min_vox.x);
    vox_coord.y = static_cast<int>(std::floor(pnts[i].y * inv_leaf_size.y) - min_vox.y);
    vox_coord.z = static_cast<int>(std::floor(pnts[i].z * inv_leaf_size.z) - min_vox.z);
    //compute the centroid leaf index
    const int leaf_idx = vox_coord.x * vox_div_mul.x + vox_coord.y * vox_div_mul.y + vox_coord.z * vox_div_mul.z;

    Leaf& leaf = m_leaf_map[leaf_idx];
    //accumulate a centroid for each leaf
    leaf.centroid.x += pnts[i].x;
    leaf.centroid.y += pnts[i].y;
    leaf.centroid.z += pnts[i].z;
    //record an index that refers to pnts[] in the leaf
    leaf.indices.push_back(i);
  }

  uint32_t num_voxel = 0;
  dst_cloud.Clear();
  dst_cloud.SetData(src_cloud.pnts);
  dst_cloud.indices.reserve( m_leaf_map.size() );
  
  for(std::unordered_map<int, Leaf>::iterator it = m_leaf_map.begin(); it != m_leaf_map.end(); ++it){
    Leaf &leaf = it->second; //'second' gives the value (which is the Leaf struct)
    const uint32_t num_leaf_pnt = leaf.indices.size();
    
    //compute the leafs centroid
    const float num_leaf_pnt_inv = 1.0f / static_cast<float>(num_leaf_pnt);
    leaf.centroid.x *= num_leaf_pnt_inv;
    leaf.centroid.y *= num_leaf_pnt_inv;
    leaf.centroid.z *= num_leaf_pnt_inv;;
    
    //find the pnts within the leaf (voxel) that is closest to the centroid
    if(true){
      float centroid_nn_dist = sm::SqrVerDist3( leaf.centroid, pnts[ leaf.indices[0] ] ), temp;
      leaf.centroid_index = leaf.indices[0];
      for(uint32_t i = 1; i < num_leaf_pnt; ++i)
        if( ( temp = sm::SqrVerDist3(leaf.centroid, pnts[ leaf.indices[i] ]) ) < centroid_nn_dist){
          centroid_nn_dist = temp;
          leaf.centroid_index = leaf.indices[i]; //leaf.centroid_index is given an index refering to pnts
        }
      dst_cloud.indices.push_back(leaf.centroid_index); //load the centroid nearest neighbor into the dst_cloud
    }
    else
      dst_cloud.indices.push_back(leaf.indices[0]);

    num_voxel++;
  }

  ruTime.stop(print_info);
}

} //namespace pclx

#endif //__PCLX_VOXEL_GRID_H__

