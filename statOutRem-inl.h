#include "pclx/core.h"

namespace pclx{

template <typename CloudType>
void StatOutRem(const CloudType &src_cloud, CloudType &dst_cloud,
                const std::vector< std::vector<float> > &nn_table_dist,
                const double std_dev_coef, const bool print_stats = false){
  STD_INVALID_ARG_E(src_cloud.pnts != nullptr && src_cloud.indices.size() > 0)

  mio::CRUTimer ruTime("pclx::StatOutRem()");
  ruTime.start();

  const uint32_t src_indices_size = src_cloud.indices.size();
  std::vector<float> dists(src_indices_size);
  dst_cloud.pnts = src_cloud.pnts;
  dst_cloud.colors = src_cloud.colors;
  
  //this is assuming src_cloud.indices and nn_table_dist are parallel
  //for each pnt (nn_table_dist[i]), we compute the average euclidean distance to it's nearest neighbors
  for(uint32_t i = 0; i < src_indices_size; ++i){
    double dist_sum = 0;
    for(int j = 1; j < nn_table_dist[i].size(); ++j)
      dist_sum += std::sqrt(nn_table_dist[i][j]);
    const double temp = nn_table_dist[i].size() > 1 ? nn_table_dist[i].size() - 1 : 1;
    dists[i] = dist_sum / temp;
  }

  cv::Scalar mean, stddev;
  cv::meanStdDev(dists, mean, stddev);
  
  if(print_stats)
    printf("Mean: %f, StdDev: %f\n", mean.val[0], stddev.val[0]);
  
  //a distance that is bigger than this signals an outlier
  const double dist_thresh = mean.val[0] + std_dev_coef * stddev.val[0];
  
  dst_cloud.indices.clear();
  for(uint32_t i = 0; i < src_indices_size; ++i)
    if(dists[i] <= dist_thresh)
      dst_cloud.indices.push_back( src_cloud.indices[i] );
  
  ruTime.stop(print_stats);
}

} //namespace pclx

