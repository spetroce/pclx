#ifndef __PCL_PCA_INL_H__
#define __PCL_PCA_INL_H__

#include <pthread.h>
#include <iostream>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include "mio/altro/thread.hpp"
#include "mio/math/math.hpp"
#include "pclx/core.h"

namespace pclx{

template<typename CloudType>
void ComputePntPCA(const CloudType &cloud, const table32S_t &nn_table,
                   std::vector<cv::Mat> &eigenvalues, std::vector<cv::Mat> &eigenvectors, 
                   const bool print_info = false, 
                   const float radius = -1, const table32F_t *nn_table_dist = nullptr){
  typedef typename CloudType::PointType PointType;
  typedef typename CloudType::PointType::DataType PointDataType;
  STD_INVALID_ARG_E(sizeof(PointType) / sizeof(PointDataType) == 3)
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  mio::CSysTimer sysTimer("ComputePntPCA()");
  sysTimer.start();

  const auto &pnts = *cloud.pnts;
  const uint32_t num_indices = cloud.indices.size();
  eigenvalues.resize(num_indices);
  eigenvectors.resize(num_indices);
  bool dist_check = (radius > 0 && nn_table_dist != nullptr);

#if ICV_OPENCV_VERSION_MAJOR < 3
  const int PCA_FLAGS = CV_PCA_DATA_AS_ROW;
#else
  const int PCA_FLAGS = cv::PCA::DATA_AS_ROW;
#endif

#define COMPUTE_SET_PCA_MACRO(num_pnt_)\
  /*check distance and increment number of points*/\
  pnt_mat_sub = pnt_mat( cv::Range(0, num_pnt_), cv::Range::all() );\
  /*The computed eigenvalues are sorted from the largest to the smallest (stored as a column vector)*/\
  /*and the corresponding eigenvectors are stored as vector rows accessed using PCA_object.eigenvectors.*/\
  pca(pnt_mat_sub, cv::Mat(), PCA_FLAGS);\
  eigenvalues[idx] = pca.eigenvalues.clone();\
  eigenvectors[idx] =  pca.eigenvectors.clone();\

  auto tbb_func = [&](const tbb::blocked_range<uint32_t> &br){
    const uint32_t br_end = br.end();
    uint32_t num_pnt;
    cv::PCA pca;
    const int cv_pnt_data_type = cv::DataType<PointDataType>::type;
    cv::Mat pnt_mat(num_indices, 3, cv_pnt_data_type), pnt_mat_sub;
    PointType *pnt_mat_data = mio::StaticCastPtr<PointType>(pnt_mat.data);

    if(dist_check){
      const float radius_sqrd = radius*radius;
      const table32F_t &nn_table_dist_ = *nn_table_dist;
      for(uint32_t idx = br.begin(); idx != br_end; ++idx){
        const auto &nn_vec = nn_table[idx];
        const auto &nn_dist_vec = nn_table_dist_[idx];
        const uint32_t num_nn = nn_vec.size();
        if(num_nn > 2){
          uint32_t pnt_count = 0;
          for(uint32_t i = 0; i < num_nn; ++i){
            if(nn_dist_vec[i] < radius_sqrd){
              pnt_mat_data[i] = pnts[ nn_vec[i] ];
              pnt_count++;
            }
            else //we can break because the distances are sorted (smallest to largest)
              break;
          }
          if(pnt_count < 3)
            continue;
          COMPUTE_SET_PCA_MACRO(pnt_count)
        }
        else{
          eigenvalues[idx] = cv::Mat::zeros(3, 1, cv_pnt_data_type);
          eigenvectors[idx] = cv::Mat::zeros(3, 3, cv_pnt_data_type);
        }
      }
    }
    else{
      for(uint32_t idx = br.begin(); idx != br_end; ++idx){
        const auto &nn_vec = nn_table[idx];
        const uint32_t num_nn = nn_vec.size();
        if(num_nn > 2){
          for(uint32_t i = 0; i < num_nn; ++i)
            pnt_mat_data[i] = pnts[ nn_vec[i] ];
          COMPUTE_SET_PCA_MACRO(num_nn)
        }
        else{
          eigenvalues[idx] = cv::Mat::zeros(3, 1, cv_pnt_data_type);
          eigenvectors[idx] = cv::Mat::zeros(3, 3, cv_pnt_data_type);
        }
      }
    }
  };

#undef COMPUTE_SET_PCA_MACRO

  tbb::parallel_for(tbb::blocked_range<uint32_t>(0, num_indices), tbb_func);
  sysTimer.stop(print_info);
}


template <typename CloudType>
int ComputeCloudPCA(CloudType &cloud, const bool sqrt_of_eigenvalues = false, const bool print_info = false){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  mio::CRUTimer ruTimer("ComputeCloudPCA()");
  ruTimer.start();

  //clear eigenvalue and eigenvector mats
  cloud.eigenvalues.release();
  cloud.eigenvectors.release();

  cv::PCA pca;
  const auto &pnts = *cloud.pnts;
  const uint32_t num_indices = cloud.indices.size();
  cv::Mat pnt_mat(pnts.size(), 3, cv::DataType<typename CloudType::PointType::DataType>::type);
  typename CloudType::PointType *pnt_mat_data = mio::StaticCastPtr<typename CloudType::PointType>(pnt_mat.data);

  uint32_t j = 0;
  if(cloud.indices.size() > 2){
    for(uint32_t i = 0; i < num_indices; ++i)
      pnt_mat_data[i] = pnts[i];

#if ICV_OPENCV_VERSION_MAJOR < 3
    const int PCA_FLAGS = CV_PCA_DATA_AS_ROW;
#else
    const int PCA_FLAGS = cv::PCA::DATA_AS_ROW;
#endif

    //The computed eigenvalues are sorted from the largest to the smallest (stored as a column vector)
    //and the corresponding eigenvectors are stored as vector rows accessed using PCA_object.eigenvectors.
    pca(pnt_mat, cv::Mat(), PCA_FLAGS);

    // use clone() to physically copy the data and so we can modify pca.eigen...
    cloud.eigenvectors = pca.eigenvectors.clone();
    cloud.eigenvalues = pca.eigenvalues.clone();
  }
  else{
    cloud.eigenvectors = cv::Mat::zeros(3, 3, CV_32F);
    cloud.eigenvalues = cv::Mat::zeros(3, 1, CV_32F);
    return -1;
  }

  if(sqrt_of_eigenvalues)
    cv::sqrt(cloud.eigenvalues, cloud.eigenvalues);

  ruTimer.stop(print_info);
  return 0;
}


template <typename CloudType>
void ComputeCloudPCA(std::vector<CloudType> &cloud_vec, const bool print_info = false){
  for(auto &cloud : cloud_vec)
    ComputeCloudPCA(cloud, print_info);
}


template<typename PointType>
void OrientNormal(const PointType &point, PointType view_point, cv::Mat &eigenvectors){
  view_point.x -= point.x;
  view_point.y -= point.y;
  view_point.z -= point.z;
  //dot product between the (view_point - point) and the normal vector
  const float cos_theta = ( view_point.x * eigenvectors.at<float>(2, 0) + 
                            view_point.y * eigenvectors.at<float>(2, 1) + 
                            view_point.z * eigenvectors.at<float>(2, 2) );
  //if the angle is more than 90 degrees, flip
  if(cos_theta < 0.0f)
    eigenvectors.row(2) *= -1.0f;
}


template <typename CloudType>
void OrientNormals(const CloudType &cloud, const typename CloudType::PointType &view_point,
                   std::vector<cv::Mat> &eigenvectors){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_E( cloud.indices.size() == eigenvectors.size() )

  if(cloud.indices.size() > 0){
    STD_INVALID_ARG_E(eigenvectors.front().depth() == CV_32F)
    const auto &pnts = *cloud.pnts;
    const uint32_t num_indices = cloud.indices.size();
    for(uint32_t i = 0; i < num_indices; ++i){
      STD_INVALID_ARG_E(eigenvectors[i].rows == 3 && eigenvectors[i].cols == 3)
      const auto ray = view_point - pnts[i];
      float *ev_data = mio::StaticCastPtr<float>(eigenvectors[i].data);
      const float cos_theta = ray.x*ev_data[6] + ray.y*ev_data[7] + ray.z*ev_data[8];
      if(cos_theta < 0.0f){
        ev_data[6] = -ev_data[6];
        ev_data[7] = -ev_data[7];
        ev_data[8] = -ev_data[8];
      }
    }
  }
}


//smush point cloud along greatest principle axis
template <typename CloudType>
void ProjectCloudPCA(const CloudType &cloud, CloudType &proj_cloud, const uint8_t proj_axis){
  STD_INVALID_ARG_EM(proj_axis >= 0 && proj_axis <= 2, std::string("proj_axis: ") + std::to_string(proj_axis))
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_E(proj_cloud.pnts != nullptr)
  STD_INVALID_ARG_EM(cloud.eigenvectors.rows == 3 && cloud.eigenvectors.cols == 3,
                     std::string("size: ") + mio::SizeStr(cloud.eigenvectors))
  STD_INVALID_ARG_EM(cloud.eigenvectors.depth() == CV_32F, std::string("depth: ") + mio::DepthStr(cloud.eigenvectors))
  STD_INVALID_ARG_E(cloud.eigenvectors.isContinuous())

  const uint32_t num_indices = cloud.indices->size();
  if(num_indices == 0)
    return;

  const auto &cloud_pnts = *cloud.pnts;
  auto &proj_cloud_pnts = *proj_cloud.pnts;
  auto &proj_cloud_indices = proj_cloud.indices;
  proj_cloud_pnts.resize(num_indices);
  proj_cloud_indices.resize(num_indices);
  typename CloudType::PointType p;

  const float *ev0 = cloud.eigenvectors.ptr(0),
              *ev1 = cloud.eigenvectors.ptr(1),
              *ev2 = cloud.eigenvectors.ptr(2);
  if(proj_axis == 0){ //eigenvector associated with largest eigenvalue
    for(uint32_t i = 0; i < num_indices; ++i){
      const uint32_t idx = cloud.indices[i];
      const auto &cloud_pnt = cloud_pnts[idx];
      proj_cloud_pnts[i].x = cloud_pnt.x * ev1[0] + cloud_pnt.y * ev1[1] + cloud_pnt.z * ev1[2];
      proj_cloud_pnts[i].y = cloud_pnt.x * ev2[0] + cloud_pnt.y * ev2[1] + cloud_pnt.z * ev2[2];
      proj_cloud_pnts[i].z = 0.0f;
      proj_cloud_indices[i] = i;
    }
  }
  if(proj_axis == 1){
    for(uint32_t i = 0; i < num_indices; ++i){
      const uint32_t idx = cloud.indices[i];
      const auto &cloud_pnt = cloud_pnts[idx];
      proj_cloud_pnts[i].x = cloud_pnt.x * ev0[0] + cloud_pnt.y * ev0[1] + cloud_pnt.z * ev0[2];
      proj_cloud_pnts[i].y = cloud_pnt.x * ev2[0] + cloud_pnt.y * ev2[1] + cloud_pnt.z * ev2[2];
      proj_cloud_pnts[i].z = 0.0f;
      proj_cloud_indices[i] = i;
    }
  }
  if(proj_axis == 2){ //eigenvector associated with smallest eigenvalue
    for(uint32_t i = 0; i < num_indices; ++i){
      const uint32_t idx = cloud.indices[i];
      const auto &cloud_pnt = cloud_pnts[idx];
      proj_cloud_pnts[i].x = cloud_pnt.x * ev0[0] + cloud_pnt.y * ev0[1] + cloud_pnt.z * ev0[2];
      proj_cloud_pnts[i].y = cloud_pnt.x * ev1[0] + cloud_pnt.y * ev1[1] + cloud_pnt.z * ev1[2];
      proj_cloud_pnts[i].z = 0.0f;
      proj_cloud_indices[i] = i;
    }
  }
}


template <typename CloudType>
void GetBoundingPlanePCA(CloudType &cloud, typename CloudType::PointType (&plane_vertices)[4],
                         const bool compute_centroid){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_EM(cloud.eigenvectors.rows == 3 && cloud.eigenvectors.cols == 3,
                     std::string("size: ") + mio::SizeStr(cloud.eigenvectors))
  typedef typename CloudType::PointType PointType;
  if(cloud.indices.size() < 4){
    for(auto &vert : plane_vertices)
      vert = PointType(0, 0, 0);
    return;
  }

  STD_INVALID_ARG_EM(cloud.eigenvectors.depth() == CV_32F, std::string("depth: ") + mio::DepthStr(cloud.eigenvectors))
  STD_INVALID_ARG_E(cloud.eigenvectors.isContinuous())
  const float *ev0 = mio::StaticCastPtr<float>(cloud.eigenvectors.ptr(0)),
              *ev1 = mio::StaticCastPtr<float>(cloud.eigenvectors.ptr(1));
  float ev0_max_dp = FLT_MIN, ev1_max_dp = FLT_MIN,
        ev0_min_dp = FLT_MAX, ev1_min_dp = FLT_MAX, ev0_dp, ev1_dp;
  const auto &cloud_pnts = *cloud.pnts;
  if(compute_centroid)
    cloud.ComputeCentroid();
  const PointType &cent = cloud.centroid;

  for(const auto &idx : cloud.indices){
    const auto &cloud_pnt = cloud_pnts[idx];
    ev0_dp = (cloud_pnt.x-cent.x) * ev0[0] + (cloud_pnt.y-cent.y) * ev0[1] + (cloud_pnt.z-cent.z) * ev0[2];
    ev1_dp = (cloud_pnt.x-cent.x) * ev1[0] + (cloud_pnt.y-cent.y) * ev1[1] + (cloud_pnt.z-cent.z) * ev1[2];
    if(ev0_dp > ev0_max_dp)
      ev0_max_dp = ev0_dp;
    else if(ev0_dp < ev0_min_dp)
      ev0_min_dp = ev0_dp;
    if(ev1_dp > ev1_max_dp)
      ev1_max_dp = ev1_dp;
    else if(ev1_dp < ev1_min_dp)
      ev1_min_dp = ev1_dp;
  }

  //compute plane edge midpoints
  PointType ev0_max_mid(ev0_max_dp*ev0[0], ev0_max_dp*ev0[1], ev0_max_dp*ev0[2]),
            ev0_min_mid(ev0_min_dp*ev0[0], ev0_min_dp*ev0[1], ev0_min_dp*ev0[2]),
            ev1_max_mid(ev1_max_dp*ev1[0], ev1_max_dp*ev1[1], ev1_max_dp*ev1[2]),
            ev1_min_mid(ev1_min_dp*ev1[0], ev1_min_dp*ev1[1], ev1_min_dp*ev1[2]);
  //compute corners of plane and add to the centroid
  plane_vertices[0] = cent + ev0_max_mid + ev1_max_mid;
  plane_vertices[1] = cent + ev0_max_mid + ev1_min_mid;
  plane_vertices[2] = cent + ev0_min_mid + ev1_min_mid;
  plane_vertices[3] = cent + ev0_min_mid + ev1_max_mid;
}


template <typename CloudType>
void GetBoundingCirclePCA(CloudType &cloud, std::vector<typename CloudType::PointType> &circ_vertices,
                          const int num_circ_vertices, const float radius, const bool compute_centroid){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  STD_INVALID_ARG_EM(cloud.eigenvectors.rows == 3 && cloud.eigenvectors.cols == 3,
                     std::string("size: ") + mio::SizeStr(cloud.eigenvectors))
  typedef typename CloudType::PointType PointType;

  circ_vertices.resize(num_circ_vertices);
  const float *ev0 = mio::StaticCastPtr<float>(cloud.eigenvectors.ptr(0)),
              *ev1 = mio::StaticCastPtr<float>(cloud.eigenvectors.ptr(1));
  if(compute_centroid)
    cloud.ComputeCentroid();
  const PointType &cent = cloud.centroid;
  const float theta_step = 2.0*M_PI/static_cast<float>(num_circ_vertices);
  float theta = 0;

  // first vertice is ev0*radius
  circ_vertices[0] = PointType(cent.x + ev0[0]*radius, cent.y + ev0[1]*radius, cent.z + ev0[2]*radius);
  for(size_t i = 1; i < num_circ_vertices; ++i){
    theta += theta_step;
    const float cos_theta = std::cos(theta);
    const float sin_theta = std::sin(theta);
    circ_vertices[i] = PointType(cent.x + (ev0[0]*cos_theta + ev1[0]*sin_theta)*radius,
                                 cent.y + (ev0[1]*cos_theta + ev1[1]*sin_theta)*radius,
                                 cent.z + (ev0[2]*cos_theta + ev1[2]*sin_theta)*radius);
  }
}

} //namespace pc

#endif //__PCL_PCA_INL_H__

