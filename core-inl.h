#ifndef __PCL_CORE_INL_H__
#define __PCL_CORE_INL_H__

#if ICV_OPENCV_VERSION_MAJOR < 3
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#else
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#endif

#include <vector>
#include <float.h>
#include "mio/altro/types.hpp"
#include "mio/altro/error.hpp"
#include "mio/altro/timers.hpp"


#define PCL_LCL_PNTS_COLORS(cloud) std::vector<typename CloudType::PointType> &pnts = *cloud.pnts; \
                                   std::vector<typename CloudType::ColorType> &colors = *cloud.colors;

namespace pclx{
  typedef std::vector< std::vector<uint32_t> > table32U_t;
  typedef std::vector< std::vector<float> > table32F_t;
  typedef std::vector< std::vector<int> > table32S_t;

  static color3ub_t red3ub(255, 0, 0);
  static color3f_t red3f(1.0f, 0, 0);
  static color4ub_t red4ub(255, 0, 0, 255);
  static color4f_t red4f(1.0f, 0, 0, 1.0f);

  static color3ub_t white3ub(255, 255, 255);
  static color3f_t white3f(1.0f, 1.0f, 1.0f);
  static color4ub_t white4ub(255, 255, 255, 255);
  static color4f_t white4f(1.0f, 1.0f, 1.0f, 1.0f);

  template <typename PNT_T = vertex3f_t, typename COLOR_T = color3ub_t>
  class CPointCloud{
    public:
      typedef PNT_T PointType;
      typedef COLOR_T ColorType;

      PNT_T centroid;
      float origin_radius, centroid_radius; //radius variable kept for legacy reasons, don't use
      cv::Mat eigenvalues,
              eigenvectors;
      vertex3f_t dim, min, max;
      bool draw_cloud;

      std::vector<PNT_T> *pnts;
      std::vector<COLOR_T> *colors;
      std::vector<uint32_t> indices; //the values of this vector are indices for 'pnts' (.eg pnts[ indices[i] ])
      //local_idx_lu satisfies the following: local_idx_lu[ indices[i] ] == i
      //that is, given an index 'i' for a pnt 'p' in the vector 'pnts', it returns an
      //index 'j' for the vector 'indices' whose value at 'j' is 'i'
      std::vector<int> local_idx_lu;

      CPointCloud() : pnts(nullptr), colors(nullptr){ Clear(); };
      CPointCloud(std::vector<PNT_T> *pnts_, std::vector<COLOR_T> *colors_) : pnts(pnts_), colors(colors_){};

      void SetData(std::vector<PNT_T> *pnts_, std::vector<COLOR_T> *colors_ = nullptr, const bool clear = true);
      bool CheckData(bool check_all = false) const;
      void Clear(const bool clear_all = true);
      void DefaultIndices();
      void ComputeCentroid();
      float ComputeCentroidRadius(const bool compute_centroid);
      float ComputeOriginRadius();
      void ComputeDimensions();
      void CenterCloud(const bool pre_compute_centroid, const bool post_compute_centroid);
      void CreateContinuousPntArray(PNT_T *&pnts_array_) const;
      void CreateContinuousPntArray(std::vector<PNT_T> &pnts_vec) const;
      void CreateContinuousColorArray(COLOR_T *&colors_array_) const;
      void CreateContinuousColorArray(std::vector<COLOR_T> &colors_vec) const;
      void BuildLocalIdxLU();
      void ScaleCloud(const float cloud_radius, const float desired_cloud_radius);

    private:
      static uint32_t id_index;
  };


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::SetData(std::vector<PNT_T> *pnts_, std::vector<COLOR_T> *colors_, const bool clear){
    if(clear)
      Clear();
    pnts = pnts_;
    colors = colors_;
  }


  template <typename PNT_T, typename COLOR_T>
  bool CPointCloud<PNT_T, COLOR_T>::CheckData(bool check_all) const{
    if(pnts == nullptr)
      return(false);
    if(pnts->size() == nullptr)
      return(false);
    if(check_all){
      if(colors == nullptr)
        return(false);
      if(colors->size() != indices.size() )
        return(false);
    }
    return(true);
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::Clear(const bool clear_all){
    if(clear_all){
      pnts = nullptr;
      colors = nullptr;
    }
    indices.clear();
    local_idx_lu.clear();
    centroid = PNT_T(0.0f, 0.0f, 0.0f);
    centroid_radius = origin_radius = 0.0f;
    dim = min = max = vertex3f_t();
    eigenvalues.release();
    eigenvectors.release();
    draw_cloud = true;
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::DefaultIndices(){
    STD_INVALID_ARG_E(pnts != nullptr)
    const uint32_t pnts_size = pnts->size();
    indices.resize(pnts_size);
    for(uint32_t i = 0; i < pnts_size; ++i)
      indices[i] = i;
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::ComputeCentroid(){
    STD_INVALID_ARG_E(pnts != nullptr)
    double x = 0, y = 0, z = 0;
    std::vector<PNT_T> &pnts_ = *pnts;
    for(auto &i : indices){
      x += pnts_[i].x;
      y += pnts_[i].y;
      z += pnts_[i].z;
    }
    const float indices_size = indices.size();
    centroid.x = x / indices_size;
    centroid.y = y / indices_size;
    centroid.z = z / indices_size;
  }


  template <typename PNT_T, typename COLOR_T>
  float CPointCloud<PNT_T, COLOR_T>::ComputeCentroidRadius(const bool compute_centroid){
    STD_INVALID_ARG_E(pnts != nullptr)
    double dist, max = DBL_MIN;
    const std::vector<PNT_T> &pnts_ = *pnts;
    float x, y, z;
    if(compute_centroid)
      ComputeCentroid();
    for(auto &i : indices){
      x = pnts_[i].x - centroid.x;
      y = pnts_[i].y - centroid.y;
      z = pnts_[i].z - centroid.z;
      dist = x*x + y*y + z*z;
      if(dist > max)
        max = dist;
    }
    centroid_radius = std::sqrt(max);
    return centroid_radius;
  }


  template <typename PNT_T, typename COLOR_T>
  float CPointCloud<PNT_T, COLOR_T>::ComputeOriginRadius(){
    STD_INVALID_ARG_E(pnts != nullptr)
    double dist, max = DBL_MIN;
    const std::vector<PNT_T> &pnts_ = *pnts;
    for(auto &i : indices){
      dist = pnts_[i].x * pnts_[i].x + pnts_[i].y * pnts_[i].y + pnts_[i].z * pnts_[i].z;
      if(dist > max)
        max = dist;
    }
    origin_radius = std::sqrt(max);
    return origin_radius;
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::ComputeDimensions(){
    STD_INVALID_ARG_E(pnts != nullptr)
    if(pnts->size() > 0){
      std::vector<PNT_T> &pnts_ = *pnts;
      min.x = max.x = pnts_[ indices[0] ].x;
      min.y = max.y = pnts_[ indices[0] ].y;
      min.z = max.z = pnts_[ indices[0] ].z;
      for(auto &i : indices){
        if(min.x > pnts_[i].x) min.x = pnts_[i].x;
        if(max.x < pnts_[i].x) max.x = pnts_[i].x;
        if(min.y > pnts_[i].y) min.y = pnts_[i].y;
        if(max.y < pnts_[i].y) max.y = pnts_[i].y;
        if(min.z > pnts_[i].z) min.z = pnts_[i].z;
        if(max.z < pnts_[i].z) max.z = pnts_[i].z;
      }
      dim.x = std::sqrt( (max.x - min.x) * (max.x - min.x) );
      dim.y = std::sqrt( (max.y - min.y) * (max.y - min.y) );
      dim.z = std::sqrt( (max.z - min.z) * (max.z - min.z) );
    }
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::CenterCloud(const bool pre_compute_centroid, const bool post_compute_centroid){
    STD_INVALID_ARG_E(pnts != nullptr)
    if(pre_compute_centroid)
      ComputeCentroid();

    std::vector<PNT_T> &pnts_ = *pnts;
    for(auto &pnt : pnts_){
      pnt.x -= centroid.x;
      pnt.y -= centroid.y;
      pnt.z -= centroid.z;
    }

    if(post_compute_centroid)
      ComputeCentroid();
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::CreateContinuousPntArray(PNT_T *&pnts_array) const{
    STD_INVALID_ARG_E(pnts != nullptr)
    const uint32_t indices_size = indices.size();
    const std::vector<PNT_T> &pnts_ = *pnts;
    if(indices_size > 0){
      if(pnts_array != nullptr)
        delete [] pnts_array;
      pnts_array = new PNT_T[indices_size];
      for(uint32_t i = 0; i < indices_size; ++i)
        pnts_array[i] = pnts_[ indices[i] ];
    }
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::CreateContinuousPntArray(std::vector<PNT_T> &pnts_vec) const{
    STD_INVALID_ARG_E(pnts != nullptr)
    const uint32_t indices_size = indices.size();
    const std::vector<PNT_T> &pnts_ = *pnts;
    if(indices_size > 0){
      pnts_vec.resize(indices_size);
      for(uint32_t i = 0; i < indices_size; ++i)
        pnts_vec[i] = pnts_[ indices[i] ];
    }
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::CreateContinuousColorArray(COLOR_T *&colors_array) const{
    STD_INVALID_ARG_E(colors != nullptr)
    const uint32_t indices_size = indices.size();
    const std::vector<COLOR_T> &colors_ = *colors;
    if(indices_size > 0){
      if(colors_array)
        delete [] colors_array;
      colors_array = new COLOR_T[indices_size];
      for(uint32_t i = 0; i < indices_size; ++i)
        colors_array[i] = colors_[ indices[i] ];
    }
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::CreateContinuousColorArray(std::vector<COLOR_T> &colors_vec) const{
    STD_INVALID_ARG_E(colors != nullptr)
    const uint32_t indices_size = indices.size();
    const std::vector<COLOR_T> &colors_ = *colors;
    if(indices_size > 0){
      colors_vec.resize(indices_size);
      for(uint32_t i = 0; i < indices_size; ++i)
        colors_vec[i] = colors_[ indices[i] ];
    }
  }


  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::BuildLocalIdxLU(){
    STD_INVALID_ARG_E(pnts != nullptr)
    STD_INVALID_ARG_E(indices.size() > 0)
    uint32_t max_idx = 0;
    for(auto &idx : indices)
      if(idx > max_idx)
        max_idx = idx;
    local_idx_lu.resize(max_idx + 1, -1);
    uint32_t i = 0;
    for(auto &idx : indices)
      local_idx_lu[idx] = i++;
  }


  //NOTE: this will only scale values that are in the indices vec
  template <typename PNT_T, typename COLOR_T>
  void CPointCloud<PNT_T, COLOR_T>::ScaleCloud(const float cloud_radius, const float desired_cloud_radius){
    STD_INVALID_ARG_E(pnts != nullptr)
    std::vector<PNT_T> &pnts_ = *pnts;
    const float scale = desired_cloud_radius / cloud_radius;
    for(const auto &idx : indices){
      pnts_[idx].x *= scale;
      pnts_[idx].y *= scale;
      pnts_[idx].z *= scale;
    }
  }


  template <typename CloudType>
  inline void ComputeCentroid(std::vector<CloudType> &cloud_vec){
    for(auto &cloud : cloud_vec)
      ComputeCentroid(cloud);
  }


  template <typename CloudType>
  inline void ArrayToStdVec(const vertex3f_t *pnt, const uint32_t numPnt,
                            std::vector<typename CloudType::PointType> vecPnt){
    vecPnt = std::vector<typename CloudType::PointType>(pnt, pnt + numPnt);
  }


  template <typename CloudType>
  inline void RemoveOutliers(CloudType &cloud, const float max_dist){
    STD_INVALID_ARG_E(cloud.pnts != nullptr)
    const std::vector<typename CloudType::PointType> &pnts_ = *cloud.pnts;
    const float sqrd_max_dist = max_dist * max_dist;
    const uint32_t indices_size = cloud.indices.size();
    std::vector<uint32_t> indices_new;
    indices_new.reserve(indices_size);
    
    for(uint32_t i = 0; i < indices_size; ++i){
      const uint32_t idx = cloud.indices[i];
      const float dist = pnts_[idx].x * pnts_[idx].x +
                         pnts_[idx].y * pnts_[idx].y +
                         pnts_[idx].z * pnts_[idx].z;
      if(dist < sqrd_max_dist)
        indices_new.push_back(idx);
    }

    cloud.indices = indices_new; //TODO - don't copy the whole vector, maybe indices should be a pointer instead
  }


  template <typename CloudType>
  inline void CloudsFromIndices(std::vector<typename CloudType::PointType> &src_pnts,
                                std::vector< std::vector<uint32_t> > indices_dvec,
                                std::vector<CloudType> &cloud_vec, const bool copy_data = true){
    const uint32_t num_indice_vec = indices_dvec.size();
    cloud_vec.resize(num_indice_vec);
    for(uint32_t i = 0; i < num_indice_vec; ++i){
      cloud_vec[i].SetData(&src_pnts);
      if(copy_data)
        cloud_vec[i].indices.assign(indices_dvec[i].begin(), indices_dvec[i].end());
      else
        cloud_vec[i].indices.swap(indices_dvec[i]);
    }
  }


  //copy/swap all the contents of cld_a into cld_b EXCEPT the 'pnts' and 'colors' std::vector pointers
  template <typename CloudType>
  inline void CopySwapCloud(CloudType &cld_a, CloudType &cld_b, const bool vec_swap = true){
    cld_b.centroid = cld_a.centroid;
    cld_b.centroid_radius = cld_a.centroid_radius;
    cld_b.origin_radius = cld_a.origin_radius;
    if(!cld_a.eigenvalues.empty())
      cld_b.eigenvalues = cld_a.eigenvalues.clone();
    if(!cld_a.eigenvectors.empty())
      cld_b.eigenvectors = cld_a.eigenvectors.clone();
    cld_b.dim = cld_a.dim;
    cld_b.min = cld_a.min;
    cld_b.max = cld_a.max;
    cld_b.draw_cloud = cld_a.draw_cloud;
    cld_b.pnts = cld_a.pnts;
    cld_b.colors = cld_a.colors;
    std::swap(cld_a.indices, cld_b.indices);
    std::swap(cld_a.local_idx_lu, cld_b.local_idx_lu);
  }

}

#endif //__PCL_CORE_INL_H__

