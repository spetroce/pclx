#ifndef __PCL_IO_INL_H__
#define __PCL_IO_INL_H__

#include <sys/types.h>
#include <fstream>
#include <iostream>
#include "mio/altro/io.hpp"
#include "mio/altro/casting.hpp"
#include "pclx/core.h"


namespace pclx{

  template <typename CloudType>
  void PrintPointCloud(const CloudType &cloud){
    STD_INVALID_ARG_E(cloud.pnts != nullptr)
    PCL_LCL_PNTS_COLORS(cloud)
    for(auto &idx : cloud.indices)
      std::cout << pnts[idx] << " ";
    std::cout << "\n";
  }


  template <typename CloudType>
  void ImportBPT(const std::string &file_full, CloudType &cloud){
    STD_INVALID_ARG_E(cloud.pnts != nullptr)
    PCL_LCL_PNTS_COLORS(cloud)

    std::ifstream file;
    file.open(file_full, std::ios::binary);
    if( !file.is_open() ){
      printf( "ImportBPT() - file %s could not be opened\n", file_full.c_str() );
      return;
    }

    // get length of file:
    file.seekg (0, std::ios::end);
    int num_byte = file.tellg();
    file.seekg (0, std::ios::beg);

    int rgba_size = sizeof(typename CloudType::ColorType), pnt_size = sizeof(typename CloudType::PointType);
    if( (num_byte % (rgba_size + pnt_size)) != 0 ){
      printf("ImportBPT() - point list error, (num_bytes %% scan_point_size) != 0\n");
      return;
    }

    const uint32_t num_pnt = num_byte/(rgba_size + pnt_size);
    pnts = std::vector< typename CloudType::PointType >(num_pnt);
    if(cloud.colors != nullptr)
      colors = std::vector< typename CloudType::ColorType >(num_pnt);

    typename CloudType::ColorType color_temp;
    for(uint32_t i = 0; i < num_pnt; ++i){
      file.read( (char *)&color_temp, rgba_size );
      file.read( (char *)&pnts[i], pnt_size );
      if(cloud.colors != nullptr){
        colors[i].r = color_temp.r;
        colors[i].g = color_temp.g;
        colors[i].b = color_temp.b;
      }
    }

    cloud.DefaultIndices();
    file.close();
  }


  template <typename CloudType>
  void SaveBPT(const std::string &file_full, const CloudType &cloud){
    STD_INVALID_ARG_E(cloud.pnts != nullptr)
    PCL_LCL_PNTS_COLORS(cloud);

    std::ofstream file;
    file.open(file_full, std::ios::binary);
    if( !file.is_open() ){
      printf( "SaveBPT() - file %s could not be opened\n", file_full.c_str() );
      return;
    }

    typedef typename CloudType::PointType PointType;
    typedef typename CloudType::ColorType ColorType;

    typedef struct SBPT{
      PointType pnt;
      ColorType color;
      SBPT() {};
      SBPT(PointType pnt_, ColorType color_) : pnt(pnt_), color(color_) {};
    } __attribute__((packed)) bpt_t;

    const uint32_t indices_size = cloud.indices.size();
    std::vector<bpt_t> bpt_vec(indices_size);

    ColorType default_color;
    if(cloud.colors)
      for(auto &i : cloud.indices)
        bpt_vec[i] = bpt_t(pnts[i], colors[i]);
    else
      for(auto &i : cloud.indices)
        bpt_vec[i] = bpt_t(pnts[i], default_color);

    file.write(mio::StaticCastPtr<char>( bpt_vec.data() ), sizeof(bpt_t) * indices_size);
    file.close();
  }


  template <typename CloudType>
  void ImportAsciiXYZ(const std::string &file_full, CloudType &cloud){
    STD_INVALID_ARG_E(cloud.pnts != nullptr)
    PCL_LCL_PNTS_COLORS(cloud);

    std::ifstream file;
    file.open(file_full);
    if( !file.is_open() ){
      printf( "ImportAsciiXYZ() - file %s could not be opened\n", file_full.c_str() );
      return;
    }

    std::string line;
    float x, y, z;
    while( std::getline(file, line) ){
      std::istringstream iss(line);
      if(iss >> x >> y >> z)
        pnts.push_back( typename CloudType::PointType(x, y, z) ); 
    }

    cloud.DefaultIndices();
    file.close();
  }


  template <typename CloudType>
  void SaveAsciiXYZ(const std::string &file_full, const CloudType &cloud){
    STD_INVALID_ARG_E(cloud.pnts != nullptr)
    PCL_LCL_PNTS_COLORS(cloud);

    std::ofstream file;
    file.open(file_full);
    if( !file.is_open() ){
      printf( "ImportAsciiXYZ() - file %s could not be opened\n", file_full.c_str() );
      return;
    }

    for(auto &i : cloud.indices)
      file << pnts[i].x << " " << pnts[i].y << " " << pnts[i].z << " " << std::endl;

    file.close();
  }


//  int loadCloud(const std::string fileFull, CPointCloud &srcCloud, urgentData_t *urgentData = NULL){
//    ERR_CHK( !srcCloud.pnt || !srcCloud.color, "loadCloud() - !srcCloud.pnt || !srcCloud.color", return(-1) );
//    std::string fileExt;
//    int nNumSrcPnt = 0;
//    FileNameExpand(fileFull, ".", NULL, NULL, NULL, &fileExt);
//    
//    if(fileExt == "bpt"){
//      ERR_CHK( ( nNumSrcPnt = loadBPT(fileFull.c_str(), srcCloud) ) == -1, 
//                                         "nLoadCloud() - nLoadBPT() - error", return(-1) );
//    }
//    else if(fileExt == "upc"){
//      ERR_CHK(loadOttawaTile(fileFull.c_str(), urgentData, false) == -1, 
//                              "nLoadCloud() - nLoadOttawaTile() - error", return(-1) );
//      ERR_CHK( ( nNumSrcPnt = extractOttawaVertices(srcCloud, urgentData) ) == -1,
//                                                     "nLoadCloud() - nLoadOttawaTile() - error", return(-1) );
//    }
//    else{
//      printf( "on file: %s - unrecognized file extension: %s\n", fileFull.c_str(), fileExt.c_str() );
//      return(-1);
//    }
//    printf( "Opened cloud: %s\n", fileFull.c_str() );
//    
//    ERR_CHK(nNumSrcPnt <= 0, "nLoadCloud() - invalid cloud size", return(-1) );
//    srcCloud.clear();
//    srcCloud.srcIdx.resize(nNumSrcPnt);
//    for(uint32_t i = 0; i < nNumSrcPnt; ++i)
//      srcCloud.srcIdx[i] = i;
//    
//    return(nNumSrcPnt);
//  }


//  #include <QDebug>
//  int saveCameraState(const char *dir_path, const char *fileName, qglviewer::Camera *myCamera){
//    char fileName_ext[128], fileName_num[128];
//    int i = 1;
//    
//    sprintf(fileName_num, "%s", fileName);
//    sprintf(fileName_ext, "%s%s.xml", dir_path, fileName_num);

//    while( bFileExists(fileName_ext) ){
//      sprintf(fileName_num, "%s%d", fileName, i);
//      sprintf(fileName_ext, "%s%s.xml", dir_path, fileName_num);
//      ++i;
//    }
//    
//    QDomDocument document(fileName_num);
//    document.appendChild( myCamera->domElement("Camera", document) );

//    QFile f(fileName_ext);
//    #if QT_VERSION < 0x040000
//    if( f.open(IO_ReadWrite) ){
//    #else
//    if( f.open(QIODevice::ReadWrite) ){
//    #endif
//      QTextStream out(&f);
//      document.save(out, 2);
//    }
//    else{
//      qDebug() << "cloudViewer: nSaveCameraState() - open() QFile::error() == " << f.error() << "-" << f.errorString();
//      return(-1);
//    }
//    
//    return(0);
//  }

} //namespace pc;

#endif //__PCL_IO_INL_H__
