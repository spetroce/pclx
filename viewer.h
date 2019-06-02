#ifndef __PCLX_VIEWER_H__
#define __PCLX_VIEWER_H__

// in order to get function prototypes from glext.h, define GL_GLEXT_PROTOTYPES before including glext.h
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>
#include <GL/glext.h>
#include <QGLViewer/qglviewer.h>
#include <qapplication.h>
#include <qtextstream.h>
#include <QFile>
#include <pthread.h>
#include "pclx/core.h"
#include "mio/altro/casting.h"
#include "mio/math/math.h"

static color3ub_t pclx_axis_colors[3] = { color3ub_t(255, 100, 160),
                                          color3ub_t(160, 100, 255),
                                          color3ub_t(100, 255, 160) };

#define GL_FLOAT_CAST(float_data) (GLfloat *)(float_data)
#define GL_UBYTE_CAST(ubyte_data) (GLubyte *)(ubyte_data)

#define GL_PUSH_DUMMY_NAME\
  glPushName(-1);\
  glBegin(GL_POINTS);\
  glEnd();\
  glPopName();

#define GL_PUSH_PNT_NAME(id_, draw_code_)\
  glPushName(id_);\
  glBegin(GL_POINTS);\
  draw_code_;\
  glEnd();\
  glPopName();


namespace pclx{

/*
class Viewer : public QGLViewer{
  public:
    void init();
    //void keyPressEvent(QKeyEvent *e);
    void draw();
    void drawWithNames();
    void postSelection(const QPoint& point);
    QString helpString() const{
      QString text("<h2>pclx::Viewer</h2>");
      text += "Use the mouse to move the camera around the object. ";
      text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
      text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
      text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
      text += "Simply press the function key again to restore it. Several keyFrames define a ";
      text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
      text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
      text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
      text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
      text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
      text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
      text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
      text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
      text += "Press <b>Escape</b> to exit the viewer.";
      return text;
    }

  private:
    qglviewer::Vec selected_point;
};

void pclx::Viewer::init(){

}

void pclx::Viewer::draw(){

}

void pclx::Viewer::drawWithNames(){

}

void pclx::Viewer::postSelection(const QPoint& point){

}
*/


//draws the column vectors of a rotation matrix
inline void DrawRotationMatrix(const cv::Mat &R, const cv::Mat &T, const double length = 10){
  STD_INVALID_ARG_E(R.rows == 3 && R.cols == 3)
  STD_INVALID_ARG_E(T.rows == 1 && T.cols == 3 || T.rows == 3 && T.cols == 1)
  cv::Mat T_, vec(3, 1, CV_32F);
  if(T.rows == 3)
    T.convertTo(T_, CV_32F);
  else
    cv::Mat( T.t() ).convertTo(T_, CV_32F);
  const GLfloat *vec_data = mio::StaticCastPtr<GLfloat>(vec.data),
                *T_data = mio::StaticCastPtr<GLfloat>(T_.data);
  glBegin(GL_LINES);
  for(uint32_t i = 0; i < 3; ++i){
    cv::normalize(R.col(i), vec, length, 0, cv::NORM_L2, CV_32F);
    vec += T_;
    glColor3ubv( (GLubyte *)(pclx_axis_colors+i) );
    glVertex3fv(T_data);
    glVertex3fv(vec_data);
  }
  glEnd();
}


inline void GenerateColorMap(std::vector<color4f_t> &color_vec, const uint32_t num_colors){
  STD_INVALID_ARG_E(num_colors > 0)
  const float step_size = 255.0f / num_colors;

  cv::Mat gray_gradient(num_colors, 1, CV_8U), color_map;
  const uint32_t temp = num_colors-1;
  gray_gradient.at<uint8_t>(0) = 0;
  for(uint32_t i = 1; i < temp; ++i)
    gray_gradient.at<uint8_t>(i) = step_size*i;
  gray_gradient.at<uint8_t>(temp) = 255;
  cv::applyColorMap(gray_gradient, color_map, cv::COLORMAP_JET);

  color_vec.resize(num_colors);
  const color3ub_t *data = mio::StaticCastPtr<color3ub_t>(color_map.data);
  const float scale = 1.0f / 255.0f;
  for(uint32_t i = 0; i < num_colors; ++i)
    color_vec[i] = color4f_t(data[i].r * scale, data[i].g * scale, data[i].b * scale, 1);
}


template <typename PointType>
inline void DrawPlane(const PointType (&plane_vertices)[4]){
  glColor4ub(255, 100, 160, 127);
  glBegin(GL_QUADS);
  for(uint8_t i = 0; i < 4; ++i)
    glVertex3f(plane_vertices[i].x, plane_vertices[i].y, plane_vertices[i].z);
  glEnd();
}


template <typename PointType>
inline void DrawCircle(const PointType &cent, const std::vector<PointType> &circ_vertices,
                       const bool repeat_first = true){
  if(circ_vertices.size() > 2){
    glColor4ub(255, 100, 160, 127);
	  glBegin(GL_TRIANGLE_FAN);
		  glVertex3f(cent.x, cent.y, cent.z);
      for(const auto &vert : circ_vertices)
        glVertex3f(vert.x, vert.y, vert.z);
      if(repeat_first)
        glVertex3f(circ_vertices[0].x, circ_vertices[0].y, circ_vertices[0].z);
	  glEnd();
  }
}


template <typename PointType>
inline void GeneratePointSphere(std::vector<PointType> &pnt, const float radius, const PointType &center){
  pnt.clear();
  for(int i = -90; i < 90; i+=3){
    const float theta = sm::DegToRad( static_cast<float>(i) );
    for(int j = 0; j < 360; j+=3){
      const float phi = sm::DegToRad( static_cast<float>(j) );
      pnt.push_back( PointType(radius * std::cos(theta) * std::cos(phi) + center.x,
                               radius * std::cos(theta) * std::sin(phi) + center.y,
                               radius * std::sin(theta) + center.z) );
    }
  }
}


template <typename PointType>
void DrawMarker(const PointType &pnt, const float line_length, const color4f_t *line_color){
  glColor3fv( GL_FLOAT_CAST(line_color) );
  const float l = line_length * 0.5f;

  glBegin(GL_LINES);
    glVertex3f(pnt.x - l, pnt.y, pnt.z);
    glVertex3f(pnt.x + l, pnt.y, pnt.z);
    glVertex3f(pnt.x, pnt.y - l, pnt.z);
    glVertex3f(pnt.x, pnt.y + l, pnt.z);
    glVertex3f(pnt.x, pnt.y, pnt.z - l);
    glVertex3f(pnt.x, pnt.y, pnt.z + l);
  glEnd();
}


template <typename CloudType>
void DrawCloud(const CloudType &cloud, color4f_t *static_color = nullptr, 
               const int pix_per_pnt = 1, const bool draw_with_names = false){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  if(cloud.pnts->size() < 1)
    return;
  if(!cloud.draw_cloud || cloud.indices.size() < 1)
    return;

  typename CloudType::PointType *pnts_data = cloud.pnts->data();
  typename CloudType::ColorType *colors_data;

  color4f_t error_color = red4f,
            default_color = white4f;
  if(static_color == nullptr)
    if(cloud.colors != nullptr){
      if( cloud.colors->size() != cloud.pnts->size() ){
        static_color = &error_color;
        colors_data = cloud.colors->data();
      }
    }
    else
      static_color = &default_color;

  glPointSize(pix_per_pnt);
  
  if(draw_with_names){ //draw with names (names are the values from cloud.indices)
    GL_PUSH_DUMMY_NAME
    if(static_color){ //draw all points with same color
      glColor4fv( GL_FLOAT_CAST(static_color) );
      for(auto &i : cloud.indices){
        GL_PUSH_PNT_NAME( i, glVertex3fv( GL_FLOAT_CAST(pnts_data+i) ) )
      }
    }
    else{ //draw points with supplied color
      for(auto &i : cloud.indices){
        glColor3ubv( GL_UBYTE_CAST(colors_data+i) );
        GL_PUSH_PNT_NAME( i, glVertex3fv( GL_FLOAT_CAST(pnts_data+i) ) )
      }
    }
  }
  else{ //draw without names
    glBegin(GL_POINTS);
      if(static_color){ //draw with static color
        glColor4fv( GL_FLOAT_CAST(static_color) );
        for(auto &i : cloud.indices)
          glVertex3fv( GL_FLOAT_CAST(pnts_data+i) );
      }
      else //draw with supplied color
        for(auto &i : cloud.indices){
          glColor3ubv( GL_UBYTE_CAST(colors_data+i) );
          glVertex3fv( GL_FLOAT_CAST(pnts_data+i) );
        }
    glEnd();
  }
}


//use this for drawing clusters with names
template <typename CloudType>
void DrawClusters(const std::vector<CloudType> &cloud_vec,
                  const bool draw_with_names, const int selected_name_id,
                  const color4f_t *select_color, const std::vector<color4f_t> &color_vec, 
                  const int pix_per_pnt, const int min_pnt_threshold){
  const uint32_t color_vec_size = color_vec.size();
  uint32_t i = 0;
  if(draw_with_names){
    GL_PUSH_DUMMY_NAME
    for(auto &cloud : cloud_vec){
      if(cloud.indices.size() > min_pnt_threshold){
        glPushName(i);
        DrawCloud(cloud, &color_vec[i % color_vec_size], pix_per_pnt, false);
        glPopName();
      }
      ++i;
    }
  }
  else
    for(auto &cloud : cloud_vec)
      if(cloud.indices.size() > min_pnt_threshold)
        if(i == selected_name_id) //draw cloud with select_color
          DrawCloud(cloud, select_color, pix_per_pnt, false);
        else
          DrawCloud(cloud, &color_vec[i % color_vec_size], pix_per_pnt, false);
}


//simpler version of the above func
template <typename CloudType>
void DrawClusters(const std::vector<CloudType> &cloud_vec, const std::vector<color4f_t> &color_vec){
  const uint32_t color_vec_size = color_vec.size();
  uint32_t i = 0;
  for(auto &cloud : cloud_vec){
    color4f_t color = color_vec[i % color_vec_size];
    DrawCloud(cloud, &color, 2);
    ++i;
  }
}


template <typename CloudType>
void DrawNN(CloudType &cloud, const pclx::table32S_t &nn_table, 
            const int selected_name_id, const float marker_line_length, const color4f_t *marker_color){
  const auto &pnts = *cloud.pnts;
  if(selected_name_id < 0)
    return;

  if(cloud.local_idx_lu.size() == 0)
    cloud.BuildLocalIdxLU();
  const auto &nn_vec = nn_table[ cloud.local_idx_lu[selected_name_id] ];
  for(auto &idx : nn_vec)
    DrawMarker(pnts[idx], marker_line_length, marker_color);
}


template <typename PointType>
inline void DrawPntEigenvector(const PointType &pnt, const cv::Mat &eigenvectors, const float &scale,
                               const bool &draw_norm_only = false, const bool &call_gl_begin_end = false){
  STD_INVALID_ARG_E( (std::is_same<typename PointType::DataType, GLfloat>::value) )
  STD_INVALID_ARG_E(eigenvectors.rows == 3 && eigenvectors.cols == 3)
  STD_INVALID_ARG_E(eigenvectors.depth() == CV_32F)

  if(call_gl_begin_end)
    glBegin(GL_LINES);

  const GLfloat *ev_data = mio::StaticCastPtr<GLfloat>(eigenvectors.data);

  if(!draw_norm_only){
    glColor3ubv( (GLubyte *)(pclx_axis_colors) );
    glVertex3f(pnt.x, pnt.y, pnt.z);
    glVertex3f(pnt.x + ev_data[0] * scale,
               pnt.y + ev_data[1] * scale,
               pnt.z + ev_data[2] * scale);
    glColor3ubv( (GLubyte *)(pclx_axis_colors+1) );
    glVertex3f(pnt.x, pnt.y, pnt.z);
    glVertex3f(pnt.x + ev_data[3] * scale,
               pnt.y + ev_data[4] * scale,
               pnt.z + ev_data[5] * scale);
  }
  glColor3ubv( (GLubyte *)(pclx_axis_colors+2) );
  glVertex3f(pnt.x, pnt.y, pnt.z);
  glVertex3f(pnt.x + ev_data[6] * scale,
             pnt.y + ev_data[7] * scale,
             pnt.z + ev_data[8] * scale);

  if(call_gl_begin_end)
    glEnd();
}


template <typename CloudType>
void DrawPntEigenvectors(const CloudType &cloud, const std::vector<cv::Mat> &eigenvalues, 
                         const std::vector<cv::Mat> &eigenvectors, const float scale,
                         const bool draw_norm_only = false){
  STD_INVALID_ARG_E(cloud.pnts != nullptr)
  const auto &pnts = *cloud.pnts;
  glBegin(GL_LINES);
  const uint32_t num_pnts = cloud.indices.size();
  for(uint32_t i = 0; i < num_pnts; ++i)
    DrawPntEigenvector(pnts[ cloud.indices[i] ], eigenvectors[i], scale, draw_norm_only);
  glEnd();
}


template <typename CloudType>
void DrawCloudEigenvectors(CloudType &cloud, const float scale_coeff = 0.1f){
  STD_INVALID_ARG_E(cloud.eigenvectors.depth() == CV_32F && cloud.eigenvalues.depth() == CV_32F)
  STD_INVALID_ARG_EM(cloud.eigenvectors.rows == 3 && cloud.eigenvectors.cols == 3,
                     std::string("size: ") + mio::SizeStr(cloud.eigenvectors))
  STD_INVALID_ARG_EM(cloud.eigenvalues.rows == 3 && cloud.eigenvalues.cols == 1,
                     std::string("size: ") + mio::SizeStr(cloud.eigenvalues))

  cv::Mat vec1, vec2(1, 3, CV_32F),
          centroid( 1, 3, CV_32F, mio::StaticCastPtr<uint8_t>(&cloud.centroid) ),
          scale_vec = cloud.eigenvalues * scale_coeff * 0.5f;
  STD_INVALID_ARG_E(scale_vec.depth() == CV_32F)
  const GLfloat *vec2_data = mio::StaticCastPtr<GLfloat>(vec2.data);

  glBegin(GL_LINES);
  for(uint32_t i = 0; i < 3; ++i){
    glColor3ubv( (GLubyte *)(pclx_axis_colors+i) );
    vec1 = cloud.eigenvectors.row(i) * scale_vec.at<float>(i);
    vec2 = centroid + vec1;
    glVertex3fv(vec2_data);
    vec2 = centroid - vec1;
    glVertex3fv(vec2_data);
  }
  glEnd();
}


template <typename CloudType>
void DrawCloudEigenvectors(const std::vector<CloudType> &cloud_vec){
  for(auto &cloud : cloud_vec)
    DrawCloudEigenvectors(cloud);
}


template <typename CloudType>
void DrawCentroid(CloudType &cloud, int pix_per_pnt, color4f_t color){
  glPointSize(pix_per_pnt);
  glColor3f(color.r, color.g, color.b);
  glBegin(GL_POINTS);
    glVertex3fv( (GLfloat *)&cloud.centroid );
  glEnd();
}

} //namespace pclx

#endif //__PCLX_VIEWER_H__
