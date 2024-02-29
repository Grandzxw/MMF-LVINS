#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <pcl/io/pcd_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include "Spatialcontext/Spatialcontext.h"
#include "omp.h"

using namespace std;


// 08 -- 4071
//05 -- 2760
//00 -- 4540

// number of sequence
const int N = 4071;

// kitti sequence
const string seq = "08";


pcl::PointCloud<pcl::PointXYZI>::Ptr getCloud(std::string file){
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud;
    int32_t num = 1000000;
    float *data = (float*)malloc(num*sizeof(float));
    float *px = data+0;
    float *py = data+1;
    float *pz = data+2;
    float *pr = data+3;
    FILE *stream;
    stream = fopen (file.c_str(),"rb");
    if(stream==NULL){
        std::cerr<<"Stream is NULL!"<<std::endl;
        return NULL;
    }
    cloud.reset(new pcl::PointCloud<pcl::PointXYZI>());
    num = fread(data,sizeof(float),num,stream)/4;
    for (int32_t i=0; i<num; i++) {
        pcl::PointXYZI p;
        p.x=(*px);
        p.y=(*py);
        p.z=(*pz);
        p.intensity=(*pr);
        cloud->points.push_back(p);
        px+=4; py+=4; pz+=4; pr+=4;
    }
    cloud->height=1;
    cloud->width=cloud->points.size();
    fclose(stream);
    free(data);
    return cloud;
}

int main(int argc, char *argv[])
{

    std::cout << "number of available processors: " << omp_get_num_procs()<< std::endl;
    std::cout << "number of threads: " << omp_get_max_threads() << std::endl;
    omp_set_num_threads(4);//设置线程数，一般设置的线程数不超过CPU核心数
    std::ofstream ofs("/home/grand/place_recognizate/Spatial_Context/results/" + seq+".txt");
    SCManager sc_manager;

    //#pragma omp paralell for
    for(int i =0; i <=N-1 ;i++)
    {
        // Eigen::Vector3d a;
        // a << 1.0, 2.0, 3.0;
        // cout << a.transpose() << endl;
        // a.normalize();
        // cout << a.transpose() << endl;


        std::stringstream ss;
        ss << setw(6) << setfill('0') << i;
        //cout << ss.str()+".bin" << std::endl;

        // kitti velodyne bins
        std::string filename = "/home/grand/place_recognizate/KITTI_DATA/" + seq + "/velodyne/" + ss.str() + ".bin";
        // cout<<"filename: "<<filename<<endl;
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud0(new pcl::PointCloud<pcl::PointXYZI>);
        cloud0 = getCloud(filename);
        // std::fstream input(filename, std::ios::in | std::ios::binary);
        // input.seekg(0, std::ios::beg);
        // for (int ii=0; input.good() && !input.eof(); ii++) {
        //     pcl::PointXYZI point;
        //     input.read((char *) &point.x, 3*sizeof(float));
        //     float intensity;
        //     input.read((char *) &intensity, sizeof(float));
        //     cloud0->push_back(point);
        // }
        //std::cout<<"point size: "<<cloud0->size()<<std::endl;
        sc_manager.makeAndSaveSpatialcontextAndKeys(*cloud0);
      
        float mindis = 1000;
        int loop_id = -1;

        auto sc_res = sc_manager.detectLoopClosureID();
        if(i>30)
        {
            loop_id = std::get<0>(sc_res);
        }
        if(loop_id != -1)
        {
            ofs << i <<" " << loop_id << " " << std::get<2>(sc_res)  << " " << 1 << std::endl;
            //cout << i << " " << loop_id << " " << std::get<2>(sc_res) << " " << 1 << std::endl;
        }
        else
        {
            ofs<< i << " " << loop_id << " " << std::get<2>(sc_res) << " " << 0 << std::endl;
            //cout << i << " " << loop_id << " " << std::get<2>(sc_res) << " " << 0 << std::endl; 
        }
    }

    return 0;
}