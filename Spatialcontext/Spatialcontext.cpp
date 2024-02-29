#include "Spatialcontext.h"


//输入
void coreImportTest (void)
{
    cout << "scancontext lib is successfully imported." << endl;
} // coreImportTest

//转换为角度值
float rad2deg(float radians)
{
    return radians * 180.0 / M_PI;
}

//转换为弧度值
float deg2rad(float degrees)
{
    return degrees * M_PI / 180.0;
}

//求x,y的tan角度，放在 0度到360度
float xy2theta( const float & _x, const float & _y )
{
    if ( (_x >= 0) & (_y >= 0)) 
        return (180/M_PI) * atan(_y / _x);

    if ( (_x < 0) & (_y >= 0)) 
        return 180 - ( (180/M_PI) * atan(_y / (-_x)) );

    if ( (_x < 0) & (_y < 0)) 
        return 180 + ( (180/M_PI) * atan(_y / _x) );

    if ( (_x >= 0) & (_y < 0))
        return 360 - ( (180/M_PI) * atan((-_y) / _x) );
} // xy2theta


//将_mat的列平移_num_shift形成新的MatrixXd
MatrixXd circshift( MatrixXd &_mat, int _num_shift )
{
    // shift columns to right direction 
    assert(_num_shift >= 0);
    if( _num_shift == 0 )
    {
        MatrixXd shifted_mat( _mat );
        return shifted_mat; // Early return 
    }
    MatrixXd shifted_mat = MatrixXd::Zero( _mat.rows(), _mat.cols() );
    for ( int col_idx = 0; col_idx < _mat.cols(); col_idx++ )
    {
        int new_location = (col_idx + _num_shift) % _mat.cols();
        shifted_mat.col(new_location) = _mat.col(col_idx);
    }
    return shifted_mat;

} // circshift


//这段代码的作用是将一个 Eigen 矩阵转换为一个包含同样元素的std::vector<float>向量
std::vector<float> eig2stdvec( MatrixXd _eigmat )
{
    std::vector<float> vec( _eigmat.data(), _eigmat.data() + _eigmat.size() );
    return vec;
} // eig2stdvec


//cosdistance的距离
double cosdistance ( MatrixXd &_dsc1, MatrixXd &_dsc2 )
{
    int num_eff_cols = 0;
    double sum_sector_similarity = 0;

    for ( int col_idx = 0; col_idx < _dsc1.cols(); col_idx++ )
    {
        VectorXd col_dsc1 = _dsc1.col(col_idx);
        VectorXd col_dsc2 = _dsc2.col(col_idx);

        if( col_dsc1.norm() == 0 | col_dsc2.norm() == 0 )
            continue;
        double sector_similarity = col_dsc1.dot(col_dsc2) / (col_dsc1.norm() * col_dsc2.norm());

        sum_sector_similarity = sum_sector_similarity + sector_similarity;
        num_eff_cols = num_eff_cols + 1;
    }

    double sc_sim = sum_sector_similarity / num_eff_cols;
    return 1.0 - sc_sim;
}


//相关性系数 pc_sim越接近于1相似度越高， 1.0 - pc_sim越小相似度越高
double Pearson_Correlation(MatrixXd &_dsc1, MatrixXd &_dsc2){
    int num_eff_cols = 0;
    double sum_sector_similarity = 0;
    for (int col_idx = 0; col_idx < _dsc1.cols(); col_idx++) {
        Eigen::VectorXd col_dsc1 = _dsc1.col(col_idx);
        Eigen::VectorXd col_dsc2 = _dsc2.col(col_idx);
        // if (col_dsc1.norm() == 0 || col_dsc2.norm() == 0) {
        //     continue;
        // }
        double mean_dsc1 = col_dsc1.mean();
        double mean_dsc2 = col_dsc2.mean();
        double cov = ((col_dsc1.array() - mean_dsc1) * (col_dsc2.array() - mean_dsc2)).sum();
        double var_dsc1 = ((col_dsc1.array() - mean_dsc1) * (col_dsc1.array() - mean_dsc1)).sum();
        double var_dsc2 = ((col_dsc2.array() - mean_dsc2) * (col_dsc2.array() - mean_dsc2)).sum();
        double sector_similarity = cov / (std::sqrt(var_dsc1) * std::sqrt(var_dsc2));
        sum_sector_similarity += sector_similarity;
        num_eff_cols++;
    }
    double pc_sim = std::abs(sum_sector_similarity / num_eff_cols);

    return 1.0 - pc_sim;
}



//wasserstein_dist越小，相似性较高，wasserstein_sim越小相似度越高
double wasserstein_distance (MatrixXd &_dsc1, MatrixXd &_dsc2){
    double wasserstein_dist = 0.0;
    int num_eff_cols = 0;
    for ( int col_idx = 0; col_idx < _dsc1.cols(); col_idx++ ){
        VectorXd col_dsc1 = _dsc1.col(col_idx);
        VectorXd col_dsc2 = _dsc2.col(col_idx);
        double diff = (col_dsc1 - col_dsc2).array().abs().sum();
        wasserstein_dist += diff;
        num_eff_cols++;
    }
    double wasserstein_dist_mean = wasserstein_dist / num_eff_cols;
    //double wasserstein_sim = 1 - (1.0 / (1.0 + wasserstein_dist_mean)); 
    return wasserstein_dist_mean;
}



// //得到每个bin的状态两种状态值
std::tuple<vector<vector<vector<int>>>,vector<vector<vector<double>>>, vector<vector<vector<double>>>> SCManager::makeSpatialState(pcl::PointCloud<SCPointType>& _scan_down){
    float azim_angle, azim_range, azim_vert;
    int ring_idx, sctor_idx,vert_idx;
    int num_pts_scan_down = _scan_down.points.size();
    vector<vector<vector<int>>> density_state(PC_NUM_VERTICAL,vector<vector<int>>(PC_NUM_RING,vector<int>(PC_NUM_SECTOR,0)));
    vector<vector<vector<double>>> intensity_state(PC_NUM_VERTICAL,vector<vector<double>>(PC_NUM_RING,vector<double>(PC_NUM_SECTOR,0)));
    vector<vector<vector<double>>> height_state(PC_NUM_VERTICAL,vector<vector<double>>(PC_NUM_RING,vector<double>(PC_NUM_SECTOR,0)));

    //vector<vector<vector<vector<SCPointType>>>> point_state(PC_NUM_RING, vector<vector<vector<SCPointType>>>(PC_NUM_SECTOR, vector<vector<SCPointType>>(PC_NUM_VERTICAL)));
    //视场角的角度大小是[-24,2]
    for(int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)
    {
        SCPointType pt;
        pt.x = _scan_down.points[pt_idx].x;
        pt.y = _scan_down.points[pt_idx].y;
        pt.z = _scan_down.points[pt_idx].z; //+ 1.73; // naive adding is ok (all points should be > 0). ?   先不加看看问题
        pt.intensity = _scan_down.points[pt_idx].intensity;
        //x,y,z to ring sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
        azim_angle = xy2theta(pt.x, pt.y);  //(atan2(pt.y, pt.x) * 180.0f / M_PI) + 180; //xy2theta(pt.x, pt.y); 
        azim_vert = (atan2(pt.z, azim_range) * 180.0f / M_PI) + 24;
        
        ring_idx = std::max( std::min( PC_NUM_RING-1, int(floor( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 0 );
        sctor_idx = std::max( std::min( PC_NUM_SECTOR-1, int(floor( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 0 );
        vert_idx = std::max( std::min( PC_NUM_VERTICAL-1, int(floor( (azim_vert / PC_MAX_HEIGHT) * PC_NUM_VERTICAL )) ), 0 );

        density_state[vert_idx][ring_idx][sctor_idx] = density_state[vert_idx][ring_idx][sctor_idx] + 1;

        double Q_intensity = pt.intensity;
        intensity_state[vert_idx][ring_idx][sctor_idx] = std::max(intensity_state[vert_idx][ring_idx][sctor_idx], Q_intensity);

        double max_height = pt.z;
        height_state[vert_idx][ring_idx][sctor_idx] = std::max(height_state[vert_idx][ring_idx][sctor_idx], max_height);
        //point_state[ring_idx][sctor_idx][vert_idx].push_back(pt);
    }
    return std::make_tuple(density_state,intensity_state,height_state);
}



// //得到描述子
// MatrixXd SCManager::getSpatialDesc( std::pair<vector<vector<vector<int>>>,vector<vector<vector<double>>>>&_istate ) {
//     //定义两个不同
//     MatrixXd scan_descs(PC_NUM_RING,PC_NUM_SECTOR);
//     auto density_state = _istate.first;
//     auto intensity_state = _istate.second;
//     vector<vector<vector<double>>> density_weight(PC_NUM_VERTICAL,vector<vector<double>>(PC_NUM_RING,vector<double>(PC_NUM_SECTOR)));
//     vector<vector<vector<double>>> intensity_weight(PC_NUM_VERTICAL,vector<vector<double>>(PC_NUM_RING,vector<double>(PC_NUM_SECTOR)));

//     //计算点云密度权重 and 强度权重
//     for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
//         for(int ring_idx = 0; ring_idx < PC_NUM_RING; ring_idx++){
//             vector<int> cnt_points;
//             double mean_intensity = 0;
//             for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
//                 cnt_points.push_back(density_state[vert_idx][ring_idx][sect_idx]);
//                 mean_intensity += intensity_state[vert_idx][ring_idx][sect_idx];
//             }
//             mean_intensity /= PC_NUM_SECTOR;
//             std::sort(cnt_points.begin(),cnt_points.end());
//             int median_density = cnt_points[20];

//             for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
//                 //density for point
//                 if (median_density ==0 || density_state[vert_idx][ring_idx][sect_idx] > 2*median_density){
//                     density_weight[vert_idx][ring_idx][sect_idx] = 1;
//                 }else{
//                     density_weight[vert_idx][ring_idx][sect_idx] = (float)density_state[vert_idx][ring_idx][sect_idx]/(2*median_density);
//                 }
//                 //intensity for point
//                 if (mean_intensity ==0 || intensity_state[vert_idx][ring_idx][sect_idx] > 2*mean_intensity){
//                     intensity_weight[vert_idx][ring_idx][sect_idx] = 1;
//                 }else{
//                     intensity_weight[vert_idx][ring_idx][sect_idx] = intensity_state[vert_idx][ring_idx][sect_idx]/(2*mean_intensity);
//                 }
//             }
//         }
//     }
//     //for the bins height
//     for(int ring_idx = 0; ring_idx < PC_NUM_RING; ring_idx++) {
//         for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
//             float element =0;
//             float intensity_element =0;
//             for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
//                 if (density_state[vert_idx][ring_idx][sect_idx]>1){
//                     element+= density_weight[vert_idx][ring_idx][sect_idx] * pow(2,vert_idx);
//                     intensity_element += (density_weight[vert_idx][ring_idx][sect_idx] * pow(2,vert_idx))*intensity_weight[vert_idx][ring_idx][sect_idx];
//                 }
//             }
//             scan_descs(ring_idx,sect_idx) = intensity_element;
//         }
//     }
//     return scan_descs;
// }


//得到每个点的个数bin
// std::pair<vector<vector<vector<int>>>,vector<vector<vector<double>>>> SCManager::makeSpatialState(pcl::PointCloud<SCPointType>  &cloud) {
//     vector<vector<vector<int>>> state(i_rows_,vector<vector<int>>(i_cols_,vector<int>(i_hor_)));
//     vector<vector<vector<double>>> intensity_map(i_rows_,vector<vector<double>>(i_cols_,vector<double>(i_hor_,0)));
// //  64-line LiDAR vertical FOV [-24,2]
//     float h_gap = (float)360/i_hor_;
//     float d_gap = (float)80/i_cols_;
//     float v_gap = (float)32/i_rows_;
//     for (pcl::PointXYZI p : cloud.points)
//     {
//         float dis = sqrt(p.data[0] * p.data[0] + p.data[1] * p.data[1]);
//         float arc = (atan2(p.data[2], dis) * 180.0f / M_PI) + 24;
//         float yaw = (atan2(p.data[1], p.data[0]) * 180.0f / M_PI) + 180;

//         int Q_dis = std::min(std::max((int)floor(dis/d_gap), 0), i_cols_-1);
//         int Q_arc = std::min(std::max((int)floor(arc /v_gap), 0), i_rows_-1);
//         int Q_yaw = std::min(std::max((int)floor(yaw /h_gap), 0), i_hor_-1);
//         state[Q_arc][Q_dis][Q_yaw] =  state[Q_arc][Q_dis][Q_yaw]+1;
//         double Q_intensity = p.intensity;
//         intensity_map[Q_arc][Q_dis][Q_yaw] = std::max(intensity_map[Q_arc][Q_dis][Q_yaw], Q_intensity);
//         //std::cout << intensity_map[Q_arc][Q_dis][Q_yaw] << std::endl;
//     }
//     return make_pair(state,intensity_map);
// }



//空间分布特性
MatrixXd SCManager::getSpatialDesc( std::tuple<vector<vector<vector<int>>>,vector<vector<vector<double>>>,vector<vector<vector<double>>>>&_istate ) {
    
    MatrixXd m(i_cols_,i_hor_);
    auto _state = std::get<0>(_istate);
    auto _intensity = std::get<1>(_istate);
    auto _height = std::get<2>(_istate);

    //for point density
    // vector<vector<vector<float>>> weight(i_rows_,vector<vector<float>>(i_cols_,vector<float>(i_hor_)));
    // //calculate the point density weight
    // for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
    //     for (int col_idx = 0; col_idx < i_cols_; col_idx++) {
    //         vector<int> cnts;
    //         double mean_cnt;
    //         for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
    //             cnts.push_back(_state[row_idx][col_idx][hor_idx]);
    //             mean_cnt += _state[row_idx][col_idx][hor_idx];
    //         }
    //         mean_cnt /= i_hor_;
    //         std::sort(cnts.begin(),cnts.end());
    //         int median = cnts[20];
    //         //cout << "median" << median << endl;
    //         for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
    //             if (median ==0 || _state[row_idx][col_idx][hor_idx] > 2*median || _state[row_idx][col_idx][hor_idx] == 0){
    //                 weight[row_idx][col_idx][hor_idx] = 1;
    //             }
    //             else
    //                 weight[row_idx][col_idx][hor_idx] = (float)_state[row_idx][col_idx][hor_idx]/(2*median);
    //             //std::cout << "_state is : " << _state[row_idx][col_idx][hor_idx] << std::endl;
    //         }
    //     }
    // }


    //for point density
    vector<vector<vector<float>>> weight(i_rows_,vector<vector<float>>(i_cols_,vector<float>(i_hor_)));
    //calculate the point density weight
    for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
        for (int col_idx = 0; col_idx < i_cols_; col_idx++) {
            vector<int> cnts;
            double mean_cnt;
            for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
                cnts.push_back(_state[row_idx][col_idx][hor_idx]);
                mean_cnt += pow(_state[row_idx][col_idx][hor_idx],2);
            }
            mean_cnt = sqrt(mean_cnt);
            std::sort(cnts.begin(),cnts.end());
            int median = cnts[20];
            //cout << "median" << median << endl;
            for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
                if (mean_cnt ==0 || _state[row_idx][col_idx][hor_idx] <= mean_cnt / 20.0){
                    weight[row_idx][col_idx][hor_idx] = 1;
                }
                else
                    weight[row_idx][col_idx][hor_idx] = (float)_state[row_idx][col_idx][hor_idx]/(mean_cnt);
                //std::cout << "_state is : " << _state[row_idx][col_idx][hor_idx] << std::endl;
            }
        }
    }


    //calculate the point intensity weight  row行中
    // vector<vector<vector<double>>> intensity_weight(i_rows_,vector<vector<double>>(i_cols_,vector<double>(i_hor_)));
    // for ( int col_idx = 0; col_idx < i_cols_; col_idx++ ) {
    //     for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
    //         vector<double> h_intensity;
    //         double mean_intensity = 0;
    //         for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
    //             //cnts.push_back(_intensity[row_idx][col_idx][hor_idx]);
    //             mean_intensity += _intensity[row_idx][col_idx][hor_idx];
    //         }
    //         mean_intensity /= i_rows_;
    //         //std::cout << "median is : " << mean_intensity<< std::endl;
    //         for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
    //             //std::cout << "intensity_element is : " << _intensity[row_idx][col_idx][hor_idx] << std::endl;
    //             if (mean_intensity ==0 || _intensity[row_idx][col_idx][hor_idx] > 2*mean_intensity || _intensity[row_idx][col_idx][hor_idx] == 0){
    //                 intensity_weight[row_idx][col_idx][hor_idx] = 1;
    //             }
    //             else
    //             {
    //                 intensity_weight[row_idx][col_idx][hor_idx] = _intensity[row_idx][col_idx][hor_idx]/(2*mean_intensity);
    //                 //std::cout << "intensity_weight is : " << intensity_weight[row_idx][col_idx][hor_idx] << std::endl;
    //                 //std::cout << "_intensity is : " << _intensity[row_idx][col_idx][hor_idx] << std::endl;
    //             }
    //         }
    //     }
    // }


    vector<vector<vector<double>>> intensity_weight(i_rows_,vector<vector<double>>(i_cols_,vector<double>(i_hor_,0)));
    for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
        for (int col_idx = 0; col_idx < i_cols_; col_idx++) {
            double mean_intensity = 0;
            vector<double> median_intensity;
            for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
                mean_intensity += _intensity[row_idx][col_idx][hor_idx];
                median_intensity.push_back(_intensity[row_idx][col_idx][hor_idx]);
            }
            mean_intensity /= i_hor_;
            double median = median_intensity[20];
            for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
                if (mean_intensity ==0 || _intensity[row_idx][col_idx][hor_idx] > 2*mean_intensity || _intensity[row_idx][col_idx][hor_idx] == 0){
                    intensity_weight[row_idx][col_idx][hor_idx] = 1;
                }
                else
                {
                    intensity_weight[row_idx][col_idx][hor_idx] = _intensity[row_idx][col_idx][hor_idx]/(2*mean_intensity);
                    //std::cout << "intensity_weight is : " << intensity_weight[row_idx][col_idx][hor_idx] << std::endl;
                    //std::cout << "_intensity is : " << _intensity[row_idx][col_idx][hor_idx] << std::endl;
                }
            }
        }
    }


    // vector<vector<vector<double>>> intensity_weight(i_rows_,vector<vector<double>>(i_cols_,vector<double>(i_hor_,0)));
    // for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
    //     for (int col_idx = 0; col_idx < i_cols_; col_idx++) {
    //         double mean_intensity = 0;
    //         double norm_intensity = 0;
    //         for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
    //             mean_intensity += _intensity[row_idx][col_idx][hor_idx];
    //             norm_intensity += pow(_intensity[row_idx][col_idx][hor_idx],2);
    //         }
    //         mean_intensity /= i_hor_;
    //         norm_intensity = sqrt(norm_intensity);
    //         for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
    //             if (norm_intensity == 0 || _intensity[row_idx][col_idx][hor_idx] < norm_intensity / 20){
    //                 intensity_weight[row_idx][col_idx][hor_idx] = 1;
    //             }
    //             else
    //             {
    //                 intensity_weight[row_idx][col_idx][hor_idx] = _intensity[row_idx][col_idx][hor_idx]/(norm_intensity);
    //                 //std::cout << "intensity_weight is : " << intensity_weight[row_idx][col_idx][hor_idx] << std::endl;
    //                 //std::cout << "_intensity is : " << _intensity[row_idx][col_idx][hor_idx] << std::endl;
    //             }
    //         }
    //     }
    // }


    //for bins height
    for ( int col_idx = 0; col_idx < i_cols_; col_idx++ ) {
        for (int hor_idx = 0; hor_idx < i_hor_; ++hor_idx) {
            float element =0;
            float intensity_element =0;
            for (int row_idx = 0; row_idx < i_rows_; ++row_idx) {
                
                if (_state[row_idx][col_idx][hor_idx]>0)
                {
                    //each element equals to products point density and elevation weight
                    //element+= weight[row_idx][col_idx][hor_idx] * pow(2,row_idx);
                    //intensity_element += _state[row_idx][col_idx][hor_idx] * _intensity[row_idx][col_idx][hor_idx] * _height[row_idx][col_idx][hor_idx];
                    //intensity_element += weight[row_idx][col_idx][hor_idx] * pow(2,row_idx) * intensity_weight[row_idx][col_idx][hor_idx]; //weight[row_idx][col_idx][hor_idx] *pow(2,row_idx) * 
                    //intensity_element += weight[row_idx][col_idx][hor_idx] * intensity_weight[row_idx][col_idx][hor_idx] * _height[row_idx][col_idx][hor_idx];
                    intensity_element += (weight[row_idx][col_idx][hor_idx] * 0.8 + intensity_weight[row_idx][col_idx][hor_idx] * 0.2) * _height[row_idx][col_idx][hor_idx];// * intensity_weight[row_idx][col_idx][hor_idx];
                    
                    //intensity_element += weight[row_idx][col_idx][hor_idx] * _height[row_idx][col_idx][hor_idx];
                    //intensity_element += intensity_weight[row_idx][col_idx][hor_idx];
                    //std::cout << "intensity_element is : " << intensity_element << std::endl;
                    //std::cout << "height is : " << pow(2,row_idx) << std::endl;
                    //std::cout << "density weight: " <<  weight[row_idx][col_idx][hor_idx] << std::endl;
                    //std::cout << "intensity weight: " <<  intensity_weight[row_idx][col_idx][hor_idx] << std::endl;
                    //std::cout << "intensity: " <<  intensity_weight[row_idx][col_idx][hor_idx] << std::endl;
                    //std::cout << "density: " <<  _state[row_idx][col_idx][hor_idx] << std::endl;
                    //std::cout << "max_density: " <<  _height[row_idx][col_idx][hor_idx] << std::endl;
                }
            }
            m(col_idx,hor_idx) = intensity_element;
            //std::cout << "intensity_element is : " << intensity_element << std::endl;
        }
    }
    return m;
}



double SCManager::distDirectSC ( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    int num_eff_cols = 0; // i.e., to exclude all-nonzero sector
    double sum_sector_similarity = 0;
    for ( int col_idx = 0; col_idx < _sc1.cols(); col_idx++ )
    {
        VectorXd col_sc1 = _sc1.col(col_idx);
        VectorXd col_sc2 = _sc2.col(col_idx);
        
        if( (col_sc1.norm() == 0) | (col_sc2.norm() == 0) )
            continue; // don't count this sector pair. 

        double sector_similarity = col_sc1.dot(col_sc2) / (col_sc1.norm() * col_sc2.norm());

        sum_sector_similarity = sum_sector_similarity + sector_similarity;
        num_eff_cols = num_eff_cols + 1;
    }
    double sc_sim = sum_sector_similarity / num_eff_cols;
    return 1.0 - sc_sim;

} // distDirectSC


//一直平移找到最小的平移
int SCManager::fastAlignUsingVkey( MatrixXd & _vkey1, MatrixXd & _vkey2)
{
    int argmin_vkey_shift = 0;
    double min_veky_diff_norm = 1000000000;
    for ( int shift_idx = 0; shift_idx < _vkey1.cols(); shift_idx++ )
    {
        MatrixXd vkey2_shifted = circshift(_vkey2, shift_idx);

        MatrixXd vkey_diff = _vkey1 - vkey2_shifted;

        double cur_diff_norm = vkey_diff.norm();
        if( cur_diff_norm < min_veky_diff_norm )
        {
            argmin_vkey_shift = shift_idx;
            min_veky_diff_norm = cur_diff_norm;
        }
    }
    return argmin_vkey_shift;
} // fastAlignUsingVkey




std::pair<double, int> SCManager::distanceBtnScanContext( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    // 1. fast align using variant key (not in original IROS18)
    MatrixXd vkey_sc1 = makeSectorkeyFromScancontext( _sc1 );
    MatrixXd vkey_sc2 = makeSectorkeyFromScancontext( _sc2 );
    int argmin_vkey_shift = fastAlignUsingVkey( vkey_sc1, vkey_sc2 );

    const int SEARCH_RADIUS = round( 0.5 * SEARCH_RATIO * _sc1.cols() ); // a half of search range 
    std::vector<int> shift_idx_search_space { argmin_vkey_shift };
    for ( int ii = 1; ii < SEARCH_RADIUS + 1; ii++ )
    {
        shift_idx_search_space.push_back( (argmin_vkey_shift + ii + _sc1.cols()) % _sc1.cols() );
        shift_idx_search_space.push_back( (argmin_vkey_shift - ii + _sc1.cols()) % _sc1.cols() );
    }
    std::sort(shift_idx_search_space.begin(), shift_idx_search_space.end());

    // 2. fast columnwise diff 
    int argmin_shift = 0;
    double min_sc_dist = 10000000;
    for ( int num_shift: shift_idx_search_space )
    {
        MatrixXd sc2_shifted = circshift(_sc2, num_shift);
        double cur_sc_dist = distDirectSC( _sc1, sc2_shifted );
        if( cur_sc_dist < min_sc_dist )
        {
            argmin_shift = num_shift;
            min_sc_dist = cur_sc_dist;
        }
    }

    return make_pair(min_sc_dist, argmin_shift);

} // distanceBtnScanContext


MatrixXd SCManager::makeScancontext( pcl::PointCloud<SCPointType> & _scan_down )
{
    TicToc t_making_desc;

    int num_pts_scan_down = _scan_down.points.size();

    // main
    const int NO_POINT = -1000;
    MatrixXd desc = NO_POINT * MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);

    SCPointType pt;
    float azim_angle, azim_range; // wihtin 2d plane
    int ring_idx, sctor_idx;
    for (int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)
    {
        pt.x = _scan_down.points[pt_idx].x; 
        pt.y = _scan_down.points[pt_idx].y;
        pt.z = _scan_down.points[pt_idx].z + LIDAR_HEIGHT; // naive adding is ok (all points should be > 0).

        // xyz to ring, sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);
        azim_angle = xy2theta(pt.x, pt.y);

        // if range is out of roi, pass
        if( azim_range > PC_MAX_RADIUS )
            continue;

        ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );
        sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );

        // taking maximum z 
        if ( desc(ring_idx-1, sctor_idx-1) < pt.z ) // -1 means cpp starts from 0
            desc(ring_idx-1, sctor_idx-1) = pt.z; // update for taking maximum value at that bin
    }

    // reset no points to zero (for cosine dist later)
    for ( int row_idx = 0; row_idx < desc.rows(); row_idx++ )
        for ( int col_idx = 0; col_idx < desc.cols(); col_idx++ )
            if( desc(row_idx, col_idx) == NO_POINT )
                desc(row_idx, col_idx) = 0;

    t_making_desc.toc("PolarContext making");

    return desc;
} // SCManager::makeScancontext


// MatrixXd SCManager::makeRingkeyFromScancontext( Eigen::MatrixXd &_desc )
// {
//     /* 
//      * summary: rowwise mean vector
//     */
//     Eigen::MatrixXd invariant_key(_desc.rows(), 1);
//     for ( int row_idx = 0; row_idx < _desc.rows(); row_idx++ )
//     {
//         Eigen::MatrixXd curr_row = _desc.row(row_idx);
//         invariant_key(row_idx, 0) = curr_row.mean();
//     }

//     return invariant_key;
// } // SCManager::makeRingkeyFromScancontext



MatrixXd SCManager::makeRingkeyFromScancontext( Eigen::MatrixXd &_desc )
{
    /* 
     * summary: rowwise mean vector
    */
    Eigen::MatrixXd invariant_key(_desc.rows(), 1);
    Eigen::VectorXd entropy_vect = Eigen::VectorXd::Zero(_desc.rows());
    Eigen::VectorXd variation_vect = Eigen::VectorXd::Zero(_desc.rows());
    Eigen::VectorXd mean_vect = Eigen::VectorXd::Zero(_desc.rows());
    for ( int row_idx = 0; row_idx < _desc.rows(); row_idx++ )
    {
        // for(int ){
        // }
        Eigen::MatrixXd curr_row = _desc.row(row_idx);
        //double max = curr_row.maxCoeff(); // 最大值
        double var = (curr_row.array() - curr_row.mean()).square().sum() / curr_row.size(); // 方差  有可能方差为0

        //double skewness = 0;
        //double kurtosis = 0;

        // if(var == 0 )
        // {
        //     //cout <<"var is zero" << endl;
        //     skewness =0;
        //     kurtosis = 0;
        // }else{
        //     skewness = ((curr_row.array() - curr_row.mean()) / std::sqrt(var)).cube().sum() / curr_row.size();
        //     kurtosis = ((curr_row.array() - curr_row.mean()) / std::sqrt(var)).pow(4).sum() / curr_row.size();
        // }

        double std_vect = std::sqrt(var);
        double entropy = 0;
        double variation = 0;
        if(curr_row.sum() > 5){
            // Compute the probability distribution
            Eigen::MatrixXd prob = curr_row / curr_row.sum();
            // Check if the probability matrix has any zero
            // If it does, we will replace zero with a very small positive number to prevent taking log of zero
            prob = prob.unaryExpr([](double v) { return v == 0 ? std::numeric_limits<double>::min() : v; });
            // Compute the entropy
            entropy = -(prob.array() * prob.array().log()).sum() / std::log(PC_NUM_SECTOR);
            variation = std_vect / curr_row.mean(); //curr_row.mean();
            entropy_vect(row_idx) = entropy;
            variation_vect(row_idx) = variation;
            mean_vect(row_idx) = curr_row.mean();
        }
        //cout << "sum: "<<curr_row.sum() << endl;
        // Eigen::MatrixXd prob = curr_row;
        // prob = prob.array() - prob.minCoeff(); // Subtract the minimum value.
        // prob = prob / prob.sum(); // Normalize to make it a probability distribution.
        // double entropy = -(prob.array() * (prob.array().log().unaryExpr([](double v) { return v > 0 ? std::log(v) : 0; }))).sum(); // Calculate entropy
        // cout << "skewness" <<  skewness <<endl;
        //cout << "kurtosis" <<  kurtosis <<endl;
        //cout << "std" <<  std <<endl;
        // cout << "mean" <<  curr_row.mean() <<endl;
        //cout << "entropy" <<  entropy <<endl;
        //cout << "variation" <<  variation <<endl;
        //invariant_key(row_idx, 0) = curr_row.mean() * skewness; //(entropy / std::log(PC_NUM_SECTOR)) + curr_row.mean() * skewness;//skewness; //curr_row.mean()  * entropy;//max;// + skewness + curr_row.mean(); //std + curr_row.mean();// + var; //max ;//curr_row.mean() ; //skewness; //curr_row.mean() + skewness + kurtosis;// entropy; // + skewness; //+ entropy  + kurtosis;
        //invariant_key(row_idx, 0) = curr_row.mean() * entropy;//skewness;//curr_row.mean();// * skewness;
    }


    for ( int row_idx = 0; row_idx < _desc.rows(); row_idx++ ){
        Eigen::MatrixXd curr_row = _desc.row(row_idx);
        //Eigen::MatrixXd prob_row = curr_row / curr_row.sum();
        double weight_entropy = (1 - entropy_vect(row_idx)) / (Eigen::VectorXd::Ones(entropy_vect.size()) - entropy_vect).sum();
        double weight_variation = variation_vect(row_idx) / (variation_vect.sum());
        invariant_key(row_idx, 0) = curr_row.sum();// *  (weight_variation * 0.9 + weight_entropy * 0.1); //* weight_variation; //(weight_variation * 0.9 + weight_entropy * 0.1);
    }

    return invariant_key;
} // SCManager::makeRingkeyFromScancontext



MatrixXd SCManager::makeSectorkeyFromScancontext( Eigen::MatrixXd &_desc )
{
    /* 
     * summary: columnwise mean vector
    */
    Eigen::MatrixXd variant_key(1, _desc.cols());
    for ( int col_idx = 0; col_idx < _desc.cols(); col_idx++ )
    {
        Eigen::MatrixXd curr_col = _desc.col(col_idx);
        variant_key(0, col_idx) = curr_col.mean();
    }

    return variant_key;
} // SCManager::makeSectorkeyFromScancontext


// MatrixXd SCManager::makeRingFingerprint(vector<vector<vector<int>>>& _density_state){
//     Eigen::MatrixXd variant_ring_key(1, PC_NUM_VERTICAL);
//     for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
//         vector<int> cnt;
//         double mean_cnt = 0;
//         int max_cnt = 0;
//         double std_cnt = 0;

//         double sum = 0;
//         double sum_sq_diff = 0.0;
//         for(int ring_idx = 0; ring_idx < PC_NUM_RING; ring_idx++){
//             int count_bin = 0;
//             for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
//                 if(_density_state[vert_idx][ring_idx][sect_idx] > 0){
//                     count_bin++;
//                     //count_bin += _density_state[vert_idx][ring_idx][sect_idx];
//                     //count_bin *= sect_idx;  否定这个方案
//                 }
//             }
//             cnt.push_back(count_bin);
//             max_cnt = std::max(max_cnt,count_bin);
//             sum += count_bin;
//         }
//         mean_cnt = sum / PC_NUM_RING;
//         for (int i = 0; i < PC_NUM_RING; ++i) {
//             sum_sq_diff += (cnt[i] - mean_cnt) * (cnt[i] - mean_cnt);
//         }
//         std_cnt = std::sqrt(sum_sq_diff / PC_NUM_RING);
//         variant_ring_key(0, vert_idx) =  mean_cnt + max_cnt + std_cnt;
//     }
//     return variant_ring_key;
// }



MatrixXd SCManager::makeRingFingerprint(vector<vector<vector<double>>>& _intensity_state){
    Eigen::MatrixXd variant_ring_key(1, PC_NUM_RING);

    // for(int ring_idx = 0; ring_idx < PC_NUM_RING; ring_idx++){
    //     vector<int> intensity;
    //     double mean_intensity = 0;
    //     int max_intensity = 0;
    //     double std_intensity = 0;

    //     double sum = 0;
    //     double sum_sq_diff = 0.0;
    //     for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
    //         int all_intensity = 0;
    //         for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
    //             if(_intensity_state[vert_idx][ring_idx][sect_idx] > 0){
    //                 //all_intensity += 1 * pow(2,vert_idx);
    //                 all_intensity++;
    //             }
    //         }
    //         intensity.push_back(all_intensity);
    //         max_intensity = std::max(max_intensity,all_intensity);
    //         sum += all_intensity;
    //     }

    //     mean_intensity = sum / PC_NUM_SECTOR;
    //     for (int i = 0; i < PC_NUM_SECTOR; ++i) {
    //         sum_sq_diff += (intensity[i] - mean_intensity) * (intensity[i] - mean_intensity);
    //     }
    //     std_intensity = std::sqrt(sum_sq_diff / PC_NUM_SECTOR);
    //     variant_ring_key(0, ring_idx) =  mean_intensity + max_intensity + std_intensity;
    // }
    // return variant_ring_key;

    for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
        vector<int> intensity;
        double mean_intensity = 0;
        int max_intensity = 0;
        double std_intensity = 0;

        double sum = 0;
        double sum_sq_diff = 0.0;

        for(int ring_idx = 0; ring_idx < PC_NUM_RING; ring_idx++){
            int all_intensity = 0;
            for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
                if(_intensity_state[vert_idx][ring_idx][sect_idx] > 0){
                    //all_intensity = std::max(_intensity_state[vert_idx][ring_idx][sect_idx],all_intensity);
                    all_intensity++;
                    //all_intensity += 100;
                    //all_intensity += _intensity_state[vert_idx][ring_idx][sect_idx];
                    //all_intensity += 1 * pow(2,vert_idx);
                }
            }
            intensity.push_back(all_intensity);
            max_intensity = std::max(max_intensity,all_intensity);
            sum += all_intensity;
        }
        mean_intensity = sum / PC_NUM_RING;
        for (int i = 0; i < PC_NUM_RING; ++i) {
            sum_sq_diff += (intensity[i] - mean_intensity) * (intensity[i] - mean_intensity);
        }
        std_intensity = std::sqrt(sum_sq_diff / PC_NUM_RING);
        variant_ring_key(0, vert_idx) =  mean_intensity + max_intensity + std_intensity;
    }
    return variant_ring_key;
}




//make fingerprint  ring-key
// MatrixXd SCManager::makeRingFingerprint(vector<vector<vector<double>>>& _intensity_state){
    // Eigen::MatrixXd variant_ring_key(1, PC_NUM_RING);
    // for(int ring_idx = 0; ring_idx < PC_NUM_RING; ring_idx++){

    // //     //使用sector
    // //     // double mean_vert_intensity = 0;
    // //     // for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
    // //     //     double mean_vert_intensity = 0;
    // //     //     for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
    // //     //         mean_vert_intensity += _intensity_state[vert_idx][ring_idx][sect_idx];
    // //     //     }
    // //     //     mean_vert_intensity /= PC_NUM_VERTICAL;
    // //     //     double variance_vert_intensity = 0;
    // //     //     double max_vert_intensity = 0;
    // //     //     double skewness_vert_intensity = 0;
    // //     //     double kurtosis_vert_intensity = 0;
    // //     //     double entropy_vert_intensity = 0;

    // //     //     for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
    // //     //         double diff = _intensity_state[vert_idx][ring_idx][sect_idx] - mean_vert_intensity;
    // //     //         variance_vert_intensity += diff * diff;
    // //     //         max_vert_intensity = std::max(max_vert_intensity, _intensity_state[vert_idx][ring_idx][sect_idx]);
    // //     //         skewness_vert_intensity += std::pow(diff, 3.0);
    // //     //         kurtosis_vert_intensity += std::pow(diff, 4.0);
    // //     //     }
    // //     //     variance_vert_intensity /= PC_NUM_VERTICAL;
    // //     //     //求标准差
    // //     //     double stdev_vert_intensity = std::sqrt(variance_vert_intensity);
    // //     //     if(stdev_vert_intensity)
    // //     //     {
    // //     //         skewness_vert_intensity = (skewness_vert_intensity/PC_NUM_VERTICAL) / std::pow(stdev_vert_intensity, 3.0);
    // //     //         kurtosis_vert_intensity = ((kurtosis_vert_intensity/PC_NUM_VERTICAL) / std::pow(stdev_vert_intensity, 4.0)) - 3.0;
    // //     //         entropy_vert_intensity = 0.5 * std::log(2.0 * M_PI * std::exp(1) * variance_vert_intensity);
    // //     //     }else{
    // //     //         skewness_vert_intensity = 0;
    // //     //         kurtosis_vert_intensity = 0;
    // //     //         entropy_vert_intensity = 0;
    // //     //     }
    // //     //     cout << "variance_vert_intensity: " << variance_vert_intensity << endl;
    // //     //     cout << "skewness_vert_intensity: " << skewness_vert_intensity << endl;
    // //     //     cout << "kurtosis_vert_intensity: " << kurtosis_vert_intensity << endl;
    // //     //     cout << "entropy_vert_intensity: " << entropy_vert_intensity << endl;

    // //     //     //求信息熵的分布
    // //     //     //double entropy_vert_intensity = 0.5 * std::log(2.0 * M_PI * std::exp(1) * variance_vert_intensity);
            
    // //     //     //不太确定是否是这样
    // //     //     //double vert_intensity_state = entropy_vert_intensity * max_vert_intensity;
    // //     //     double vert_intensity_state = max_vert_intensity + kurtosis_vert_intensity;
    // //     //     mean_vert_intensity += vert_intensity_state;
    // //     // }
    // //     // variant_ring_key(0, ring_idx) = mean_vert_intensity / PC_NUM_SECTOR;  //求每一个列的intensity的最大值以及分布情况去对应  感觉会比单一的z最大值更好


    // //     //使用vertial
    //     double mean_vert_intensity = 0;
    //     for(int vert_idx = 0; vert_idx < PC_NUM_VERTICAL; vert_idx++){
    //         double mean_sect_intensity = 0;
    //         for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
    //             mean_sect_intensity += _intensity_state[vert_idx][ring_idx][sect_idx];
    //         }
    //         mean_sect_intensity /= PC_NUM_SECTOR;
    //         double variance_sect_intensity = 0;
    //         double max_sect_intensity = 0;
    //         double skewness_sect_intensity = 0;
    //         double kurtosis_sect_intensity = 0;
    //         double entropy_sect_intensity = 0;

    //         for(int sect_idx = 0; sect_idx < PC_NUM_SECTOR; sect_idx++){
    //             double diff = _intensity_state[vert_idx][ring_idx][sect_idx] - mean_sect_intensity;
    //             variance_sect_intensity += diff * diff;
    //             max_sect_intensity = std::max(max_sect_intensity, _intensity_state[vert_idx][ring_idx][sect_idx]);
    //             skewness_sect_intensity += std::pow(diff, 3.0);
    //             kurtosis_sect_intensity += std::pow(diff, 4.0);
    //         }

    //         variance_sect_intensity /= PC_NUM_SECTOR;
    //         double stdev_vert_intensity = std::sqrt(variance_sect_intensity);
    //         if(stdev_vert_intensity)
    //         {
    //             skewness_sect_intensity = (skewness_sect_intensity/PC_NUM_SECTOR) / std::pow(stdev_vert_intensity, 3.0);
    //             kurtosis_sect_intensity = ((kurtosis_sect_intensity/PC_NUM_SECTOR) / std::pow(stdev_vert_intensity, 4.0)) - 3.0;
    //             entropy_sect_intensity = 0.5 * std::log(2.0 * M_PI * std::exp(1) * variance_sect_intensity);
    //         }else{
    //             skewness_sect_intensity = 0;
    //             kurtosis_sect_intensity = 0;
    //             entropy_sect_intensity = 0;
    //         }
    //         cout << "variance_sect_intensity: " << variance_sect_intensity << endl;
    //         cout << "skewness_sect_intensity: " << skewness_sect_intensity << endl;
    //         cout << "kurtosis_sect_intensity: " << kurtosis_sect_intensity << endl;
    //         cout << "entropy_sect_intensity: " << entropy_sect_intensity << endl;
    //         cout << "max_sect_intensity: " << max_sect_intensity << endl;
            
    //         //不太确定是否是这样
    //         double sect_intensity_state = max_sect_intensity + skewness_sect_intensity + variance_sect_intensity + kurtosis_sect_intensity;
    //         mean_vert_intensity += sect_intensity_state;
    //     }
    //     variant_ring_key(0, ring_idx) = mean_vert_intensity / PC_NUM_VERTICAL;
    // }
    // return variant_ring_key;
// }



const Eigen::MatrixXd& SCManager::getConstRefRecentSCD(void)
{
    return polarcontexts_.back();
}

double SCManager::simDirect ( MatrixXd &_vec1, MatrixXd &_vec2 )
{
    VectorXd vec1 = _vec1.row(0);
    VectorXd vec2 = _vec2.row(0);
    assert(vec1.size() == vec2.size());  // Ensure same size
    double dot_product = vec1.dot(vec2);
    double norm_product = vec1.norm() * vec2.norm();
    double similarity = dot_product / norm_product;
    return similarity;
} 




void SCManager::makeAndSaveSpatialcontextAndKeys( pcl::PointCloud<SCPointType> & _scan_down )
{
    auto Spatialcontext = makeSpatialState(_scan_down);
    vector<vector<vector<int>>> density_state = std::get<0>(Spatialcontext);
    vector<vector<vector<double>>> intensity_state = std::get<1>(Spatialcontext);
    vector<vector<vector<double>>> height_state = std::get<2>(Spatialcontext);
    //Eigen::MatrixXd ringFingerkey = makeRingFingerprint(intensity_state); // or density state
    //Eigen::MatrixXd ringFingerkey = makeRingFingerprint(density_state); // or density state
    MatrixXd m_spatial = getSpatialDesc(Spatialcontext);
    Eigen::MatrixXd ringFingerkey = makeRingkeyFromScancontext(m_spatial);


    Eigen::MatrixXd sc = makeScancontext(_scan_down); // v1 
    //Eigen::MatrixXd ringkey = makeRingkeyFromScancontext( sc );
    Eigen::MatrixXd sectorkey = makeSectorkeyFromScancontext( sc );

    std::vector<float> polarcontext_invkey_vec = eig2stdvec( ringFingerkey );  //本质上是相同的
    polarcontexts_.push_back( sc ); 
    polarcontext_invkeys_.push_back( ringFingerkey );
    polarcontext_vkeys_.push_back( sectorkey );
    polarcontext_invkeys_mat_.push_back( polarcontext_invkey_vec );

    spatial_desc.push_back(m_spatial);
    // cout <<polarcontext_vkeys_.size() << endl;
} // SCManager::makeAndSaveSpatialcontextAndKeys


std::tuple<int, float, float>  SCManager::detectLoopClosureID ( void )
{
    int loop_id { -1 }; // init with -1, -1 means no loop (== LeGO-LOAM's variable "closestHistoryFrameID")

    auto curr_key = polarcontext_invkeys_mat_.back(); // current observation (query)
    auto curr_desc = polarcontexts_.back(); // current observation (query)

    //desc for spatial
    auto curr_spatial_desc = spatial_desc.back();
    /* 
     * step 1: candidates from ringkey tree_  NUM_EXCLUDE_RECENT = 30 代表
     */
    if( (int)polarcontext_invkeys_mat_.size() < NUM_EXCLUDE_RECENT + 1)
    {
        std::tuple<int, float, float> result {loop_id, 0.0, 0.0};
        return result; // Early return 
    }
    // tree_ reconstruction (not mandatory to make everytime)  TREE_MAKING_PERIOD_
    if( tree_making_period_conter % TREE_MAKING_PERIOD_ == 0) // to save computation cost
    {
        TicToc t_tree_construction;
        polarcontext_invkeys_to_search_.clear();
        polarcontext_invkeys_to_search_.assign( polarcontext_invkeys_mat_.begin(), polarcontext_invkeys_mat_.end() - NUM_EXCLUDE_RECENT ) ;
        polarcontext_tree_.reset(); 
        polarcontext_tree_ = std::make_unique<InvKeyTree>(PC_NUM_RING /* dim */, polarcontext_invkeys_to_search_, 10 /* max leaf */ );
        // tree_ptr_->index->buildIndex(); // inernally called in the constructor of InvKeyTree (for detail, refer the nanoflann and KDtreeVectorOfVectorsAdaptor)
        t_tree_construction.toc("Tree construction");
    }
    tree_making_period_conter = tree_making_period_conter + 1;
    double min_dist = 10000000; // init with somthing large
    int nn_align = 0;
    int nn_idx = 0;

    double person_min_dist = 10000000;
    int person_nn_align = 0;
    int person_nn_idx = 0;

    double wasserstein_min_dist = 10000000;
    int wasserstein_nn_align = 0;
    int wasserstein_nn_idx = 0;

    // knn search  NUM_CANDIDATES_FROM_TREE
    std::vector<size_t> candidate_indexes( NUM_CANDIDATES_FROM_TREE ); 
    std::vector<float> out_dists_sqr( NUM_CANDIDATES_FROM_TREE );

    TicToc t_tree_search;
    nanoflann::KNNResultSet<float> knnsearch_result( NUM_CANDIDATES_FROM_TREE );
    knnsearch_result.init( &candidate_indexes[0], &out_dists_sqr[0] );
    polarcontext_tree_->index->findNeighbors( knnsearch_result, &curr_key[0] /* query */, nanoflann::SearchParams(10) ); 
    t_tree_search.toc("Tree search");

    /* 
     *  step 2: pairwise distance (find optimal columnwise best-fit using cosine distance)
     */
    TicToc t_calc_dist;   
    for ( int candidate_iter_idx = 0; candidate_iter_idx < NUM_CANDIDATES_FROM_TREE; candidate_iter_idx++ )
    {
        MatrixXd polarcontext_candidate = polarcontexts_[ candidate_indexes[candidate_iter_idx] ];
        std::pair<double, int> sc_dist_result = distanceBtnScanContext( curr_desc, polarcontext_candidate ); 
        
        double candidate_dist = sc_dist_result.first;
        int candidate_align = sc_dist_result.second;

        //还是用的余弦距离
        auto spatial_feature_candidate = spatial_desc[ candidate_indexes[candidate_iter_idx] ];
        float min_distance = 1;
        float min_person_distance = 1;
        float min_wasserstein_distance = 100000;
        int min_wasserstein = 0;
        int  min_sect = 0;

        double max_simlarity = 0;
        int max_sect = 0;
        MatrixXd vectorM1 = makeSectorkeyFromScancontext(curr_spatial_desc);
        for (int i = 0; i < i_hor_; ++i) {
        MatrixXd transM2 = circshift(spatial_feature_candidate,i);
        Eigen::MatrixXd vectorM2 = makeSectorkeyFromScancontext(transM2);

        double simlarity = simDirect(vectorM1,vectorM2);
        if(simlarity > max_simlarity)
        {
            max_simlarity = simlarity;
            max_sect = i;
        }

        float dis = cosdistance(curr_spatial_desc,transM2);
        // float person_dis = Pearson_Correlation(curr_spatial_desc,transM2);
        //float wasserstein_dis = wasserstein_distance(curr_spatial_desc,transM2);
        float wasserstein_dis = wasserstein_distance(vectorM1,vectorM2);
        //cout << "dis: " << dis <<endl;
        //cout << "person_dis: " << person_dis <<endl;
        //cout << "wasserstein_dis: " << wasserstein_dis <<endl;
        if (dis<min_distance){
            min_distance = dis;
            min_sect = i;
        }
        // if (person_dis<min_person_distance){
        //     min_person_distance = person_dis;
        // }
        if (wasserstein_dis<min_wasserstein_distance){
            min_wasserstein_distance = wasserstein_dis;
            min_wasserstein = i;
        }

        }

        cout << "max_sect " << max_sect << " min_sect " << min_sect << " min_wasserstein " << min_wasserstein <<endl;
        MatrixXd transM3 = circshift(spatial_feature_candidate,min_wasserstein);
        min_distance = cosdistance(curr_spatial_desc,transM3);
        // std::cout << "Spatial desc is : " <<  min_distance << std::endl;
        // std::cout << "SC desc is : " <<  candidate_dist << std::endl;
        // if( candidate_dist < min_dist )
        // {
        //     min_dist = candidate_dist;
        //     nn_align = candidate_align;

        //     nn_idx = candidate_indexes[candidate_iter_idx];
        // }

        if( min_distance < min_dist )
        {
            min_dist = min_distance;
            nn_align = candidate_align;
            nn_idx = candidate_indexes[candidate_iter_idx];
        }

        // if( min_person_distance < person_min_dist )
        // {
        //     person_min_dist = min_person_distance;
        //     person_nn_align = candidate_align;

        //     person_nn_idx = candidate_indexes[candidate_iter_idx];
        // }

        // if( min_wasserstein_distance < wasserstein_min_dist )
        // {
        //     wasserstein_min_dist = min_wasserstein_distance;
        //     wasserstein_nn_align = candidate_align;

        //     wasserstein_nn_idx = candidate_indexes[candidate_iter_idx];
        // }
    }

    t_calc_dist.toc("Distance calc");

    cout << "min_dist: " << min_dist << " nn_align: " << nn_align << " nn_idx: " <<nn_idx << endl;
    //cout << "person_min_dist: " << person_min_dist << " person_nn_align: " << person_nn_align << " person_nn_idx: " <<person_nn_idx << endl;
    //cout << "wasserstein_min_dist: " << wasserstein_min_dist << " wasserstein_nn_align: " << wasserstein_nn_align << " wasserstein_nn_idx: " <<wasserstein_nn_idx << endl;

    /* 
     * loop threshold check  SC_DIST_THRES  SC_DIST_THRES
     */
    if( min_dist < 0.6 )
    {
        loop_id = nn_idx; 
    
        // std::cout.precision(3); 
        //cout << "[Loop found] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        //cout << "[Loop found] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }
    else
    {
        //std::cout.precision(3); 
        //cout << "[Not loop] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        //cout << "[Not loop] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }
    // To do: return also nn_align (i.e., yaw diff)
    float yaw_diff_rad = deg2rad(nn_align * PC_UNIT_SECTORANGLE);
    std::tuple<int, float, float> result(loop_id, yaw_diff_rad, min_dist);
    return result;
} // SCManager::detectLoopClosureID

// } // namespace SC2