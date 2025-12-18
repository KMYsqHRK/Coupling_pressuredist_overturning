#pragma once
#if !defined(PW_SDK_EXAMPLE_POINT_PRESSURE_H_INCLUDED)
#define PW_SDK_EXAMPLE_POINT_PRESSURE_H_INCLUDED

#include "example_module.hpp"
#include <string>
#include <vector>
#include <eigen/Eigen/Sparse>

namespace pw { namespace example { namespace point_pressure {

// 境界点の設定
struct BoundaryPoint {
    std::string name;
    double x, y, z;
};

// 設定データ
struct Settings {
    // for pressure distribution
    std::string output_file;
    int log_interval;
    int DF_index; 
    int division_number;
    std::vector<BoundaryPoint> boundary_points;
    double measurement_radius;
    bool exclude_ghost_particles;
    double calculation_start_time;
    double hydraulic_conductivity;//ダルシー則のための透水係数
    double kinematic_viscosity_coefficient;//動粘性係数
    
    // for moment and overturning calculation
    int object_id;
    double FSI_enabled;
    int window_size;
    double mass;
    double gravity[3];
    double rotation_center[3];
    double center_of_gravity[3];
    double moment_of_inertia;  // Moment of inertia around rotation center
    double resistance_moment;   // Resistance moment by piles
    
    // for force history
    std::string force_history_output_file;
    int force_history_save_interval;

    // for pressure grid export
    std::string pressure_grid_output_file;
    int pressure_grid_save_interval;
    bool pressure_grid_include_wet_dry;

    // for units
    std::string length_in_PW;
};

// モジュール継承
class Module : public pw::example::Module {
private:
    Settings m_settings;
    // Iterative solver settings (no pre-factorization needed)
    // Default: max_iterations=1000, tolerance=1e-6
    // These defaults work well for typical grid sizes (5x5 to 50x50)
    // Convergence typically occurs in 10-100 iterations with IncompleteLUT preconditioner
    int m_max_iterations;
    double m_solver_tolerance;
    int m_cached_n;

    Eigen::MatrixXd m_grid_coords; // Grid coordinates for pressure point
    bool m_grid_coords_initialized;

    std::vector<std::array<double, 3>> m_original_boundary_coords; // Boundary point coordinates (x,y)
    bool m_original_boundary_coords_initialized;

    // Distance matrix from rotation center to pressure points
    Eigen::MatrixXd m_distance_matrix;
    bool m_distance_matrix_initialized;

    // Wet-dry state matrix
    Eigen::MatrixXd m_wet_dry_matrix;
    bool m_wet_dry_matrix_initialized;

    Eigen::MatrixXd m_grid_infil_distance;
    bool m_grid_infil_distance_initialized;

    bool pile_is_broken = false;

    int tid = 0;

	std::vector<double>  Externalforce_x;
    std::vector<double>  Upliftforce;
    std::vector<double>  UpliftTorque;
    std::vector<double>  HorizontalTorque;
    std::vector<double>  TotalTorque;
	std::vector<double>  Angle;
	std::vector<double>  AngularVelocity;
	std::vector<double>  Time;

public:
    // Constructor
    Module() : m_max_iterations(1000), m_solver_tolerance(1e-6), m_cached_n(0), m_grid_coords_initialized(false), m_original_boundary_coords_initialized(false), m_distance_matrix_initialized(false),
               m_wet_dry_matrix_initialized(false), m_grid_infil_distance_initialized(false) {
        Externalforce_x.push_back(0.0);
        Upliftforce.push_back(0.0);
        UpliftTorque.push_back(0.0);
        HorizontalTorque.push_back(0.0);
        TotalTorque.push_back(0.0);
        Angle.push_back(0.0);
        AngularVelocity.push_back(0.0);
        Time.push_back(0.0);
    }
    
    // members
    Eigen::MatrixXd boundary_conditions;
    
    // methods
    virtual void initialize();

    // 設定ファイルの読み込み(only once)
    void load_settings(const std::string& filename = "settings.json");

    // ポアソン方程式の係数マトリクス作成(for iterative solver)
    // Compressed system: builds matrix only for wet points
    Eigen::SparseMatrix<double> create_coefficient_matrix(int n, const Eigen::MatrixXd& wet_dry_matrix,
                                                          const Eigen::MatrixXi& grid_to_compressed);
    // Distance matrix calculation (called once at initialization)
    void initialize_distance_matrix(int n);
    void initialize_wet_dry_matrix(int n);
    void initialize_grid_coordinates(int n, double lattice_x, double lattice_y);
    void initialize_infil_distance_matrix(int n);

    // main loop(every loop)
    void main(const pw::api::Session& session);

    // Step 1 in main loop :境界条件作成
    std::vector<std::array<double, 3>> generate_boundary_points(int n, double dx, double dy, double dz, double dTheta_y);
    Eigen::MatrixXd create_pressure_matrix(const std::vector<double>& boundary_pressures, int n);
    // Compressed RHS: constructs RHS for wet points only with P=0 at wet-dry interface
    Eigen::VectorXd construct_rhs_vector(const Eigen::MatrixXd& boundary_matrix, int n,
                                         const Eigen::MatrixXd& wet_dry_matrix,
                                         const Eigen::MatrixXi& grid_to_compressed);

    // Step 2 in main loop :ポアソン方程式求解
    Eigen::MatrixXd solve_pressure_distribution(const Eigen::MatrixXd& boundary_matrix, int n);

    // step 2 in main loop :圧力分布行列を用いた浸潤距離行列の更新
    void update_infil_distance_matrix(const pw::api::Session & session, const Eigen::MatrixXd& pressure_matrix);

    // step 3 : Update wet-dry matrix
    void update_wet_dry_matrix(const pw::api::Session & session, const std::vector<double>& boundary_pressures);

    // Step 5 in main loop :圧力分布からモーメントを計算
    void calculate_moment(const pw::api::Session & session, const Eigen::MatrixXd& pressure_distribution, const Eigen::MatrixXd& wet_dry_matrix);
    void calculate_new_position(const pw::api::Session & session, const Eigen::MatrixXd& pressure_distribution, const Eigen::MatrixXd& wet_dry_matrix);

    // Step 6 in main loop :運動方程式より与えられる速度からノードの位置を更新
    void set_node_position(const pw::api::Session & session);

    // Pretime process: handles initialization before calculation start time
    void pretime_process(const pw::api::Session & session);

    // Force history output
    void save_force_history() const;

    // Pressure grid export
    void save_pressure_grid(const Eigen::MatrixXd& pressure_distribution,
                           const Eigen::MatrixXd& wet_dry_matrix,
                           double current_time) const;

    // Unit conversion helper: returns conversion factor from PW units to meters
    // For mm: returns 1000.0 (1 m = 1000 mm)
    // For cm: returns 100.0 (1 m = 100 cm)
    // For m: returns 1.0 (1 m = 1 m)
    double get_pw_to_meters_conversion() const {
        if (m_settings.length_in_PW == "mm") {
            return 1000.0;
        } else if (m_settings.length_in_PW == "cm") {
            return 100.0;
        } else if (m_settings.length_in_PW == "m") {
            return 1.0;
        } else {
            std::cerr << "Warning: Unknown length unit '" << m_settings.length_in_PW
                      << "'. Defaulting to mm." << std::endl;
            return 1000.0;
        }
    }
};

// 単一の点で圧力測定関数(structure for single point pressure measurement)
double calculate_pressure_at_point(
    const pw::api::Session& session,
    double x, double y, double z,
    double radius,
    int DF_index
);

// 重み関数
inline double weight(double distance, double radius) {
    if (distance > radius) return 0.0;
    return (radius / distance) - 1;
}

//移動平均を取得する関数
inline double calculateWindowAverage(const std::vector<double>& vec, const size_t windowSize) {
    size_t actualSize = vec.size();
    double sum = 0.0;
    // 指定されたwindowSize（または配列サイズ分）の要素を合計
    for (size_t i = 0; i < std::min(windowSize, actualSize); ++i) {
        sum += vec[actualSize - 1 - i];
    }
    // windowSizeで割って平均を計算（不足分は0として扱われる）
    return sum / windowSize;
}

}}}
#endif