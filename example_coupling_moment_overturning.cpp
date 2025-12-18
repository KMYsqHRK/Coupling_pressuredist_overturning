#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include "nlohmann/single_include/nlohmann/json.hpp"
#include "example_coupling_moment_overturning.hpp"

using json = nlohmann::json;

namespace pw { namespace example { namespace point_pressure {

// 前方宣言
class PressureDistributionSolver;

// DF周辺粒子について特定の点近傍の粒子を探しEigen/SparceのMatrixに格納するアルゴリズム
// パーセンタイル法
double calculate_pressure_at_point(
    const pw::api::Session& session,
    double x, double y, double z,
    double radius,
    int DF_index)
{
    std::vector<pw::api::DistanceField> distance_fields = session.distance_fields();
    auto particles = session.mps_particles();
    auto arrays = particles.arrays();
    const auto positions = arrays.vec3_array("position", PW_ACCESS_read_only_c);
    const auto pressures = arrays.scalar_array("pressure", PW_ACCESS_read_only_c);

    // 距離場が存在しない場合はエラーを返す
    if (distance_fields.empty()) {
        std::cerr << "Error: No distance fields available." << std::endl;
        return 0.0;
    }

    pw::api::DistanceField& distance_field = distance_fields[DF_index];
    pw::api::Arrays distance_field_arrays = distance_field.arrays();

    // 距離場の粒子情報を取得
    auto interaction = distance_field.interaction(particles);

    // 近傍粒子がいない場合は0.0を返す
    if (interaction.size == 0)
    {
        return 0.0;
    }

    // 有効な圧力値を収集
    std::vector<double> valid_pressures;
    valid_pressures.reserve(interaction.size);

    for (size_t i = 0; i < interaction.size; i++) {
        // 距離チェック
        auto pos = positions[interaction.index[i]];
        double dx = pos[0] - x;
        double dy = pos[1] - y;
        double dz = pos[2] - z;
        double dist_sq = dx * dx + dy * dy + dz * dz;

        if (dist_sq <= radius * radius) {
            double pressure = pressures[interaction.index[i]];
            // Exclude free surface particles (pressure < 0.001)
            if (pressure >= 0.001) {
                valid_pressures.push_back(pressure);
            }
        }
    }

    // 有効な粒子がない場合
    if (valid_pressures.empty()) {
        return 0.0;
    }

    size_t n = valid_pressures.size();

    // 粒子数が少ない場合は全て使用（外れ値除外の意味がない）
    if (n < 4) {
        double sum = 0.0;
        for (double p : valid_pressures) {
            sum += p;
        }
        return sum / n;
    }

    // Q1 (25パーセンタイル) を取得
    size_t q1_idx = n / 4;
    std::nth_element(valid_pressures.begin(),
                     valid_pressures.begin() + q1_idx,
                     valid_pressures.end());

    // Q3 (75パーセンタイル) を取得
    size_t q3_idx = 3 * n / 4;
    std::nth_element(valid_pressures.begin() + q1_idx,
                     valid_pressures.begin() + q3_idx,
                     valid_pressures.end());

    // Q1～Q3の範囲（中央50%）の平均を計算
    double sum = 0.0;
    for (size_t i = q1_idx; i < q3_idx; i++) {
        sum += valid_pressures[i];
    }

    return sum / (q3_idx - q1_idx);
}

// 設定読み込み
void Module::load_settings(const std::string& filename)
{
    std::cout << "Loading settings from: " << filename << std::endl;
    
    // Initialize default settings
    m_settings.output_file = "pressure_distribution_results.csv";
    m_settings.log_interval = 10;
    m_settings.DF_index = 0;
    m_settings.division_number = 5;
    m_settings.measurement_radius = 0.01;
    m_settings.exclude_ghost_particles = true;
    m_settings.calculation_start_time = 0.0;
    
    // Moment and overturning calculation defaults
    m_settings.object_id = 0;
    m_settings.FSI_enabled = 1.0;
    m_settings.window_size = 1;
    m_settings.mass = 1.0;

    // Gravity defaults (standard Earth gravity in z-direction)
    m_settings.gravity[0] = 0.0;
    m_settings.gravity[1] = 0.0;
    m_settings.gravity[2] = 9.8;

    // Rotation center defaults
    m_settings.rotation_center[0] = 0.0; // x
    m_settings.rotation_center[1] = 0.0; // y
    m_settings.rotation_center[2] = 0.0; // z

    // Center of gravity defaults
    m_settings.center_of_gravity[0] = 0.0; // x
    m_settings.center_of_gravity[1] = 0.0; // y
    m_settings.center_of_gravity[2] = 0.0; // z

    // Moment of inertia and resistance moment defaults
    m_settings.moment_of_inertia = 1.0;    // kg·m^2
    m_settings.resistance_moment = 0.0;    // N·m
    
    // Force history defaults
    m_settings.force_history_output_file = "force_history.csv";
    m_settings.force_history_save_interval = 10;

    // Pressure grid export defaults
    m_settings.pressure_grid_output_file = "pressure_grid.csv";
    m_settings.pressure_grid_save_interval = 100;
    m_settings.pressure_grid_include_wet_dry = true;

    m_settings.length_in_PW = "mm";
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open settings file " << filename 
                  << ". Using default settings." << std::endl;
        return;
    }
    
    try {
        json config;
        file >> config;
        
        // 基本設定
        if (config.contains("output_file")) {
            m_settings.output_file = config["output_file"];
        }
        
        if (config.contains("log_interval")) {
            m_settings.log_interval = config["log_interval"];
        }
        
        // フィルタリング設定
        if (config.contains("filtering") && config["filtering"].contains("exclude_ghost_particles")) {
            m_settings.exclude_ghost_particles = config["filtering"]["exclude_ghost_particles"];
        }
        
        // 分割数設定
        if (config.contains("division_number")) {
            m_settings.division_number = config["division_number"];
        }
        
        // 測定半径設定
        if (config.contains("measurement_radius")) {
            m_settings.measurement_radius = config["measurement_radius"];
        }
        
        // 境界点設定
        if (config.contains("boundary_points")) {
            m_settings.boundary_points.clear();
            for (const auto& point_config : config["boundary_points"]) {
                BoundaryPoint point;
                point.name = point_config.value("name", "Point");
                
                if (point_config.contains("position")) {
                    point.x = point_config["position"].value("x", 0.0);
                    point.y = point_config["position"].value("y", 0.0);
                    point.z = point_config["position"].value("z", 0.0);
                }
                
                m_settings.boundary_points.push_back(point);
            }
        }

        if (config.contains("df_index")) {
            int df_index = config["df_index"];
            if (df_index < 0) {
                throw std::out_of_range("DF_index must be non-negative");
            }
            m_settings.DF_index = df_index;
        }

        // シミュレーション設定
        if (config.contains("simulation")) {
            const auto& sim_config = config["simulation"];
            
            if (sim_config.contains("object_id")) {
                m_settings.object_id = sim_config["object_id"];
            }
            
            if (sim_config.contains("fsi_enabled")) {
                m_settings.FSI_enabled = sim_config["fsi_enabled"];
            }
            
            if (sim_config.contains("windowSize")) {
                m_settings.window_size = sim_config["windowSize"];
            }

            if (sim_config.contains("calculation_start_time")){
                m_settings.calculation_start_time = sim_config["calculation_start_time"];
            }
        }
        
        // 物理パラメータ設定
        if (config.contains("physics")) {
            const auto& physics_config = config["physics"];

            if (physics_config.contains("mass")) {
                m_settings.mass = physics_config["mass"];
            }

            if (physics_config.contains("moment_of_inertia")) {
                m_settings.moment_of_inertia = physics_config["moment_of_inertia"];
            }

            if (physics_config.contains("resistance_moment")) {
                m_settings.resistance_moment = physics_config["resistance_moment"];
            }

            // Center of gravity
            if (physics_config.contains("center_of_gravity")) {
                const auto& cog_config = physics_config["center_of_gravity"];

                if (cog_config.contains("x")) {
                    m_settings.center_of_gravity[0] = cog_config["x"];
                }

                if (cog_config.contains("y")) {
                    m_settings.center_of_gravity[1] = cog_config["y"];
                }

                if (cog_config.contains("z")) {
                    m_settings.center_of_gravity[2] = cog_config["z"];
                }
            }
        }
        
        // 重力設定
        if (config.contains("gravity")) {
            const auto& gravity_config = config["gravity"];
            
            if (gravity_config.contains("x")) {
                m_settings.gravity[0] = gravity_config["x"];
            }
            
            if (gravity_config.contains("y")) {
                m_settings.gravity[1] = gravity_config["y"];
            }
            
            if (gravity_config.contains("z")) {
                m_settings.gravity[2] = gravity_config["z"];
            }
        }

        // Rotation center settings
        if (config.contains("rotation_center")) {
            const auto& rotation_center_config = config["rotation_center"];
            
            if (rotation_center_config.contains("x")) {
                m_settings.rotation_center[0] = rotation_center_config["x"];
            }
            
            if (rotation_center_config.contains("y")) {
                m_settings.rotation_center[1] = rotation_center_config["y"];
            }
            
            if (rotation_center_config.contains("z")) {
                m_settings.rotation_center[2] = rotation_center_config["z"];
            }
        }

        //infiltration settings
        if (config.contains("infiltration")){
            const auto& infiltration_config = config["infiltration"];
            if (infiltration_config.contains("hydraulic_conductivity")){
                m_settings.hydraulic_conductivity = infiltration_config["hydraulic_conductivity"];
            }
            if (infiltration_config.contains("kinematic_viscosity_coefficient")){
                m_settings.kinematic_viscosity_coefficient = infiltration_config["kinematic_viscosity_coefficient"];
            }
        }
        
        // Force history settings
        if (config.contains("force_history")) {
            const auto& force_history_config = config["force_history"];

            if (force_history_config.contains("output_file")) {
                m_settings.force_history_output_file = force_history_config["output_file"];
            }

            if (force_history_config.contains("save_interval")) {
                m_settings.force_history_save_interval = force_history_config["save_interval"];
            }
        }

        // Pressure grid export settings
        if (config.contains("pressure_grid_export")) {
            const auto& pressure_grid_config = config["pressure_grid_export"];

            if (pressure_grid_config.contains("output_file")) {
                m_settings.pressure_grid_output_file = pressure_grid_config["output_file"];
            }

            if (pressure_grid_config.contains("save_interval")) {
                m_settings.pressure_grid_save_interval = pressure_grid_config["save_interval"];
            }

            if (pressure_grid_config.contains("include_wet_dry")) {
                m_settings.pressure_grid_include_wet_dry = pressure_grid_config["include_wet_dry"];
            }
        }

        if (config.contains("units")){
            const auto& units_config = config["units"];
            if (units_config.contains("length_in_PW")){
                m_settings.length_in_PW = units_config["length_in_PW"];
            }
        }

        // Iterative solver settings with validation
        if (config.contains("solver")) {
            const auto& solver_config = config["solver"];

            if (solver_config.contains("max_iterations")) {
                int config_max_iter = solver_config["max_iterations"];
                if (config_max_iter > 0 && config_max_iter <= 100000) {
                    m_max_iterations = config_max_iter;
                    std::cout << "Using solver max_iterations from config: " << m_max_iterations << std::endl;
                } else {
                    std::cerr << "Warning: Invalid max_iterations (" << config_max_iter
                              << "). Must be in range [1, 100000]. Using default: " << m_max_iterations << std::endl;
                }
            } else {
                std::cout << "Using default solver max_iterations: " << m_max_iterations << std::endl;
            }

            if (solver_config.contains("tolerance")) {
                double config_tol = solver_config["tolerance"];
                if (config_tol > 0.0 && config_tol < 1.0) {
                    m_solver_tolerance = config_tol;
                    std::cout << "Using solver tolerance from config: " << m_solver_tolerance << std::endl;
                } else {
                    std::cerr << "Warning: Invalid tolerance (" << config_tol
                              << "). Must be in range (0, 1). Using default: " << m_solver_tolerance << std::endl;
                }
            } else {
                std::cout << "Using default solver tolerance: " << m_solver_tolerance << std::endl;
            }
        } else {
            // No solver config found - using defaults
            std::cout << "No 'solver' section found in config - using default iterative solver settings:" << std::endl;
            std::cout << "  Default max_iterations: " << m_max_iterations << std::endl;
            std::cout << "  Default tolerance: " << m_solver_tolerance << std::endl;
        }

        std::cout << "Settings loaded successfully:" << std::endl;
        std::cout << "  Output file: " << m_settings.output_file << std::endl;
        std::cout << "  Log interval: " << m_settings.log_interval << std::endl;
        std::cout << "  DF index: " << m_settings.DF_index << std::endl;
        std::cout << "  Division number: " << m_settings.division_number << std::endl;
        std::cout << "  Boundary points: " << m_settings.boundary_points.size() << std::endl;
        std::cout << "  Measurement radius: " << m_settings.measurement_radius << std::endl;
        std::cout << "  Object ID: " << m_settings.object_id << std::endl;
        std::cout << "  FSI enabled: " << m_settings.FSI_enabled << std::endl;
        std::cout << "  Window size: " << m_settings.window_size << std::endl;
        std::cout << "  Mass: " << m_settings.mass << " kg" << std::endl;
        std::cout << "  Moment of inertia: " << m_settings.moment_of_inertia << " kg·m^2" << std::endl;
        std::cout << "  Resistance moment: " << m_settings.resistance_moment << " N·m" << std::endl;
        std::cout << "  Gravity: (" << m_settings.gravity[0] << ", " << m_settings.gravity[1] << ", " << m_settings.gravity[2] << ") m/s^2" << std::endl;
        std::cout << "  Rotation center: (" << m_settings.rotation_center[0] << ", " << m_settings.rotation_center[1] << ", " << m_settings.rotation_center[2] << ")" << std::endl;
        std::cout << "  Center of gravity: (" << m_settings.center_of_gravity[0] << ", " << m_settings.center_of_gravity[1] << ", " << m_settings.center_of_gravity[2] << ")" << std::endl;
        std::cout << "  Force history output file: " << m_settings.force_history_output_file << std::endl;
        std::cout << "  Force history save interval: " << m_settings.force_history_save_interval << std::endl;
        std::cout << "  Pressure grid output file: " << m_settings.pressure_grid_output_file << std::endl;
        std::cout << "  Pressure grid save interval: " << m_settings.pressure_grid_save_interval << std::endl;
        std::cout << "  Pressure grid include wet/dry: " << (m_settings.pressure_grid_include_wet_dry ? "true" : "false") << std::endl;
        std::cout << "  Iterative solver max iterations: " << m_max_iterations << std::endl;
        std::cout << "  Iterative solver tolerance: " << m_solver_tolerance << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing settings file: " << e.what() << std::endl;
        std::cerr << "Using default settings." << std::endl;
    }
}

// 境界上の4n個の点を生成する関数
// 回転移動の考慮はしていない
std::vector<std::pair<double, double>> Module::generate_boundary_points(int n, double deformation_x, double deformation_y)
{
    std::vector<std::pair<double, double>> boundary_points;
    
    if (m_settings.boundary_points.size() < 3) {
        std::cerr << "Error: Need at least 3 boundary points to define parallelogram" << std::endl;
        return boundary_points;
    }
    
    // 3つの境界点から矩形（平行四辺形）を定義
    double x1 = m_settings.boundary_points[0].x + deformation_x, y1 = m_settings.boundary_points[0].y + deformation_y;
    double x2 = m_settings.boundary_points[1].x + deformation_x, y2 = m_settings.boundary_points[1].y + deformation_y; 
    double x3 = m_settings.boundary_points[2].x + deformation_x, y3 = m_settings.boundary_points[2].y + deformation_y;
    
    // 4番目の点を計算（平行四辺形の対角点）
    double x4 = x2 + x3 - x1;
    double y4 = y2 + y3 - y1;
    
    // 各辺でn個の点を生成（4辺で4n個）
    // 辺1: (x1,y1) -> (x2,y2)
    for (int i = 0; i < n; i++) {
        double t = (double)i / n;
        double x = x1 + t * (x2 - x1);
        double y = y1 + t * (y2 - y1);
        boundary_points.push_back({x, y});
    }
    
    // 辺2: (x2,y2) -> (x4,y4)
    for (int i = 0; i < n; i++) {
        double t = (double)i / n;
        double x = x2 + t * (x4 - x2);
        double y = y2 + t * (y4 - y2);
        boundary_points.push_back({x, y});
    }
    
    // 辺3: (x4,y4) -> (x3,y3)
    for (int i = 0; i < n; i++) {
        double t = (double)i / n;
        double x = x4 + t * (x3 - x4);
        double y = y4 + t * (y3 - y4);
        boundary_points.push_back({x, y});
    }
    
    // 辺4: (x3,y3) -> (x1,y1)
    for (int i = 0; i < n; i++) {
        double t = (double)i / n;
        double x = x3 + t * (x1 - x3);
        double y = y3 + t * (y1 - y3);
        boundary_points.push_back({x, y});
    }
    
    return boundary_points;
}

// 境界の圧力値だけが格納された(n+1)×(n+1)の行列を作成
Eigen::MatrixXd Module::create_pressure_matrix(const std::vector<double>& boundary_pressures, int n)
{
    Eigen::MatrixXd pressure_matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
    
    if (boundary_pressures.size() != 4 * n) {
        std::cerr << "Error: Expected " << 4 * n << " boundary pressure values, got " 
                  << boundary_pressures.size() << std::endl;
        return pressure_matrix;
    }
    
    // 境界値を行列の境界要素に配置
    int idx = 0;
    
    // 下辺 (i=0, j=0 to n)
    for (int j = 0; j < n; j++) {
        pressure_matrix(0, j) = boundary_pressures[idx++];
    }
    
    // 右辺 (i=0 to n, j=n)
    for (int i = 0; i < n; i++) {
        pressure_matrix(i, n) = boundary_pressures[idx++];
    }
    
    // 上辺 (i=n, j=n to 0)
    for (int j = n; j > 0; j--) {
        pressure_matrix(n, j) = boundary_pressures[idx++];
    }
    
    // 左辺 (i=n to 0, j=0)
    for (int i = n; i > 0; i--) {
        pressure_matrix(i, 0) = boundary_pressures[idx++];
    }
    
    return pressure_matrix;
}

// Update wet-dry matrix based on infiltration spreading (Darcy's law)
// Models infiltration spreading from each cell similar to Conway's Game of Life
void Module::update_wet_dry_matrix(const pw::api::Session & session, const std::vector<double>& boundary_pressures)
{
    std::cout << "Updating wet-dry matrix based on infiltration spreading (Darcy's law)" << std::endl;

    int n = m_settings.division_number;
    double lattice_x = std::abs(m_settings.boundary_points[1].x - m_settings.boundary_points[0].x);
    double lattice_y = std::abs(m_settings.boundary_points[2].y - m_settings.boundary_points[0].y);

    // Calculate cell dimensions
    double dx = lattice_x / n;  // Cell size in X direction
    double dy = lattice_y / n;  // Cell size in Y direction

    // Define threshold values
    const double eps = 1e-6;  // Minimum infiltration threshold
    double threshold_x = dx * dx;  // Threshold for X-direction spreading
    double threshold_y = dy * dy;  // Threshold for Y-direction spreading
    double threshold_diagonal = dx * dx + dy * dy;  // Threshold for diagonal spreading

    // Create a temporary matrix to store new wet cells (to avoid modifying during iteration)
    Eigen::MatrixXd new_wet_dry_matrix = m_wet_dry_matrix;

    // Process each grid cell
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            double infil_dist = m_grid_infil_distance(i, j);

            // Step 1: If infiltration distance exceeds eps, mark current cell as wet
            if (infil_dist > eps) {
                new_wet_dry_matrix(i, j) = 1.0;

                // Step 2: Check if we should propagate to X-direction neighbors
                if (infil_dist > threshold_x) {
                    // Set adjacent cells in X direction (i-1, i+1) to wet
                    if (i > 0) {
                        new_wet_dry_matrix(i - 1, j) = 1.0;  // South neighbor
                    }
                    if (i < n) {
                        new_wet_dry_matrix(i + 1, j) = 1.0;  // North neighbor
                    }
                }

                // Step 3: Check if we should propagate to Y-direction neighbors
                if (infil_dist > threshold_y) {
                    // Set adjacent cells in Y direction (j-1, j+1) to wet
                    if (j > 0) {
                        new_wet_dry_matrix(i, j - 1) = 1.0;  // West neighbor
                    }
                    if (j < n) {
                        new_wet_dry_matrix(i, j + 1) = 1.0;  // East neighbor
                    }
                }

                // Step 4: Check if we should propagate to diagonal neighbors
                if (infil_dist > threshold_diagonal) {
                    // Set diagonal neighbors to wet
                    if (i > 0 && j > 0) {
                        new_wet_dry_matrix(i - 1, j - 1) = 1.0;  // Southwest
                    }
                    if (i > 0 && j < n) {
                        new_wet_dry_matrix(i - 1, j + 1) = 1.0;  // Southeast
                    }
                    if (i < n && j > 0) {
                        new_wet_dry_matrix(i + 1, j - 1) = 1.0;  // Northwest
                    }
                    if (i < n && j < n) {
                        new_wet_dry_matrix(i + 1, j + 1) = 1.0;  // Northeast
                    }
                }
                // Once diagonal threshold is exceeded, we don't need to track further neighbors
                // (already handled above)
            }
        }
    }

    // Update the wet-dry matrix
    m_wet_dry_matrix = new_wet_dry_matrix;

    // Count wet points for logging
    int wet_count = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (m_wet_dry_matrix(i, j) > 0.5) {
                wet_count++;
            }
        }
    }

    std::cout << "Wet-dry matrix updated: " << wet_count << "/" << ((n+1)*(n+1)) << " points are wet" << std::endl;
    std::cout << "  Cell size: dx=" << dx << ", dy=" << dy << std::endl;
    std::cout << "  Thresholds: X=" << threshold_x << ", Y=" << threshold_y << ", Diagonal=" << threshold_diagonal << std::endl;
}

void Module::update_infil_distance_matrix(const pw::api::Session & session, const Eigen::MatrixXd& pressure_matrix)
{   
    int n = m_settings.division_number;
    double dt;
	PW_SOLVER_get_delta_time(session.solver().handle(), &dt);
    for (int i=0; i <= n; ++i){
        for (int j=0; j <= n; ++j){
            m_grid_infil_distance(i,j) += m_settings.hydraulic_conductivity * pressure_matrix(i,j) * dt / m_settings.kinematic_viscosity_coefficient;
        }
    }
}

// メイン圧力分布計算関数(this function is called from the session)
void Module::main(const pw::api::Session& session)
{
    // Get current simulation time
    double ctime;
    PW_SOLVER_get_current_time(session.solver().handle(), &ctime);

    int n = m_settings.division_number;
    double lattice_x = std::abs(-m_settings.boundary_points[0].x + m_settings.boundary_points[1].x);
    double lattice_y = std::abs(-m_settings.boundary_points[0].y + m_settings.boundary_points[2].y);

    if (!m_grid_coords_initialized) {
        initialize_grid_coordinates(n, lattice_x, lattice_y);
    }

    if (!m_original_boundary_coords_initialized) {
        m_original_boundary_coords = generate_boundary_points(n, 0.0, 0.0);
        m_original_boundary_coords_initialized = true;
    }

    // Initialize distance matrix once
    if (!m_distance_matrix_initialized) {
        initialize_distance_matrix(n);
    }

    // Initialize wet-dry matrix once
    if (!m_wet_dry_matrix_initialized) {
        initialize_wet_dry_matrix(n);
    }

    // Initialize infil distance matrix
    if (!m_grid_infil_distance_initialized){
        initialize_infil_distance_matrix(n);
    }

    // Check if current time is before calculation start time
    if (ctime < m_settings.calculation_start_time) {
        std::cout << "Skipping pressure distribution calculation (current time: " << ctime
                  << " < start time: " << m_settings.calculation_start_time << ")" << std::endl;

        // Still need to update position/velocity for next timestep
        // Use zero forces before calculation starts
        pretime_process(session);
        set_node_position(session);

        // Save force history periodically even before calculation starts
        if (tid % m_settings.force_history_save_interval == 0) {
            save_force_history();
        }
        return;
    }

    std::cout << "Calculating pressure distribution with n=" << n << " at time=" << ctime << std::endl;

    // 物体の移動を考慮して境界上の4n個の点を生成
    // Position is in PW units (mm/cm/m), convert to meters for pressure calculation
    double pw_to_m = get_pw_to_meters_conversion();
    auto boundary_coords = generate_boundary_points(n, Position.back() / pw_to_m, 0.0);
    std::cout << "Generated " << boundary_coords.size() << " boundary points" << std::endl;
    
    // 各境界点での圧力を計算
    std::vector<double> boundary_pressures;
    boundary_pressures.reserve(boundary_coords.size());//メモリを予約
    
    double z = m_settings.boundary_points[0].z; // Z座標は最初の境界点から取得
    
    for (const auto& point : boundary_coords) {
        double pressure = calculate_pressure_at_point(
            session, point.first, point.second, z,
            m_settings.measurement_radius, m_settings.DF_index);
        boundary_pressures.push_back(pressure);
    }
    
    // (n+1)×(n+1)行列を作成
    boundary_conditions = create_pressure_matrix(boundary_pressures, n);

    // 圧力分布,上向き合力を計算
    Eigen::MatrixXd pressure_distribution = solve_pressure_distribution(boundary_conditions, n);

    // 浸潤半径行列を更新
    update_infil_distance_matrix(session, pressure_distribution);

    // 圧力分布が全て0でなく、dry点が存在しているときのみwet-dry行列の更新
    int total_points = (n+1) * (n+1);
    int count_wet = (m_wet_dry_matrix.array() > 0.9).count();
    bool has_dry_points = (count_wet != total_points);
    bool has_pressure = !pressure_distribution.isZero(1e-6);

    if (has_pressure && has_dry_points) {
        update_wet_dry_matrix(session, boundary_pressures);
    }

    calculate_new_position(session, pressure_distribution, m_wet_dry_matrix);
    set_node_position(session);

    // Save force history periodically
    if (tid % m_settings.force_history_save_interval == 0) {
        save_force_history();
    }

    // Save pressure grid periodically
    if (tid % m_settings.pressure_grid_save_interval == 0) {
        save_pressure_grid(pressure_distribution, m_wet_dry_matrix, ctime);
    }
}

// Creates coefficient matrix for solving Poisson equation with iterative solver
// Compressed system: builds matrix only for wet points with P=0 at wet-dry interface
Eigen::SparseMatrix<double> Module::create_coefficient_matrix(int n, const Eigen::MatrixXd& wet_dry_matrix,
                                                               const Eigen::MatrixXi& grid_to_compressed)
{
    // Count wet points to determine matrix size
    int wet_point_count = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (grid_to_compressed(i, j) >= 0) {
                wet_point_count++;
            }
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(9 * wet_point_count);  // Estimate for 9-point stencil

    std::cout << "Building compressed coefficient matrix for " << wet_point_count << " wet points (9-point stencil)..." << std::endl;

    // Create coefficient matrix for 2D Poisson equation using finite differences
    // For each wet interior point, apply standard 9-point stencil
    // At wet-dry interfaces, apply P=0 Dirichlet boundary condition
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            // Skip dry points and original boundary points
            if (grid_to_compressed(i, j) < 0) continue;

            // Skip original boundary points (i=0, i=n, j=0, j=n)
            if (i == 0 || i == n || j == 0 || j == n) continue;

            int row = grid_to_compressed(i, j);
            double diagonal_value = 0.0;

            // 9-point stencil: check all 8 neighbors (4 direct + 4 diagonal)
            // Direct neighbors (N, S, E, W) and diagonal neighbors (NE, NW, SE, SW)

            // South neighbor (i-1, j)
            if (i > 0) {
                if (wet_dry_matrix(i-1, j) > 0.5 && grid_to_compressed(i-1, j) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i-1, j);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // North neighbor (i+1, j)
            if (i < n) {
                if (wet_dry_matrix(i+1, j) > 0.5 && grid_to_compressed(i+1, j) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i+1, j);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // West neighbor (i, j-1)
            if (j > 0) {
                if (wet_dry_matrix(i, j-1) > 0.5 && grid_to_compressed(i, j-1) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i, j-1);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // East neighbor (i, j+1)
            if (j < n) {
                if (wet_dry_matrix(i, j+1) > 0.5 && grid_to_compressed(i, j+1) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i, j+1);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // Southwest neighbor (i-1, j-1)
            if (i > 0 && j > 0) {
                if (wet_dry_matrix(i-1, j-1) > 0.5 && grid_to_compressed(i-1, j-1) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i-1, j-1);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // Southeast neighbor (i-1, j+1)
            if (i > 0 && j < n) {
                if (wet_dry_matrix(i-1, j+1) > 0.5 && grid_to_compressed(i-1, j+1) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i-1, j+1);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // Northwest neighbor (i+1, j-1)
            if (i < n && j > 0) {
                if (wet_dry_matrix(i+1, j-1) > 0.5 && grid_to_compressed(i+1, j-1) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i+1, j-1);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // Northeast neighbor (i+1, j+1)
            if (i < n && j < n) {
                if (wet_dry_matrix(i+1, j+1) > 0.5 && grid_to_compressed(i+1, j+1) >= 0) {
                    // Wet neighbor - add to matrix
                    int col = grid_to_compressed(i+1, j+1);
                    triplets.emplace_back(row, col, 1.0);
                    diagonal_value -= 1.0;
                } else {
                    // Dry neighbor or original boundary - P=0 Dirichlet BC (no contribution)
                    diagonal_value -= 1.0;
                }
            }

            // Add diagonal entry
            triplets.emplace_back(row, row, diagonal_value);
        }
    }

    // Build sparse matrix
    Eigen::SparseMatrix<double> A(wet_point_count, wet_point_count);
    A.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "Compressed coefficient matrix created: " << wet_point_count << "x" << wet_point_count
              << " system with " << triplets.size() << " non-zero entries" << std::endl;

    return A;
}

// Constructs the right-hand side vector from boundary conditions
// Compressed system: constructs RHS for wet points only with P=0 at wet-dry interface
Eigen::VectorXd Module::construct_rhs_vector(const Eigen::MatrixXd& boundary_matrix, int n,
                                             const Eigen::MatrixXd& wet_dry_matrix,
                                             const Eigen::MatrixXi& grid_to_compressed)
{
    // Count wet points to determine vector size
    int wet_point_count = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (grid_to_compressed(i, j) >= 0) {
                wet_point_count++;
            }
        }
    }

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(wet_point_count);

    std::cout << "Building compressed RHS vector for " << wet_point_count << " wet points (9-point stencil)..." << std::endl;

    // Fill RHS vector with boundary condition contributions
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            // Skip dry points and original boundary points
            if (grid_to_compressed(i, j) < 0) continue;
            if (i == 0 || i == n || j == 0 || j == n) continue;

            int idx = grid_to_compressed(i, j);

            // Add boundary contributions (negative because moved to RHS)
            // The boundary can be either:
            // 1. Original boundary (i=0, i=n, j=0, j=n) with measured pressure
            // 2. Wet-dry interface with P=0 (Dirichlet BC)

            // Direct neighbors (4-point contribution)

            // South neighbor (i-1, j)
            if (i > 0) {
                if (i-1 == 0) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(0, j);
                } else if (wet_dry_matrix(i-1, j) < 0.5 || grid_to_compressed(i-1, j) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // North neighbor (i+1, j)
            if (i < n) {
                if (i+1 == n) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(n, j);
                } else if (wet_dry_matrix(i+1, j) < 0.5 || grid_to_compressed(i+1, j) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // West neighbor (i, j-1)
            if (j > 0) {
                if (j-1 == 0) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(i, 0);
                } else if (wet_dry_matrix(i, j-1) < 0.5 || grid_to_compressed(i, j-1) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // East neighbor (i, j+1)
            if (j < n) {
                if (j+1 == n) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(i, n);
                } else if (wet_dry_matrix(i, j+1) < 0.5 || grid_to_compressed(i, j+1) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // Diagonal neighbors (4-point contribution for 9-point stencil)

            // Southwest neighbor (i-1, j-1)
            if (i > 0 && j > 0) {
                if ((i-1 == 0 || j-1 == 0)) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(i-1, j-1);
                } else if (wet_dry_matrix(i-1, j-1) < 0.5 || grid_to_compressed(i-1, j-1) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // Southeast neighbor (i-1, j+1)
            if (i > 0 && j < n) {
                if ((i-1 == 0 || j+1 == n)) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(i-1, j+1);
                } else if (wet_dry_matrix(i-1, j+1) < 0.5 || grid_to_compressed(i-1, j+1) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // Northwest neighbor (i+1, j-1)
            if (i < n && j > 0) {
                if ((i+1 == n || j-1 == 0)) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(i+1, j-1);
                } else if (wet_dry_matrix(i+1, j-1) < 0.5 || grid_to_compressed(i+1, j-1) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }

            // Northeast neighbor (i+1, j+1)
            if (i < n && j < n) {
                if ((i+1 == n || j+1 == n)) {
                    // Original boundary
                    rhs(idx) -= boundary_matrix(i+1, j+1);
                } else if (wet_dry_matrix(i+1, j+1) < 0.5 || grid_to_compressed(i+1, j+1) < 0) {
                    // Wet-dry interface: P=0 (no contribution to RHS)
                    rhs(idx) -= 0.0;
                }
            }
        }
    }

    std::cout << "Compressed RHS vector constructed with " << wet_point_count << " elements" << std::endl;
    return rhs;
}

// Solves the Poisson equation to obtain pressure distribution of inner points using iterative solver
// Uses compressed system that only includes wet points with P=0 at wet-dry interface
Eigen::MatrixXd Module::solve_pressure_distribution(const Eigen::MatrixXd& boundary_matrix, int n)
{
    std::cout << "========== Solving pressure distribution with compressed wet-only system ==========" << std::endl;

    // Step 1: Create renumbering map from (i,j) grid coordinates to compressed indices
    // -1 means point is excluded from compressed system (dry or original boundary)
    Eigen::MatrixXi grid_to_compressed = Eigen::MatrixXi::Constant(n + 1, n + 1, -1);

    int compressed_index = 0;

    // First pass: assign indices to wet interior points
    for (int i = 1; i < n; ++i) {  // Skip original boundaries (i=0, i=n)
        for (int j = 1; j < n; ++j) {  // Skip original boundaries (j=0, j=n)
            if (m_wet_dry_matrix(i, j) > 0.5) {  // Wet point
                grid_to_compressed(i, j) = compressed_index++;
            }
        }
    }

    int total_wet_points = compressed_index;
    std::cout << "Grid renumbering complete: " << total_wet_points << " wet interior points" << std::endl;
    std::cout << "Total grid points (including boundaries): " << (n+1)*(n+1) << std::endl;
    std::cout << "Compression ratio: " << (double)total_wet_points / ((n+1)*(n+1)) * 100.0 << "%" << std::endl;

    // Handle edge case: if no wet points, return zero matrix
    if (total_wet_points == 0) {
        std::cout << "Warning: No wet interior points found. Returning zero pressure matrix." << std::endl;
        Eigen::MatrixXd pressure_matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
        // Keep original boundary values
        pressure_matrix = boundary_matrix;
        // Set all interior points to 0
        for (int i = 1; i < n; ++i) {
            for (int j = 1; j < n; ++j) {
                pressure_matrix(i, j) = 0.0;
            }
        }
        return pressure_matrix;
    }

    // Step 2: Create compressed coefficient matrix for wet points only
    std::cout << "Creating compressed coefficient matrix for iterative solver..." << std::endl;
    Eigen::SparseMatrix<double> A = create_coefficient_matrix(n, m_wet_dry_matrix, grid_to_compressed);

    // Step 3: Construct compressed RHS vector with P=0 at wet-dry interface
    Eigen::VectorXd rhs = construct_rhs_vector(boundary_matrix, n, m_wet_dry_matrix, grid_to_compressed);

    // Step 4: Configure and solve iterative system
    std::cout << "Configuring BiCGSTAB iterative solver..." << std::endl;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;

    solver.setMaxIterations(m_max_iterations);
    solver.setTolerance(m_solver_tolerance);

    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Error: Failed to initialize iterative solver" << std::endl;
        return boundary_matrix;  // Return boundary matrix as fallback
    }

    std::cout << "Solving compressed linear system (" << total_wet_points << " unknowns)..." << std::endl;
    Eigen::VectorXd solution = solver.solve(rhs);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Error: Iterative solver failed to converge" << std::endl;
        std::cerr << "  Iterations: " << solver.iterations() << std::endl;
        std::cerr << "  Estimated error: " << solver.error() << std::endl;
        return boundary_matrix;  // Return boundary matrix as fallback
    }

    std::cout << "Iterative solver converged successfully:" << std::endl;
    std::cout << "  Iterations: " << solver.iterations() << "/" << m_max_iterations << std::endl;
    std::cout << "  Estimated error: " << solver.error() << " (tolerance: " << m_solver_tolerance << ")" << std::endl;

    // Step 5: Map compressed solution back to full grid
    std::cout << "Mapping compressed solution back to full grid..." << std::endl;
    Eigen::MatrixXd pressure_matrix = boundary_matrix;  // Start with boundary values

    // Fill wet interior points with solution
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < n; ++j) {
            int compressed_idx = grid_to_compressed(i, j);
            if (compressed_idx >= 0) {
                // Wet point - use solved pressure
                pressure_matrix(i, j) = solution(compressed_idx);
            } else {
                // Dry point - P=0 (atmospheric pressure)
                pressure_matrix(i, j) = 0.0;
            }
        }
    }

    std::cout << "Pressure distribution solved successfully with physically correct wet-dry treatment" << std::endl;
    std::cout << "==============================================================================\n" << std::endl;

    return pressure_matrix;
}

void Module::calculate_moment(const pw::api::Session & session, const Eigen::MatrixXd& pressure_distribution, const Eigen::MatrixXd& wet_dry_matrix){
    // Calculate moments from pressure distribution for overturning analysis
    double aload[3];
    double atorque[3];

    PW_DF_get_load(session.distance_fields()[m_settings.DF_index].handle(),PW_LOAD_sum_c,aload,atorque);

    Externalforce_x.emplace_back(aload[0]);
    HorizontalTorque.emplace_back(atorque[1]);

    // Grid parameters
    int n = m_settings.division_number;
    double lattice_x = std::abs(m_settings.boundary_points[1].x - m_settings.boundary_points[0].x);
    double lattice_y = std::abs(m_settings.boundary_points[2].y - m_settings.boundary_points[0].y);
    double dx = lattice_x / n;
    double dy = lattice_y / n;
    double area_element = dx * dy;

    int rows = pressure_distribution.rows();
    int cols = pressure_distribution.cols();

    // Calculate total uplift force and moment from pressure distribution
    double total_pressure_force = 0.0;
    double pressure_moment = 0.0;  // Moment from pressure around rotation center

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double weight = 1.0;

            // Apply trapezoidal rule weights for numerical integration
            if ((i == 0 || i == rows - 1) && (j == 0 || j == cols - 1)) {
                weight = 0.25;  // Corner points
            } else if (i == 0 || i == rows - 1 || j == 0 || j == cols - 1) {
                weight = 0.5;   // Edge points
            }

            // Calculate pressure force at this point
            double local_pressure_force = pressure_distribution(i, j) * area_element * weight;
            total_pressure_force += local_pressure_force;

            // Calculate moment arm (distance from rotation center to this point)
            // Using the distance matrix that was calculated during initialization
            double moment_arm = m_distance_matrix(i, j);

            // Pressure creates upward force, moment is force × distance
            pressure_moment += local_pressure_force * moment_arm;
        }
    }

    // Calculate gravity moment (weight acting at center of gravity)
    double total_weight = m_settings.mass * m_settings.gravity[2];
    double cog_x = m_settings.center_of_gravity[0];
    double cog_y = m_settings.center_of_gravity[1];
    double rotation_x = m_settings.rotation_center[0];
    double rotation_y = m_settings.rotation_center[1];

    // Moment arm from rotation center to center of gravity
    double gravity_moment_arm = std::abs(cog_x - rotation_x);
    double gravity_moment = total_weight * gravity_moment_arm;

    // Net moment (pressure moment - gravity moment + resistance moment)
    // Positive moment tends to overturn, negative moment stabilizes
    double net_moment = pressure_moment + atorque[1]  - gravity_moment;

    // Store uplift force for history
    Upliftforce.emplace_back(total_pressure_force);
    UpliftTorque.emplace_back(pressure_moment);
    TotalTorque.emplace_back(net_moment);

    std::cout << "=== Moment and Overturning Analysis ===" << std::endl;
    std::cout << "Total pressure force (uplift): " << total_pressure_force << " N" << std::endl;
    std::cout << "Pressure moment (overturning): " << pressure_moment << " N·m" << std::endl;
    std::cout << "Gravity moment (stabilizing): " << gravity_moment << " N·m" << std::endl;
    std::cout << "Resistance moment (piles): " << m_settings.resistance_moment << " N·m" << std::endl;
    std::cout << "Net moment: " << net_moment << " N·m" << std::endl;
    if (net_moment > 0) {
        std::cout << "WARNING: Positive net moment - building may overturn!" << std::endl;
    } else {
        std::cout << "Building is stable (negative net moment)" << std::endl;
    }
    std::cout << "========================================" << std::endl;
}


void Module::calculate_new_position(const api::Session & session, const Eigen::MatrixXd& pressure_distribution, const Eigen::MatrixXd& wet_dry_matrix)
{
    double dt, ctime;
	PW_SOLVER_get_delta_time(session.solver().handle(), &dt);
    PW_SOLVER_get_current_time(session.solver().handle(), &ctime);

    calculate_moment(session, pressure_distribution, wet_dry_matrix);

	// -- External Force parameters
    // Force is in N (kg·m/s²), mass is in kg, dt is in s
    // Acceleration is in m/s², multiply by pw_to_m to convert to PW units/s²
    // Velocity is stored in PW units/s (e.g., mm/s if PW uses mm)
    double rad_to_deg = 57.2958;
    AngularVelocity.emplace_back(AngularVelocity.back() + (TotalTorque.back() / m_settings.moment_of_inertia) * dt * rad_to_deg);
	double new_Angle =  Angle.back() +  AngularVelocity.back() * dt;
    Angle.emplace_back(new_Angle);
    Time.emplace_back(ctime);
}

void
Module::set_node_position(const api::Session & session)
{
    double ctime,dtime;
	PW_SOLVER_get_current_time(session.solver().handle(), &ctime);
	PW_SOLVER_get_delta_time(session.solver().handle(), &dtime);
        for (int i=0;i<session.nodes().size();i++) {
            if (session.nodes()[i].object_id()==m_settings.object_id) {
         	    pw::api::Animation ani=session.nodes()[i].animation("transform.location",0); //移動についてX軸方向についてのアニメーション設定を取得
                ani.x.emplace_back(ctime); // xに時間
                ani.y.emplace_back(Position.back()); // yに座標
   	            session.nodes()[i].animation("transform.location",0,ani);//位置を変更、X軸、移動の設定を入力
		        session.nodes()[i].update_motion(ctime);
	        }
        }
    tid++;
}

void
Module::pretime_process(const api::Session & session)
{
    double dt, ctime;
    PW_SOLVER_get_delta_time(session.solver().handle(), &dt);
    PW_SOLVER_get_current_time(session.solver().handle(), &ctime);

    // 揚圧力が0の状態で運動を計算
    // Create zero pressure distribution and dry matrix for pretime processing
    int n = m_settings.division_number;
    Eigen::MatrixXd zero_pressure = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::MatrixXd dry_matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);  // All dry initially

    calculate_new_position(session, zero_pressure, dry_matrix);

    std::cout << "Pretime process: No forces applied, Position=0.0, Velocity=0.0" << std::endl;
}

void Module::save_force_history() const
{
    std::ofstream file(m_settings.force_history_output_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << m_settings.force_history_output_file << " for writing" << std::endl;
        return;
    }

    // Write header(Upliftforce is the total uplift force)
    file << "Time,Externalforce_x,Upliftforce,Force,Position,Velocity" << std::endl;

    // Write data
    for (size_t i = 0; i < Time.size(); ++i) {
        file << Time[i] << ","
             << (i < Externalforce_x.size() ? Externalforce_x[i] : 0.0) << ","
             << (i < Upliftforce.size() ? Upliftforce[i] : 0.0) << ","
             << (i < Force.size() ? Force[i] : 0.0) << ","
             << (i < Position.size() ? Position[i] : 0.0) << ","
             << (i < Velocity.size() ? Velocity[i] : 0.0) << std::endl;
    }

    file.close();
    std::cout << "Force history saved to " << m_settings.force_history_output_file << " with " << Time.size() << " data points" << std::endl;
}

void Module::save_pressure_grid(const Eigen::MatrixXd& pressure_distribution,
                                const Eigen::MatrixXd& wet_dry_matrix,
                                double current_time) const
{
    // Check if file exists to determine if we need to write header
    bool file_exists = std::ifstream(m_settings.pressure_grid_output_file).good();

    // Open file in append mode
    std::ofstream file(m_settings.pressure_grid_output_file, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << m_settings.pressure_grid_output_file << " for writing" << std::endl;
        return;
    }

    // Write header if file is new
    if (!file_exists) {
        file << "Time,i,j,x,y,z,Pressure";
        if (m_settings.pressure_grid_include_wet_dry) {
            file << ",WetDry";
        }
        file << std::endl;
    }

    // Get grid parameters
    int n = m_settings.division_number;
    int grid_size = n + 1;  // (n+1) x (n+1) grid

    // Get boundary coordinates to calculate grid coordinates
    const auto& boundary_points = m_settings.boundary_points;
    if (boundary_points.size() < 3) {
        std::cerr << "Error: Not enough boundary points to calculate grid coordinates" << std::endl;
        return;
    }

    // Calculate lattice dimensions from boundary points
    // Point1 and Point2 define x-direction, Point1 and Point3 define y-direction
    double origin_x = boundary_points[0].x;
    double origin_y = boundary_points[0].y;
    double origin_z = boundary_points[0].z;

    double lattice_x = std::sqrt(std::pow(boundary_points[1].x - boundary_points[0].x, 2) +
                                  std::pow(boundary_points[1].y - boundary_points[0].y, 2));
    double lattice_y = std::sqrt(std::pow(boundary_points[2].x - boundary_points[0].x, 2) +
                                  std::pow(boundary_points[2].y - boundary_points[0].y, 2));

    double dx = lattice_x / n;
    double dy = lattice_y / n;

    // Write grid data
    int points_written = 0;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            // Calculate grid coordinates
            double x = origin_x + j * dx;
            double y = origin_y + i * dy;
            double z = origin_z;

            // Get pressure value
            double pressure = pressure_distribution(i, j);

            // Write data line
            file << current_time << ","
                 << i << ","
                 << j << ","
                 << x << ","
                 << y << ","
                 << z << ","
                 << pressure;

            // Write wet/dry status if requested
            if (m_settings.pressure_grid_include_wet_dry) {
                double wet_dry = wet_dry_matrix(i, j);
                file << "," << wet_dry;
            }

            file << std::endl;
            points_written++;
        }
    }

    file.close();
    std::cout << "Pressure grid saved to " << m_settings.pressure_grid_output_file
              << " at time " << current_time << " (" << points_written << " points)" << std::endl;
}

// Calculate (n+1)×(n+1) grid coordinates representing pressure distribution points
void Module::initialize_grid_coordinates(int n, double lattice_x, double lattice_y)
{
    // Create 3D coordinate matrix: (n+1)×(n+1)×3 where last dimension is [x,y,z]
    m_grid_coords.resize(n + 1, (n + 1) * 3);

    if (m_settings.boundary_points.size() < 3) {
        std::cerr << "Error: Need at least 3 boundary points to calculate grid coordinates" << std::endl;
    }

    // Use the first boundary point as the origin reference
    double origin_x = m_settings.boundary_points[0].x;
    double origin_y = m_settings.boundary_points[0].y;
    double origin_z = m_settings.boundary_points[0].z;

    // Calculate grid spacing
    double dx = lattice_x / n;
    double dy = lattice_y / n;

    // Generate grid coordinates
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            double x = origin_x + j * dx;
            double y = origin_y + i * dy;
            double z = origin_z;  // Assuming constant z-level

            // Store coordinates in flattened format: [x1,y1,z1,x2,y2,z2,...]
            int col_base = j * 3;
            m_grid_coords(i, col_base) = x;
            m_grid_coords(i, col_base + 1) = y;
            m_grid_coords(i, col_base + 2) = z;
        }
    }
    m_grid_coords_initialized = true;
    std::cout << "Grid coordinates calculated for " << (n+1) << "×" << (n+1) << " points" << std::endl;
}

// Initialize distance matrix from rotation center to each pressure measurement point
// calculate_moment_arm
void Module::initialize_distance_matrix(int n)
{
    if (m_distance_matrix_initialized && m_cached_n == n) {
        std::cout << "Distance matrix already initialized for n=" << n << std::endl;
        return;
    }

    double lattice_x = std::abs(m_settings.boundary_points[1].x - m_settings.boundary_points[0].x);
    double lattice_y = std::abs(m_settings.boundary_points[2].y - m_settings.boundary_points[0].y);

    // Initialize distance matrix
    m_distance_matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);

    double rotation_center_x = m_settings.rotation_center[0];
    double rotation_center_y = m_settings.rotation_center[1];
    double rotation_center_z = m_settings.rotation_center[2];

    // Calculate distances from rotation center to each grid point
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            int col_base = j * 3;
            double point_x = m_grid_coords(i, col_base);
            double point_y = m_grid_coords(i, col_base + 1);
            double point_z = m_grid_coords(i, col_base + 2);

            // Calculate 3D Euclidean distance
            double dx = point_x - rotation_center_x;
            double distance = std::abs(dx);

            m_distance_matrix(i, j) = distance;
        }
    }

    m_distance_matrix_initialized = true;
    std::cout << "Distance matrix initialized successfully for " << (n+1) << "×" << (n+1) << " grid" << std::endl;
    std::cout << "Rotation center: (" << rotation_center_x << ", " << rotation_center_y << ", " << rotation_center_z << ")" << std::endl;
}

// Initialize wet-dry matrix (1 = wet, 0 = dry)
void Module::initialize_wet_dry_matrix(int n)
{
    if (m_wet_dry_matrix_initialized && m_cached_n == n) {
        std::cout << "Wet-dry matrix already initialized for n=" << n << std::endl;
        return;
    }

    // Initialize with all points as dry (value = 0.0)
    m_wet_dry_matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
    m_wet_dry_matrix_initialized = true;
    std::cout << "Wet-dry matrix initialized successfully for " << (n+1) << "×" << (n+1) << " grid (all dry)" << std::endl;
}

void Module::initialize_infil_distance_matrix(int n)
{
    if (m_grid_infil_distance_initialized && m_cached_n == n){
        std::cout << "grid infil distance matrix already initialized for n=" << n << std::endl;
        return;
    }
    m_grid_infil_distance = Eigen::MatrixXd::Zero(n + 1, n + 1);
    m_grid_infil_distance_initialized = true;
    std::cout << "grid infil distance matrix initialized successfully for " << (n+1) << "×" << (n+1) << " grid (all dry)" << std::endl;
}

// ユーザー関数クラス実装
class PressureDistributionSolver : public api::UserFunction {
protected:
    api::Session& m_session;
    Module& m_module;
    
public:
    PressureDistributionSolver(api::Session& session, Module& module) 
        : m_session(session), m_module(module) {}
    
    virtual std::string name() {
        return "point_pressure_measurement";
    }
    
    virtual Arrays input_arrays() {
        return Arrays(PW_ARRAY_OWNER_mps_c, {"position", "pressure"});
    }
    
    virtual Arrays output_arrays() {
        return Arrays();
    }
    
    virtual void execute(PW_CALL_PHASE_t phase) {
        if (phase != PW_CALL_PHASE_execute_c) {
            return;
        }
        m_module.main(m_session);
    }
};

// Module初期化
void Module::initialize()
{
    std::cout << "Initializing Point Pressure Module" << std::endl;
    
    // 設定読み込み
    load_settings("settings.json");
        
    // ユーザー関数登録 - 各ステップ後に圧力測定＆圧力分布計算を実行
    add_user_function(
        "update_position_and_velocity_from_force",
        PW_CALL_POINT_post_c,
        std::make_shared<PressureDistributionSolver>(m_session, *this));

    std::cout << "Point Pressure Module initialized successfully" << std::endl;
}

}}}