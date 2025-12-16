#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "sdk/particleworks_api.hpp"
#include "example_module.hpp"
#include "example_coupling_friction_sliding.hpp"

// よくわからない事、ループ処理されていることは明示的にしなくてもよいのはわかるが、何秒毎に処理を行うかはどこで決しているのか
int main(int argc, char * argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();
    // ここでキャッシュ行列とキャッシュソルバーを作っておく
    try {
        // Particleworks セッション開始
        pw::api::Session session("point_pressure", argc, argv);
        
        // 点圧力測定モジュールの作成と初期化
        // ここの点圧力測定と同時に右辺を作成するように変更
        std::shared_ptr<pw::example::Module> pressure_module;
        pressure_module = std::make_shared<pw::example::point_pressure::Module>();
        pressure_module->set_session(session.handle());
        pressure_module->initialize();
        
        auto solver_start = std::chrono::high_resolution_clock::now();
        session.solver().step(-1);  // 全ステップ実行
        auto solver_end = std::chrono::high_resolution_clock::now();
        
        auto solver_duration = std::chrono::duration_cast<std::chrono::seconds>(solver_end - solver_start);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred." << std::endl;
        return 1;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n=== Execution Summary ===" << std::endl;
    std::cout << "Total execution time: " << total_duration.count() << " seconds" << std::endl;
    std::cout << "=========================" << std::endl;
    
    return 0;
}