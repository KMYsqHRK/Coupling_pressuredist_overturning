# Coupling_pressuredist_sliding
ParticleWorks上でランタイム中に同時に圧力マッピングを行い、ポアソン方程式を解き、転倒モーメントを考慮して転倒計算を行うコード

次にコンパイル、実行時に注意するべきことを記載する。
なお、Lliux版で実行する前提である。

2025/12/16より開発を開始、ポリゴンによって構成されたポリゴンに対するシミュレーションを行う。
1次元移動を考慮。
衝突についても簡易的に解く。

## 0. Make , OneAPI libraries LD_LIBRARY_PATHS set.

   Before make, you have to set environmental vars to your session.
   
   `source /opt/intel/oneapi/setvars.sh` 
   
   In addition, make sure you have PW libraries(/opt/pw800/lib) set.
   
   `export PW_SDK_PATH=/home/labkamiya/prometech/sdk`
   
   If starvar.sh is not found, Please for installing icp(Intel Compiler)
   ~~~bash
   wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt update

sudo apt install intel-oneapi-compiler-dpcpp-cpp
   ~~~
   After install icp, you have to do this section. 

## 1. In this directory , first do 
~~~bash
make clean
make all
~~~
to compile.

## 2. Files are
   - example.cpp : main program with point pressure measurement
   - example_module.cpp : base module class
   - example_point_pressure.cpp : point pressure measurement module
   - example_point_pressure.hpp : module header
## 3. To run the program go to the project directory where you see,
   - df
   - pre
   - result
## 4. Point Pressure Measurement Configuration

The point pressure measurement module requires a `settings.json` file in the execution directory.

### 4.1 settings.json Format
```json
{
    "description": "Point pressure measurement configuration",
    "output_file": "pressure_results.csv",
    "log_interval": 10,
    "measurement_points": [
        {
            "name": "Point1",
            "initial_position": {
                "x": 0.0,
                "y": 0.0,
                "z": 1.0
            },
            "radius": 2.0
        },
        {
            "name": "Point2",
            "initial_position": {
                "x": 5.0,
                "y": 0.0,
                "z": 1.0
            },
            "radius": 1.5
        }
    ],
    "filtering": {
        "exclude_ghost_particles": true
    },
    "units": {
        "length": "m",
        "time": "s",
        "pressure": "Pa"
    }
}
```

### 4.2 Configuration Parameters
- `output_file`: CSV file name for pressure results
- `log_interval`: Measurement frequency (every N simulation steps)
- `measurement_points`: Array of measurement points with:
  - `name`: Point identifier
  - `initial_position`: x, y, z coordinates
  - `radius`: Measurement radius around the point
- `filtering.exclude_ghost_particles`: Whether to exclude ghost particles
- `units`: Unit descriptions for documentation

### 4.3 Output Format
Results are saved to CSV with columns:
- Time, Point_Name, X, Y, Z, Radius, Average_Pressure, Max_Pressure, Min_Pressure, Particle_Count
## 5. Running the Simulation

### 5.1 Prerequisites
Make sure you have:
- `settings.json` file in the execution directory
- CUDA libraries, PW libraries, OneAPI libraries in LD_LIBRARY_PATHS

### 5.2 Execution Commands

For preprocessing:
```bash
/home/labkamiya/prometech/bin/app.pre.double -p . -f sub.json -n 4
```

For simulation with point pressure measurement:
```bash
/home/labkamiya/repos/Coupling_pressuredist_overturning/bin/pw.executable.example -p . -f sub.json -n 2 -k cuda -g 0
```

### 5.3 Performance Features
- **OpenMP Parallelization**: The pressure calculation loop is parallelized for better performance
- **Optimized Compilation**: Built with `-O3` optimization and Intel compiler
- **Efficient Distance Field Interaction**: Uses ParticleWorks distance field API for neighbor finding

## 6. Result Analysis

### 6.1 Pressure Results
After simulation, check the output CSV file (default: `pressure_results.csv`) for:
- Time series pressure data at each measurement point
- Statistical information (average, min, max pressure)
- Particle count validation

### 6.2 Visualization in ParticleWorks
To visualize the moving polygon:
1. Copy contents of `result/***.motion.csv`
2. Import into `Input/node/***.obj keyframe` in ParticleWorks
3. Use the CSV import feature

### 6.3 Data Processing Tips
- Plot Time vs Average_Pressure for each measurement point
- Monitor Particle_Count to ensure sufficient particles in measurement radius
- Use Excel or similar tools to analyze the CSV output

## 7. Somethings to keep in mind.
   The steps in analysis should probably be,
   - a. Make case in particleworks.
   - b. Make user ***.obj is the parent of the domain
   - c. Open the result/sub.json file to find the object_id of the polygon, so that
      you can write the number in the settings_sdof.csv file.
   - d. Make sure you load the earthquake acceleration file into the 
      -Input/gravity keyframe
       as an external acceleration term to model the inertial effect on the fluid.
       (Should eventually make this automatic).
   - e. When in need look at the API headers in pw800/sdk/include/sdk and examples in
      pw800/sdk/share/src

Currently working on case,

/home/tkoyama/ParticleworksProjects/pw202308071750/scene
