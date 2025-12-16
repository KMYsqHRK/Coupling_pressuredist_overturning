# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a C++ project for fluid simulation using ParticleWorks SDK. The project simulates friction on structures in fluid environments and provides point pressure measurement capabilities. It uses the Intel OneAPI compiler and requires specific ParticleWorks SDK installation.

## Build Commands

### Prerequisites Setup
Before building, set up the environment:
```bash
source /opt/intel/oneapi/setvars.sh
export PW_SDK_PATH=/home/labkamiya/prometech/sdk
```

### Build Process
```bash
make clean
make all
```

This creates the executable at `bin/pw.executable.example`.

## Running the Application

### Preprocessing
```bash
/home/labkamiya/prometech/bin/app.pre.double -p . -f sub.json -n 4
```

### Main Execution
```bash
/home/labkamiya/repos/mapping/bin/pw.executable.example -p . -f sub.json -n 1 -k cuda -g 0
```

## Architecture

### Core Components

- **example.cpp**: Main program entry point that initializes ParticleWorks session and coordinates pressure measurements
- **example_module.hpp/cpp**: Base module class providing user function management and session handling
- **example_point_pressure.hpp/cpp**: Pressure measurement module with configurable measurement points

### Key Design Patterns

- **Module System**: Uses polymorphic modules inherited from `pw::example::Module` base class
- **Session Management**: ParticleWorks API session (`pw::api::Session`) coordinates solver execution
- **User Functions**: Custom functions are registered with the solver through `add_user_function()`
- **Configuration**: JSON-based configuration in `settings.json` for measurement points

### Data Flow

1. Main program initializes ParticleWorks session
2. Point pressure module loads configuration from `settings.json`
3. Module registers user functions with the solver
4. Solver executes simulation steps while measuring pressure at configured points
5. Results are saved to CSV files for analysis

## Configuration

### settings.json Structure
- `measurement_points`: Array of points with position (x,y,z) and radius
- `output_file`: CSV output filename
- `log_interval`: Measurement frequency
- `filtering.exclude_ghost_particles`: Boolean for particle filtering

## Dependencies

- ParticleWorks SDK (PW_SDK_PATH environment variable)
- Intel OneAPI compiler (icpx)
- Boost libraries (filesystem, program_options, serialization)
- nlohmann/json header (included in third_party/)

## File Organization

- `bin/`: Compiled executables
- `obj/`: Object files
- `third_party/`: External dependencies
- Main source files in root directory