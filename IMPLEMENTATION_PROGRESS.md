# DREAM Quasilinear Diffusion Implementation Progress

## Overview
Implementation of quasilinear diffusion for runaway electron control using whistler waves, based on Zehua Guo's theory.

**Start Date:** 2026-05-01  
**Status:** Phase 1 - Infrastructure Setup (In Progress)

---

## ✅ Completed Work

### Phase 1: Infrastructure Setup (Week 1-2)

#### 1. Extended OptionConstants Enum ✓
**File:** `include/DREAM/Settings/OptionConstants.enum.hpp`

Added three new enumerations:
- `wave_spectrum_type`: WAVE_SPECTRUM_UNIFORM, GAUSSIAN, CUSTOM
- `ql_diffusion_mode`: NEGLECT, INCLUDE
- `ql_harmonic_mode`: N_MINUS_1, N_PLUS_1, BOTH

**Lines added:** 25 lines

#### 2. Created WaveSpectrum C++ Class ✓
**Files:**
- `include/DREAM/Settings/WaveSpectrum.hpp` (147 lines)
- `src/Settings/WaveSpectrum.cpp` (247 lines)

**Features implemented:**
- Uniform spectrum generation in (k, θ_k) space
- Gaussian spectrum with configurable center and spread
- Custom spectrum loading from file
- Amplitude management (uniform or per-mode)
- Grid initialization and accessors
- Information printing

**Key methods:**
```cpp
WaveSpectrum(len_t nk=100, len_t nkt=20);
void setUniformSpectrum(real_t k_min, real_t k_max, 
                       real_t ktheta_min, real_t ktheta_max);
void setGaussianSpectrum(real_t k0, real_t ktheta0, 
                        real_t delta_k, real_t delta_ktheta);
void setCustomSpectrum(const std::string& filename);
void setAmplitude(real_t amp);
```

#### 3. Created WaveSpectrum Python Interface ✓
**File:** `py/DREAM/Settings/WaveSpectrum.py` (172 lines)

**Features:**
- Mirror of C++ functionality in Python
- Dictionary conversion for HDF5 output
- NumPy integration for amplitude arrays
- Comprehensive documentation

**Usage example:**
```python
ws = WaveSpectrum(num_k=100, num_ktheta=20)
ws.setUniformSpectrum(35, 45, 0.1, 0.3)
ws.setAmplitude(1e-10)
```

#### 4. Updated Build System ✓
**File:** `src/CMakeLists.txt`

Added WaveSpectrum.cpp to the build configuration.

#### 5. Testing ✓
**File:** `test_wave_spectrum.py` (167 lines)

Created comprehensive test suite with 5 tests:
1. Uniform spectrum creation
2. Gaussian spectrum creation
3. Custom spectrum setup
4. Individual amplitudes
5. Dictionary conversion

**Test Results:** ✅ ALL TESTS PASSED

---

## 📋 Next Steps

### Immediate Tasks (This Week)

#### 1. Create WhistlerDispersion Class
**Priority:** HIGH  
**Estimated time:** 2-3 days

Implement dispersion relation solver for whistler waves:
- Solve full cold plasma dispersion relation
- Calculate ω(k, θ_k) for given parameters
- Compute group velocity dω/dk
- Calculate polarization vectors (Ex, Ey, Ez)

**Files to create:**
- `include/DREAM/Equations/Kinetic/WhistlerDispersion.hpp`
- `src/Equations/Kinetic/WhistlerDispersion.cpp`

**Key challenge:** Polynomial root finding (11th order)

#### 2. Create ResonanceSolver Class
**Priority:** HIGH  
**Estimated time:** 2-3 days

Implement resonance condition solver:
- Find resonant k values for given (p, ξ, n)
- Handle multiple resonance solutions
- Filter non-physical roots
- Validate against dispersion relation

**Files to create:**
- `include/DREAM/Equations/Kinetic/ResonanceSolver.hpp`
- `src/Equations/Kinetic/ResonanceSolver.cpp`

**Reference:** QUADRE's `calculate_exel_k()` function

#### 3. Create WaveParticleCoupling Class
**Priority:** MEDIUM  
**Estimated time:** 2 days

Calculate coupling strength |Ψ_{n,k}|²:
- Bessel function evaluation
- Polarization vector dot products
- Geometric factors

**Files to create:**
- `include/DREAM/Equations/Kinetic/WaveParticleCoupling.hpp`
- `src/Equations/Kinetic/WaveParticleCoupling.cpp`

**Dependencies:** GSL or Boost for Bessel functions

### Medium-term Tasks (Next 2 Weeks)

#### 4. Create QuasilinearDiffusionTerm Class
**Priority:** HIGH  
**Estimated time:** 1 week

Main diffusion operator implementation:
- Inherit from FVM::DiffusionTerm
- Implement Rebuild() method
- Integrate over wave spectrum
- Assemble diffusion matrix

**Files to create:**
- `include/DREAM/Equations/Kinetic/QuasilinearDiffusionTerm.hpp`
- `src/Equations/Kinetic/QuasilinearDiffusionTerm.cpp`

**Key algorithm:**
```cpp
for each mode m in spectrum:
    find resonant k values
    calculate omega, polarization
    calculate group velocity factor
    calculate coupling strength
    interpolate to grid
    accumulate to D_pp, D_pxi, D_xixi
```

#### 5. Extend DistributionFunction Settings
**Priority:** MEDIUM  
**Estimated time:** 2-3 days

Add quasilinear diffusion settings to distribution function:
- Parse quasilinear diffusion options
- Store WaveSpectrum object
- Enable/disable diffusion term

**Files to modify:**
- `include/DREAM/Settings/Equations/DistributionFunction.hpp`
- `src/Settings/Equations/distribution.cpp`
- `py/DREAM/Settings/Equations/DistributionFunction.py`

#### 6. Integrate into Equation System
**Priority:** HIGH  
**Estimated time:** 2-3 days

Hook up the diffusion term to the kinetic equation:
- Add term construction logic
- Register with UnknownQuantityHandler
- Ensure proper matrix assembly

**Files to modify:**
- `src/EquationSystem/UnknownQuantityHandler.cpp`
- `src/EqsysInitializer.cpp`

### Long-term Tasks (Month 2-3)

#### 7. Performance Optimization
**Priority:** MEDIUM  
**Estimated time:** 1-2 weeks

- Parallelize outer loops (OpenMP)
- Optimize polynomial root finding
- Cache frequently computed values
- Profile and optimize bottlenecks

#### 8. Validation and Testing
**Priority:** HIGH  
**Estimated time:** 2-3 weeks

- Compare with QUADRE results
- Test energy control scenarios
- Verify conservation properties
- Benchmark performance

#### 9. Documentation
**Priority:** LOW  
**Estimated time:** 1 week

- User guide for quasilinear diffusion
- Example scripts
- API documentation
- Theory-to-code mapping document

---

## 📊 Progress Summary

| Task | Status | Completion |
|------|--------|------------|
| **Phase 1: Infrastructure** | **IN PROGRESS** | **~30%** |
| - OptionConstants enum | ✅ Complete | 100% |
| - WaveSpectrum C++ class | ✅ Complete | 100% |
| - WaveSpectrum Python interface | ✅ Complete | 100% |
| - Build system update | ✅ Complete | 100% |
| - Basic testing | ✅ Complete | 100% |
| **Phase 2: Physics Models** | **PENDING** | **0%** |
| - WhistlerDispersion | ⏳ Not started | 0% |
| - ResonanceSolver | ⏳ Not started | 0% |
| - WaveParticleCoupling | ⏳ Not started | 0% |
| **Phase 3: Diffusion Operator** | **PENDING** | **0%** |
| - QuasilinearDiffusionTerm | ⏳ Not started | 0% |
| **Phase 4: Integration** | **PENDING** | **0%** |
| - Settings extension | ⏳ Not started | 0% |
| - Equation system integration | ⏳ Not started | 0% |
| **Phase 5: Testing & Validation** | **PENDING** | **0%** |

**Overall Progress:** ~5% complete

---

## 🔧 Technical Notes

### Design Decisions Made

1. **Wave Grid Strategy:** Use uniform grid in (k, θ_k) space to match QUADRE's approach
   - Allows flexible spectrum specification
   - Compatible with pre-allocated grid systems
   - Efficient numerical integration

2. **Class Architecture:** Separate concerns into distinct classes
   - WaveSpectrum: Manages spectrum parameters
   - WhistlerDispersion: Handles dispersion physics
   - ResonanceSolver: Finds resonance conditions
   - WaveParticleCoupling: Calculates interaction strength
   - QuasilinearDiffusionTerm: Assembles diffusion operator

3. **Python Interface:** Mirror C++ functionality for ease of use
   - Consistent API between languages
   - Easy parameter exploration
   - Seamless HDF5 output

### Challenges Identified

1. **Polynomial Root Finding:** 11th order polynomial may have numerical stability issues
   - Solution: Use robust library (Eigen/GSL) with proper filtering

2. **Performance:** Triple nested loops (radial × momentum × spectrum) will be slow
   - Solution: OpenMP parallelization, caching, vectorization

3. **Integration Complexity:** Need to hook into DREAM's equation system properly
   - Solution: Follow existing patterns (PitchScatterTerm, SynchrotronTerm)

---

## 📝 Files Created/Modified

### New Files (7)
1. `include/DREAM/Settings/OptionConstants.enum.hpp` (modified, +25 lines)
2. `include/DREAM/Settings/WaveSpectrum.hpp` (new, 147 lines)
3. `src/Settings/WaveSpectrum.cpp` (new, 247 lines)
4. `py/DREAM/Settings/WaveSpectrum.py` (new, 172 lines)
5. `src/CMakeLists.txt` (modified, +1 line)
6. `test_wave_spectrum.py` (new, 167 lines)
7. `IMPLEMENTATION_PROGRESS.md` (this file)

**Total lines added:** ~759 lines

### Files to Create Next (6)
1. `include/DREAM/Equations/Kinetic/WhistlerDispersion.hpp`
2. `src/Equations/Kinetic/WhistlerDispersion.cpp`
3. `include/DREAM/Equations/Kinetic/ResonanceSolver.hpp`
4. `src/Equations/Kinetic/ResonanceSolver.cpp`
5. `include/DREAM/Equations/Kinetic/WaveParticleCoupling.hpp`
6. `src/Equations/Kinetic/WaveParticleCoupling.cpp`

---

## 🎯 Milestones

- [x] **Milestone 1:** WaveSpectrum infrastructure (COMPLETED)
- [ ] **Milestone 2:** Physics model classes (Target: Week 2)
- [ ] **Milestone 3:** Diffusion operator implementation (Target: Week 3-4)
- [ ] **Milestone 4:** Full integration and testing (Target: Week 5-6)
- [ ] **Milestone 5:** Validation against QUADRE (Target: Week 7-8)

---

## 💡 Lessons Learned

1. **Start with infrastructure:** Building solid foundation makes later work easier
2. **Test early and often:** Catching bugs early saves time
3. **Follow existing patterns:** DREAM has well-established architecture
4. **Document decisions:** Future maintainers will thank you

---

**Last Updated:** 2026-05-01  
**Next Review:** After completing WhistlerDispersion class
