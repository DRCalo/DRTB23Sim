
#ifndef TB_GEOM_CONSTANTS
#define TB_GEOM_CONSTANTS

// Inlcudes G4Types, G4PhysicalConstants and G4SystemOfUnits
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <array>
#include <cstdint>
#include <math.h>
#include <stdint.h>

class TestBeamGeometryConstants {
  private:
    static constexpr uint8_t m_NoModulesX = 3;
    static constexpr uint8_t m_NoModulesY = 3;
    static constexpr uint8_t m_NoSiPMModules = 1;
public:
  static TestBeamGeometryConstants *GetInstance() {
    // In C++11 static function members are thread safe
    // Read access to this class is thread safe, write is undefined behaivor!
    static TestBeamGeometryConstants *instance = new TestBeamGeometryConstants();
    return instance;
  }

  constexpr uint8_t nModulesOnX() const { return m_NoModulesX; }
  constexpr uint8_t nModulesOnY() const { return m_NoModulesY; }
  constexpr uint16_t nModules() const { return m_NoModulesY * m_NoModulesX; }

  constexpr double tubeOuterRadius() const { return m_TubeOuterRadius; }
  constexpr double tubeInnerRadius() const { return m_TubeInnerRadius; }
  constexpr double tubeOuterDiameter() const { return 2 * m_TubeOuterRadius; }
  constexpr double tubeInnerDiameter() const { return 2 * m_TubeInnerRadius; }
  constexpr double tubeLength() const { return m_TubeLength; }

  constexpr double fiberCladdingOuterRadius() const {
    return m_FiberCladdingOuterRadius;
  }
  constexpr double fiberCladdingInnerRadius() const {
    return m_FiberCladdingInnerRadius;
  }
  constexpr double fiberCladdingOuterDiameter() const {
    return 2 * m_FiberCladdingOuterRadius;
  }
  constexpr double fiberCladdingInnerDiameter() const {
    return 2 * m_FiberCladdingInnerRadius;
  }
  constexpr double fiberCoreOuterRadius() const {
    return m_FiberCoregOuterRadius;
  }
  constexpr double fiberCoreOuterDiameter() const {
    return 2 * m_FiberCoregOuterRadius;
  }

  constexpr uint8_t nTubesInModuleX() const { return m_NoTubesInModuleOnX; }
  constexpr uint8_t nTubesInModuleY() const { return m_NoTubesInModuleOnY; }

  constexpr double sipmSizeX() const { return m_SiPMCoatingSizeXY; }
  constexpr double sipmSizeY() const { return m_SiPMCoatingSizeXY; }
  constexpr double sipmSizeZ() const { return m_SIPMCoatingSizeZ; }
  constexpr double sipmSiliconSizeX() const { return m_SiPMSiliconSizeXY; }
  constexpr double sipmSiliconSizeY() const { return m_SiPMSiliconSizeXY; }
  constexpr double sipmSiliconSizeZ() const { return m_SIPMSiliconSizeZ; }

  constexpr double moduleSizeX() const {
    return 2 * m_TubeOuterRadius * m_NoTubesInModuleOnX;
  }
  constexpr double moduleSizeY() const {
    return m_TubeOuterRadius * m_NoTubesInModuleOnY * std::sqrt(3);
  }
  constexpr double moduleSizeZ() const {
    return m_TubeLength + m_SIPMCoatingSizeZ;
  }

  constexpr double calorimeterSizeX() const {
    return m_NoModulesX * moduleSizeX();
  }
  constexpr double calorimeterSizeY() const {
    return m_NoModulesY * moduleSizeY();
  }
  constexpr double calorimeterSizeZ() const { return moduleSizeZ(); }

  constexpr double rotationV() const { return m_RotationVertical; }
  constexpr double rotationH() const { return m_RotationHorizontal; }
  
  constexpr std::array<int32_t, m_NoModulesX * m_NoModulesY> moduleFlag() const { 
    return m_ModuleFlag;
  }
  constexpr uint8_t nActiveModules() const { return m_NoActiveModules; }
  constexpr uint8_t nSiPMModules() const { return m_NoSiPMModules; }
  constexpr std::array<int32_t, m_NoSiPMModules> idxSiPMModules() const { 
    return m_SiPMModuleIdx;
  }

private:
  static constexpr double m_TubeOuterRadius = 1.0 * mm;
  static constexpr double m_TubeInnerRadius = 0.55 * mm;
  static constexpr double m_TubeLength = 1000.0 * mm;

  static constexpr double m_FiberCladdingInnerRadius = 0.4925 * mm;
  static constexpr double m_FiberCladdingOuterRadius = 0.5 * mm;
  static constexpr double m_FiberCoregOuterRadius = m_FiberCladdingInnerRadius;

  static constexpr uint8_t m_NoTubesInModuleOnX = 16;
  static constexpr uint8_t m_NoTubesInModuleOnY = 20;

  static constexpr double m_SiPMCoatingSizeXY = 1.3 * mm;
  static constexpr double m_SiPMSiliconSizeXY = 1.0 * mm;
  static constexpr double m_SIPMCoatingSizeZ = 500 * um;
  static constexpr double m_SIPMSiliconSizeZ = 300 * um;

  static constexpr double m_ModuleSizeZ = m_TubeLength + m_SIPMCoatingSizeZ;

  static constexpr double m_RotationVertical = 0.0 * deg;
  static constexpr double m_RotationHorizontal = 0.0 * deg;

  static constexpr std::array<int32_t, m_NoModulesX * m_NoModulesY> m_ModuleFlag = {3,2,1,5,0,4,8,7,6};
  static constexpr uint8_t m_NoActiveModules = 9;
  static constexpr std::array<int32_t, m_NoSiPMModules> m_SiPMModuleIdx = {0};
  // A singleton class cannot be copied moved or assigned
  // Delete all "rule of 5" methods and declare ctor as private
  //
  // Private ctor
  constexpr TestBeamGeometryConstants() = default;
  // Deleted copy ctor
  TestBeamGeometryConstants(const TestBeamGeometryConstants &) = delete;
  // Deleted copy assign
  TestBeamGeometryConstants &operator=(const TestBeamGeometryConstants &) = delete;
  // Deleted move ctor
  TestBeamGeometryConstants(TestBeamGeometryConstants &&) = delete;
  // Deleted move assign
  TestBeamGeometryConstants &operator=(TestBeamGeometryConstants &&) = delete;
};

#endif // !TB_GEOM_CONSTANTS
