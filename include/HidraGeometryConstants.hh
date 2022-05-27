#ifndef HIDRA_GEOM_CONSTANTS
#define HIDRA_GEOM_CONSTANTS

// Inlcudes G4Types, G4PhysicalConstants and G4SystemOfUnits
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <array>
#include <math.h>
#include <stdint.h>

class HidraGeometryConstants {
public:
  static HidraGeometryConstants *GetInstance() {
    // In C++11 static function members are thread safe
    // Read access to this class is thread safe, write is undefined behaivor!
    static HidraGeometryConstants *instance = new HidraGeometryConstants();
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
  
  // TODO: avoid hard coded "24*5". It should be m_NoModulesX * m_NoModulesY
  // but those variables are defined as private members later due to ODR
  constexpr std::array<int32_t, 24*5> moduleFlag() const { 
    return m_ModuleFlag;
  }
  constexpr uint8_t nActiveModules() const { return m_NoActiveModules; }
  constexpr uint8_t nSiPMModules() const { return m_NoSiPMModules; }
  // TODO: avoid hard coded "10". It shoul be m_NoSiPMModules
  // but this variable is defined as private member later due to ODR
  constexpr std::array<int32_t, 10> idxSiPMModules() const { 
    return m_SiPMModuleIdx;
  }

private:
  static constexpr uint8_t m_NoModulesX = 24;
  static constexpr uint8_t m_NoModulesY = 5;

  static constexpr double m_TubeOuterRadius = 1.0 * mm;
  static constexpr double m_TubeInnerRadius = 0.55 * mm;
  static constexpr double m_TubeLength = 2500.0 * mm;

  static constexpr double m_FiberCladdingInnerRadius = 0.4925 * mm;
  static constexpr double m_FiberCladdingOuterRadius = 0.5 * mm;
  static constexpr double m_FiberCoregOuterRadius = m_FiberCladdingInnerRadius;

  static constexpr uint8_t m_NoTubesInModuleOnX = 64;
  static constexpr uint8_t m_NoTubesInModuleOnY = 16;

  static constexpr double m_SiPMCoatingSizeXY = 1.3 * mm;
  static constexpr double m_SiPMSiliconSizeXY = 1.0 * mm;
  static constexpr double m_SIPMCoatingSizeZ = 500 * um;
  static constexpr double m_SIPMSiliconSizeZ = 300 * um;

  static constexpr double m_ModuleSizeZ = m_TubeLength + m_SIPMCoatingSizeZ;

  static constexpr double m_RotationVertical = 0.0 * deg;
  static constexpr double m_RotationHorizontal = 0.0 * deg;

  static constexpr std::array<int32_t, m_NoModulesX * m_NoModulesY>
      m_ModuleFlag = {-1, -1, -1, -1, -1, -1, -1, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  -1, -1, -1, -1, -1, -1, -1,
                      -1, -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, -1, -1,
                      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
                      -1, -1, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, -1, -1,
                      -1, -1, -1, -1, -1, -1, -1, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, -1, -1, -1, -1, -1, -1, -1};
  static constexpr uint8_t m_NoActiveModules = 84;
  static constexpr uint8_t m_NoSiPMModules = 10;
  static constexpr std::array<int32_t, m_NoSiPMModules> m_SiPMModuleIdx = {
      37, 38, 39, 40, 41, 42, 43, 44, 45, 46};
  // A singleton class cannot be copied moved or assigned
  // Delete all "rule of 5" methods and declare ctor as private
  //
  // Private ctor
  constexpr HidraGeometryConstants() = default;
  // Deleted copy ctor
  HidraGeometryConstants(const HidraGeometryConstants &) = delete;
  // Deleted copy assign
  HidraGeometryConstants &operator=(const HidraGeometryConstants &) = delete;
  // Deleted move ctor
  HidraGeometryConstants(HidraGeometryConstants &&) = delete;
  // Deleted move assign
  HidraGeometryConstants &operator=(HidraGeometryConstants &&) = delete;
};

#endif // !HIDRA_GEOM_CONSTANTS
