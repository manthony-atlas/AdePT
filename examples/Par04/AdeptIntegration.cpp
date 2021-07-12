// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#include "AdeptIntegration.h"
#include "AdeptIntegration.cuh"

#include "G4RunManager.hh"

#include "VecGeom/management/GeoManager.h"

void AdeptIntegration::AddTrack(int tid, int pdg, double energy, double x, double y, double z, double dirx, double diry,
                                double dirz)
{
  fBuffer[tid].toDevice.emplace_back(pdg, energy, x, y, z, dirx, diry, dirz);
  if (pdg == 11)
    fBuffer[tid].nelectrons++;
  else if (pdg == -11)
    fBuffer[tid].npositrons++;
  else if (pdg == 22)
    fBuffer[tid].ngammas++;
}

void AdeptIntegration::Initialize()
{
  if (fInit) return;
  assert(fMaxBatch > 0 && "AdeptImtegration::Initialize - Maximum batch size not set.");

  assert(vecgeom::GeoManager::Instance().IsClosed() && "VecGeom geometry not closed.");
  const vecgeom::cxx::VPlacedVolume *world = vecgeom::GeoManager::Instance().GetWorld();
  fNthreads                                = G4RunManager::GetRunManager()->GetNumberOfThreads();
  std::cout << "=== AdeptIntegration: Number of threads: " << fNthreads << std::endl;

  AdeptIntegration::InitializeGPU(world, fMaxBatch);
  fInit = true;
}

void AdeptIntegration::Cleanup()
{
  if (!fCleanup.test_and_set()) {
    std::cout << "=== Cleaning up GPU ...\n";
    AdeptIntegration::FreeGPU();
  }
}

void AdeptIntegration::Shower(int event, int tid)
{
  AdeptIntegration::ShowerGPU(event, tid, fBuffer[tid]);
  fBuffer[tid].Clear();
}
