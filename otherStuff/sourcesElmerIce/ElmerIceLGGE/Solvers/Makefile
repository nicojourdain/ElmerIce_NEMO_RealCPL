SOLVERSOURCE := \
Solvers/ComputeDevStressSSA.F90\
Solvers/Flotation.F90\
Solvers/SSABulkVelocities.F90\
Solvers/SSAStar.F90\
Solvers/SSAStrainRate.F90\
Solvers/ComputeChi.F90 \
Solvers/CostSolver_Robin_Stresses.F90\
Solvers/DensitySolve.F90\
Solvers/DJDE_Adjoint.F90\
Solvers/DJDE_Adjoint_Z.F90\
Solvers/DJDE_Robin.F90\
Solvers/DJDE_Robin_Z.F90\
Solvers/DJDmu_Adjoint_Z.F90\
Solvers/DJDmu_Robin_Z.F90\
Solvers/FrontDisplacement.F90\
Solvers/GetArea.F90\
Solvers/GetLateralFrictionForce.F90\
Solvers/GroundedSolver_withSave.F90\
Solvers/IntegrateViscositySSA.F90\
Solvers/InterpBedrock.F90\
Solvers/InterpolateMeshToMesh.F90\
Solvers/OutPutCavity.F90\
Solvers/SeaDamagePressure.F90\
Solvers/SIAFlux.F90\
Solvers/SmoothFS.F90\
Solvers/SSABasalMelt.F90\
Solvers/ThicknessSolverDG.F90\
Solvers/MISMIP_Melt_Param.F90\
Solvers/MISMIP_Melt_Param_ForIsmip6.F90\
Solvers/Optimize_m1qn3Parallel.Multi.F90\
Solvers/EffectivePressure.F90
#Solvers/ToLonLat.F90
#### Exclus à cause de pb à la compilation ou executionn
#Probleme compilation :
#-- Solvers/CalvingSolver.F90\ -> L.98 LoadMesh No implicit type
#-- Solvers/CELSolver.F90:741:41: Pas a jour pour l'appel à DiffuseConvectiveBoundary
#-- Solvers/HydroSolver.F90:1289:40: Pas a jour pour l'appel à DiffuseConvectiveBoundary
#-- Solvers/SedSolver.F90:706:41:Pas a jour pour l'appel à DiffuseConvectiveBoundary
#-- Solvers/TwoMeshes.F90:108:14 -> L. 108 LoadMesh No implicit type
#Probleme execution :
#-- Solvers/AIFlowSolve_nlD2.F90\ r2ro_
#-- Solvers/AIFlowSolve_nlS2.F90\ r2ro_
#-- Solvers/ComputeDevStress_General.F90\ !!Necessite  GolfLaw.F90 dans CSC_SVN/elmerice/Solvers
#-- Solvers/LumiSolve.F90\ !! Pb fns parameterA and B not included in the librairy
#-- Solvers/GradTransFSolve.F90\ !! Manque fn SocleB
#-- Solvers/hRefinementSolver.F90 !!cputime_
#-- Solvers/Optimize_m1qn3Parallel_Z.F90 !!ctonbe_
#-- Solvers/Optimize_m1qn3Serial.F90
#-- Solvers/Optimize_m1qn3Serial_Z.F90
