include(MiniFEMAddTest)

# general
minifem_add_test(hdf5write)
minifem_add_test(io)
minifem_add_test(iodisc)
minifem_add_test(iotime)
minifem_add_test(timeder)


# geometry
minifem_add_test(mesh_refine)
minifem_add_test(normal)


# finite elements
minifem_add_test(assemblyvector)
minifem_add_test(fevar)
minifem_add_test(inside_test)
minifem_add_test(inverse_mapping)
minifem_add_test(qr)
minifem_add_test(reffe_test)
minifem_add_test(vectorfe1d)


# linear algebra
minifem_add_test(blockmatrix)
minifem_add_test(prolongation)
minifem_add_test(solver_test)


# bc
minifem_add_test(fourier_nb)
minifem_add_test(neumann)
minifem_add_test(robin_bc)
minifem_add_test(periodic)
minifem_add_test(projection_p1top2)
minifem_add_test(projection_p2top1)
minifem_add_test(projection_rt0)


# pde 1d
minifem_add_test(advection1d DATA advection1d.yaml)
minifem_add_test(poisson1d)
minifem_add_test(poisson1d_p2)
minifem_add_test(supg1d)


# pde 2d tri
minifem_add_test(advection2dtri)
minifem_add_test(dual)
minifem_add_test(poisson2dcr)
minifem_add_test(poisson2dmixedtri)
minifem_add_test(poisson2dtri)
minifem_add_test(poisson2dtri_p2)


# pde 2d quad
minifem_add_test(advection2dquad DATA square_q.msh)
minifem_add_test(dr2dquad)
minifem_add_test(poisson2dmixedquad)
minifem_add_test(poisson2dquad)
minifem_add_test(poisson2dquad_q2)


# pde 3d tet
minifem_add_test(poisson3dtet DATA poisson3dtet.yaml cube_uns.msh)
minifem_add_test(poisson3dtet_p2)


# pde 3d hex
minifem_add_test(poisson3dhex)
minifem_add_test(poisson3dhex_q2)


# stokes monolithic 2d
minifem_add_test(cavity)
minifem_add_test(cavitytime DATA cavitytime.yaml)
minifem_add_test(ns_test)
minifem_add_test(stokes2dquad)
minifem_add_test(slip)
minifem_add_test(steptime DATA steptime.yaml)
minifem_add_test(stokes2dtri DATA square_uns.msh)
minifem_add_test(stokes2dcircle)

# stokes split 2d
minifem_add_test(ipcs)
minifem_add_test(ipcs_rt)
minifem_add_test(nipcs)


# stokes monolithic 3d
minifem_add_test(stokes3dtet)
minifem_add_test(stokes3dhex)
