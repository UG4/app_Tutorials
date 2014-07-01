--------------------------------------------------------------------------------
--	tut11_geometry_generation.lua
--
--	This tutorial requires the 'promesh' plugin. If you're building ug from
--	source, make sure to call cmake with -DProMesh=ON
--
--	The tutorial will show how a circular geometry filled with nicely shaped
--	triangles can be generated. Also subset assignment is performed.
--	The resulting geometry is then ready for simulation.
--------------------------------------------------------------------------------

--	Create a MeshObject first. All promesh-methods operate on a MeshObject
obj = PM_MeshObject()


--	This method creates a circle with center (0, 0, 0), radius 2,
--	64 rim-vertices, subset-index 0 and fills the interior with triangles.
PM_CreateCircle(obj, MakeVec(0, 0, 0), 2, 64, 0, true)


--	Retriangulation leads to a minimal inner angle of 30°. Smoothing
--	will move the vertices a little so that the grid appears 'smoother'.
PM_SelectAll(obj)
PM_Retriangulate(obj, 30)
PM_ClearSelection(obj)
PM_SelectInnerVertices(obj)
PM_LaplacianSmooth(obj, 0.1, 10)


--	Now we'll assign subsets so that the resulting grid is ready for simulation
--	with ug4.
PM_SelectAll(obj)
PM_AssignSubset(obj, 0)
PM_SetSubsetName(obj, 0, "inner")
PM_ClearSelection(obj)
PM_SelectBoundaryEdges(obj)
PM_SelectAssociatedVertices(obj)
PM_AssignSubset(obj, 1)
PM_SetSubsetName(obj, 1, "boundary")
PM_AssignSubsetColors(obj)


--	Finally save the resulting mesh to a file
PM_SaveMesh(obj, "tut11-circle.ugx")
