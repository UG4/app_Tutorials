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

--	Create a Mesh first. All promesh-methods operate on a Mesh
obj = Mesh()


--	This method creates a circle with center (0, 0, 0), radius 2,
--	64 rim-vertices, subset-index 0 and fills the interior with triangles.
CreateCircle(obj, MakeVec(0, 0, 0), 2, 64, 0, true)


--	Retriangulation leads to a minimal inner angle of 30°. Smoothing
--	will move the vertices a little so that the grid appears 'smoother'.
SelectAll(obj)
Retriangulate(obj, 30)
ClearSelection(obj)
SelectInnerVertices(obj)
LaplacianSmooth(obj, 0.1, 10)


--	Now we'll assign subsets so that the resulting grid is ready for simulation
--	with ug4.
SelectAll(obj)
AssignSubset(obj, 0)
SetSubsetName(obj, 0, "inner")
ClearSelection(obj)
SelectBoundaryEdges(obj)
SelectAssociatedVertices(obj)
AssignSubset(obj, 1)
SetSubsetName(obj, 1, "boundary")
AssignSubsetColors(obj)


--	Finally save the resulting mesh to a file
SaveMesh(obj, "tut11-circle.ugx")
