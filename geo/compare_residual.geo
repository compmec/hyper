// Given mesh files (reference solution)
//Merge "../msh/triangle-tri56-ref.msh";
Merge "../msh/triangle-quad44-ref.msh";
/*
Solution file first:
view[0] = ref nodal displacement field (U)
view[1] = ref nodal residual field (R)
view[2] = ref element def. gradient field (F)
view[3] = ref element eng. stress field (P)
*/

// Your mesh files
//Merge "../msh/triangle-tri56-val.msh";
Merge "../msh/triangle-quad44-val.msh";
/*
Your solution second:
view[4] = nodal displacement field (U)
view[5] = nodal residual field (R)
view[6] = element def. gradient field (F)
view[7] = element eng. stress field (P)
*/

// Apply the view as the current background mesh
Background Mesh View[0];
For numview In {0:PostProcessing.NbViews-1}
  View[numview].Visible = 0;
EndFor
// Compare views
// vector fields : 2D vectors are stored in v0 and v1 components of the view (vi = 0, i>=2)
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).OtherTimeStep = -1; //If `TimeStep' < 0, the plugin extracts data from all the time steps in the view.
Plugin(MathEval).ForceInterpolation = 0;
Plugin(MathEval).PhysicalRegion = -1; //If `PhysicalRegion' < 0, the plugin is run on all physical regions.
Plugin(MathEval).View = 1;
Plugin(MathEval).OtherView = 5;
Plugin(MathEval).Expression0 = " Fabs(v0 - w0)";
Plugin(MathEval).Expression1 = " Fabs(v1 - w1)";
Plugin(MathEval).Run;
View[PostProcessing.NbViews-1].Name = "Absolute error on the residual field";
// tensor fields : 2D tensors are stored in v0, v1, v3 and v4 components of the view (_11, _12, _21, _22 components respectively)
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).OtherTimeStep = -1; //If `TimeStep' < 0, the plugin extracts data from all the time steps in the view.
Plugin(MathEval).ForceInterpolation = 0;
Plugin(MathEval).PhysicalRegion = -1; //If `PhysicalRegion' < 0, the plugin is run on all physical regions.
Plugin(MathEval).View = 2;
Plugin(MathEval).OtherView = 6;
Plugin(MathEval).Expression0 = " Fabs(v0 - w0)";
Plugin(MathEval).Expression1 = " Fabs(v1 - w1)";
Plugin(MathEval).Expression3 = " Fabs(v3 - w3)";
Plugin(MathEval).Expression4 = " Fabs(v4 - w4)";
Plugin(MathEval).Run;
View[PostProcessing.NbViews-1].Name = "Absolute error on the def. gradient field";
// tensor fields : 2D tensors are stored in v0, v1, v3 and v4 components of the view (_11, _12, _21, _22 components respectively)
Plugin(MathEval).TimeStep = -1;
Plugin(MathEval).OtherTimeStep = -1; //If `TimeStep' < 0, the plugin extracts data from all the time steps in the view.
Plugin(MathEval).ForceInterpolation = 0;
Plugin(MathEval).PhysicalRegion = -1; //If `PhysicalRegion' < 0, the plugin is run on all physical regions.
Plugin(MathEval).View = 3;
Plugin(MathEval).OtherView = 7;
Plugin(MathEval).Expression0 = " Fabs(v0 - w0)";
Plugin(MathEval).Expression1 = " Fabs(v1 - w1)";
Plugin(MathEval).Expression3 = " Fabs(v3 - w3)";
Plugin(MathEval).Expression4 = " Fabs(v4 - w4)";
Plugin(MathEval).Run;
View[PostProcessing.NbViews-1].Name = "Absolute error on the eng. stress field";
