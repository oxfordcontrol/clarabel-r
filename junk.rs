"max_iter" => 
   match typed_value {
       TypedSexp::Integer(i) => settings.max_iter = i.as_slice()[0],
_ => settings.max_iter = settings.max_iter,
},
"time_limit" => 
   match typed_value {
       TypedSexp::Real(f) => settings.time_limit = f.as_slice()[0],
_ => settings.time_limit = settings.time_limit,
},
"verbose" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.verbose = b.as_slice_raw()[0] != 0,
_ => settings.verbose = settings.verbose,
},
"max_step_fraction" => 
   match typed_value {
       TypedSexp::Real(f) => settings.max_step_fraction = f.as_slice()[0],
_ => settings.max_step_fraction = settings.max_step_fraction,
},
"tol_gap_abs" => 
   match typed_value {
       TypedSexp::Real(f) => settings.tol_gap_abs = f.as_slice()[0],
_ => settings.tol_gap_abs = settings.tol_gap_abs,
},
"tol_gap_rel" => 
   match typed_value {
       TypedSexp::Real(f) => settings.tol_gap_rel = f.as_slice()[0],
_ => settings.tol_gap_rel = settings.tol_gap_rel,
},
"tol_feas" => 
   match typed_value {
       TypedSexp::Real(f) => settings.tol_feas = f.as_slice()[0],
_ => settings.tol_feas = settings.tol_feas,
},
"tol_infeas_abs" => 
   match typed_value {
       TypedSexp::Real(f) => settings.tol_infeas_abs = f.as_slice()[0],
_ => settings.tol_infeas_abs = settings.tol_infeas_abs,
},
"tol_infeas_rel" => 
   match typed_value {
       TypedSexp::Real(f) => settings.tol_infeas_rel = f.as_slice()[0],
_ => settings.tol_infeas_rel = settings.tol_infeas_rel,
},
"tol_ktratio" => 
   match typed_value {
       TypedSexp::Real(f) => settings.tol_ktratio = f.as_slice()[0],
_ => settings.tol_ktratio = settings.tol_ktratio,
},
"reduced_tol_gap_abs" => 
   match typed_value {
       TypedSexp::Real(f) => settings.reduced_tol_gap_abs = f.as_slice()[0],
_ => settings.reduced_tol_gap_abs = settings.reduced_tol_gap_abs,
},
"reduced_tol_gap_rel" => 
   match typed_value {
       TypedSexp::Real(f) => settings.reduced_tol_gap_rel = f.as_slice()[0],
_ => settings.reduced_tol_gap_rel = settings.reduced_tol_gap_rel,
},
"reduced_tol_feas" => 
   match typed_value {
       TypedSexp::Real(f) => settings.reduced_tol_feas = f.as_slice()[0],
_ => settings.reduced_tol_feas = settings.reduced_tol_feas,
},
"reduced_tol_infeas_abs" => 
   match typed_value {
       TypedSexp::Real(f) => settings.reduced_tol_infeas_abs = f.as_slice()[0],
_ => settings.reduced_tol_infeas_abs = settings.reduced_tol_infeas_abs,
},
"reduced_tol_infeas_rel" => 
   match typed_value {
       TypedSexp::Real(f) => settings.reduced_tol_infeas_rel = f.as_slice()[0],
_ => settings.reduced_tol_infeas_rel = settings.reduced_tol_infeas_rel,
},
"reduced_tol_ktratio" => 
   match typed_value {
       TypedSexp::Real(f) => settings.reduced_tol_ktratio = f.as_slice()[0],
_ => settings.reduced_tol_ktratio = settings.reduced_tol_ktratio,
},
"equilibrate_enable" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.equilibrate_enable = b.as_slice_raw()[0] != 0,
_ => settings.equilibrate_enable = settings.equilibrate_enable,
},
"equilibrate_max_iter" => 
   match typed_value {
       TypedSexp::Integer(i) => settings.equilibrate_max_iter = i.as_slice()[0],
_ => settings.equilibrate_max_iter = settings.equilibrate_max_iter,
},
"equilibrate_min_scaling" => 
   match typed_value {
       TypedSexp::Real(f) => settings.equilibrate_min_scaling = f.as_slice()[0],
_ => settings.equilibrate_min_scaling = settings.equilibrate_min_scaling,
},
"equilibrate_max_scaling" => 
   match typed_value {
       TypedSexp::Real(f) => settings.equilibrate_max_scaling = f.as_slice()[0],
_ => settings.equilibrate_max_scaling = settings.equilibrate_max_scaling,
},
"linesearch_backtrack_step" => 
   match typed_value {
       TypedSexp::Real(f) => settings.linesearch_backtrack_step = f.as_slice()[0],
_ => settings.linesearch_backtrack_step = settings.linesearch_backtrack_step,
},
"min_switch_step_length" => 
   match typed_value {
       TypedSexp::Real(f) => settings.min_switch_step_length = f.as_slice()[0],
_ => settings.min_switch_step_length = settings.min_switch_step_length,
},
"min_terminate_step_length" => 
   match typed_value {
       TypedSexp::Real(f) => settings.min_terminate_step_length = f.as_slice()[0],
_ => settings.min_terminate_step_length = settings.min_terminate_step_length,
},
"direct_kkt_solver" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.direct_kkt_solver = b.as_slice_raw()[0] != 0,
_ => settings.direct_kkt_solver = settings.direct_kkt_solver,
},
"direct_solve_method" => 
   match typed_value {
       TypedSexp::String(s) => settings.direct_solve_method = String::from(s.to_vec().get(0)),
_ => settings.direct_solve_method = settings.direct_solve_method,
},
"static_regularization_enable" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.static_regularization_enable = b.as_slice_raw()[0] != 0,
_ => settings.static_regularization_enable = settings.static_regularization_enable,
},
"static_regularization_constant" => 
   match typed_value {
       TypedSexp::Real(f) => settings.static_regularization_constant = f.as_slice()[0],
_ => settings.static_regularization_constant = settings.static_regularization_constant,
},
"static_regularization_proportional" => 
   match typed_value {
       TypedSexp::Real(f) => settings.static_regularization_proportional = f.as_slice()[0],
_ => settings.static_regularization_proportional = settings.static_regularization_proportional,
},
"dynamic_regularization_enable" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.dynamic_regularization_enable = b.as_slice_raw()[0] != 0,
_ => settings.dynamic_regularization_enable = settings.dynamic_regularization_enable,
},
"dynamic_regularization_eps" => 
   match typed_value {
       TypedSexp::Real(f) => settings.dynamic_regularization_eps = f.as_slice()[0],
_ => settings.dynamic_regularization_eps = settings.dynamic_regularization_eps,
},
"dynamic_regularization_delta" => 
   match typed_value {
       TypedSexp::Real(f) => settings.dynamic_regularization_delta = f.as_slice()[0],
_ => settings.dynamic_regularization_delta = settings.dynamic_regularization_delta,
},
"iterative_refinement_enable" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.iterative_refinement_enable = b.as_slice_raw()[0] != 0,
_ => settings.iterative_refinement_enable = settings.iterative_refinement_enable,
},
"iterative_refinement_reltol" => 
   match typed_value {
       TypedSexp::Real(f) => settings.iterative_refinement_reltol = f.as_slice()[0],
_ => settings.iterative_refinement_reltol = settings.iterative_refinement_reltol,
},
"iterative_refinement_abstol" => 
   match typed_value {
       TypedSexp::Real(f) => settings.iterative_refinement_abstol = f.as_slice()[0],
_ => settings.iterative_refinement_abstol = settings.iterative_refinement_abstol,
},
"iterative_refinement_max_iter" => 
   match typed_value {
       TypedSexp::Integer(i) => settings.iterative_refinement_max_iter = i.as_slice()[0],
_ => settings.iterative_refinement_max_iter = settings.iterative_refinement_max_iter,
},
"iterative_refinement_stop_ratio" => 
   match typed_value {
       TypedSexp::Real(f) => settings.iterative_refinement_stop_ratio = f.as_slice()[0],
_ => settings.iterative_refinement_stop_ratio = settings.iterative_refinement_stop_ratio,
},
"presolve_enable" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.presolve_enable = b.as_slice_raw()[0] != 0,
_ => settings.presolve_enable = settings.presolve_enable,
},
"chordal_decomposition_enable" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.chordal_decomposition_enable = b.as_slice_raw()[0] != 0,
_ => settings.chordal_decomposition_enable = settings.chordal_decomposition_enable,
},
"chordal_decomposition_merge_method" => 
   match typed_value {
       TypedSexp::String(s) => settings.chordal_decomposition_merge_method = String::from(s.to_vec().get(0)),
_ => settings.chordal_decomposition_merge_method = settings.chordal_decomposition_merge_method,
},
"chordal_decomposition_compact" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.chordal_decomposition_compact = b.as_slice_raw()[0] != 0,
_ => settings.chordal_decomposition_compact = settings.chordal_decomposition_compact,
},
"chordal_decomposition_complete_dual" => 
   match typed_value {
       TypedSexp::Logical(b) => settings.chordal_decomposition_complete_dual = b.as_slice_raw()[0] != 0,
_ => settings.chordal_decomposition_complete_dual = settings.chordal_decomposition_complete_dual,
},
