use extendr_api::prelude::*;
use clarabel::algebra::*;
use clarabel::solver::*;
use lazy_static::lazy_static;
use regex::Regex;

// Solve using Clarabel Rust solver.
#[extendr]
fn clarabel_solve(m: i32, n: i32, Ai: &[i32], Ap: &[i32], Ax: &[f64], b: &[f64], q: &[f64],
		  Pi: &[i32], Pp: &[i32], Px: &[f64], cone_spec: List, r_settings: List) -> Robj {
    #![allow(non_snake_case)]
    // QP Example
    let P = if Px.len() > 0 {
	let mut SPp = Vec::<usize>::with_capacity(Pp.len());
	for i in 0..Pp.len() {
	    SPp.push(Pp[i] as usize);
	}
	let mut SPi = Vec::<usize>::with_capacity(Pi.len());
	for i in 0..Pi.len() {
	    SPi.push(Pi[i] as usize);
	}
	
	CscMatrix::new(
            n as usize,             // m
            n as usize,             // n
            SPp,           // colptr
            SPi,           // rowval
            Px.to_vec(),   // nzval
	)
    } else {
	CscMatrix::zeros((n as usize, n as usize)) // For P = 0
    };

    // println!("P: {:?}", P);
    
    // assert!(P.check_format().is_ok());    

    let A = {
	let mut SAp = Vec::<usize>::with_capacity(Ap.len());
	for i in 0..Ap.len() {
	    SAp.push(Ap[i] as usize);
	}
	let mut SAi = Vec::<usize>::with_capacity(Ai.len());
	for i in 0..Ai.len() {
	    SAi.push(Ai[i] as usize);
	}

	CscMatrix::new(
            m as usize,             // m
            n as usize,             // n
            SAp,           // colptr
            SAi,           // rowval
            Ax.to_vec(),   // nzval
	)
    };

    // println!("A: {:?}", A);
    // assert!(A.check_format().is_ok());
    
    // println!("b: {:?}", b);
    // println!("q: {:?}", q);    
    
    // Handle cones
    lazy_static! {
	static ref ZC: Regex = Regex::new("^z").unwrap();  // ZeroCone
	static ref NNC: Regex = Regex::new("^l").unwrap(); // NonNegativeCone
	static ref SOC: Regex = Regex::new("^q").unwrap(); // Second Order Cone
	static ref EPC: Regex = Regex::new("^ep").unwrap(); // Exponential Cone
	static ref PC: Regex = Regex::new("^p").unwrap();  // Power Cone
	static ref PSDTC: Regex = Regex::new("^s").unwrap();  // PSD Triangle Cone
    }

    let mut cones: Vec::<SupportedConeT<f64>> = Vec::new();
    for (key, value) in cone_spec.iter() {
	if ZC.is_match(key.as_ref()) {
	    cones.push(ZeroConeT(value.as_integer().expect("Positive integer expected") as usize));
	} else if NNC.is_match(key.as_ref()) {
	    cones.push(NonnegativeConeT(value.as_integer().expect("Positive integer expected") as usize));	    
	} else if SOC.is_match(key.as_ref()) {
	    cones.push(SecondOrderConeT(value.as_integer().expect("Positive integer expected") as usize));
	} else if EPC.is_match(key.as_ref()) {
	    for _i in 0..value.as_integer().expect("Positive integer expected") {
		cones.push(ExponentialConeT());
	    }
	} else if PSDTC.is_match(key.as_ref()) {
	    cones.push(PSDTriangleConeT(value.as_integer().expect("Positive integer expected") as usize));	    
	} else if PC.is_match(key.as_ref()) {
	    cones.push(PowerConeT(value.as_real().expect("Positive real value expected")));	
        } else {
	    println!("Unknown cone {}; ignoring", &key);
	}
    }

    // println!("cones: {:?}", cones);    
    // Update default settings with specified R settings for use below
    let settings = update_settings(r_settings);
    
    let mut solver = DefaultSolver::new(&P, &q, &A, &b, &cones, settings);
    solver.solve();
    r!(list!(x = solver.solution.x.iter().collect_robj(),
	     z = solver.solution.z.iter().collect_robj(),
	     s = solver.solution.s.iter().collect_robj(),
	     obj_val = solver.solution.obj_val,
	     status = solver.solution.status as i32 + 1, // R's 1-based index
	     solve_time = solver.solution.solve_time,
	     iterations = solver.solution.iterations,
	     r_prim = solver.solution.r_prim,
	     r_dual = solver.solution.r_dual))
}

fn update_settings(r_settings: List) -> DefaultSettings<f64> {
    let mut settings = DefaultSettings::default();
    // Should really be using something like structmap (https://github.com/ex0dus-0x/structmap) but
    // unsuccessful so far, so a directly generated implementation
    for (key, value) in r_settings.iter() {
	match key.as_ref() {
	    "max_iter" => settings.max_iter = value.as_integer().expect("max_iter should be scalar integer") as u32,
	    "time_limit" => settings.time_limit = value.as_real().expect("should be scalar double"),
	    "verbose" => settings.verbose = value.as_bool().expect("verbose should be scalar boolean"),
	    "max_step_fraction" => settings.max_step_fraction = value.as_real().expect("should be scalar double"),
	    "tol_gap_abs" => settings.tol_gap_abs = value.as_real().expect("should be scalar double"),
	    "tol_gap_rel" => settings.tol_gap_rel = value.as_real().expect("should be scalar double"),
	    "tol_feas" => settings.tol_feas = value.as_real().expect("should be scalar double"),
	    "tol_infeas_abs" => settings.tol_infeas_abs = value.as_real().expect("should be scalar double"),
	    "tol_infeas_rel" => settings.tol_infeas_rel = value.as_real().expect("should be scalar double"),
	    "tol_ktratio" => settings.tol_ktratio = value.as_real().expect("should be scalar double"),
	    "reduced_tol_gap_abs" => settings.reduced_tol_gap_abs = value.as_real().expect("should be scalar double"),
	    "reduced_tol_gap_rel" => settings.reduced_tol_gap_rel = value.as_real().expect("should be scalar double"),
	    "reduced_tol_feas" => settings.reduced_tol_feas = value.as_real().expect("should be scalar double"),
	    "reduced_tol_infeas_abs" => settings.reduced_tol_infeas_abs = value.as_real().expect("should be scalar double"),
	    "reduced_tol_infeas_rel" => settings.reduced_tol_infeas_rel = value.as_real().expect("should be scalar double"),
	    "reduced_tol_ktratio" => settings.reduced_tol_ktratio = value.as_real().expect("should be scalar double"),
	    "equilibrate_enable" => settings.equilibrate_enable = value.as_bool().expect("equilibrate_enable should be scalar boolean"),
	    "equilibrate_max_iter" => settings.equilibrate_max_iter = value.as_integer().expect("Should be scalar integer") as u32,
	    "equilibrate_min_scaling" => settings.equilibrate_min_scaling = value.as_real().expect("should be scalar double"),
	    "equilibrate_max_scaling" => settings.equilibrate_max_scaling = value.as_real().expect("should be scalar double"),
	    "linesearch_backtrack_step" => settings.linesearch_backtrack_step = value.as_real().expect("should be scalar double"),
	    "min_switch_step_length" => settings.min_switch_step_length = value.as_real().expect("should be scalar double"),
	    "min_terminate_step_length" => settings.min_terminate_step_length = value.as_real().expect("should be scalar double"),
	    "direct_kkt_solver" => settings.direct_kkt_solver = value.as_bool().expect("direct_kkt_solver should be scalar boolean"),
	    "direct_solve_method" => settings.direct_solve_method = String::from(value.as_str().expect("direct_solve_method should be qdldl")),
	    "static_regularization_enable" => settings.static_regularization_enable = value.as_bool().expect("static_regularization_enable should be scalar boolean"),
	    "static_regularization_constant" => settings.static_regularization_constant = value.as_real().expect("should be scalar double"),
	    "static_regularization_proportional" => settings.static_regularization_proportional = value.as_real().expect("should be scalar double"),
	    "dynamic_regularization_enable" => settings.dynamic_regularization_enable = value.as_bool().expect("dynamic_regularization_enable should be scalar boolean"),
	    "dynamic_regularization_eps" => settings.dynamic_regularization_eps = value.as_real().expect("should be scalar double"),
	    "dynamic_regularization_delta" => settings.dynamic_regularization_delta = value.as_real().expect("should be scalar double"),
	    "iterative_refinement_enable" => settings.iterative_refinement_enable = value.as_bool().expect("iterative_refinement_enable should be scalar boolean"),
	    "iterative_refinement_reltol" => settings.iterative_refinement_reltol = value.as_real().expect("should be scalar double"),
	    "iterative_refinement_abstol" => settings.iterative_refinement_abstol = value.as_real().expect("should be scalar double"),
	    "iterative_refinement_max_iter" => settings.iterative_refinement_max_iter = value.as_integer().expect("Should be scalar integer") as u32,
	    "iterative_refinement_stop_ratio" => settings.iterative_refinement_stop_ratio = value.as_real().expect("should be scalar double"),
	    "presolve_enable" => settings.presolve_enable = value.as_bool().expect("presolve_enable should be scalar boolean"),
	     _ => (),
	    
	}
    }
    settings
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod clarabel;
    fn clarabel_solve;
}
