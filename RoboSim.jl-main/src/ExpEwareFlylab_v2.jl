module EwareExpDIFL

"""
Eware implementation for the the experiments in the Flylab
- Robot Model: Double integerator model
- Planer/High-level controller: Convex MPC
- Controller: Double integerator trajectories (x,y,z, dx, dy, dz, u1. u2, u3) are tracked by the geometric controller
Author: Kaleb Ben Naveed
"""

using DifferentialEquations
using LinearAlgebra
using BlockDiagonals
using ForwardDiff
using ControlSystems
import Convex as cvx
using ECOS
import DifferentialEquations: solve as DEsolve
import DifferentialEquations: ODEProblem, Tsit5, PeriodicCallback, CallbackSet
using ComponentArrays
using Parameters:@unpack

include("ExpDynamicsLibrary.jl")
include("ExpControllerLibrary.jl")
include("RefTrajectoryLibrary.jl")
include("utils.jl")


NomTrajectory = ExpController.NomTrajectory


"""
I have two days to finish this within 2 days:
1. Generate the reference trajectory using the convex mpc
2. Use differential flatness property to convert the (pos, vel, acc) to quadrotor reference states
"""

# For Multiple quads, this quad defines whats associated with each quad
mutable struct QuadObj{F}
    id::Int64                                 # robot number
    u::Vector{F}                              # Control input
    nom_traj::NomTrajectory{Float64}          # Nominal trajectory
    can_traj::NomTrajectory{Float64}          # Candidate trajectory
    com_traj::NomTrajectory{Float64}          # Committed trajectory
    cand_traj_time::Vector{F}                 # candidate trajectory end_time
    can_soc_end::Vector{F}                    # Soc at the end of the candidate trajectory
    charge_pos::Vector{F}                     # Charging position
    return_pos::Vector{F}                      # The return spot before quad lands
    init_soc::Float64                         # Initial soc of the robot
    discharge_rate::Float64                   # Discharge rate of the robot
    charge_start_time::Float64                # Charge start of the robot
    T_char::Float64                           # charge_time for the robot; same for all
    flight_mode::String                       # flight mode pretty descriptive 
    landed::Int64                             # Landed or not
    compute_time::Vector{F}                   # Compute_time for something # TODO: check this
    return_mode::Bool                         # Return mode for the quad
end

QuadObj() = QuadObj(
    Int64, # id 
    Float64[], # u
    NomTrajectory(), # nom traj 
    NomTrajectory(), # can traj
    NomTrajectory(), #cand_traj
    Float64[], # cand_traj_time
    Float64[], # can_soc_end
    Float64[], # charge_pos
    Float64,   # init soc
    Float64,   # discharge rate
    Float64,   # charge start time
    Float64,   # charging duration
    String,    # flight mode
    Int64,     # flight status: 1 = landed, 0: flying
    Float64[], # Compute time for candidate trajectory
    Bool)
    


# Differenet Battery dynamics
battery_dynamics_control(u, p, t) = -0.1*norm(p, 2)^2
battery_dynamics_linear(u,p,t) = -(p)

"""
Checks whether the quad is within a norm ball around the charging station 
"""
function reached_charging_station(comb_state, flight_mode, charging_pos)
    # Here comb state is of only single quad + its battery level
    robot_pos = comb_state[1:3]
    dist = norm(robot_pos - charging_pos) <= 0.05
    if (dist) && (flight_mode == "Landing")
        return true, dist
    else
        return false, dist
    end
end

function recompute_u_track_cand!(integrator)
    
    # The reason I made this as a seperate function was because i wanted one function just
    # for the forward integration of the candidate trajectory and the other one for strictly 
    # following the candidate trajectories for all of the quadrotors
    
    t = integrator.t
    qi_state = integrator.u # Get the quad state

    at_charging, dist = reached_charging_station(integrator.u.x, integrator.p[1].flight_mode, integrator.p[1].charge_pos)

    # if (integrator.p[1].landed == 0)
    #     println("Distance to the charging station: $(dist)")
    # end

    # Check is it the first time quad has been reported to be at the charging station since it last charged
    if (integrator.p[1].landed == 0) && (at_charging == true)
        
        # Store the current time so that charging time can be tracked
        integrator.p[1].charge_start_time = t

        # Refill the battery
        integrator.u.soc = integrator.p[1].init_soc
        printstyled("-----------------------\n"; color = :blue)
        printstyled("charge_start_time: $(t)\n"; color = :blue)
        printstyled("cand_end_time: $(integrator.p[1].com_traj.ts[end])\n"; color = :blue)
        # end-1 because soc at the end is of the candidate trajectory which was not committed 
        printstyled("soc_at_end: $(integrator.p[1].can_soc_end[end-1])\n"; color = :blue)
        printstyled("-----------------------\n")
        
    end

    # Control for Quad 
    if (at_charging)
        # If at charging station then keep flag the landed status
        integrator.p[1].landed = 1
    end 

    ui = ExpController.geometric_π_tracker(t, qi_state, integrator.p[1].com_traj; landing_status = integrator.p[1].landed, return_mode = integrator.p[1].return_mode, ren_point_1 = integrator.p[1].return_pos)
    
    integrator.p[1].u = ui
end

function simulate_u_track_cand!(integrator)
    # This is only for 1 quad
    t = integrator.t
    quad_state = integrator.u
    
    # params = [control_input, candidadate trajectory]
    u = ExpController.geometric_π_tracker(t, quad_state, integrator.p[2])
    integrator.p[1] = u
end

function simulate_tracking_candidate(t0, initial_state, cand_traj::NomTrajectory, discharge_rate)    
    
    tspan = (t0, cand_traj.ts[end])
    params = [[0,0.,0,0], cand_traj, discharge_rate]
    
    prob = ODEProblem(ExpDyn.closed_loop_simulate!, initial_state, tspan, params)
    cb_u = PeriodicCallback(simulate_u_track_cand!, 0.01; initial_affect=true)
    sol = solve(prob, Tsit5(), callback=cb_u)
    return sol
end


function get_candidate_traj!(integrator, t_k, state_k, traj_nom_k, discharge_rate)

    ## Generate the b2b Trajectory
    # tk_n  = end of time for nominal trajectory
    t_kN = traj_nom_k.ts[end]
    
    # state_kn = end state at nominal trajectory
    state_kN_end = traj_nom_k.xs[end]
    
    ## Get the backup trajectory  ###(x0; t0 = 0.0, N_mpc = 40, T_H_sec = 2.0, TRAJ = "NM")
    # println("here")   
    # Construct the first candidate trajectory
    comp_time_b2b = @elapsed traj_b2b_kn = ExpController.DI_get_nominal_mpc_π(state_kN_end, 1; t0 = t_kN, T_H_sec = 9.0, N_mpc = 180, TRAJ = "b2b_mid", ren_point_1 = integrator.p[1].return_pos, ren_point_2 = integrator.p[1].charge_pos)
    # println("here2")
    # Construct the nom + b2b_1
    traj_cand = NomTrajectory()
    traj_cand = deepcopy(traj_nom_k)       
    
    # Pop the last element in the candidate trajectory because that corresponds to the last state of the nominal trajectory which will be replaced  by the b2b_1
    pop!(traj_cand.xs)
    pop!(traj_cand.ts)

    ##  Then append b2b to nominal trajectory
    # index starts from 2 since candidate trajectory already contains the start state of the b2b trajectory as nom traj last state
    for i = 1:length(traj_b2b_kn.ts)
        push!(traj_cand.ts, traj_b2b_kn.ts[i])
        push!(traj_cand.xs, traj_b2b_kn.xs[i])
    end

    # Remember u has 1 lentgh less than x
    for i = 1:length(traj_b2b_kn.us)
        push!(traj_cand.us, traj_b2b_kn.us[i])
    end

    # Record the compute time
    push!(integrator.p[1].compute_time, comp_time_b2b)

    # println("Candidate trajectory time: $(traj_cand.ts)")
    # println("Candidate trajectory time: $(traj_cand.ts)")
    
    
    ##  Evaluate the candidate trajectory
    
    # First forward propogate the trajectory
    # Here the quad state has to be given i.e. state_k
    traj_cand_sol_k = simulate_tracking_candidate(t_k, state_k, traj_cand, discharge_rate)

    # Time taken by  the candidate trajectory (end of candidate tracjectory - start state) by each quad
    cand_traj_duration = traj_cand_sol_k.t[end] - t_k
    push!(integrator.p[1].cand_traj_time, cand_traj_duration)

    # Remember here traj_sol_cand_k and store the SoC at the end
    combined_state_sol = traj_cand_sol_k.u
    push!(integrator.p[1].can_soc_end, combined_state_sol[end].soc)

    # Store the candidate trajectory
    integrator.p[1].can_traj = traj_cand 

end



function recompute_committed_timed_quads!(integrator)
    """
    Recompute the committed trajectory
    """
    println("=========== In recompute_committed ==============")

    # current time and current state
    t_k = integrator.t

    #Unpack the states
    state_k = integrator.u
    println("Current time: $(t_k)")

    # ------------ TBDel ------------
    # state_quad = state_k[1:end-1]
    # ------------ TBDel ------------

    # If the quad is done charging station, make sure to change the status of flight mode and landed status
    charge_time_dur = 5.0
    return_time_dur = 2.0
    if (t_k >= integrator.p[1].charge_start_time + charge_time_dur) && (integrator.p[1].landed  == 1)

        integrator.p[1].return_mode = true
    end

    if (t_k >= integrator.p[1].charge_start_time + charge_time_dur + return_time_dur) && (integrator.p[1].landed  == 1)
        # Change the flight mode to Mission
        integrator.p[1].flight_mode = "Mission"

        # Change the landed mode to negative
        integrator.p[1].landed = 0

        # set the return mode to false
        integrator.p[1].return_mode = false
    end

    # If Quad is not in mission mode, no need to plan more; The last committed trajectory is being followed right now or 
    # it is at the charging station 
    if integrator.p[1].flight_mode != "Mission"
        return
    end

    ## get nominal trajectory for the quad
    # Send the state of the robot 
    di_state = ExpDyn.get_DI_state_from_quad(state_k)
    # If TRAJ given: It has three options: SP (For fixed setpoint), NM (Default one; for tracking the nominal trajectory), B2B (Back-to-base traj)
    comp_time_nom = @elapsed traj_nom_k = ExpController.DI_get_nominal_mpc_π(di_state, 1; t0 = t_k) ##TRAJ = "SP"
    integrator.p[1].nom_traj = traj_nom_k

    ## Generate the candidate trajectory
    get_candidate_traj!(integrator, integrator.t, integrator.u, integrator.p[1].nom_traj, integrator.p[1].discharge_rate)

    
    # ====== Eware check ========
    # Now check the battery SoC condition
    if integrator.p[1].can_soc_end[end] > 0
        # candidate become committed {done}
        integrator.p[1].com_traj = integrator.p[1].can_traj
    else
        println("Trajectory not committed due to eware; Quad returning")
        integrator.p[1].flight_mode = "Landing"
    end

    println("=========== Out recompute_committed ==============")
    return

end


### ========================================================================

function simulate(quad_init_state, t0, t_max, quad_id, discharge_rate)
    tspan = (t0, t_max)

    # Let this code by used by multiple robots 
    # Eware algorithm requires robot to have a separate charging base 
    # For now their initial position will be considered as the place where they return for charging

    # Here initial state consists of Quad_state and state of charge of the battery
    # so it is (14x1) , where (13x1) is Quad state and 1 is SOC
    
    # Now the states of multiple of quadrotors are concatenated
    # So, for the case of two quadrotors it is (28x1), where (14x1) is the augmentated state for quad 1
    
    # Quad Object
    # id::Int64                                 # robot number
    # u::Vector{F}                              # Control input
    # nom_traj::NomTrajectory{Float64}          # Nominal trajectory
    # can_traj::NomTrajectory{Float64}          # Candidate trajectory
    # com_traj::NomTrajectory{Float64}          # Committed trajectory
    # cand_traj_time::Vector{F}                 # candidate trajectory end_time
    # can_soc_end::Vector{F}                    # Soc at the end of the candidate trajectory
    # charge_pos::Vector{F}                     # Charging position
    # init_soc::Float64                         # Initial soc of the robot
    # discharge_rate::Float64                   # Discharge rate of the robot
    # charge_start_time::Float64                # Charge start of the robot
    # T_char::Float64                           # charge_time for the robot; same for all
    # flight_mode::String                       # flight mode pretty descriptive 
    # landed::Int64                             # Landed or not
    # compute_time::Vector{F}                   # Compute_time for something # TODO: check this
    
    # Create N quad objects
    u0 = [0., 0, 0, 0] # same for the other quads
    init_flight_mode = "Mission"
    @unpack x, soc = quad_init_state

    quad_charge_pos = [0.0; 0.0; 3.0]
    quad_return_pos = [0.0; 0.0; 6.0]
    
    params = QuadObj[] # Initialing a vector of type QuadObj()
    quad = QuadObj(quad_id,
                u0,
                NomTrajectory(), 
                NomTrajectory(),
                NomTrajectory(),
                Float64[],
                Float64[],
                quad_charge_pos,
                quad_return_pos, 
                soc, # initial soc
                discharge_rate,
                0.0,
                15.0,
                init_flight_mode,
                0,
                Float64[],
                false)

    push!(params, quad)

    prob = ODEProblem(ExpDyn.closed_loop_quad!, quad_init_state, tspan, params) 
    
    # Callbacks
    
    # Committed trajectories for the Quad
    cb_traj_com_quads = PeriodicCallback(recompute_committed_timed_quads!, 1.0; initial_affect=true)

    # Control_input
    cb_u = PeriodicCallback(recompute_u_track_cand!,0.01; initial_affect=true, save_positions=(false, false) )
    
    # Combine the set of callbacks
    cbs = CallbackSet(cb_traj_com_quads, cb_u)

    sol = DEsolve(prob, Tsit5(), callback=cbs)
    
    return sol, params
end

end