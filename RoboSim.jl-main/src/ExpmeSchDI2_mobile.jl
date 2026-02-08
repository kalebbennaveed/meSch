module meSchExpDImobile

"""
Eware implementation for the the experiments
- Robot Model: Double integrator model
- Planer/High-level controller: Convex MPC
- Controller: Double integrator trajectories (x,y,z, dx, dy, dz, u1. u2, u3) are tracked by the geometric controller
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
    init_soc::Float64                         # Initial soc of the robot
    discharge_rate::Float64                   # Discharge rate of the robot
    charge_start_time::Float64                # Charge start of the robot
    T_char::Float64                           # charge_time for the robot; same for all
    flight_mode::String                       # flight mode pretty descriptive 
    landed::Int64                             # Landed or not
    compute_time::Vector{F}                   # Compute_time for something # TODO: check this
    flight_time::Vector{F}
    charge_dist::Vector{F}
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
    Float64[],
    Float64[])

    # For Multiple quads, this quad defines whats associated with each quad
mutable struct RoverObj{F}
    u::Vector{F} # Control input
end

RoverObj() = RoverObj(Float64[])


# Differenet Battery dynamics
battery_dynamics_control(u, p, t) = -0.1*norm(p, 2)^2
battery_dynamics_linear(u,p,t) = -(p)

"""
Checks whether the quad is within a norm ball around the charging station 
"""
function reached_charging_station(comb_state, flight_mode, charging_pos)
    # Here comb state is of only single quad + its battery level
    robot_pos = comb_state[1:3]
    dist = norm(robot_pos - charging_pos)
    # if (flight_mode == "Landing")
    #     println("flight_mode: $(flight_mode)")
    #     println("robot pos: $(robot_pos)")
    #     println("charging pos: $(charging_pos)")
    #     println("dist:  $(dist)")
    # end
    if (0.00 <= norm(robot_pos - charging_pos) <= 0.05) && (flight_mode == "Landing")
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
    # Get the quad state
    quad_keys = keys(integrator.u)[1:2]
    for i = 1:length(integrator.p[1])
        qi_state = integrator.u[quad_keys[i]] 
        
        at_charging, dist = reached_charging_station(integrator.u[quad_keys[i]].x, integrator.p[1][i].flight_mode, integrator.p[1][i].charge_pos)
        dist_to_rover = norm(integrator.u[quad_keys[i]].x[1:3] - integrator.u["rover"].pos)
        push!(integrator.p[1][i].charge_dist, dist_to_rover)
        # Check is it the first time quad has been reported to be at the charging station since it last charged
        if (integrator.p[1][i].landed == 0) && (at_charging == true)
            
            # Store the current time so that charging time can be tracked
            integrator.p[1][i].charge_start_time = t
            println("integrator.p[1][i].init_soc : $(integrator.p[1][i].init_soc)")
            if i == 1
                integrator.u.quad1.soc = integrator.p[1][i].init_soc
            elseif i == 2
                integrator.u.quad2.soc = integrator.p[1][i].init_soc
            end
            printstyled("-----------------------\n"; color = :blue)
            printstyled("charge_start_time: $(t)\n"; color = :blue)
            printstyled("cand_end_time: $(integrator.p[1][i].com_traj.ts[end])\n"; color = :blue)
            # end-1 because soc at the end is of the candidate trajectory which was not committed 
            printstyled("soc_at_end: $(integrator.p[1][i].can_soc_end[end-1])\n"; color = :blue)
            printstyled("-----------------------\n")
            
        end
        # Control for Quad 
        if (at_charging)
            # If at charging station then keep flag the landed status
            integrator.p[1][i].landed = 1
        end 
        # println("Going to controller")
        ui = ExpController.geometric_π_tracker(t, qi_state, integrator.p[1][i].com_traj; landing_status = integrator.p[1][i].landed)
        # println("ui : $(ui)")
        integrator.p[1][i].u = ui

    # Controlle for rover
    u = ExpController.NonLinController_MD(t, integrator.u.rover)
    integrator.p[2].u = u

    end
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
   
function get_candidate_traj!(integrator, t_k, state_k, traj_nom_k, quad_key, q_num, discharge_rate, ren_state)
    
    # println("IN cand")
    # println("quad_key: $(quad_key)")
    # println("q_num: $(q_num)")
    ## Generate the b2b Trajectory
    # tk_n  = end of time for nominal trajectory
    t_kN = traj_nom_k.ts[end]
    
    # state_kn = end state at nominal trajectory
    state_kN_end = traj_nom_k.xs[end]
    
    ## Get the backup trajectory  ###(x0; t0 = 0.0, N_mpc = 40, T_H_sec = 2.0, TRAJ = "NM")
    # println("Going to b2b")
    comp_time_b2b = @elapsed traj_b2b_kn = ExpController.DI_get_nominal_mpc_π(state_kN_end, integrator.p[1][q_num]; t0 = t_kN, T_H_sec = 8.0, N_mpc = 160, TRAJ = "B2B", ren_point = ren_state)
    push!(integrator.p[1][q_num].compute_time, comp_time_b2b)
    # println("Back from b2b")
    ## Combine the nominal trajectory and backup trajectory
    traj_cand = NomTrajectory()
    
    # First candidate == nominal trajectory
    traj_cand = deepcopy(traj_nom_k)
    
    # pop the last state from the traj_cand
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
    
    ##  Evaluate the candidate trajectory
    # println("Reached here in cand")
    # First forward propogate the trajectory
    # Here the quad state has to be given i.e. state_k
    traj_cand_sol_k = simulate_tracking_candidate(t_k, state_k, traj_cand, discharge_rate)

    # Time taken by  the candidate trajectory (end of candidate tracjectory - start state) by each quad
    cand_traj_duration = traj_cand_sol_k.t[end] - t_k
    push!(integrator.p[1][q_num].cand_traj_time, cand_traj_duration)

    # Remember here traj_sol_cand_k and store the SoC at the end
    combined_state_sol = traj_cand_sol_k.u
    push!(integrator.p[1][q_num].can_soc_end, combined_state_sol[end].soc)

    # Store the candidate trajectory
    integrator.p[1][q_num].can_traj = traj_cand 

    # println("Leaving cand")

end

function flight_time_estimation_cf_lin(state_k, discharge_rate)
    """
    At each time step, compute the remaining flight time of the quad
    Forward propogate through the quadrotor battery model
    Assuming the battery model to be linear right now
    """
    current_soc = state_k.soc
    println("current_soc: $(current_soc)")
    T_RF = current_soc / discharge_rate
    return T_RF   
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
    quad_keys = keys(state_k)[1:2]
    println("Current time: $(t_k)")

    # Get the state of  the charging station
    ren_x, ren_y, ren_z, ϕd = RefTrajLib.figureEight(t_k + 10.0)
    ren_state  = [ren_x[1];ren_y[1];0.40]

    for i = 1:length(integrator.p[1])
        # ------------ TBDel ------------
        # state_quad = state_k[1:end-1]
        # ------------ TBDel ------------

        # If the quad is done charging station, make sure to change the status of flight mode and landed status
        charge_time_dur = 5.0
        println("integrator.p[1][i].charge_start_time: $(integrator.p[1][i].charge_start_time)")
        println("charge_time_dur: $(charge_time_dur)")
        println("integrator.p[1][i].landed : $(integrator.p[1][i].landed )")

        if (t_k >= integrator.p[1][i].charge_start_time + charge_time_dur) && (integrator.p[1][i].landed  == 1)
            println("Here: Quad $(i) returning to mission")

            # Change the flight mode to Mission
            integrator.p[1][i].flight_mode = "Mission"

            # Change the landed mode to negative
            integrator.p[1][i].landed = 0

            # # Refill the battery
            # # TODO: Fix this stupid hack
            # println("integrator.p[i].init_soc : $(integrator.p[i].init_soc)")
            # if i == 1
            #     integrator.u.quad1.soc = integrator.p[i].init_soc
            # elseif i == 2
            #     integrator.u.quad2.soc = integrator.p[i].init_soc
            # elseif i == 3
            #     integrator.u.quad3.soc = integrator.p[i].init_soc
            # end
            # println("integrator.u[quad_keys[i]].soc : $(integrator.u[quad_keys[i]].soc)")
        end
        println("integrator.u[quad_keys[i]].soc : $(integrator.u[quad_keys[i]].soc)")

        # Update flight time
        flight_time = flight_time_estimation_cf_lin(integrator.u[quad_keys[i]], integrator.p[1][i].discharge_rate)
        push!(integrator.p[1][i].flight_time, flight_time)       

        # If Quad is not in mission mode, no need to plan more; The last committed trajectory is being followed right now or 
        # it is at the charging station 
        if integrator.p[1][i].flight_mode != "Mission"
            continue
        end

        ## get nominal trajectory for the quad
        # Send the state of the robot 
        di_state = ExpDyn.get_DI_state_from_quad(integrator.u[quad_keys[i]])
        # If TRAJ given: It has three options: SP (For fixed setpoint), NM (Default one; for tracking the nominal trajectory), B2B (Back-to-base traj)
        comp_time_nom = @elapsed traj_nom_k = ExpController.DI_get_nominal_mpc_π(di_state, i; t0 = t_k) ##TRAJ = "SP"
        integrator.p[1][i].nom_traj = traj_nom_k

        ## Generate the candidate trajectory
        get_candidate_traj!(integrator, integrator.t, integrator.u[quad_keys[i]], integrator.p[1][i].nom_traj, quad_keys[i], i, integrator.p[1][i].discharge_rate, ren_state)
    end


    ## meSch check ##
    quads_list = [quad for quad in integrator.p[1]]
    deleteat!(quads_list, findall(x -> x.flight_mode == "Landing", quads_list))

    if length(quads_list) > 1
        sort!(quads_list, by = x -> x.flight_time[end])
        gap_flags = Bool[] # initialize the gap flags
        req_gap = 15.0 
        println("Remaining flight time 1st quad: $(quads_list[1].flight_time[end])")
        for i = length(quads_list):-1:2
            gap_distance_i = quads_list[i].flight_time[end] - (quads_list[1].cand_traj_time[end] - 1.0)
            # println("quads_list[i].flight_time[end]: $(quads_list[i].flight_time[end])")
            # println("quads_list[1].cand_traj_time[end]: $(quads_list[1].cand_traj_time[end])")
            gap_flag_i = gap_distance_i > (((i - 1) * req_gap) + 1.0)
            println("Remaining flight time: $(quads_list[i].flight_time[end])")
            println("TC: $((quads_list[1].cand_traj_time[end] - 1.0))")
            println("TF - TC > k(T_c + T_delta): $(gap_distance_i) > $((((i - 1) * req_gap) + 1.0))")
            println("Expected gap: $(gap_distance_i - (((i - 1) * req_gap) + 1.0))")
            if gap_flag_i == true
                printstyled("Gap satisfied: $(gap_flag_i)\n"; color = :green)
            else
                printstyled("Gap satisfied: $(gap_flag_i)\n"; color = :red)
            end
            
            push!(gap_flags, gap_flag_i)
        end   
        
        # [Gap Violation] Check if any of the gap flags is false or 0, which means the gap requirement is not met
        gap_violation = false in gap_flags
        println("current time: $(integrator.t)")
        if gap_violation == true
            printstyled("GAP FLAG: $(gap_violation)\n"; color = :red)

            # Change the mode of the quad with Lower flight time to "Landing"
            quad_LFT_id = quads_list[1].id
            integrator.p[1][quad_LFT_id].flight_mode = "Landing"
            for i = length(quads_list):-1:2
                quad_id = quads_list[i].id
                integrator.p[1][quad_id].com_traj = integrator.p[1][quad_id].can_traj
                integrator.p[1][quad_id].charge_pos = ren_state
            end               
        else
            printstyled("GAP FLAG: $(gap_violation)\n"; color = :green)
            # Being here means gaps were fine 
            for i = 1:length(quads_list)
                quad_id = quads_list[i].id
                if integrator.p[1][quad_id].can_soc_end[end] > 0. 
                    integrator.p[1][quad_id].com_traj = integrator.p[1][quad_id].can_traj
                    integrator.p[1][quad_id].charge_pos = ren_state
                else
                    integrator.p[1][quad_id].flight_mode = "Landing"
                end 
            end           
        end           

    elseif length(quads_list) == 1

        # get the id of the robot still in the air
        quad_air_id = quads_list[1].id
        # Now check the battery SoC condition
        if integrator.p[1][quad_air_id].can_soc_end[end] > 0
            # candidate become committed {done}
            integrator.p[1][quad_air_id].com_traj = integrator.p[1][quad_air_id].can_traj
            integrator.p[1][quad_air_id].charge_pos = ren_state
        else
            integrator.p[1][quad_air_id].flight_mode = "Landing"
        end
    else
        return 
    end
    println("=========== Out recompute_committed ==============")
    return
end





### ========================================================================

function simulate(sys_ic, t0, t_max, discharge_rate)
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
    u0_rover = [0., 0]
    init_flight_mode = "Mission"
    quad_charge_pos = [0.0; 0.5; 0.11]
    quad_keys = keys(sys_ic)[1:2]

    q_params = QuadObj[] # Initialing a vector of type QuadObj()
    quad1 = QuadObj(1,
                u0,
                NomTrajectory(), 
                NomTrajectory(),
                NomTrajectory(),
                Float64[],
                Float64[],
                quad_charge_pos, 
                sys_ic[quad_keys[1]].soc, # initial soc
                discharge_rate,
                0.0,
                15.0,
                init_flight_mode,
                0,
                Float64[],
                Float64[],
                Float64[])


    quad2 = QuadObj(2,
                u0,
                NomTrajectory(), 
                NomTrajectory(),
                NomTrajectory(),
                Float64[],
                Float64[],
                quad_charge_pos, 
                sys_ic[quad_keys[2]].soc, # initial soc
                discharge_rate,
                0.0,
                15.0,
                init_flight_mode,
                0,
                Float64[],
                Float64[],
                Float64[])


    push!(q_params, quad1)
    push!(q_params, quad2)

    r_param = RoverObj(u0_rover)
    params = [q_params, r_param]

    prob = ODEProblem(ExpDyn.N2_quad_rover_closed_loop!, sys_ic, tspan, params) 
    
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