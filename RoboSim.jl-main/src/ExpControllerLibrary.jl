module ExpController

"""
This module is created for the Double Integrator Convex MPC. 
This is more of a branch code which will be merged with the Controller module after development
1. Set the equilibrium state with control inout equal to g and the state equal some height. Maybe have the same setting for the
   
"""


## TODO Remove the g acceleration from the double integrator model
## TODO For the landing state, the DI u_ref should be [0, 0, 0]
## TODO Probably remove the control input bounds or have something like 


using LinearAlgebra
using BlockDiagonals
using ForwardDiff
using ControlSystems
import Convex as cvx
using ECOS
import DifferentialEquations: solve as DEsolve
import DifferentialEquations: ODEProblem, Tsit5, PeriodicCallback, CallbackSet
using StaticArrays
using ComponentArrays
using Parameters: @unpack


include("ExpDynamicsLibrary.jl")
include("RefTrajectoryLibrary.jl")
include("utils.jl")

# Some utils

SV3{F} = SVector{3,F}
SV4{F} = SVector{4,F}
SM3{F} = SMatrix{3,3,F,9}
SM4{F} = SMatrix{4,4,F,16}


# Number of states and control inputs
DINx = ExpDyn.DINx
DINu = ExpDyn.DINu

# Equilibrium states 
DIuhover = ExpDyn.DIuhover
DIxeq = ExpDyn.DIx0


# DI Control input limits
DIumin = ExpDyn.DIumin
DIumax = ExpDyn.DIumax

# Get the A and B matrices  and the cost matrices Q and R for the DI model
A, B = ExpDyn.get_DI_3Dmodel()
DI_Q, DI_R = ExpDyn.get_DI_cost_matrices()


# Geometric Controller Parameters
# Geometric Controller parameters

quadrotor_parameters = ExpDyn.quadrotor_parameters
controller_params = ComponentArray(
    kx = 3.0, 
    kv = 5.0,
    kR = 0.90,
    kΩ = 0.120, 
    m = 0.680,
    J = quadrotor_parameters.J, # this is cheating, but is convenient,
    invG = quadrotor_parameters.invG, # this is also cheating...
    g = 9.81
)


QuadStateAxes = Axis(
    x = 1:3,
    v = 4:6,
    R = ViewAxis(7:15, ShapedAxis((3, 3), NamedTuple())),
    Ω = 16:18,
    ω = 19:22,
)
QuadTargetStateAxes =
    Axis(xd = 1:3, vd = 4:6, ad = 7:9, b1d = 10:12, Ωd = 13:15, αd = 16:18)
QuadFlatStateAxes =
    Axis(x = 1:3, v = 4:6, a = 7:9, j = 10:12, s = 13:15, ψ = 16, ψd = 17, ψdd = 18)


# Nominal Trajectory from MPC
struct NomTrajectory{F}
    ts::Vector{F}
    xs::Vector{Vector{F}}
    us::Vector{Vector{F}}
end

# The type of the nominal trajectory
NomTrajectory() = NomTrajectory(Float64[],Vector{Float64}[],Vector{Float64}[])

# Use the same mpc to create a backpack andno
function DI_get_nominal_mpc_π(x0, q_n; t0 = 0.0, N_mpc = 80, T_H_sec = 4.0, TRAJ = "NM", ren_point_1 = [0.0; 0.0; 0.11], ren_point_2 = [0.0; 0.0; 0.11], ω_val = 1.0)
    
    # NM: Nominal Trajectory
    # B2B: Back to Base Trajectory

    # Her x0 is the initial state for the DI model; Make sure the double integrator initial state is sent
    # x0 = [pos_x, pos_y, vel_x, vel_y]
    
    # define variables
    dt = T_H_sec / N_mpc
    
    # =================================================
    
    # Get the DI Desired trajectories for state and control
    traj_type = 1 
    if TRAJ == "NM"
        # println("here")
        DIX_ref, DIU_ref = RefTrajLib.DI_3D_get_desired_trajectory(t0, traj_type, N_mpc, dt, q_n)
    elseif TRAJ == "SP"
        # println("")
        # Get the setpoint which does not correspond to the landing
        DIX_ref, DIU_ref = RefTrajLib.DI_3D_get_desired_nonlanding_point(t0, N_mpc, dt, q_n)
    elseif TRAJ == "b2b_mid"
        DIX_ref, DIU_ref = RefTrajLib.DI_3D_get_desired_point_mid_stage(t0, N_mpc, dt, q_n, ren_point_1, ren_point_2)
    else
        # Get the B2B trajectory
        DIX_ref, DIU_ref = RefTrajLib.DI_3D_get_desired_point(t0, N_mpc, dt, q_n, ren_point_1)
    end
    

#     X_traj = cvx.Variable(Nx̅, N_mpc)...
    DIX_traj = [cvx.Variable(DINx) for i=1:N_mpc] 
#     U_traj = cvx.Variable(Nu, N_mpc - 1)...
    DIU_traj = [cvx.Variable(DINu) for i=1:(N_mpc-1)];
    
    # Cost function 
    cost = 0.0
    
    # state cost
    for i = 1:N_mpc
        cost += 0.5*cvx.quadform(DIX_traj[i] - DIX_ref[i], DI_Q)
    end
    
    # control input cost
    for i = 1:(N_mpc - 1)
        cost += 0.5*cvx.quadform(DIU_traj[i] , DI_R) # - DIU_ref[i]
    end   
    
    prob = cvx.minimize(cost)
    
    # initial condition constraint
    prob.constraints += DIX_traj[1] == x0


    # Dynamics constraint
    for i = 1:(N_mpc - 1)
        prob.constraints += DIX_traj[i+1] == ExpDyn.DI_3D_dynamics_rk4(DIX_traj[i], DIU_traj[i], dt)        
    end  
    
    # Control input constraints # For now, no control input constarints
    # DIumin = [-4.; -4; 0.8*controller_params.m*controller_params.g]
    # DIumax = [4.; 4; 2.0*controller_params.m*controller_params.g]
    #     umin = [0.2*m*g; 0.2*m*g; 0.2*m*g; 0.2*m*g]
    # umax = [0.6*m*g; 0.6*m*g; 0.6*m*g; 0.6*m*g]
    # for i = 1:(N_mpc - 1)
    #     prob.constraints += DIU_traj[i] <= DIumax
    #     prob.constraints += DIU_traj[i] >= DIumin
    # end


    cvx.solve!(prob, ECOS.Optimizer; silent_solver = true)
    
    # Get the result
    DIX_traj_sol = cvx.evaluate.(DIX_traj)
    DIU_traj_sol = cvx.evaluate.(DIU_traj)
    
    return NomTrajectory(
        [t0 + (i-1) * dt for i=1:N_mpc],
        DIX_traj_sol,
        DIU_traj_sol
    )

    # TODO: Potential dynamics problem here
    # vectors = [0.0, 0.0, 3.0, 0, 0, 0]
    # vectoru =  [0.0, 0.0, 0.0]
    
    # return NomTrajectory(
    #     [t0 + (i-1) * dt for i=1:N_mpc],
    #     [vectors for i=1:N_mpc],
    #     [vectoru for i=1:N_mpc-1]
    # )
end

# Interpolation for trajectories
function DI_interpolate_traj(t, traj::NomTrajectory)
    @assert t >= traj.ts[1]
       
    if t >= traj.ts[end]
        # beyond the horizon - extrapolate
        h = t - traj.ts[end]
        # di_state_ref = A * traj.xs[end] + B * traj.us[end]
        di_state_ref = ExpDyn.DI_3D_dynamics_rk4(traj.xs[end],traj.us[end],h)
        di_u_ref = traj.us[:,end]
    else

        N = length(traj.ts)
        ind = findfirst(traj.ts .> t) - 1 # cant be nothing
        @assert !isnothing(ind)
        h = t - traj.ts[ind]
        di_state_ref = ExpDyn.DI_3D_dynamics_rk4(traj.xs[ind],traj.us[ind],h)
        # di_state_ref = A * traj.xs[ind] + B * traj.us[ind]
        di_u_ref = traj.us[ind]
    end
    
    return di_state_ref, di_u_ref
end




# function flat_state_to_quad_state(flat_states::ComponentVector{F}, g::F = F(9.81)) where {F}

#     @unpack x, v, a, j, s, ψ, ψd, ψdd = flat_states

#     # First lets see what does it return
#     out = flat_state_to_quad_state(x, v, a, j, s, ψ, ψd, ψdd, g)

#     return ComponentArray(vcat(out...), QuadTargetStateAxes)

# end



# """
# Based on "Accurate Tracking of Aggressive Quadrotor Trajectories using Incremental Nonlinear Dynamic Inversion and Differential Flatness" by Ezra Tal, Sertac Karaman
# """
function flat_state_to_quad_state(t, di_state, di_u, landing_status)

    # Define the flat states (pos, vel, acc, jerk and snap)
    x = di_state[1:3]
    ## define the deired acceleration
    ## define the desired yaw, yaw rate, and yaw acceleration at current time
    if landing_status == 1
        v = zeros(3)
        a = zeros(3)
        ψ = 0.0
        ψd = 0.0
        ψdd = 0.0
    else
        v = di_state[4:6]
        a = di_u
        ψ = 3.0*t/8
        ψd = 3.0/8
        ψdd = 0.0
    end

    ## define desried jerk and snap
    j = zeros(3)
    s = zeros(3)  

    g = controller_params.g

    # desired rotation matrix
    force = SA[a[1], a[2], a[3]+g]
    zb = normalize(force)
    xc = SA[cos(ψ), sin(ψ), 0]
    yb = cross(zb, xc) |> normalize
    xb = cross(yb, zb) |> normalize

    R = SMatrix{3,3}([xb yb zb])
    b1d = R * SA[1, 0, 0]

    # construct τ
    τ = norm(force)

    # construct S matrix
    bx1 = R[1, 1]
    bx2 = R[2, 1]
    bx3 = R[3, 1]
    by1 = R[1, 2]
    by2 = R[2, 2]
    by3 = R[3, 2]
    bz1 = R[1, 3]
    bz2 = R[2, 3]
    bz3 = R[3, 3]

    # 1x3 matrix
    S = SA[
        0;;
        (bx2 * bz1 - bx1 * bz2) / (bx1^2 + bx2^2);;
        (-bx2 * by1 + bx1 * by2) / (bx1^2 + bx2^2)
    ]

    # solve for Ω, τdot
    iz = SA[0, 0, 1]
    hatizT = QuaternionUtils.hat(iz)

    bz = R[:, 3]

    M = zero(MMatrix{4,4,Float64,16})
    M[1:3, 1:3] = τ * R * hatizT
    M[1:3, 4] = bz
    M[4, 1:3] = S
    invM = inv(SMatrix{4,4,Float64,16}(M))

    Ωτd = invM * [SVector{3}(j); ψd]

    Ω = SV3{Float64}(Ωτd[1], Ωτd[2], Ωτd[3])
    τd = Ωτd[4]

    # construct Sdot matrix
    # expression derived using mathematica
    w1, w2, w3 = Ω[1], Ω[2], Ω[3]

    Sd = SA[
        0;;
        (bx1 * Ω[1]) / (bx1^2 + bx2^2) +
        (bx2 * Ω[2]) / (bx1^2 + bx2^2) +
        ((bx1^2 * bz1 - bx2^2 * bz1 + 2 * bx1 * bx2 * bz2) * Ω[3]) / (bx1^2 + bx2^2)^2;;
        ((bx1^2 * bx2 + bx2^3 - bx1^2 * by1 + bx2^2 * by1 - 2 * bx1 * bx2 * by2) * Ω[3]) / (bx1^2 + bx2^2)^2
    ]

    # solve for α, τdd
    B1 = R * (2 * τd * hatizT + τ * QuaternionUtils.hat(Ω) * hatizT) * Ω
    B2 = (Sd * Ω)
    B = [B1; B2]
    ατdd = invM * [SVector{3}(s - B1); (ψdd - B2[1])]

    α = SV3{Float64}(ατdd[1], ατdd[2], ατdd[3])
    #     τdd = ατdd[4]
    
    return x, v, a, b1d, Ω, α

end


function geometric_controller(state, xd, vd, ad, b1d, Ωd, αd)

    @unpack x, v, R, Ω = state
    @unpack kx, kv, kR, kΩ, m, J, g = controller_params
    e3 = [0, 0, 1]
    ex = x - xd
    ev = v - vd
    
    # construct desired rotation matrix
    b3d = (-kx * ex - kv * ev + m * g * e3 + m * ad) |> normalize
    b2d = cross(b3d, normalize(b1d)) |> normalize
    b1d_n = cross(b2d, b3d) |> normalize
    
    Rd = [b1d_n b2d b3d]

    eR = 0.5 * QuaternionUtils.vee(Rd' * R - R' * Rd)
    eΩ = Ω - R' * Rd * Ωd

    f = dot(-kx * ex - kv * ev + m * g * e3 + m * ad, R * e3)
    M = -kR * eR - kΩ * eΩ + cross(Ω, (J * Ω)) - J * (QuaternionUtils.hat(Ω) * R' * Rd * Ωd - R' * Rd * αd)
    
    return [f, M...]

end

# Controller for the Quad 
function fM_to_ω(fM)
    @unpack invG = quadrotor_parameters
    ωsq = SMatrix{4,4}(invG) * SVector{4}(fM)
    return SVector{4}(QuaternionUtils.signedsqrt(w) for w in ωsq)
end


function geometric_π_tracker(t, state, traj::NomTrajectory; landing_status = 0,return_mode = false,ren_point_1 = [0.0; 0.0; 6.0])

    @unpack x, v, R, Ω = state

    # Interpolate the trajectory
    if (landing_status == 1) && (return_mode == false)
        # If landed then stay there 
        di_state_ref = x
        di_u_ref = zeros(3)
    elseif (landing_status == 1) && (return_mode = true)
        # If returning first return to the fixed point in the middle to avoid collisiion with the control room
        di_state_ref = ren_point_1
        di_u_ref = zeros(3)
    else
        di_state_ref, di_u_ref = DI_interpolate_traj(t, traj::NomTrajectory)
    end

    # Convert into the target quadrotor state\
    xd, vd, ad, b1d, Ωd, αd = flat_state_to_quad_state(t, di_state_ref, di_u_ref, landing_status)

    # determine required thrust and moments using the geometric controller
    fM = geometric_controller(state, xd, vd, ad, b1d, Ωd, αd)
    # Determine the desired control input
    control = fM_to_ω(fM)
    # println("Control input: $(control)")
    return control

end

function NonLinController_MD(t, state)
    # println("Here is nonlinear controller")
    @unpack pos, phi = state
    x = pos[1]
    y = pos[2]
    z = pos[3]
    ϕ = phi
    # println("Here2")
    # Get desired state (x_des and y_des ∈ (3x1) with pos, vel and acceleration)
    x_des, y_des, z_des, ϕd = RefTrajLib.figureEight(t)
    
    # Compute desired linear and angular velocity
    vd = sqrt(x_des[2]^2 + y_des[2]^2)
    ωd = (y_des[3] * x_des[2] - x_des[3] * y_des[2])/(x_des[2]^2 + y_des[2]^2)
    
    # Gains 
    k1 = 2.
    k2 = 1.
    k3 = 1.
    
    # Error coordinates
    xe = cos(ϕd) * (x - x_des[1]) + sin(ϕd) * (y - y_des[1])
    ye = -sin(ϕd) * (x - x_des[1]) + cos(ϕd) * (y - y_des[1])
    ϕe = ϕ - ϕd
    
    # Nonlinear Tracking controller
    v = (vd - k1 * abs(vd) * (xe + ye * tan(ϕe)))/cos(ϕe)
    ω = ωd - (k2 * vd * ye + k3 * abs(vd) * tan(ϕe)) * cos(ϕe)^2
    u = [v ; ω]
    
    return u
     
end

end