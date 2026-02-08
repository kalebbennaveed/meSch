module ExpDyn

using LinearAlgebra
using BlockDiagonals
using ForwardDiff
using ControlSystems
using ComponentArrays
using Rotations
using StaticArrays
using Parameters:@unpack

include("utils.jl")

# Quad Parameters

SV3{F} = SVector{3,F}
SV4{F} = SVector{4,F}
SM3{F} = SMatrix{3,3,F,9}
SM4{F} = SMatrix{4,4,F,16}

function initialize_quad_params(quad_p)

    @unpack J, k_f, k_μ, motor_pos, motor_dir = quad_p

    invJ = inv(J)

    G = zeros(4, 4)

    Gf = k_f * ones(1, 4)

    GM = zeros(3, 4)
    for i = 1:4
        G[1, i] = k_f
        G[2:4, i] = cross(motor_pos[:, i], [0, 0, k_f])
        G[4, i] = -motor_dir[i] * k_μ
    end

    invG = inv(G)

    return ComponentArray(quad_p; invJ = invJ, invG = invG)
end


# Quad Parameters
# m = 0.5
# ℓ = 0.1750
# J = Diagonal([0.0023, 0.0023, 0.004])
# J = diagm([0.0049, 0.0049, 0.0069])
# g = 9.81
# kt=1.0
# km=0.0245

quad_LX = 0.08 # 0.08
quad_LY = 0.08

const quadrotor_parameters =
    ComponentArray(
        mass = 0.680,            # mass of quadrotor [kg]
        J = diagm([0.005, 0.005, 0.003]) |> collect,                # moment of inertia of quadrotor in body axes [kg/m^2]
        τ_m = 0.02,            # time constant for quadrotor motors [s]
        J_m = 6.62e-6,         # moment of inertia of each motor+propeller [kg/m^2]
        g = 9.81,              # acceleration due to gravity [m/s^2]
        k_f = 1.91e-6,         # thrust coefficient, [N / (rad/s)^2]
        k_μ = 2.7e-7,          # torque coefficient, [Nm / (rad/s)^2]
        k_drag_f = 0.1,        # linear drag coefficient, [N / (m/s)^2]
        k_drag_μ = 0.003 * I(3), # angular drag coefficient, [Nm / (rad/s)^2]
        ω_max = 2520.0,          # maximum rotation rate of motors [rad/s]
        motor_pos = hcat(
            [
                [quad_LX, -quad_LY, 0.0],
                [-quad_LX, quad_LY, 0.0],
                [quad_LX, quad_LY, 0.0],
                [-quad_LX, -quad_LY, 0.0],
            ]...,
        ),                        # motor positions in body-fixed frame (i-th column contains i-th motor pos)
        motor_dir = [1, 1, -1, -1], # motor directions in body-fixed frame
    ) |> initialize_quad_params



m = quadrotor_parameters.mass
g = quadrotor_parameters.g
h = 0.01 #20 Hz

####  ========== DI Model [S] =========
# Quad Hover Condition (Equilibrium point)
DINx = 6 # (pos_x, pos_y, pos_z, vel_x, vel_y, vel_z)
DINu = 3 # (acc_x, acc_y, acc_z)
DIuhover = [0; 0; 0.0] #mg
DIx0 = [0.; 0; 1.0; 0; 0; 0] 
DIumin = [-2.; -2; 0.8*m*g]
DIumax = [2.; 2; 2.0*m*g]

function get_DI_cost_matrices()
    Q = Array(3.0*I(DINx ))
    R = Array(.1*I(DINu));
    return Q, R
end

function get_DI_3Dmodel()
    A = [0. 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1; 
        0 0 0 0 0 0; 
        0 0 0 0 0 0; 
        0 0 0 0 0 0]
    B = [0. 0 0; 
        0 0 0;
        0 0 0;
        1 0 0;
        0 1 0; 
        0 0 1]
    return A, B
end

function DI_3D_dynamics(x,u)
    r = x[1:3] #(3d position)
    v = x[4:6] #(velocity)
    ṙ = v
    v̇ = u
    return [ṙ; v̇]
end

#RK4 integration with zero-order hold on u
function DI_3D_dynamics_rk4(x,u,h)
    f1 = DI_3D_dynamics(x, u)
    f2 = DI_3D_dynamics(x + 0.5*h*f1, u)
    f3 = DI_3D_dynamics(x + 0.5*h*f2, u)
    f4 = DI_3D_dynamics(x + h*f3, u)
    xn = x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
    
    return xn
end


####  ========== DI Model [E] =========


#### =========== Quadrotor Model [S]==================

# Dynamics 
# ω : Angular acceleration in the body fixed frame
# Ω : Angular velocity in the body fixed frame
# signedsq(.) preserves the sign of the term

function motor_thrust(ω, k_f)

    return k_f * QuaternionUtils.signedsq(ω)
end

function hover_ω(params)
    @unpack k_f, mass, g = params
    return SVector{4}(sqrt(mass * g / (4 * k_f)) for i = 1:4)
end

function f_grav(state, params)
    @unpack mass, g = params
    return SA[0.0, 0, -mass*g]
end

function f_motors(state, params)
    @unpack ω = state
    @unpack k_f = params
    # in body frame
    return (SA[0.0, 0, motor_thrust(ω[i], k_f)] for i = 1:4) |> collect
end

function f_thrust(state, params)
    @unpack R = state
    # R is body-> world
    return SMatrix{3,3}(R) * sum(f_motors(state, params))
end

function f_drag(state, params, wind_ext)
    @unpack v = state
    @unpack k_drag_f = params

    δv = SVector{3}(v - wind_ext)

    return -k_drag_f * norm(δv) * δv
end

function μ_thrust(Dω, state, params)
    @unpack Ω, ω = state
    @unpack J_m, k_μ, motor_pos, motor_dir = params

    # contribution from thrust of each motor
    
    f_motors_ = f_motors(state, params)
    μ_1 = (sum(cross(motor_pos[:, i], f_motors_[i]) for i = 1:4))

    # contribution from acceleration of each motor
    μ_2 = SA[0, 0, sum(-motor_dir[i] * (k_μ * QuaternionUtils.signedsq(ω[i]) + J_m * (Dω)[i]) for i = 1:4)]
    μ_thrust = μ_1 + μ_2
    return μ_thrust
end

function μ_drag(state, params)
    @unpack Ω = state
    @unpack k_drag_μ = params
    return -k_drag_μ * norm(Ω) * SVector{3}(Ω)
end

function angular_momentum(state, params)
    @unpack Ω, ω = state
    @unpack J, motor_dir, J_m = params

    H = SVector{3}(J * Ω) - SA[0, 0, sum(motor_dir[i] * J_m * ω[i] for i = 1:4)]

    return H
end






# ω_des is the desired control input
#ds is the next state or the derivative of the state
function quadrotor3D!(ds, state, ω_des, discharge_rate;
    wind_ext = SA[0, 0, 0.0],
    f_ext = SA[0, 0, 0.0],
    μ_ext = SA[0, 0, 0.0],
)

    @unpack x, v, R, Ω, ω, soc = state
    @unpack mass,
    J,
    τ_m,
    J_m,
    g,
    k_f,
    k_μ,
    k_drag_f,
    k_drag_μ,
    motor_pos,
    motor_dir,
    ω_max = quadrotor_parameters

    ## Motor Dynamics
    ds.ω = (1 / τ_m) .* (clamp.(ω_des, 0, ω_max) - ω)

    ## Position Dynamics
    f_total =
        f_grav(state, quadrotor_parameters) +
        f_thrust(state, quadrotor_parameters) +
        f_drag(state, quadrotor_parameters, wind_ext) +
        f_ext
    ds.x = v
    ds.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.ω, state, quadrotor_parameters) + μ_drag(state, quadrotor_parameters) + μ_ext
    H = angular_momentum(state, quadrotor_parameters)

    ds.R = R * QuaternionUtils.hat(Ω)
    ds.Ω = J \ (μ_total - cross(Ω, H))
    ds.soc = -discharge_rate

    return

end





# ω_des is the desired control input
#ds is the next state or the derivative of the state
function N_quadrotor3D!(ds, quad_states, ω_list, discharge_rates;
    wind_ext = SA[0, 0, 0.0],
    f_ext = SA[0, 0, 0.0],
    μ_ext = SA[0, 0, 0.0],
)
    
    @unpack mass,
    J,
    τ_m,
    J_m,
    g,
    k_f,
    k_μ,
    k_drag_f,
    k_drag_μ,
    motor_pos,
    motor_dir,
    ω_max = quadrotor_parameters

    # Get the field names dynamically
    # For first quad 

    @unpack x, v, R, Ω, ω, soc = quad_states.quad1

    ## Motor Dynamics
    ds.quad1.ω = (1 / τ_m) .* (clamp.(ω_list[1], 0, ω_max) - ω)
    # println("ds.quad1.ω: $(ds.quad1.ω)")
    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad1, quadrotor_parameters) +
        f_thrust(quad_states.quad1, quadrotor_parameters) +
        f_drag(quad_states.quad1, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad1.x = v
    ds.quad1.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad1.ω, quad_states.quad1, quadrotor_parameters) + μ_drag(quad_states.quad1, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad1, quadrotor_parameters)

    ds.quad1.R = R * QuaternionUtils.hat(Ω)
    ds.quad1.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad1.soc = -(discharge_rates[1])

    @unpack x, v, R, Ω, ω, soc = quad_states.quad2

    ## Motor Dynamics
    
    ds.quad2.ω = (1 / τ_m) .* (clamp.(ω_list[2], 0, ω_max) - ω)

    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad2, quadrotor_parameters) +
        f_thrust(quad_states.quad2, quadrotor_parameters) +
        f_drag(quad_states.quad2, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad2.x = v
    ds.quad2.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad2.ω, quad_states.quad2, quadrotor_parameters) + μ_drag(quad_states.quad2, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad2, quadrotor_parameters)

    ds.quad2.R = R * QuaternionUtils.hat(Ω)
    ds.quad2.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad2.soc = -(discharge_rates[2])
    
    return 
end


function N3_quadrotor3D!(ds, quad_states, ω_list, discharge_rates;
    wind_ext = SA[0, 0, 0.0],
    f_ext = SA[0, 0, 0.0],
    μ_ext = SA[0, 0, 0.0],
)
    
    @unpack mass,
    J,
    τ_m,
    J_m,
    g,
    k_f,
    k_μ,
    k_drag_f,
    k_drag_μ,
    motor_pos,
    motor_dir,
    ω_max = quadrotor_parameters

    # Get the field names dynamically
    # For first quad 

    @unpack x, v, R, Ω, ω, soc = quad_states.quad1

    ## Motor Dynamics
    ds.quad1.ω = (1 / τ_m) .* (clamp.(ω_list[1], 0, ω_max) - ω)
    # println("ds.quad1.ω: $(ds.quad1.ω)")
    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad1, quadrotor_parameters) +
        f_thrust(quad_states.quad1, quadrotor_parameters) +
        f_drag(quad_states.quad1, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad1.x = v
    ds.quad1.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad1.ω, quad_states.quad1, quadrotor_parameters) + μ_drag(quad_states.quad1, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad1, quadrotor_parameters)

    ds.quad1.R = R * QuaternionUtils.hat(Ω)
    ds.quad1.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad1.soc = -(discharge_rates[1])

    @unpack x, v, R, Ω, ω, soc = quad_states.quad2

    ## Motor Dynamics
    
    ds.quad2.ω = (1 / τ_m) .* (clamp.(ω_list[2], 0, ω_max) - ω)

    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad2, quadrotor_parameters) +
        f_thrust(quad_states.quad2, quadrotor_parameters) +
        f_drag(quad_states.quad2, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad2.x = v
    ds.quad2.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad2.ω, quad_states.quad2, quadrotor_parameters) + μ_drag(quad_states.quad2, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad2, quadrotor_parameters)

    ds.quad2.R = R * QuaternionUtils.hat(Ω)
    ds.quad2.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad2.soc = -(discharge_rates[2])

    @unpack x, v, R, Ω, ω, soc = quad_states.quad3

    ## Motor Dynamics
    
    ds.quad3.ω = (1 / τ_m) .* (clamp.(ω_list[3], 0, ω_max) - ω)

    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad3, quadrotor_parameters) +
        f_thrust(quad_states.quad3, quadrotor_parameters) +
        f_drag(quad_states.quad3, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad3.x = v
    ds.quad3.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad3.ω, quad_states.quad3, quadrotor_parameters) + μ_drag(quad_states.quad3, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad3, quadrotor_parameters)

    ds.quad3.R = R * QuaternionUtils.hat(Ω)
    ds.quad3.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad3.soc = -(discharge_rates[2])    
    
    return 
end



function N2_quadrotor_rover_3D!(ds, quad_states, ω_list, u_rover, discharge_rates;
    wind_ext = SA[0, 0, 0.0],
    f_ext = SA[0, 0, 0.0],
    μ_ext = SA[0, 0, 0.0],
)
    
    @unpack mass,
    J,
    τ_m,
    J_m,
    g,
    k_f,
    k_μ,
    k_drag_f,
    k_drag_μ,
    motor_pos,
    motor_dir,
    ω_max = quadrotor_parameters

    # Get the field names dynamically
    # For first quad 

    @unpack x, v, R, Ω, ω, soc = quad_states.quad1

    ## Motor Dynamics
    ds.quad1.ω = (1 / τ_m) .* (clamp.(ω_list[1], 0, ω_max) - ω)
    # println("ds.quad1.ω: $(ds.quad1.ω)")
    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad1, quadrotor_parameters) +
        f_thrust(quad_states.quad1, quadrotor_parameters) +
        f_drag(quad_states.quad1, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad1.x = v
    ds.quad1.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad1.ω, quad_states.quad1, quadrotor_parameters) + μ_drag(quad_states.quad1, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad1, quadrotor_parameters)

    ds.quad1.R = R * QuaternionUtils.hat(Ω)
    ds.quad1.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad1.soc = -(discharge_rates[1])

    @unpack x, v, R, Ω, ω, soc = quad_states.quad2

    ## Motor Dynamics
    
    ds.quad2.ω = (1 / τ_m) .* (clamp.(ω_list[2], 0, ω_max) - ω)

    ## Position Dynamics
    f_total =
        f_grav(quad_states.quad2, quadrotor_parameters) +
        f_thrust(quad_states.quad2, quadrotor_parameters) +
        f_drag(quad_states.quad2, quadrotor_parameters, wind_ext) +
        f_ext
    ds.quad2.x = v
    ds.quad2.v = f_total / mass

    ## Rotational Dynamics
    μ_total = μ_thrust(ds.quad2.ω, quad_states.quad2, quadrotor_parameters) + μ_drag(quad_states.quad2, quadrotor_parameters) + μ_ext
    H = angular_momentum(quad_states.quad2, quadrotor_parameters)

    ds.quad2.R = R * QuaternionUtils.hat(Ω)
    ds.quad2.Ω = J \ (μ_total - cross(Ω, H))
    ds.quad2.soc = -(discharge_rates[2])

    # Rover states 
    # println("Updating rover states")
    @unpack pos, phi = quad_states.rover

    ds.rover.pos[1] = u_rover[1] * cos(phi)
    ds.rover.pos[2] = u_rover[1] * sin(phi)
    ds.rover.pos[3] = 0.
    ds.rover.phi = u_rover[2]

    return 
end



#### =========== Quadrotor Model [E] ==================

function get_DI_state_from_quad(state)
    @unpack x, v, R, Ω, ω, soc = state
    return [x;v]
end


#### Closed lopp dynamics for a single candidate trajectory forward integration (Different params)
function closed_loop_simulate!(dx, x, params, t)
    ω_des = params[1]
    discharge_rate = params[3]
    quadrotor3D!(dx, x, ω_des, discharge_rate)
    return
end

#### Closed Loop Dynamics for a single quad for whole sim
function closed_loop_quad!(dx, x, params, t)

    # First get the control inputs and landing statuses for each robot 
    # The landing logic is handled inside the controller: hovering within certain radius of charging station is considered as landding 
    ω_des = params[1].u
    discharge_rate = params[1].discharge_rate
    quadrotor3D!(dx, x, ω_des, discharge_rate)
    return
end
#### Closed 

#### Closed Loop Dynamics for a single quad for whole sim
function N_closed_loop_quad!(dx, x, params, t)

    # First get the control inputs and landing statuses for each robot 
    # The landing logic is handled inside the controller: hovering within certain radius of charging station is considered as landding 

    ω_list = Vector{Float64}[]
    discharge_rates = Float64[]
    for i = 1:length(params)
        push!(ω_list, params[i].u)
        push!(discharge_rates, params[i].discharge_rate)
    end
    # println("updated wlist : $(ω_list)")
    N_quadrotor3D!(dx, x, ω_list, discharge_rates)
    return
end


function N3_closed_loop_quad!(dx, x, params, t)

    # First get the control inputs and landing statuses for each robot 
    # The landing logic is handled inside the controller: hovering within certain radius of charging station is considered as landding 

    ω_list = Vector{Float64}[]
    discharge_rates = Float64[]
    for i = 1:length(params)
        push!(ω_list, params[i].u)
        push!(discharge_rates, params[i].discharge_rate)
    end
    # println("updated wlist : $(ω_list)")
    N3_quadrotor3D!(dx, x, ω_list, discharge_rates)
    return
end

function N2_quad_rover_closed_loop!(dx, x, params, t)

    # First get the control inputs and landing statuses for each robot 
    # The landing logic is handled inside the controller: hovering within certain radius of charging station is considered as landding 
    ω_list = Vector{Float64}[]
    discharge_rates = Float64[]
    for i = 1:length(params[1])
        push!(ω_list, params[1][i].u)
        push!(discharge_rates, params[1][i].discharge_rate)
    end
    u_rover = params[2].u
    # println("updated wlist : $(ω_list)")
    N2_quadrotor_rover_3D!(dx, x, ω_list, u_rover, discharge_rates)
    return
end

#### Closed 

end