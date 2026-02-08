module RefTrajLib

using ComponentArrays
using Parameters: @unpack
using StaticArrays, Serialization
using Interpolations

include("ExpDynamicsLibrary.jl")

m = ExpDyn.m
g = ExpDyn.g






#### DEMO 1 (2 Quads)'
# trajectory_params_0 = ComponentArray(
#     A = [0.5,1.0, 0.0],    # [5,5, 2.0]
#     ω = [0.4,0.5,0.], #ω = [0.5,0.5,0.],
#     ψ = [π/2, 0, 0],
#     off = [-1.2,0.2,1.5] # off
# )

# trajectory_params_1 = ComponentArray(
#     A = [0.5,1.0, 0.0],    # [5,5, 2.0] [2.7,2.7, 0.7]
#     ω = [0.4,0.5,0.], #
#     ψ = [π/2, 0, 0],
#     off = [0.7,0.2,1.5] # off
# )



### FLYLAB DEMO 1
# trajectory_params_0 = ComponentArray(
#     A = [2.6,2.9, 0.0],    # [5,5, 2.0]
#     ω = [5/8,4/8,6/8], #ω = [0.4,0.5,0.],
#     ψ = [π/2, 0, 0],
#     off = [-0.5,-1.4,7.0] # off
# )

### FLYLAB DEMO 2
trajectory_params_0 = ComponentArray(
    A = [1.1,1.4, 0.0],    # [5,5, 2.0]
    ω = [5/8,4/8,6/8], #ω = [0.4,0.5,0.],
    ψ = [π/2, 0, 0],
    off = [-2.0, 0.0,7.0] # off
)


trajectory_params_1 = ComponentArray(
    A = [1.1,1.4, 0.0],    # [5,5, 2.0]
    ω = [5/8,4/8,6/8], #ω = [0.4,0.5,0.],
    ψ = [π/2, 0, 0],
    off = [1.0,0.0,7.0] # off
)

trajectory_params_2 = ComponentArray(
    A = [2.5,1.2, 0.0],    # [5,5, 2.0]
    ω = [5/8,4/8,6/8], #ω = [0.4,0.5,0.],
    ψ = [π/2, 0, 0],
    off = [-0.5,-3.2,7.0] # off
)







# ====================== Lissa Jous trajectory params ===================
function lissajous(t, params, D::Integer = 0)
    @unpack A, ω, ψ, off = params
    return lissajous(t, A, ω, ψ, off, D)
end


function lissajous(t, A, ω, ψ, off, D = 0)
    N = length(A)
    if D == 0
        return SVector{N}(off[i] + A[i] * sin(ω[i] * t + ψ[i]) for i = 1:N)
    else
        return SVector{N}(A[i] * ω[i]^D * sin(ω[i] * t  + π * D /2 + ψ[i]) for i = 1:N)
    end
end
# ===============================================


# =============== Point Trajectory =================
function get_desired_trajectory_point(t0, N, dt, q_n)
    Xref = [zeros(12) for i=1:N]
    for i = 1:N
        ti = t0 + (i-1)*dt
        Xref[i][1:3] = [0.0; 0.0; 0.11]
    end
    # Control Input Trajectory
    Uref = [(g*m/4)*ones(4) for i = 1:(N-1)]

    return Xref, Uref
end
# ==================================================


# =============== Point Trajectory =================
function get_desired_trajectory_point_mobile(t0, N, dt, q_n, ren_point)
    Xref = [zeros(12) for i=1:N]
    for i = 1:N
        ti = t0 + (i-1)*dt
        Xref[i][1:3] = ren_point
    end
    # Control Input Trajectory
    Uref = [(g*m/4)*ones(4) for i = 1:(N-1)]

    return Xref, Uref
end
# ==================================================






### ------------------- Reference trajectory states for the DI Experimental code [Start] ------------------
# For 3 Quads 
# function get_desired_trajectory(t0, traj_type, N, dt, param_type)



function DI_3D_get_desired_trajectory(t0, traj_type, N, dt, param_type)
    
    if (traj_type == 1)
        DIXref = [zeros(6) for i=1:N]
        for i = 1:(N)
            ti = t0 + (i-1)*dt
            # Desired Position
            if param_type == 1
                DIXref[i][1:3] = lissajous(ti, trajectory_params_0)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_0, 1)
            elseif param_type == 2
                DIXref[i][1:3] = lissajous(ti, trajectory_params_1)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_1, 1)
            elseif param_type == 3
                DIXref[i][1:3] = lissajous(ti, trajectory_params_2)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_2, 1)
            elseif param_type == 4
                DIXref[i][1:3] = lissajous(ti, trajectory_params_4)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_4, 1)
            elseif param_type == 5
                DIXref[i][1:3] = lissajous(ti, trajectory_params_5)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_5, 1)
            end

        end

        # Control Input Trajectory
        DIUref = [[0; 0; 0.0] for i = 1:(N-1)]        
        return DIXref, DIUref
    end 
end

# This function gives a fixed setpoint as the target reference state
function DI_3D_get_desired_nonlanding_point(t0, N, dt, q_n)
    DIXref = [zeros(6) for i=1:N]
    for i = 1:N
        ti = t0 + (i-1)*dt
        DIXref[i][1:3] = [0.; 0; 1.5] # This could could be the position of the mobile charging robot or the state of the 
    end
    # Control Input Trajectory
    DIUref = [[0; 0; 0.0] for i = 1:(N-1)]
    return DIXref, DIUref 
end

function DI_3D_get_desired_point(t0, N, dt, q_n, ren_point)
    DIXref = [zeros(6) for i=1:N]
    for i = 1:N
        ti = t0 + (i-1)*dt
        DIXref[i][1:3] = ren_point # This could could be the position of the mobile charging robot or the state of the 
    end
    # Control Input Trajectory
    DIUref = [[0; 0; 0.0] for i = 1:(N-1)]
    return DIXref, DIUref 
end

function DI_3D_get_desired_point_mid_stage(t0, N, dt, q_n, ren_point_1, ren_point_2)
    DIXref = [zeros(6) for i = 1:N]
    for i = 1:N
        if i <= (N / 2)
            DIXref[i][1:3] = ren_point_1
        else
            DIXref[i][1:3] = ren_point_2
        end
    end
    # Control Input Trajectory
    DIUref = [[0; 0; 0.0] for i = 1:(N-1)]
    return DIXref, DIUref 
end

function DI_3D_get_desried_nominal_setpoint(t0, traj_type, N, dt, param_type)
    
    if (traj_type == 1)
        DIXref = [zeros(6) for i=1:N]
        for i = 1:(N)
            ti = t0 + (i-1)*dt
            # Desired Position
            if param_type == 1
                DIXref[i][1:3] = [0; 0; 3.0]
                DIXref[i][4:6] = zeros(3)
            elseif param_type == 2
                DIXref[i][1:3] = lissajous(ti, trajectory_params_2)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_2, 1)
            elseif param_type == 3
                DIXref[i][1:3] = lissajous(ti, trajectory_params_3)
                DIXref[i][4:6] = lissajous(ti, trajectory_params_3, 1)
            end

        end

        # Control Input Trajectory
        DIUref = [[0; 0; 0.0] for i = 1:(N-1)]        
        return DIXref, DIUref
    end 

end


function figureEight(t, k=6.0)  # k > 1 speeds up the trajectory

    N = 1.3
    xd = N * sin(k * t / 10) - 0.3
    yd = N * sin(k * t / 20) + 0.3
    zd = 0.
    
    # Velocity Desired
    xd_dot = k * N * cos(k * t / 10) / 10
    yd_dot = k * N * cos(k * t / 20) / 20
    zd_dot = 0.
    
    # Acceleration Desired
    xd_ddot = -k^2 * N * sin(k * t / 10) / 100
    yd_ddot = -k^2 * N * sin(k * t / 20) / 400
    zd_ddot = 0.
    
    ϕd = atan(yd_dot, xd_dot) 
    
    x_des = [xd; xd_dot; xd_ddot]
    y_des = [yd; yd_dot; yd_ddot]
    z_des = [zd; zd_dot; zd_ddot]
    
    return x_des, y_des, z_des, ϕd
    
end

function get_desired_trajectory_eight(t0, N, dt)

    # This is only for the plotting purposes thus only send the pos
    Xref = [zeros(3) for i=1:N]
    for i = 1:(N)
        ti = t0 + (i-1)*dt
        # Desired Position
        x_des, y_des, z_des, ϕd = figureEight(ti)
        Xref[i] = [x_des[1]; y_des[1]; z_des[1]]
    end

    return Xref
    
end

end