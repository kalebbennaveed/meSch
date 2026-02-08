module Visualization

import MeshCat as mc
using TrajOptPlots
using StaticArrays
using LinearAlgebra

using ColorTypes
using GeometryBasics
using ColorSchemes
using Colors
cmap = ColorSchemes.lajolla
using ComponentArrays
using Parameters: @unpack


include("utils.jl")

# To plot the reference trajectory
# Create a cylinder joining two points of the reference trajectory
# Nom: (3,3)
function vis_traj!(vis, num, name, X; R = 0.2, color = mc.RGBA(1.0, 0.0, 0.0, 2.0))
    if num == 1
        for i = 1:3:(length(X)-1)
            a = X[i][1:3]
            b = X[i+1][1:3]
            cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
            mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color) )
        end
    elseif num == 2
        for i = 1:length(X)
            a = X[i][1:3]
            sph = mc.HyperSphere(mc.Point(a...), R)
            mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color) )
        end       
    else
        for i = 1:(length(X)-1)
            a = X[i][1:3]
            b = X[i+1][1:3]
            cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
            mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color) )
        end  
    end      

    # for i = 1:3:length(X)
    #     a = X[i][1:3]
    #     sph = mc.HyperSphere(mc.Point(a...), R)
    #     mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color) )
    # end
end

# function bound!(vis, num, name, X; R = 0.2, color = mc.RGBA(1.0, 0.0, 0.0, 2.0))
#     if num == 1
#         for i = 1:(length(X)-1)
#             a = X[i][1:3]
#             b = X[i+1][1:3]
#             cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
#             mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color) )
#         end
#     elseif num == 2
#         for i = 1:length(X)
#             a = X[i][1:3]
#             sph = mc.HyperSphere(mc.Point(a...), R)
#             mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color) )
#         end       
#     else
#         for i = 1:(length(X)-1)
#             a = X[i][1:3]
#             b = X[i+1][1:3]
#             cyl = mc.Cylinder(mc.Point(a...), mc.Point(b...), R)
#             mc.setobject!(vis[name]["p"*string(i)], cyl, mc.MeshPhongMaterial(color=color) )
#         end  
#     end      

#     # for i = 1:3:length(X)
#     #     a = X[i][1:3]
#     #     sph = mc.HyperSphere(mc.Point(a...), R)
#     #     mc.setobject!(vis[name]["s"*string(i)], sph, mc.MeshPhongMaterial(color=color) )
#     # end
# end





function animate_rover!(vis, Xsim, dt; R_sphere = 0.3)
    
    # Currently only robots; But this eventually will be used to animate anyhing (not necessarily robots)
    
    # Rover object 
    rover_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "rover_v5.obj"))   
    # mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.41,0.41,0.41,1.0)))
    mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.61,0.61,0.61,1.0)))

    # # Draw a charging sphere around the rover
    # pos = Xsim[1][1:3] 
    # sphere = mc.HyperSphere(mc.Point(pos...), R_sphere)
    # mc.setobject!(vis[:charge_sphere], sphere, mc.MeshPhongMaterial(color = mc.RGBA(0.0, 1.0, 0.0, 0.3)))   

    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt))
        
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            
            # Rover 
            r = Xsim[k][1:3]
            ϕ = Xsim[k][4]
            dcm = QuaternionUtils.Rot_Yaw(ϕ)
            mc.settransform!(vis[:rover], mc.compose(mc.Translation(r),mc.LinearMap(dcm)))

            # r_sphere = [r[1:2]; 0.4]
            # mc.settransform!(vis[:charge_sphere], mc.compose(mc.Translation(r_sphere),mc.LinearMap(dcm)))
        end
    end
    mc.setanimation!(vis, anim)

#     return (mc.render(vis))
end

function vis_charging_station!(vis, name, pos; R = 0.3, color = mc.RGBA(0.0, 1.0, 0.0, 0.3))
    sphere = mc.HyperSphere(mc.Point(pos...), R)
    mc.setobject!(vis[name], sphere, mc.MeshPhongMaterial(color=color))

end

function animate_system_2!(vis, Xsim, dt)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    # Quad 1 object
    quad1_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled.obj"))
    mc.setobject!(vis[:quad1], quad1_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

    # Quad 2 object 
    # grey color: (0.41,0.41,0.41,1.0)
    quad2_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled.obj"))   
    mc.setobject!(vis[:quad2], quad2_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
    
    # Set the quad 1 to the initial position
    r_init_q1 = [4., 4., 0]
    r_init_q2 = [-4., -4., 0]
    q_init = [1.0; 0.0; 0; 0]

    dcm_init = QuaternionUtils.qtoQ(q_init) # Convert quaternion to dcm    
    mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_init_q1),mc.LinearMap(dcm_init)))
    mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_init_q2),mc.LinearMap(dcm_init)))
    
    
    dcms_q1 = [begin
            q_q1 = Xsim[k][4:7]
            dcm_q1 = QuaternionUtils.qtoQ(q_q1) 
        end
            for k = 1:length(Xsim)
        ]

    dcms_q2 = [begin
            q_q2 = Xsim[k][18:21]
            dcm_q2 = QuaternionUtils.qtoQ(q_q2) 
        end
            for k = 1:length(Xsim)
        ]    
    
    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N

            r_q1 = Xsim[k][1:3]
            r_q2 = Xsim[k][15:17]

            # Set transform for Quad1 and Quad2
            mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_q1),mc.LinearMap(dcms_q1[k]))) # Quad 1
            mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_q2),mc.LinearMap(dcms_q2[k]))) # Quad 2
            
        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end


function animate_system_3!(vis, Xsim, dt)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    # Quad 1 object
    quad1_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled.obj"))
    mc.setobject!(vis[:quad1], quad1_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

    # Quad 2 object 
    # grey color: (0.41,0.41,0.41,1.0)
    quad2_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled.obj"))   
    mc.setobject!(vis[:quad2], quad2_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

    # Quad 3 object 
    # grey color: (0.41,0.41,0.41,1.0)
    quad3_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled.obj"))   
    mc.setobject!(vis[:quad3], quad3_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
    
    # Set the quad 1 to the initial position
    r_init_q1 = [4., 4., 0]
    r_init_q2 = [-4., -4., 0]
    r_init_q3 = [-4.0, 4, 0]
    q_init = [1.0; 0.0; 0; 0]

    dcm_init = QuaternionUtils.qtoQ(q_init) # Convert quaternion to dcm    
    mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_init_q1),mc.LinearMap(dcm_init)))
    mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_init_q2),mc.LinearMap(dcm_init)))
    mc.settransform!(vis[:quad3], mc.compose(mc.Translation(r_init_q3),mc.LinearMap(dcm_init)))    
    
    dcms_q1 = [begin
            q_q1 = Xsim[k][4:7]
            dcm_q1 = QuaternionUtils.qtoQ(q_q1) 
        end
            for k = 1:length(Xsim)
        ]

    dcms_q2 = [begin
            q_q2 = Xsim[k][18:21]
            dcm_q2 = QuaternionUtils.qtoQ(q_q2) 
        end
            for k = 1:length(Xsim)
        ] 

    dcms_q3 = [begin
            q_q3 = Xsim[k][32:35]
            dcm_q3 = QuaternionUtils.qtoQ(q_q3) 
        end
            for k = 1:length(Xsim)
        ] 
    
    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N

            r_q1 = Xsim[k][1:3]
            r_q2 = Xsim[k][15:17]
            r_q3 = Xsim[k][29:31]

            # Set transform for Quad1 and Quad2
            mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_q1),mc.LinearMap(dcms_q1[k]))) # Quad 1
            mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_q2),mc.LinearMap(dcms_q2[k]))) # Quad 2
            mc.settransform!(vis[:quad3], mc.compose(mc.Translation(r_q3),mc.LinearMap(dcms_q3[k]))) # Quad 2

        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end

function animate_system_4!(vis, Xsim, dt)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    # Quad 1 object
    quad1_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    mc.setobject!(vis[:quad1], quad1_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))


    # Quad 2 object 
    # grey color: (0.41,0.41,0.41,1.0)
    quad2_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
    mc.setobject!(vis[:quad2], quad2_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))


    # Quad 3 object 
    # grey color: (0.41,0.41,0.41,1.0)
    quad3_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
    mc.setobject!(vis[:quad3], quad3_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

    # Quad 4 object 
    # grey color: (0.41,0.41,0.41,1.0)
    quad4_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
    mc.setobject!(vis[:quad4], quad4_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
    
    

    # Set the quad 1 to the initial position
    r_init_q1 = Xsim[1][1:3]
    r_init_q2 = Xsim[1][15:17]
    r_init_q3 = Xsim[1][29:31]
    r_init_q4 = Xsim[1][43:45]
    q_init = [1.0; 0.0; 0; 0]
    
    dcm_init = QuaternionUtils.qtoQ(q_init) # Convert quaternion to dcm    
    mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_init_q1),mc.LinearMap(dcm_init)))
    mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_init_q2),mc.LinearMap(dcm_init)))
    mc.settransform!(vis[:quad3], mc.compose(mc.Translation(r_init_q3),mc.LinearMap(dcm_init)))    
    mc.settransform!(vis[:quad4], mc.compose(mc.Translation(r_init_q4),mc.LinearMap(dcm_init))) 
    
    
    dcms_q1 = [begin
            q_q1 = Xsim[k][4:7]
            dcm_q1 = QuaternionUtils.qtoQ(q_q1) 
        end
            for k = 1:length(Xsim)
        ]

    dcms_q2 = [begin
            q_q2 = Xsim[k][18:21]
            dcm_q2 = QuaternionUtils.qtoQ(q_q2) 
        end
            for k = 1:length(Xsim)
        ] 

    dcms_q3 = [begin
            q_q3 = Xsim[k][32:35]
            dcm_q3 = QuaternionUtils.qtoQ(q_q3) 
        end
            for k = 1:length(Xsim)
        ] 

    dcms_q4 = [begin
        q_q4 = Xsim[k][46:49]
        dcm_q4 = QuaternionUtils.qtoQ(q_q4) 
    end
        for k = 1:length(Xsim)
    ] 
    
    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N

            r_q1 = Xsim[k][1:3]
            r_q2 = Xsim[k][15:17]
            r_q3 = Xsim[k][29:31]
            r_q4 = Xsim[k][43:45]

            # Set transform for Quad1 and Quad2
            mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_q1),mc.LinearMap(dcms_q1[k]))) # Quad 1
            mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_q2),mc.LinearMap(dcms_q2[k]))) # Quad 2
            mc.settransform!(vis[:quad3], mc.compose(mc.Translation(r_q3),mc.LinearMap(dcms_q3[k]))) # Quad 3
            mc.settransform!(vis[:quad4], mc.compose(mc.Translation(r_q4),mc.LinearMap(dcms_q4[k]))) # Quad 4

        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end



# function animate_system_N!(vis, Xsim, N, dt)
    
#     # vis is an instance of the meshcat visualizer
#     # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
#     # Quad 1 object
#     quad1_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
#     mc.setobject!(vis[:quad1], quad1_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 2 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad2_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad2], quad2_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 3 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad3_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad3], quad3_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 4 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad4_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad4], quad4_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 5 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad5_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad5], quad5_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 6 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad6_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad6], quad6_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0))) 
    
#     # Quad 7 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad7_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad7], quad7_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 8 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad8_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad8], quad8_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 9 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad9_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad9], quad9_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 10 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad10_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad10], quad10_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 11 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad11_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad11], quad11_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 12 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad12_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad12], quad12_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 12 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad13_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad13], quad13_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 12 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad14_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad14], quad14_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 12 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad15_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad15], quad15_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))

#     # Quad 12 object 
#     # grey color: (0.41,0.41,0.41,1.0)
#     quad16_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))   
#     mc.setobject!(vis[:quad16], quad16_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
    
#     # Set the quad 1 to the initial position
#     r_init_q1 = Xsim[1][1:3]
#     r_init_q2 = Xsim[1][15:17]
#     r_init_q3 = Xsim[1][29:31]
#     r_init_q4 = Xsim[1][43:45]
#     r_init_q5 = Xsim[1][57:59]
#     r_init_q6 = Xsim[1][71:73]
#     r_init_q7 = Xsim[1][85:87]
#     r_init_q8 = Xsim[1][99:101]
#     r_init_q9 = Xsim[1][113:115]
#     r_init_q10 = Xsim[1][127:129]
#     r_init_q11 = Xsim[1][141:143]
#     r_init_q12 = Xsim[1][155:157]
#     r_init_q13 = Xsim[1][169:171]
#     r_init_q14 = Xsim[1][183:185]
#     r_init_q15 = Xsim[1][197:199]
#     r_init_q16 = Xsim[1][211:213]
#     q_init = [1.0; 0.0; 0; 0]

#     dcm_init = QuaternionUtils.qtoQ(q_init) # Convert quaternion to dcm    
#     mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_init_q1),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_init_q2),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad3], mc.compose(mc.Translation(r_init_q3),mc.LinearMap(dcm_init)))    
#     mc.settransform!(vis[:quad4], mc.compose(mc.Translation(r_init_q4),mc.LinearMap(dcm_init))) 
#     mc.settransform!(vis[:quad5], mc.compose(mc.Translation(r_init_q5),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad6], mc.compose(mc.Translation(r_init_q6),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad7], mc.compose(mc.Translation(r_init_q7),mc.LinearMap(dcm_init)))    
#     mc.settransform!(vis[:quad8], mc.compose(mc.Translation(r_init_q8),mc.LinearMap(dcm_init))) 
#     mc.settransform!(vis[:quad9], mc.compose(mc.Translation(r_init_q9),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad10], mc.compose(mc.Translation(r_init_q10),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad11], mc.compose(mc.Translation(r_init_q11),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad12], mc.compose(mc.Translation(r_init_q12),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad13], mc.compose(mc.Translation(r_init_q13),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad14], mc.compose(mc.Translation(r_init_q14),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad15], mc.compose(mc.Translation(r_init_q15),mc.LinearMap(dcm_init)))
#     mc.settransform!(vis[:quad16], mc.compose(mc.Translation(r_init_q16),mc.LinearMap(dcm_init)))

#     dcms = Vector{Matrix{Float64}}[]
#     for i = 1:N
#         start_idx = 14*(i - 1) + 1
#         end_idx = 14*i
#         dcms_i = [begin
#                 q_qi = Xsim[k][(start_idx + 3):(start_idx + 6)]
#                 dcms_i = QuaternionUtils.qtoQ(q_qi) 
#             end
#                 for k = 1:length(Xsim)
#             ] 
        
#         push!(dcms, dcms_i)
#     end


#     # Set Animation 
#     anim = mc.Animation(floor(Int,1/dt)) 
    
#     for k = 1:length(Xsim)
#         mc.atframe(anim, k) do
 
#             ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N

#             r_q1 = Xsim[k][1:3]
#             r_q2 = Xsim[k][15:17]
#             r_q3 = Xsim[k][29:31]
#             r_q4 = Xsim[k][43:45]
#             r_q5 = Xsim[k][57:59]
#             r_q6 = Xsim[k][71:73]
#             r_q7 = Xsim[k][85:87]
#             r_q8 = Xsim[k][99:101]
#             r_q9 = Xsim[k][113:115]
#             r_q10 = Xsim[k][127:129]
#             r_q11 = Xsim[k][141:143]
#             r_q12 = Xsim[k][155:157]
#             r_q13 = Xsim[k][169:171]
#             r_q14 = Xsim[k][183:185]
#             r_q15 = Xsim[k][197:199]
#             r_q16 = Xsim[k][211:213]

#             # Set transform for Quad1 and Quad2
#             mc.settransform!(vis[:quad1], mc.compose(mc.Translation(r_q1),mc.LinearMap(dcms[1][k]))) # Quad 1
#             mc.settransform!(vis[:quad2], mc.compose(mc.Translation(r_q2),mc.LinearMap(dcms[2][k]))) # Quad 2
#             mc.settransform!(vis[:quad3], mc.compose(mc.Translation(r_q3),mc.LinearMap(dcms[3][k]))) # Quad 3
#             mc.settransform!(vis[:quad4], mc.compose(mc.Translation(r_q4),mc.LinearMap(dcms[4][k]))) # Quad 4
#             mc.settransform!(vis[:quad5], mc.compose(mc.Translation(r_q5),mc.LinearMap(dcms[5][k]))) # Quad 5
#             mc.settransform!(vis[:quad6], mc.compose(mc.Translation(r_q6),mc.LinearMap(dcms[6][k]))) # Quad 6
#             mc.settransform!(vis[:quad7], mc.compose(mc.Translation(r_q7),mc.LinearMap(dcms[7][k]))) # Quad 7
#             mc.settransform!(vis[:quad8], mc.compose(mc.Translation(r_q8),mc.LinearMap(dcms[8][k]))) # Quad 8
#             mc.settransform!(vis[:quad9], mc.compose(mc.Translation(r_q9),mc.LinearMap(dcms[9][k]))) # Quad 9
#             mc.settransform!(vis[:quad10], mc.compose(mc.Translation(r_q10),mc.LinearMap(dcms[10][k]))) # Quad 10
#             mc.settransform!(vis[:quad11], mc.compose(mc.Translation(r_q11),mc.LinearMap(dcms[11][k]))) # Quad 11
#             mc.settransform!(vis[:quad12], mc.compose(mc.Translation(r_q12),mc.LinearMap(dcms[12][k]))) # Quad 12

#             mc.settransform!(vis[:quad13], mc.compose(mc.Translation(r_q13),mc.LinearMap(dcms[13][k]))) # Quad 13
#             mc.settransform!(vis[:quad14], mc.compose(mc.Translation(r_q14),mc.LinearMap(dcms[14][k]))) # Quad 14
#             mc.settransform!(vis[:quad15], mc.compose(mc.Translation(r_q15),mc.LinearMap(dcms[15][k]))) # Quad 15
#             mc.settransform!(vis[:quad16], mc.compose(mc.Translation(r_q16),mc.LinearMap(dcms[16][k]))) # Quad 16

#         end
#     end
    
#     # Sets animation 
#     mc.setanimation!(vis, anim)

# end

function animate_system_N!(vis, Xsim, N, dt)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    q_init = [1.0; 0.0; 0; 0]
    dcm_init = QuaternionUtils.qtoQ(q_init)

    # Set quadrotors
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        r_init_i = Xsim[1][start_idx:(start_idx + 2)]
        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_init_i),mc.LinearMap(dcm_init)))
    end

    # Set the rover


    dcms = Vector{Matrix{Float64}}[]
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        end_idx = 14*i
        dcms_i = [begin
                q_qi = Xsim[k][(start_idx + 3):(start_idx + 6)]
                dcms_i = QuaternionUtils.qtoQ(q_qi) 
            end
                for k = 1:length(Xsim)
            ] 
        push!(dcms, dcms_i)
    end


    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N


            for i = 1:N
                start_idx = 14*(i - 1) + 1
                r_qi = Xsim[k][start_idx:(start_idx + 2)]
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_qi),mc.LinearMap(dcms[i][k]))) # Quad 1
            end
        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end


# Visualize the quadrotors with the communication range
function animate_system_geometric!(vis, Xsim, N, dt; R_sphere = 1.0)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    q_init = [1.0; 0.0; 0; 0]
    dcm_init = QuaternionUtils.qtoQ(q_init)

    # Set quadrotors
    for i = 1:N

        @unpack x, R = Xsim[1]
        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(R)))

        # Draw a charging sphere around the quads
        r_sphere = [0; 0; 0.] # Not sure why to do this
        sphere = mc.HyperSphere(mc.Point(r_sphere...), R_sphere)
        mc.setobject!(vis[:charge_sphere][string(i)], sphere, mc.MeshPhongMaterial(color = mc.RGBA(0.0, 1.0, 0.0, 0.1)))
        mc.settransform!(vis[:charge_sphere][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(R)))
        mc.setvisible!(vis[:charge_sphere][string(i)], false)
    end



    dcms = Vector{Matrix{Float64}}[]
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        end_idx = 14*i
        dcms_i = [begin
                @unpack R = Xsim[k]
                dcms_i = R
            end
                for k = 1:length(Xsim)
            ] 
        push!(dcms, dcms_i)
    end


    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N

            for i = 1:N
                @unpack x = Xsim[k]
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(dcms[i][k]))) # Quad 1
                mc.settransform!(vis[:charge_sphere][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(dcms[i][k])))

            end
        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end

function N_animate_system_geometric!(vis, Xsim, N, dt; R_sphere = 1.0)

    # vis is an instance of the meshcat visualizer
    # Xsim is a ComponentVector containing state history for multiple quadrotors
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    quad_keys = keys(Xsim[1])
    # Set quadrotors in the visualizer
    for i = 1:N
        @unpack x, R = Xsim[1][quad_keys[i]]  # Extract initial state for each quad

        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(R)))

    end

    
    # Store rotation matrices for all quads over time
    dcms = [Matrix{Float64}[] for _ in 1:N]

    for i = 1:N
        dcms[i] = [Xsim[k][quad_keys[i]].R for k in 1:length(Xsim)]
    end

    # # Set Animation
    anim = mc.Animation(floor(Int, 1/dt))

    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            for i = 1:N
                @unpack x = Xsim[k][quad_keys[i]]
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(dcms[i][k]))) # Quad 1
            end
        end
    end

    # Sets animation
    mc.setanimation!(vis, anim)
end


function N_animate_system_geometric_rover!(vis, Xsim, N, dt; R_sphere = 1.0)

    # vis is an instance of the meshcat visualizer
    # Xsim is a ComponentVector containing state history for multiple quadrotors
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    quad_keys = keys(Xsim[1])
    # Set quadrotors in the visualizer
    for i = 1:N
        @unpack x, R = Xsim[1][quad_keys[i]]  # Extract initial state for each quad

        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(R)))

    end

    # Set the rover object
    rover_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "rover_v5.obj"))   
    mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.45,0.45,0.45,1.0)))
    
    # Store rotation matrices for all quads over time
    dcms = [Matrix{Float64}[] for _ in 1:N]

    for i = 1:N
        dcms[i] = [Xsim[k][quad_keys[i]].R for k in 1:length(Xsim)]
    end

    # # Set Animation
    anim = mc.Animation(floor(Int, 1/dt))

    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
            for i = 1:N
                @unpack x = Xsim[k][quad_keys[i]]
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(x),mc.LinearMap(dcms[i][k]))) # Quad 1

                # Rover 
                @unpack pos, phi = Xsim[k][quad_keys[end]]
                dcm = QuaternionUtils.Rot_Yaw(phi)
                mc.settransform!(vis[:rover], mc.compose(mc.Translation(pos),mc.LinearMap(dcm)))
            end
        end
    end

    # Sets animation
    mc.setanimation!(vis, anim)
end




# Visualize the quadrotors with the communication range
function animate_system_N_comms_range!(vis, Xsim, N, dt; R_sphere = 1.0)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    q_init = [1.0; 0.0; 0; 0]
    dcm_init = QuaternionUtils.qtoQ(q_init)

    # Set quadrotors
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        r_init_i = Xsim[1][start_idx:(start_idx + 2)]
        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_init_i),mc.LinearMap(dcm_init)))

        # Draw a charging sphere around the quads
        r_sphere = [0; 0; 0.] # Not sure why to do this
        sphere = mc.HyperSphere(mc.Point(r_sphere...), R_sphere)
        mc.setobject!(vis[:charge_sphere][string(i)], sphere, mc.MeshPhongMaterial(color = mc.RGBA(0.0, 1.0, 0.0, 0.1)))
        mc.settransform!(vis[:charge_sphere][string(i)], mc.compose(mc.Translation(r_init_i),mc.LinearMap(dcm_init)))
        mc.setvisible!(vis[:charge_sphere][string(i)], false)
    end



    dcms = Vector{Matrix{Float64}}[]
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        end_idx = 14*i
        dcms_i = [begin
                q_qi = Xsim[k][(start_idx + 3):(start_idx + 6)]
                dcms_i = QuaternionUtils.qtoQ(q_qi) 
            end
                for k = 1:length(Xsim)
            ] 
        push!(dcms, dcms_i)
    end


    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N

            for i = 1:N
                start_idx = 14*(i - 1) + 1
                r_qi = Xsim[k][start_idx:(start_idx + 2)]
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_qi),mc.LinearMap(dcms[i][k]))) # Quad 1
                mc.settransform!(vis[:charge_sphere][string(i)], mc.compose(mc.Translation(r_qi),mc.LinearMap(dcms[i][k])))
                # Check if the robots are within the commns range
                # if i == 1
                #     r_q1 = r_qi
                #     start2_idx = 14*(2 - 1) + 1
                #     r_q2 = Xsim[k][start2_idx:(start2_idx + 2)]                    
                #     if norm(r_q1 - r_q2) <= 1.5
                #         mc.setvisible!(vis[:charge_sphere][string(1)], true)
                #         mc.setvisible!(vis[:charge_sphere][string(2)], true)
                #     else
                #         mc.setvisible!(vis[:charge_sphere][string(1)], false)
                #         mc.setvisible!(vis[:charge_sphere][string(2)], false)
                #     end                
                # end
            end
        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end




function animate_Nquads_rover!(vis, Xsim, N, dt; R_sphere = 0.3)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    q_init = [1.0; 0.0; 0; 0]
    dcm_init = QuaternionUtils.qtoQ(q_init)

    # Set quadrotors
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        r_init_i = Xsim[1][start_idx:(start_idx + 2)]
        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_init_i),mc.LinearMap(dcm_init)))
    end

    # Set the rover
    # Rover object 
    rover_start_idx = N*14 + 1
    rover_init_pos = Xsim[1][rover_start_idx : rover_start_idx + 2]

    rover_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "rover_v5.obj"))   
    # mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.41,0.41,0.41,1.0)))
    mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.45,0.45,0.45,1.0)))

    # Draw a charging sphere around the rover
    sphere = mc.HyperSphere(mc.Point(rover_init_pos...), R_sphere)
    mc.setobject!(vis[:charge_sphere], sphere, mc.MeshPhongMaterial(color = mc.RGBA(0.0, 1.0, 0.0, 0.3)))   
    
    
    dcms = Vector{Matrix{Float64}}[]
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        end_idx = 14*i
        dcms_i = [begin
                q_qi = Xsim[k][(start_idx + 3):(start_idx + 6)]
                dcms_i = QuaternionUtils.qtoQ(q_qi) 
            end
                for k = 1:length(Xsim)
            ] 
        push!(dcms, dcms_i)
    end


    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N


            for i = 1:N
                start_idx = 14*(i - 1) + 1
                r_qi = Xsim[k][start_idx:(start_idx + 2)]
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_qi),mc.LinearMap(dcms[i][k]))) # Quad 1
            end

            # Rover 
            r = Xsim[k][rover_start_idx : rover_start_idx + 2]
            ϕ = Xsim[k][rover_start_idx + 3]
            dcm = QuaternionUtils.Rot_Yaw(ϕ)
            mc.settransform!(vis[:rover], mc.compose(mc.Translation(r),mc.LinearMap(dcm)))

            r_sphere = [r[1:2]; 0.4]
            mc.settransform!(vis[:charge_sphere], mc.compose(mc.Translation(r_sphere),mc.LinearMap(dcm)))
        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end




##################3
function animate_Nquads_rover_pc!(vis, Xsim, N, dt, point_cloud, colors; R_sphere = 0.3)
    
    # vis is an instance of the meshcat visualizer
    # Xsim is the state history of the system [Quadrotor1; Battery; Quadrotor2; Battery2]
    
    quad_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "quadrotor_scaled_smaller.obj"))
    q_init = [1.0; 0.0; 0; 0]
    dcm_init = QuaternionUtils.qtoQ(q_init)


    # Cone Object 
    cone_radius = 0.45

    # Set quadrotors
    
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        r_init_i = Xsim[1][start_idx:(start_idx + 2)]
        mc.setobject!(vis[:quad][string(i)], quad_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.0,0.0,0.0,1.0)))
        mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_init_i),mc.LinearMap(dcm_init)))
        cone_i = mc.Cone(mc.Point(r_init_i[1], r_init_i[2], -0.2), mc.Point(r_init_i[1], r_init_i[2], 0.2), cone_radius)
        mc.setobject!(vis[:cone][string(i)], cone_i, mc.MeshPhongMaterial(color= mc.RGBA(1.0, 0.7431372549, 0., 0.4)))
    end


    # Set the rover
    # Rover object 
    rover_start_idx = N*14 + 1
    rover_init_pos = Xsim[1][rover_start_idx : rover_start_idx + 2]

    rover_obj = mc.MeshFileGeometry(joinpath(@__DIR__, "utils", "rover_v5.obj"))   
    # mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.41,0.41,0.41,1.0)))
    mc.setobject!(vis[:rover], rover_obj, mc.MeshPhongMaterial(color = mc.RGBA(0.45,0.45,0.45,1.0)))

    # Draw a charging sphere around the rover
    sphere = mc.HyperSphere(mc.Point(rover_init_pos...), R_sphere)
    mc.setobject!(vis[:charge_sphere], sphere, mc.MeshPhongMaterial(color = mc.RGBA(0.0, 1.0, 0.0, 0.3)))   


    # Point Cloud
    mc.setobject!(vis[:pointcloud1], mc.PointCloud(point_cloud, colors) )
    mc.setvisible!(vis[:pointcloud1], true)


    
    
    dcms = Vector{Matrix{Float64}}[]
    for i = 1:N
        start_idx = 14*(i - 1) + 1
        end_idx = 14*i
        dcms_i = [begin
                q_qi = Xsim[k][(start_idx + 3):(start_idx + 6)]
                dcms_i = QuaternionUtils.qtoQ(q_qi) 
            end
                for k = 1:length(Xsim)
            ] 
        push!(dcms, dcms_i)
    end


    # Set Animation 
    anim = mc.Animation(floor(Int,1/dt)) 
    
    for k = 1:length(Xsim)
        mc.atframe(anim, k) do
 
            ## issue: If I access the QuaternionUtils.qtoQ inside this loop it says generator function error: N


            for i = 1:N
                start_idx = 14*(i - 1) + 1
                r_qi = Xsim[k][start_idx:(start_idx + 2)]

                # For cone
                cone_pos = deepcopy(r_qi)
                # @show cone_pos

                # I have no idea why I have to do this; For now I am doing quick fix; look into this later
                if i == 1
                    cone_pos[1] -= 2.0
                    cone_pos[2] -= 2.0
                elseif i == 2
                    cone_pos[1] += 2.0
                    cone_pos[2] += 2.0
                elseif i == 3
                    cone_pos[1] += 2.0
                    cone_pos[2] -= 2.0
                else
                    cone_pos[1] -= 2.0
                    cone_pos[2] += 2.0                    
                end
                cone_pos[3] -= 0.5 * cone_pos[3]

                # Set transform for Quad 
                mc.settransform!(vis[:quad][string(i)], mc.compose(mc.Translation(r_qi),mc.LinearMap(dcms[i][k]))) # Quad 1

                # Set transform for Cone
                cone_scaling_diag = SVector(1.0, 1.0, 5.0*cone_pos[3])
                # @show cone_pos
                cone_transformation = mc.compose(mc.Translation(cone_pos), mc.LinearMap(Diagonal(cone_scaling_diag)) )
                mc.settransform!(vis[:cone][string(i)], cone_transformation)
            end

            # Rover 
            r = Xsim[k][rover_start_idx : rover_start_idx + 2]
            ϕ = Xsim[k][rover_start_idx + 3]
            dcm = QuaternionUtils.Rot_Yaw(ϕ)
            mc.settransform!(vis[:rover], mc.compose(mc.Translation(r),mc.LinearMap(dcm)))

            r_sphere = [r[1:2]; 0.4]
            mc.settransform!(vis[:charge_sphere], mc.compose(mc.Translation(r_sphere),mc.LinearMap(dcm)))
        end
    end
    
    # Sets animation 
    mc.setanimation!(vis, anim)

end
##########################3




# end of the module
end