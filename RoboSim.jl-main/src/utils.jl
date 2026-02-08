module QuaternionUtils

using LinearAlgebra
using BlockDiagonals
using StaticArrays
using ComponentArrays


SV3{F} = SVector{3,F}
SV4{F} = SVector{4,F}
SM3{F} = SMatrix{3,3,F,9}
SM4{F} = SMatrix{4,4,F,16}

# Skew symmetric matrix
function hat(v)
    return [0 -v[3] v[2];
            v[3] 0 -v[1];
            -v[2] v[1] 0]
end

function hat(v::SV3{F}) where {F}
    return @SMatrix [
        [zero(F);; -v[3];; v[2]]
        [v[3];; zero(F);; -v[1]]
        [-v[2];; v[1];; zero(F)]
    ]
end

function vee(M)
    return [
        (M[3, 2] - M[2, 3]) / 2
        (M[1, 3] - M[3, 1]) / 2
        (M[2, 1] - M[1, 2]) / 2
    ]
end

signedsq(a) = a * abs(a)
signedsqrt(a) = sign(a) * sqrt(abs(a))

function vee(M::SM3)
    return @SVector [
        (M[3, 2] - M[2, 3]) / 2
        (M[1, 3] - M[3, 1]) / 2
        (M[2, 1] - M[1, 2]) / 2
    ]
end



function L(q)
    s = q[1]
    v = q[2:4]
    L_ = [s    -v';
         v  s*I+hat(v)]
    return L_
end

const T = Diagonal([1; -ones(3)])
const H = [zeros(1,3); I]

function qtoQ(q)
    # T = Diagonal([1; -ones(3)])
    # H = [zeros(1,3); I]
    return H'*T*L(q)*T*L(q)*H
end

function G(q)
    H = [zeros(1,3); I]
    G = L(q)*H
end

function rptoq(ϕ)
    (1/sqrt(1+ϕ'*ϕ))*[1; ϕ]
end

function qtorp(q)
    q[2:4]/q[1]
end

function E(q)
    E = BlockDiagonal([1.0*I(3), G(q), 1.0*I(6)])
end


# Rotation around Yaw
function Rot_Yaw(ϕ)
    Rz = [[cos(ϕ) ;; -sin(ϕ) ;; 0];
            [sin(ϕ) ;; cos(ϕ) ;; 0];
             [0   ;;  0  ;; 1]]
    return Rz
end

end



