# functions to create spin operators


ξop=[0 0.5 ; 0.5 0]
δnop=diagm([-0.5,0.5])
idop=diagm([1.0,1.0])


function promoteOp(op,idx::Integer,num::Integer)
    ops=[idop for _ in 1:num]
    ops[idx]=op
    kron(ops...)
end

function promoteOp(op,idx::Vector,num::Integer)
    ops=[idop for _ in 1:num]
    for (pos_,idx_ ) in enumerate(idx)
        ops[idx_]=op[pos_]
    end    
    kron(ops...)
end


function spin_orbital_idx(orb_idx,spin_idx)
    (orb_idx-1)*2+spin_idx
end

# generate interactions
function cal_O1(n_orb)
    n_spin_orb=2*n_orb
    result=zeros(2^n_spin_orb,2^n_spin_orb)
    for α in 1:n_orb
        αup=spin_orbital_idx(α,1)
        αdn=spin_orbital_idx(α,2)
        result+=promoteOp([δnop,δnop],[αup,αdn],n_spin_orb)
    end
    result
end


function cal_O2(n_orb)
    n_spin_orb=2*n_orb
    result=zeros(2^n_spin_orb,2^n_spin_orb)
    for α in 1:n_orb
        for β in (α+1):n_orb
            for σ in [1,2]
                ασ=spin_orbital_idx(α,σ)
                βσbar=spin_orbital_idx(β,3-σ)
                result+=promoteOp([δnop,δnop],[ασ,βσbar],n_spin_orb)
            end            
        end        
    end
    result
end

function cal_O3(n_orb)
    n_spin_orb=2*n_orb
    result=zeros(2^n_spin_orb,2^n_spin_orb)
    for α in 1:n_orb
        for β in (α+1):n_orb
            for σ in [1,2]
                ασ=spin_orbital_idx(α,σ)
                βσ=spin_orbital_idx(β,σ)
                result+=promoteOp([δnop,δnop],[ασ,βσ],n_spin_orb)
            end            
        end        
    end
    result
end


# notice this is the summation of all ξ_ℓ
function cal_ξ(n_orb)
    n_spin_orb=2*n_orb
    result=zeros(2^n_spin_orb,2^n_spin_orb)
    for ℓ in 1:n_spin_orb
        result+=promoteOp(ξop,ℓ,n_spin_orb)
    end
    result
end


function cal_O(r,n_orb)
    cal_O1(n_orb)+(1-2*r)*cal_O2(n_orb)+(1-3*r)*cal_O3(n_orb)
end

# we use the lanczos to find the minimal ground state
# we change the original version to use lanczos
# θ=0
# n_orb=3
# r=0.1
# we also change code so the ξ is the average value
function cal_ξ_O(r,n_orb;show_progress=false,Nθ=1000)
    ξmat=cal_ξ(n_orb)
    Omat=cal_O(r,n_orb)
    θs=linspace(0,pi/2,Nθ)
    data=[]
    gs_val=0.0
    gs_vec=[1/(2)^n_orb for _ in 1:2^(2*n_orb)]
    for θ in θs
        H=-cos(θ)*ξmat+sin(θ)*Omat
        gs_val,gs_vec=eigs(H;nev=1,which=:SR,v0=reshape(gs_vec,:))
        data_=[θ,dot(gs_vec,ξmat*gs_vec)/2.0/n_orb,dot(gs_vec,Omat*gs_vec)]
        if(show_progress)
            print("r $r n_orb $(n_orb) θ $θ data_ $(data_)\n")
        end
        push!(data,data_)
    end
    data
end

function O_atomic(norb,r)
    0.25*norb*(-r*(norb-1)-1)
end

