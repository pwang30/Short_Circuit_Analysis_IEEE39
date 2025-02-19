# Author: Peng Wang       from Technical University of Madrid (UPM)
# Supervisor: Luis Badesa

# Now, this is single-level+UC+SCL TEST CASE with a modified IEEE-30 bus system
# 03.Dec.2024

import Pkg 
using JuMP,Gurobi, CSV,DataFrames,LinearAlgebra, XLSX # (here we choose Gurobi as it can solve the product of binary and continuous variables, we don't need to linearlize the nonlinear terms manually)

#----------------------------------IEEE-30 Bus System Data Introduction----------------------------------
df = DataFrame(CSV.File("/Users/ME2/Desktop/scl_ieee30/Loadcurves.csv"))        
loadcurve=df[:,:]  
df = DataFrame(CSV.File( "/Users/ME2/Desktop/scl_ieee30/Windcurves.csv") ) 
windcurves=df[:,:]
df = DataFrame(CSV.File( "/Users/ME2/Desktop/scl_ieee30/SGpara.csv") ) 
SGpara=df[:,:]
df = DataFrame(CSV.File( "/Users/ME2/Desktop/scl_ieee30/Linespara.csv" )) 
linepara=df[:,:]



#-----------------------------------Define Parameters for Calculating SCC
branch_num = size(linepara, 1)                # Number of branches in the network

branch_num= size(linepara, 1)       # number of branches
numnodes=30                         # number of nodes
Yₗᵢₙₑ= zeros(numnodes, numnodes)     # define admittance matrix of the transmission lines
Yₛ₉ = zeros(numnodes, numnodes)      # define admittance matrix of the SGs

for k in 1:branch_num               #  calculate the admittance matrix of the transmission lines
    i = linepara[k, 1]                  # bus from
    j = linepara[k, 2]                   # bus to
    Yₗᵢₙₑ[i, j] = -1/linepara[k, 4]        # off-diagonal elements
    Yₗᵢₙₑ[j, i] = Yₗᵢₙₑ[i, j]                # symmetry
end
for k in 1:numnodes
    Yₗᵢₙₑ[k, k] = -sum(Yₗᵢₙₑ[k, :])           # diagonal elements 
end  

for k in 1:size(SGpara, 1)                     # calculate the admittance matrix of the SGs
    i = SGpara[k, 1]                         # bus from
    Yₛ₉[i, i] = 1/SGpara[k, 2]                # diagonal elements 
end 

I_IBG=1  # SCC (p.u): Iₛ₉=Eₛ₉/X_dg  & I_IBG (pre-defined)
Iₗᵢₘ=2    # SCC limit
β=0.95
Eₛ₉=1
Eₛ₉=β*Eₛ₉
Iₛ₉=zeros(1,size(SGpara, 1))
for k in 1:size(SGpara, 1)
    Iₛ₉[1,k]=Eₛ₉/SGpara[k, 2]
end

buses_load_operator_1=[2 3 4 5 7 12 13 14 15 16 17 18 19 20 23]
buses_load_operator_2=[8 10 21 24 26 29]



#-----------------------------------Define Parameters for Optimization
Pˢᴳₘₐₓ=[15,14,12,13,   15,13   ]
Pˢᴳₘᵢₙ=Pˢᴳₘₐₓ- 6*ones(1,size(SGpara, 1))'
Rₘₐₓ=Pˢᴳₘᵢₙ

Kᵁ=[3.25,2.72,1.43,2.03,  3.31,2.11]
Kᴰ=[0.285,0.201,0.153,0.201,  0.312,0.189]
Cᵍᵐ=[0.9,0.6,0.5,0.7,  0.8,0.6]
Cⁿˡ=[1.2,1,0.8,1.1,  1.1,0.9]
T=24


     
#-----------------------------------Define Model
model= Model()

#-------Define Variales
# For operator 1
@variable(model, Pˢᴳ²[1:T])                # generation of SG in bus 2
@variable(model, Pˢᴳ³[1:T])                # generation of SG in bus 3
@variable(model, Pˢᴳ⁴[1:T])                # generation of SG in bus 4
@variable(model, Pˢᴳ⁵[1:T])                # generation of SG in bus 5

@variable(model, yˢᴳ²[1:T],Bin)            # on/off status of SG in bus 2
@variable(model, yˢᴳ³[1:T],Bin)            # on/off status of SG in bus 3
@variable(model, yˢᴳ⁴[1:T],Bin)            # on/off status of SG in bus 2
@variable(model, yˢᴳ⁵[1:T],Bin)            # on/off status of SG in bus 3

@variable(model, Pᴵᴮᴳ¹[1:T])               # generation of IBG (WT) in bus 1
@variable(model, Pᴵᴮᴳ²³[1:T])              # generation of IBG (WT) in bus 23

@variable(model, Cᵁ²[1:T])                 # startup costs for SG in bus 2
@variable(model, Cᴰ²[1:T])                 # shutdown costs for SG in bus 2
@variable(model, Cᵁ³[1:T])                 # startup costs for SG in bus 3
@variable(model, Cᴰ³[1:T])                 # shutdown costs for SG in bus 3
@variable(model, Cᵁ⁴[1:T])                 # startup costs for SG in bus 4
@variable(model, Cᴰ⁴[1:T])                 # shutdown costs for SG in bus 4
@variable(model, Cᵁ⁵[1:T])                 # startup costs for SG in bus 5
@variable(model, Cᴰ⁵[1:T])                 # shutdown costs for SG in bus 5

# For operator 2
@variable(model, Pˢᴳ²⁷[1:T])                # generation of SG in bus 27
@variable(model, Pˢᴳ³⁰[1:T])                # generation of SG in bus 30

@variable(model, yˢᴳ²⁷[1:T],Bin)            # on/off status of SG in bus 27
@variable(model, yˢᴳ³⁰[1:T],Bin)            # on/off status of SG in bus 30

@variable(model, Pᴵᴮᴳ²⁶[1:T])               # generation of IBG (WT) in bus 26

@variable(model, Cᵁ²⁷[1:T])                 # startup costs for SG in bus 27
@variable(model, Cᴰ²⁷[1:T])                 # shutdown costs for SG in bus 27
@variable(model, Cᵁ³⁰[1:T])                 # startup costs for SG in bus 30
@variable(model, Cᴰ³⁰[1:T])                 # shutdown costs for SG in bus 30


# @variable(model, α[1:T])                   # percentage of IBG penetration  



#-------Define Constraints
Pᴰ¹=zeros(1,T)
for t in 1:T
    for k in 1:size(buses_load_operator_1, 2)
        Pᴰ¹[t]=Pᴰ¹[t]+loadcurve[buses_load_operator_1[k],t+2]
    end
end
Pᴰ²=zeros(1,T)
for t in 1:T
    for k in 1:size(buses_load_operator_2, 2)
        Pᴰ²[t]=Pᴰ²[t]+loadcurve[buses_load_operator_2[k],t+2]
    end
end

@constraint(model, Pˢᴳ²+Pˢᴳ³+Pˢᴳ⁴+Pˢᴳ⁵+Pᴵᴮᴳ¹+Pᴵᴮᴳ²³==Pᴰ¹')     # power balance for operator 1
@constraint(model, Pˢᴳ²⁷+Pˢᴳ³⁰+Pᴵᴮᴳ²⁶==Pᴰ²')                   # power balance for operator 2

@constraint(model, Pᴵᴮᴳ¹ .<= Vector(windcurves[1, 2:end]))        # wind power limit
@constraint(model, Pᴵᴮᴳ¹>=0)
@constraint(model, Pᴵᴮᴳ²³.<= Vector(windcurves[2,2:end]))       
@constraint(model, Pᴵᴮᴳ²³>=0)
@constraint(model, Pᴵᴮᴳ²⁶.<=Vector(windcurves[3,2:end]))       
@constraint(model, Pᴵᴮᴳ²⁶>=0)

# @constraint(model, α.<=1)                   # IBG online capacity limit
# @constraint(model, α>=0)
# @constraint(model, α.*Pᴰ==Pᴵᴮᴳ)

@constraint(model, Pˢᴳ¹.<=yˢᴳ¹*Pˢᴳₘₐₓ[1])    # bounds for the output of SG in bus 1 & 2
@constraint(model, yˢᴳ¹*Pˢᴳₘᵢₙ[1].<=Pˢᴳ¹)    
@constraint(model, Pˢᴳ².<=yˢᴳ²*Pˢᴳₘₐₓ[2])    
@constraint(model, yˢᴳ²*Pˢᴳₘᵢₙ[2].<=Pˢᴳ²)    

@constraint(model, Cᵁ¹>=0)      # lower bounds for the startup costs of SG in bus 1 & 2                 
@constraint(model, Cᴰ¹>=0)                        
@constraint(model, Cᵁ²>=0)                       
@constraint(model, Cᴰ²>=0)                       
for t in 1:T-1
@constraint(model, Cᵁ¹[t]>=(yˢᴳ¹[t+1]-yˢᴳ¹[t])*Kᵁ[1])        
@constraint(model, Cᴰ¹[t]>=(yˢᴳ¹[t]-yˢᴳ¹[t+1])*Kᴰ[1])        
@constraint(model, Cᵁ²[t]>=(yˢᴳ²[t+1]-yˢᴳ²[t])*Kᵁ[2])        
@constraint(model, Cᴰ²[t]>=(yˢᴳ²[t]-yˢᴳ²[t+1])*Kᴰ[2])        
end

for t in 1:T-1       # bounds for the ramp of SGs in bus 1 & 2 
    @constraint(model, Pˢᴳ¹[t+1]-Pˢᴳ¹[t]<=Rₘₐₓ[1])        
    @constraint(model, -Rₘₐₓ[1]<=Pˢᴳ¹[t+1]-Pˢᴳ¹[t])        
    @constraint(model, Pˢᴳ²[t+1]-Pˢᴳ²[t]<=Rₘₐₓ[2])        
    @constraint(model, -Rₘₐₓ[2]<=Pˢᴳ²[t+1]-Pˢᴳ²[t])         
end
         
# bounds for the SCC of SG in bus 1 
@variable(model, Z[1:T,  1:3,  1:3])  # Z₁₁, Z₁₂, Z₁₃ 

# 定义动态导纳矩阵
@variable(model, Y_total[1:T, 1:3, 1:3])
for t in 1:T
    @constraint(model, Y_total[t,:,:] .== Y_0 + Y_g .* [yˢᴳ¹[t] 0 0; 0 yˢᴳ²[t] 0; 0 0 0])
    @constraint(model, Z[t,:,:] * Y_total[t,:,:] .== [1 0 0;0 1 0;0 0 1])  # 矩阵求逆约束
end


for t in 1:T   # bounds for the SCC of SG in bus 1
    @constraint(model, Z[t,1,1]*I_SGs[1]*yˢᴳ¹[t]+Z[t,1,2]*I_SGs[2]*yˢᴳ²[t]+Z[t,1,3]*I_IBG*α[t]>=Iₗᵢₘ*Z[t,1,1])
end

#-------Define Objective Functions
No_load_cost=sum(Cⁿˡ[1].*yˢᴳ¹)+sum(Cⁿˡ[2].*yˢᴳ²)       # no-load cost
Generation_cost=sum(Cᵍᵐ[1].*Pˢᴳ¹)+sum(Cᵍᵐ[2].*Pˢᴳ²)    # generation cost
Onoff_cost=sum(Cᵁ¹)+sum(Cᴰ¹)+sum(Cᵁ²)+sum(Cᴰ²)        # on/off cost
@objective(model, Min, No_load_cost+   Generation_cost+     Onoff_cost)  # objective function


#-----------------------------------Solve and Output Results
set_optimizer(model , Gurobi.Optimizer)
# set_attribute(model, "limits/gap", 0.0280)
# set_time_limit_sec(model, 700.0)
optimize!(model)

yˢᴳ¹=JuMP.value.(yˢᴳ¹)
yˢᴳ²=JuMP.value.(yˢᴳ²)
# α=JuMP.value.(α)
#Z=JuMP.value.(Z)
#println(value.( Z[1,:,:]))