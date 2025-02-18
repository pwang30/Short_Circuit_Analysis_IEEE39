# （1） single-level+UC+SCL TEST CASE
# （2） bilevel+UC+SCL      TEST CASE
# （3） bilevel+UC+SCL in an IEEE system  
#      Now, for a system with 3-buses
import Pkg 
using JuMP,SCIP,Gurobi,CPLEX,Ipopt

#-----------------------------------Define Parameters
# reactance for the SGs
# 	bus_SGs	  x	
data_SGs = [
	 1	    0.0697	
	 2	    0.0531
     3        Inf]

# reactance for the branches
# 	branch_num   fbus      tbus	       x	
data_branch = [
	   1           1	     2	     0.0250	
	   2           3	     2	     0.0213
	   3           1	     3	     0.0133]

numbranches= length(data_branch[:,1])       # number of branches
numnodes= Int(maximum(data_branch[:,2]))      # number of branches
Y_0= zeros(numnodes, numnodes)              # define admittance matrix of the transmission lines
Y_g= zeros(numnodes, numnodes)                            # define admittance matrix of the SGs' reactance

for k in 1:numbranches                 #  calculate the admittance matrix of the transmission lines
    i = Int(data_branch[k, 2])                                        # bus from
    j = Int(data_branch[k, 3])                                        # bus to
    Y_0[i, j] = -1/data_branch[k, 4]                                  # Off-diagonal elements
    Y_0[j, i] = Y_0[i, j]                                             # symmetry
end
for k in 1:numnodes
    Y_0[k, k] = -sum(Y_0[k, :])                                        # Diagonal elements 
end  

for k in 1:numnodes                     #   calculate the admittance matrix of the SGs
    i = Int(data_SGs[k, 1])                                        # bus from
    Y_g[k, k] = 1/data_SGs[k, 2]                                   # Diagonal elements 
end 

# SCC (p.u): I_SGs=E_g/X_dg  & I_IBG (pre-defined)
β=0.95
E_g=1
E_SGs=β*1
I_SGs=[E_SGs/data_SGs[1,2],E_SGs/data_SGs[2,2]]
I_IBG=0.5
Iₗᵢₘ=2    # SCC limit
#Zₘₐₓ=1

Pᴰ=[2.2,1.8,3,6,5.8,5.2,5.6,3.8,2.5,2.7,3,2.6,2.2,2.1,4.2,5.8,6.2,6.3,6.5,6.6,6.3,6.2,6,5.7]*2 # modify the cofficient to lead feasible solutions!!!!
P_wt=[2,1.5,1.6,1.8,1.3,0.6,2.8,3.3,3.9,4,3.3,2.9,2.7,2,0.2,3.2,5.1,3.1,1.8,2,1.3,1,2,3.8]*2
Pˢᴳₘₐₓ=[6,5]
Pˢᴳₘᵢₙ=[3,2]
Rₘₐₓ=[3,2]
     
Kᵁ=[3.25,1.42]
Kᴰ=[0.285,0.185]
Cᵍᵐ=[0.9,0.6]
Cⁿˡ=[1.2,1]
T=length(Pᴰ)


     
#-----------------------------------Define Model
model= Model()

#-------Define Variales
@variable(model, Pˢᴳ¹[1:T])                # generation of SG in bus 1
@variable(model, Pˢᴳ²[1:T])                # generation of SG in bus 2
@variable(model, Pᴵᴮᴳ[1:T])                # generation of IBG in bus 3
@variable(model, yˢᴳ¹[1:T],Bin)            # on/off status of SG in bus 1
@variable(model, yˢᴳ²[1:T],Bin)            # on/off status of SG in bus 2
@variable(model, Cᵁ¹[1:T])                 # startup costs for SG in bus 1
@variable(model, Cᴰ¹[1:T])                 # shutdown costs for SG in bus 1
@variable(model, Cᵁ²[1:T])                 # startup costs for SG in bus 2
@variable(model, Cᴰ²[1:T])                 # shutdown costs for SG in bus 2
@variable(model, α[1:T])                   # percentage of IBG online capacity  


#-------Define Constraints
@constraint(model, Pˢᴳ¹+Pˢᴳ²+Pᴵᴮᴳ==Pᴰ)     # power balance

@constraint(model, Pᴵᴮᴳ<=P_wt)             # wind power limit
@constraint(model, Pᴵᴮᴳ>=0)

@constraint(model, α.<=1)                   # IBG online capacity limit
@constraint(model, α>=0)
@constraint(model, α.*Pᴰ==Pᴵᴮᴳ)

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